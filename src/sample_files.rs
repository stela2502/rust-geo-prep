use std::collections::{BTreeMap, HashSet, BTreeSet};
use std::fs::{self, File};
use std::io::{self, Read, BufRead, BufReader, BufWriter, Write};
use std::path::{Path, PathBuf};
use std::process::Command;

#[derive(Debug, Default)]
pub struct LaneFastqs {
    pub reads: BTreeMap<String, String>,
}

impl LaneFastqs {
    /// Add a FASTQ for a lane under a specific role (R1/R2/I1/...)
    pub fn add_read(&mut self, role: &str, path: String) {
        if let Some(existing) = self.reads.get(role) {
            eprintln!(
                "Duplicate read role '{}' for lane: already have '{}', tried to add '{}' - file is ignored!",
                role, existing, path
            );
        }else {
           self.reads.insert(role.to_string(), path); 
        }
        
    }

    /// Render FASTQ cells for this lane in the provided `roles` order.
    pub fn row_cells<F>(&self, roles: &[String], fmt: &F) -> Vec<String>
    where
        F: Fn(&str) -> String,
    {
        roles
            .iter()
            .map(|role| self.reads.get(role).map(|p| fmt(p)).unwrap_or_default())
            .collect()
    }
}

#[derive(Debug, Default)]
pub struct SampleRecord {
    pub name: String,

    /// 10x bundle (zip), optional
    pub tenx: Option<String>,

    /// h5 file, optional
    pub h5_files: Option<String>,

    /// FASTQ lanes grouped by lane key, each containing role→path (R1/R2/I1/...)
    pub lanes: BTreeMap<String, LaneFastqs>,
}

impl SampleRecord {
    
    pub fn fastq_source_folders(&self) -> String
    {
        let mut folders: BTreeSet<String> = BTreeSet::new();

        for (_lane_key, lane) in &self.lanes {
            for (_role, path) in &lane.reads {
                if let Some(parent) = Path::new(path).parent() {
                    folders.insert(parent.to_string_lossy().to_string());
                }
            }
        }

        folders.into_iter().collect::<Vec<_>>().join(",")
    }

    /// Render a single flattened row for this sample: Sample + TenX + H5 + (lane blocks...)
    pub fn row_cells<F>(&self, roles: &[String], fmt: &F, max_lanes: usize) -> Vec<String>
    where
        F: Fn(&str) -> String,
    {
        let mut out = Vec::new();

        // first columns
        out.push(self.fastq_source_folders());
        out.push(self.name.clone());
        out.push(self.tenx.as_ref().map(|p| fmt(p)).unwrap_or_default());
        out.push(self.h5_files.as_ref().map(|p| fmt(p)).unwrap_or_default());

        // lane blocks (sorted by key)
        let mut lane_count = 0usize;
        for (_lane_key, lane) in &self.lanes {
            out.extend(lane.row_cells(roles, fmt));
            lane_count += 1;
        }

        // pad missing lane blocks to max_lanes
        let missing_lanes = max_lanes.saturating_sub(lane_count);
        if missing_lanes > 0 {
            out.extend(std::iter::repeat(String::new()).take(missing_lanes * roles.len()));
        }

        out
    }

    /// Iterate all file paths that belong to this sample record:
    /// - TenX bundle (if any)
    /// - H5 file (if any)
    /// - all lane read files (FASTQs)
    pub fn all_paths<'a>(&'a self) -> impl Iterator<Item = &'a str> + 'a {
        let tenx = self.tenx.as_deref().into_iter();
        let h5   = self.h5_files.as_deref().into_iter();
        let fastqs = self
            .lanes
            .values()
            .flat_map(|lane| lane.reads.values())
            .map(|s| s.as_str());

        tenx.chain(h5).chain(fastqs)
    }

    /// Number of lanes
    pub fn len(&self) -> usize {
        self.lanes.len()
    }
    pub fn total_len(&self) -> usize{
        let fastq = self.len();
        let tenx  = self.tenx.iter().count();
        let h5    = self.h5_files.iter().count();

        fastq + tenx + h5
    }
}

#[derive(Debug, Clone)]
enum ParsedKind {
    TenX,
    H5,
    Fastq { lane: String, role: String },
}

#[derive(Debug, Clone)]
struct ParsedFile {
    sample: String,
    kind: ParsedKind,
    path: String, // authoritative path (possibly transformed, e.g. 10x zip)
}

#[derive(Debug)]
pub struct SampleFiles {
    pub omit_md5: bool,
    pub samples: BTreeMap<String, SampleRecord>,

    processed_paths: HashSet<String>, // canonical paths AFTER parsing
    roles: HashSet<String>,           // FASTQ roles seen: R1, R2, I1, ...

    // kept for md5 reporting + deterministic lists
    files: Vec<(String, String, String)>, // (path, label, md5)
}

impl SampleFiles {
    /// Construct the container
    pub fn new(omit_md5: bool) -> Self {
        SampleFiles {
            omit_md5,
            samples: BTreeMap::new(),
            processed_paths: HashSet::new(),
            roles: HashSet::new(),
            files: Vec::new(),
        }
    }

    /// Total number of recorded samples
    pub fn samples(&self) -> usize {
        self.samples.len()
    }

    /// Total number of recorded files (including 10x/h5/fastq)
    pub fn len(&self) -> usize {
        self.files.len()
    }
    


    /// Write sample table (full paths)
    pub fn write_sample_files(&self, path: &str) -> io::Result<()> {
        self.write_sample_files_with(path, '\t', |p| p.to_string())
    }

    /// Write sample table (basenames only)
    pub fn write_sample_files_basename(&self, path: &str) -> io::Result<()> {
        self.write_sample_files_with(path, '\t', |p| self.extract_basename(p).unwrap_or_else(|| p.to_string()))
    }

    /// Internal sample-table writer
    fn write_sample_files_with<F>(&self, path: &str, sep: char, fmt: F) -> io::Result<()>
    where
        F: Fn(&str) -> String,
    {
        let file = File::create(path)?;
        let mut w = BufWriter::new(file);

        // FASTQ roles columns
        let mut roles: Vec<String> = self.roles.iter().cloned().collect();
        roles.sort();

        let max_lanes = self.samples.values().map(|s| s.len()).max().unwrap_or(0);

        // header
        write!(w, "Source_Path(s){sep}Sample_Lane{sep}TenX{sep}H5")?;
        for _block in 0..max_lanes {
            for role in &roles {
                write!(w, "{sep}{role}")?;
            }
        }
        writeln!(w)?;

        // rows (one per sample, flattened lane blocks)
        for sample in self.samples.values() {
            let row = sample.row_cells(&roles, &fmt, max_lanes);
            writeln!(w, "{}", row.join(&sep.to_string()))?;
        }

        Ok(())
    }

    /// Full paths + md5
    pub fn write_md5_files(&self, path: &str) -> io::Result<()> {
        self.write_md5_files_with(path, |p| p.to_string())
    }

    /// Basename only + md5
    pub fn write_md5_files_basename(&self, path: &str) -> io::Result<()> {
        self.write_md5_files_with(path, |p| self.extract_basename(p).unwrap_or_else(|| p.to_string()))
    }

    /// Internal md5 writer
    fn write_md5_files_with<F>(&self, path: &str, fmt: F) -> io::Result<()>
    where
        F: Fn(&str) -> String,
    {
        let file = File::create(path)?;
        let mut w = BufWriter::new(file);

        writeln!(w, "file_name\tmd5sum")?;

        let mut entries = self.files.clone();
        entries.sort_by(|a, b| a.0.cmp(&b.0));

        for (file_path, _label, md5sum) in entries {
            writeln!(w, "{}\t{}", fmt(&file_path), md5sum)?;
        }

        Ok(())
    }

    /// Reject public-repo downloaded / re-packed FASTQs (GEO/SRA/ENA/ArrayExpress-ish).
    fn looks_like_public_repo_dump(&self, file_path: &str) -> bool {
        let fname = match self.extract_basename(file_path) {
            Some(b) => b,
            None => return false,
        };
        let lower = fname.to_ascii_lowercase();
        let path_lower = file_path.to_ascii_lowercase();

        // GEO / SRA / ENA / DDBJ typical run accessions in filenames:
        // SRR/ERR/DRR + digits, also sometimes SRX/ERX/DRX, GSM/GSE
        let starts_with_run = |p: &str| {
            let u = fname.as_bytes();
            let pu = p.as_bytes();
            if u.len() < pu.len() + 6 { return false; } // require some digits
            if !fname.starts_with(p) { return false; }
            fname[p.len()..].chars().take(12).any(|c| c.is_ascii_digit())
        };

        // only care about fastqs
        if !(lower.ends_with(".fastq.gz") || lower.ends_with(".fq.gz") || lower.ends_with(".fastq") || lower.ends_with(".fq")) {
            return false;
        }

        

        if starts_with_run("SRR") || starts_with_run("ERR") || starts_with_run("DRR")
            || starts_with_run("SRX") || starts_with_run("ERX") || starts_with_run("DRX")
            || fname.starts_with("GSM") || fname.starts_with("GSE")
        {
            return true;
        }

        // Common “repacked from bam/sra” hints
        if lower.contains(".bam.") || lower.contains("annotated") || lower.contains("sra") {
            return true;
        }

        // Folder hints (GEO series / ArrayExpress / ENA style staging)
        if path_lower.contains("/geo/") || path_lower.contains("/gse") || path_lower.contains("/gsm")
            || path_lower.contains("/arrayexpress/") || path_lower.contains("/ena/") || path_lower.contains("/sra/")
        {
            return true;
        }

        false
    }



    /// Add a file into the NEW SampleFiles structure (samples → lanes → role→path; plus 10x/h5 on SampleRecord).
    pub fn add_file(&mut self, file_path: &str) {
        // Step 1: ignore junk FASTQs
        if let Some(basename) = self.extract_basename(file_path) {
            if basename.starts_with("Undetermined")
                || basename.starts_with("Unmapped")
                || basename.starts_with("Umapped")
            {
                return;
            }
        }

        if self.looks_like_public_repo_dump( &file_path ){
            eprintln!("Looks like public data - re-submitting is not a good idea and it will break the logics later on:\n{file_path}");
            return;
        }

        // Step 2: parse filename (may transform path e.g. 10x → zip)
        let parsed = match self.parse_file(file_path) {
            Some(v) => v,
            None => return,
        };

        // canonicalize the authoritative path (not the incoming one)
        let canonical = match std::fs::canonicalize(&parsed.path) {
            Ok(p) => p.to_string_lossy().to_string(),
            Err(_) => parsed.path.clone(),
        };

        // Deduplicate after parsing
        if !self.processed_paths.insert(canonical) {
            return;
        }

        // Step 3: md5 of the authoritative file
        let md5sum = self.get_md5sum(&parsed.path);

        // Step 4: record file for md5 report
        let label = match &parsed.kind {
            ParsedKind::TenX => "10x".to_string(),
            ParsedKind::H5 => "h5".to_string(),
            ParsedKind::Fastq { lane, role } => format!("{lane}:{role}"),
        };
        eprintln!("I obtained a usable file: {}, {}, with {}", parsed.path, label, md5sum);
        self.files.push((parsed.path.clone(), label, md5sum.clone()));

        // Step 5: update NEW structure
        let rec = self
            .samples
            .entry(parsed.sample.clone())
            .or_insert_with(|| SampleRecord {
                name: parsed.sample.clone(),
                ..Default::default()
            });

        

        match parsed.kind {
            ParsedKind::TenX => {
                if rec.tenx.is_some() {
                    // keep "true to the idea": 10x is single bundle per sample
                    panic!("Duplicate 10x bundle for sample '{}': \n{}\n{:?}", rec.name, parsed.path, rec.tenx);
                }
                rec.tenx = Some(parsed.path);
            }
            ParsedKind::H5 => {
                if let Some(existing) = rec.h5_files.as_ref() {
                    if existing == &parsed.path {
                        return; // identical duplicate → ignore
                    }
                    eprintln!(
                        "WARNING: Additional h5 for sample '{}': keeping first\n  keep: {}\n  skip: {}",
                        rec.name, existing, parsed.path
                    );
                    return; // keep first (ignore all secondaries)
                }
                rec.h5_files = Some(parsed.path);
            }
            ParsedKind::Fastq { lane, role } => {
                self.roles.insert(role.clone());
                rec.lanes.entry(lane).or_default().add_read(&role, parsed.path);
            }
        }
    }

    /// Parse a path into (sample, kind, authoritative_path). Handles 10x triplets, h5, and FASTQs.
    fn parse_file(&self, file_path: &str) -> Option<ParsedFile> {

        if should_ignore_fastq_by_name(&file_path) {
            eprintln!("WARNING: Ignoring FASTQ with non-Illumina naming: '{file_path}'");
            return None
        }

        // 10x triplets → <sample>.10x.zip
        if let Some((sample, _tag, zip_path)) = self.parse_matrix_triplets(file_path) {
            return Some(ParsedFile {
                sample,
                kind: ParsedKind::TenX,
                path: zip_path,
            });
        }

        // h5
        if file_path.ends_with(".h5") {
            let sample = self.sample_from_parent_dir_or_stem(file_path)?;
            let new_path = match self.materialize_prefixed_h5(file_path, &sample) {
                Ok(p) => p,
                Err(e) => {
                    eprintln!(
                        "WARNING: could not materialize prefixed h5 for {} ({e:?}); using original",
                        file_path
                    );
                    file_path.to_string()
                }
            };
            return Some(ParsedFile {
                sample,
                kind: ParsedKind::H5,
                path: new_path.to_string(),
            });
        }

        // FASTQs
        self.parse_fastq(file_path)
    }

    /// Parse FASTQ names into sample + lane_key + role (R1/R2/I1) + original path.
    fn parse_fastq(&self, file_path: &str) -> Option<ParsedFile> {
        let file_name = Path::new(file_path).file_name()?.to_str()?.to_string();

        if !file_name.ends_with(".fastq.gz") && ! file_name.ends_with(".fq.gz") {
            return None
        }
        if should_ignore_fastq_by_name(&file_name) {
            eprintln!("WARNING: Ignoring FASTQ with non-Illumina naming: '{file_name}'");
            return None
        }

        let stem = file_name
            .strip_suffix(".fastq.gz")
            .or_else(|| file_name.strip_suffix(".fq.gz"))
            .unwrap_or(&file_name);

        let parts: Vec<&str> = stem.split('_').collect();
        if parts.is_empty() {
            return None;
        }



        // 1) Find the read role by scanning from the BACK
        let mut role: Option<String> = None;
        let mut stop_idx = parts.len();

        for (i, p) in parts.iter().enumerate().rev() {
            if *p == "R1" || *p == "R2" || *p == "I1" || *p == "I2" {
                role = Some((*p).to_string());
                stop_idx = i;
                break;
            }
            if p.starts_with('R') && p.len() > 1 && p[1..].chars().all(|c| c.is_ascii_digit()) {
                role = Some((*p).to_string());
                stop_idx = i;
                break;
            }
            if p.starts_with('I') && p.len() > 1 && p[1..].chars().all(|c| c.is_ascii_digit()) {
                role = Some((*p).to_string());
                stop_idx = i;
                break;
            }
        }
        // Fallback ONLY if the LAST token is numeric and no explicit role was found
        if role.is_none() {
            if let Some((i, last)) = parts.iter().enumerate().last() {
                if *last == "1" {
                    role = Some("R1".to_string());
                    stop_idx = i;
                } else if *last == "2" {
                    role = Some("R2".to_string());
                    stop_idx = i;
                }
            }
        }

        // ❌ NO silent fallback anymore — missing role is a hard error
        let role = match role {
            Some(r) => r,
            None => {
                eprintln!(
                "Could not determine read role (R1/R2/I1/I2) from FASTQ name: '{}' -> assuming single end and call it R1",
                file_path
                );
                "R1".to_string()
            },
        };



        let mut sample_parts: Vec<String> = Vec::new();
        let mut s_part: Option<String> = None;
        let mut l_part: Option<String> = None;

        for p in &parts[..stop_idx] {
            // lane-ish
            if p.starts_with('S') && p.len() > 1 && p[1..].chars().all(|c| c.is_ascii_digit()) {
                s_part = Some((*p).to_string());
                continue;
            }
            if p.starts_with('L') && p.len() > 1 && p[1..].chars().all(|c| c.is_ascii_digit()) {
                l_part = Some((*p).to_string());
                continue;
            }

            // numeric token like "_1_" → lane index (NOT a read role)
            if p.chars().all(|c| c.is_ascii_digit()) {
                l_part = Some((*p).to_string());
                continue;
            }

            // otherwise part of sample name
            sample_parts.push((*p).to_string());
        }
        let sample =  if sample_parts.is_empty(){
            self.sample_from_parent_dir_or_stem(file_path)
                .unwrap_or_else(|| {
                    panic!(
                        "Could not determine sample name from path or filename: '{}'",
                        file_path
                    )
                })
        }else {
            sample_parts.join("_")
        };
        let lane = format!("{:?}_{:?}", s_part, l_part );

        Some(ParsedFile {
            sample,
            kind: ParsedKind::Fastq { lane, role },
            path: file_path.to_string(),
        })
    }


    /// Create a prefixed h5 as <sample>_<original>.h5.
    /// Accepts existing file or existing hardlink; otherwise tries hardlink, falls back to copy.
    fn materialize_prefixed_h5(&self, file_path: &str, sample: &str) -> io::Result<String> {
        let src = Path::new(file_path);

        let fname = self
            .extract_basename(file_path)
            .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidInput, "h5 has no basename"))?;

        let dst = src
            .parent()
            .map(|p| p.join(format!("{sample}_{fname}")))
            .unwrap_or_else(|| PathBuf::from(format!("{sample}_{fname}")));

        // ✅ Accept if destination already exists (idempotent)
        if dst.exists() {
            return Ok(dst.to_string_lossy().to_string());
        }

        // Try hard link first (cheap, no space)
        match fs::hard_link(src, &dst) {
            Ok(()) => return Ok(dst.to_string_lossy().to_string()),
            Err(e) if e.kind() == io::ErrorKind::AlreadyExists => {
                // race condition / parallel runs
                return Ok(dst.to_string_lossy().to_string());
            }
            Err(e) => {
                // EXDEV or permissions → fallback to copy
                eprintln!(
                    "WARNING: hard_link failed for {} -> {} ({e:?}); falling back to copy",
                    src.display(),
                    dst.display()
                );
            }
        }

        fs::copy(src, &dst)?;
        Ok(dst.to_string_lossy().to_string())
    }

    fn sample_from_parent_dir_or_stem(&self, file_path: &str) -> Option<String> {
        let p = Path::new(file_path);
        if let Some(parent) = p.parent()
            .and_then(|pp| pp.parent())
            .and_then(|gp| gp.file_name())
            .and_then(|s| s.to_str()) {
            if !parent.is_empty() {
                return Some(parent.to_string());
            }
        }
        p.file_stem().and_then(|s| s.to_str()).map(|s| s.to_string())
    }

    /// Extract basename
    fn extract_basename(&self, file_path: &str) -> Option<String> {
        Path::new(file_path)
            .file_name()
            .and_then(|name| name.to_str())
            .map(|s| s.to_string())
    }

    /// When any of the 10x triplet files is encountered, zip them into <sample>.10x.zip and return it.
    fn handle_10x_triplet(&self, file_path: &str) -> Option<(String, PathBuf)> {
        let path = Path::new(file_path);
        let fname = path.file_name()?.to_str()?;

        if !matches!(fname, "matrix.mtx.gz" | "features.tsv.gz" | "barcodes.tsv.gz") {
            return None;
        }

        // directory containing the triplet
        let level3 = path.parent()?;
        let level3_name = level3.file_name()?.to_string_lossy();

        // check triplet completeness
        let has_triplet = level3.join("matrix.mtx.gz").exists()
            && level3.join("features.tsv.gz").exists()
            && level3.join("barcodes.tsv.gz").exists();

        if !has_triplet {
            return None;
        }

        // parent dir
        let level2 = level3.parent()?;
        let level2_name = level2.file_name()?.to_string_lossy();

        // parent of parent
        let level1 = level2.parent()?;
        let level1_name = level1.file_name()?.to_string_lossy().to_string();

        let sample_name = if (level3_name.ends_with("feature_matrix") || level3_name.ends_with("feature_bc_matrix"))
            && level2_name == "outs"
        {
            level1_name
        } else if level2_name.starts_with("out") {
            level1_name
        } else if level3_name.ends_with("feature_matrix")
            || level3_name.ends_with("feature_bc_matrix")
            || level3_name.ends_with("filtered_feature_bc_matrix")
            || level3_name.ends_with("filtered_files")
        {
            // GEO-style: sample/filtered_feature_bc_matrix
            level2_name.to_string()
        } else {
            eprintln!(
                "WARNING: Found 10x triplet at '{}' but directory structure does not match any rule → ignoring.",
                level3.display()
            );
            return None;
        };

        let zip_path = level2.join(format!("{sample_name}.10x.zip"));
        if zip_path.exists() {
            return Some((sample_name, zip_path));
        }

        let status = Command::new("zip")
            .arg("-r")
            .arg(&zip_path)
            .arg(&level3)
            .status()
            .ok()?;

        if !status.success() {
            eprintln!("Failed zipping {}", level3.display());
            return None;
        }

        Some((sample_name, zip_path))
    }

    /// Detect and convert 10x triplets to a single (sample, "10x", zip_path) record.
    fn parse_matrix_triplets(&self, file_path: &str) -> Option<(String, String, String)> {
        if let Some((sample_name, zip_path)) = self.handle_10x_triplet(file_path) {
            return Some((
                sample_name,
                "10x".to_string(),
                zip_path.to_string_lossy().to_string(),
            ));
        }
        None
    }

    fn compute_file_md5_incremental(&self, file_path: &str) -> io::Result<String> {
        if self.omit_md5 {
            return Ok("omitted".to_string());
        }

        let mut f = File::open(file_path)?;
        let mut ctx = md5::Context::new();

        // HEAP allocation, not stack
        let mut buf = vec![0u8; 1024 * 1024]; // 1 MiB buffer on heap
        loop {
            let n = f.read(&mut buf)?;
            if n == 0 {
                break;
            }
            ctx.consume(&buf[..n]);
        }

        Ok(format!("{:x}", ctx.compute()))
    }

    fn get_md5sum(&self, file_path: &str) -> String {
        let path = Path::new(file_path);

        // keep your old convention
        let md5_file = if file_path.ends_with(".fastq.gz") {
            path.with_extension("fastq.gz.md5sum")
        } else if file_path.ends_with(".fq.gz") {
            path.with_extension("fq.gz.md5sum")
        } else {
            path.with_extension("md5sum")
        };

        if md5_file.exists() {
            if let Ok(file) = File::open(&md5_file) {
                let reader = BufReader::new(file);
                if let Some(Ok(line)) = reader.lines().next() {
                    return line;
                }
            }
        }

        if let Ok(md5sum) = self.compute_file_md5_incremental(file_path) {
            let _ = fs::write(&md5_file, &md5sum);
            return md5sum;
        }

        "none".to_string()
    }

    pub fn write_collect_all_files_script_sh<P: AsRef<Path>>(&self, path: P, dest:&str ) -> std::io::Result<()> {
        let f = File::create(path)?;
        let mut w = BufWriter::new(f);

        writeln!(w, "#!/usr/bin/env bash")?;
        writeln!(w, "set -euo pipefail")?;
        writeln!(w)?;
        writeln!(w, "# Collect ALL files referenced by rust-geo-prep outputs into one folder.")?;
        writeln!(w, "# Change copy tool, e.g.: COPY_CMD=(rsync --progress)   or   COPY_CMD=(cp -v)")?;
        writeln!(w, "COPY_CMD=(cp -v)")?;
        writeln!(w, "DEST=\"{}\"", dest)?;
        writeln!(w, "mkdir -p \"$DEST\"")?;
        writeln!(w)?;

        // Iterate samples in stable order (BTreeMap iteration is sorted already).
        for sample in self.samples.values() {
            writeln!(w, "###############################################################################")?;
            writeln!(w, "# Sample: {}", sample.name)?;
            writeln!(w, "###############################################################################")?;

            // De-duplicate within a sample and keep stable order
            let mut uniq: BTreeSet<String> = BTreeSet::new();
            for p in sample.all_paths() {
                uniq.insert(p.to_string());
            }

            for src in uniq {
                let p = std::path::Path::new(&src);
                let base = p
                    .file_name()
                    .map(|x| x.to_string_lossy().to_string())
                    .unwrap_or_else(|| "unknown.file".to_string());

                // Collision-resistant dest name
                let dest_name =  base.clone();

                // Escape quotes for bash
                let src_esc = src.replace('"', "\\\"");
                let dest_esc = dest_name.replace('"', "\\\"");

                writeln!(w, "\"${{COPY_CMD[@]}}\" \"{}\" \"$DEST/{}\"", src_esc, dest_esc)?;
            }

            writeln!(w)?;
        }

        Ok(())
    }

    pub fn write_collect_all_files_script_ps1<P: AsRef<Path>>(
        &self,
        path: P,
        dest: &str,
    ) -> std::io::Result<()> {
        use std::io::Write;
        let mut w = std::io::BufWriter::new(std::fs::File::create(path)?);

        writeln!(w, "# PowerShell collection script generated by rust-geo-prep")?;
        writeln!(w, "# Copies ALL referenced files into a single folder WITHOUT renaming.")?;
        writeln!(w, "# Filenames remain identical to md5/sample tables.")?;
        writeln!(w)?;
        writeln!(w, "$ErrorActionPreference = \"Stop\"")?;
        writeln!(w)?;
        writeln!(w, "$DEST = \"{}\"", dest)?;
        writeln!(w, "New-Item -ItemType Directory -Force -Path $DEST | Out-Null")?;
        writeln!(w)?;
        writeln!(w, "# Change copy tool here if needed")?;
        writeln!(w, "$COPY = {{ param($src,$dst) Copy-Item -LiteralPath $src -Destination $dst -Force }}")?;
        writeln!(w)?;

        for (sample, rec) in &self.samples {

            writeln!(w)?;
            writeln!(w, "###############################################################################")?;
            writeln!(w, "# Sample: {}", sample)?;
            writeln!(w, "###############################################################################")?;

            // collect ALL file paths belonging to this sample
            let mut files = Vec::new();

            if let Some(p) = &rec.tenx {
                files.push(p.clone());
            }
            if let Some(p) = &rec.h5_files {
                files.push(p.clone());
            }

            for lane in rec.lanes.values() {
                for p in lane.reads.values() {
                    files.push(p.clone());
                }
            }

            for src in files {
                let basename = std::path::Path::new(&src)
                    .file_name()
                    .unwrap()
                    .to_string_lossy()
                    .to_string();

                writeln!(
                    w,
                    "& $COPY \"{}\" (Join-Path $DEST \"{}\")",
                    src, basename
                )?;
            }
        }

        Ok(())
    }
}



fn looks_like_sra_accession(basename: &str) -> bool {
    // SRR/ERR/DRR + digits (common SRA run accessions)
    // Example: SRR6333601
    let b = basename.as_bytes();
    if b.len() < 4 { return false; }
    let prefix = &basename[..3];
    if prefix != "SRR" && prefix != "ERR" && prefix != "DRR"  && prefix != "GSE" { return false; }
    basename[3..].chars().all(|c| c.is_ascii_digit())
}

fn has_explicit_read_token(name: &str) -> bool {
    // covers typical conventions: _R1_, _R2_, _I1_, _I2_, _R1., _R2., _I1., _I2.
    let n = name;
    n.contains("_R1") || n.contains("_R2") || n.contains("_I1") || n.contains("_I2")
}

fn should_ignore_fastq_by_name(file_name: &str) -> bool {
    // Most direct: your exact failing pattern
    if file_name.contains(".bam.") || file_name.contains(".bam.annotated.") {
        return true;
    }

    if file_name.starts_with("test") {
        return true;
    }

    // More general: SRR/ERR/DRR-run FASTQs without explicit read role token
    // (these frequently won't match lane-style naming)
    let stem = file_name
        .trim_end_matches(".fastq.gz")
        .trim_end_matches(".fq.gz");

    looks_like_sra_accession(stem) && !has_explicit_read_token(file_name)
}

