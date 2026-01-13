// src/sample_files/parsed_file.rs
use std::fs;
use std::fs::File;
use std::io::{self, BufRead, BufReader, Read};
use std::path::{Component, Path, PathBuf};


use walkdir::WalkDir;
use std::io::Write;


#[derive(Debug, Clone)]
pub enum ParsedKind {
    TenX,
    H5,
    Fastq { lane: String, role: String },
}

#[derive(Debug, Clone)]
pub struct ParsedFile {
    pub sample: String,
    pub experiment: String,
    pub kind: ParsedKind,
    pub path: String,            // authoritative source path
    pub md5sum: Option<String>,  // computed for files; None for dirs until archived
}

impl ParsedFile {

     fn tenx_zip_path(dir: &Path) -> PathBuf {
        // put zip next to the directory, name it "<dirname>.zip"
        let parent = dir.parent().unwrap_or(dir);
        let name = dir.file_name().and_then(|s| s.to_str()).unwrap_or("tenx");
        parent.join(format!("{name}.zip"))
    }

    fn materialize_tenx_zip(dir: &Path) -> io::Result<PathBuf> {

        use zip::write::FileOptions;
        use zip::CompressionMethod;

        let opts: FileOptions<()> = FileOptions::default()
            .compression_method(CompressionMethod::Deflated)
            .unix_permissions(0o644);
        let zip_path = Self::tenx_zip_path(dir);

        // reuse if already exists and has some content
        if let Ok(md) = fs::metadata(&zip_path) {
            if md.is_file() && md.len() > 0 {
                return Ok(zip_path);
            }
        }

        // write to tmp then rename (avoid partial zips on crash)
        let tmp_path = zip_path.with_extension("zip.tmp");

        // ensure parent exists
        if let Some(par) = zip_path.parent() {
            fs::create_dir_all(par)?;
        }

        // create zip
        let f = File::create(&tmp_path)?;
        let mut zw = zip::ZipWriter::new(f);

        for entry in WalkDir::new(dir).follow_links(false).into_iter().filter_map(Result::ok) {
            let p = entry.path();

            // skip the dir itself
            if p == dir {
                continue;
            }

            let rel = p.strip_prefix(dir).unwrap_or(p);
            let rel_str = rel.to_string_lossy().replace('\\', "/"); // zip wants forward slashes

            if entry.file_type().is_dir() {
                // add directory entry (optional but fine)
                zw.add_directory(rel_str, opts)?;
            } else if entry.file_type().is_file() {
                zw.start_file(rel_str, opts)?;

                let mut rf = File::open(p)?;
                let mut buf = vec![0u8; 1024 * 1024];
                loop {
                    let n = rf.read(&mut buf)?;
                    if n == 0 { break; }
                    zw.write_all(&buf[..n])?;
                }
            }
        }

        zw.finish()?; // flush/close

        // replace old zip if present
        let _ = fs::remove_file(&zip_path);
        fs::rename(&tmp_path, &zip_path)?;

        Ok(zip_path)
    }

    fn looks_like_public_accession(fname: &str) -> bool {
        // Common run / experiment / sample / project accessions seen in public archives
        const PREFIXES: &[&str] = &[
            // SRA/ENA/DDBJ runs
            "SRR", "ERR", "DRR", "CRR",
            // experiments
            "SRX", "ERX", "DRX", "CRX",
            // samples
            "SRS", "ERS", "DRS", "CRS",
            // studies/projects
            "SRP", "ERP", "DRP", "CRP",
            "PRJNA", "PRJEB", "PRJDB",
            // biosample
            "SAMN", "SAMEA", "SAMD",
            // GEO
            "GSM", "GSE",
        ];

        let f = fname.trim();

        // Quick content markers typical for "converted" artifacts
        if f.contains(".bam.") || f.contains(".cram.") || f.contains(".sam.") || f.contains(".annotated.") {
            return true;
        }

        // Prefix + digits heuristic (avoids lots of false positives)
        for &pre in PREFIXES {
            if let Some(rest) = f.strip_prefix(pre) {
                // require at least 5 digits to avoid "SRR1" type accidental matches
                let digits: String = rest.chars().take_while(|c| c.is_ascii_digit()).collect();
                if digits.len() >= 5 {
                    return true;
                }
            }
        }

        false
    }

    /// One entrypoint: decide if path is relevant, classify, infer sample+experiment, compute md5 if file.
    pub fn from_path(scan_root: &Path, p: &Path) -> io::Result<Option<Self>> {
        let mut md = match fs::metadata(p) {
            Ok(m) => m,
            Err(e) => return Err(e),
        };

        let (effective_path ,kind) = if md.is_file() {

            let s = p.to_string_lossy();
            if Self::looks_like_public_accession( &s ) {
                // ignore public/archive-derived artifacts (SRR/ERR/DRR..., bam->fastq, annotated, etc.)
                return Ok(None);
            } else if s.ends_with(".fastq.gz") || s.ends_with(".fq.gz") {
                let (lane, role) = Self::parse_fastq_lane_role(p)?;
                ( None, ParsedKind::Fastq { lane, role })
            } else if s.ends_with(".h5") {
                (None, ParsedKind::H5)
            } else if let Some(path) = Self::tenx_triplet_dir_from_file(p) {
                if Self::looks_like_10x_triplet_dir(p)? {
                    let zip_path = Self::materialize_tenx_zip(&path)?;
                    (Some(zip_path), ParsedKind::TenX)
                }else { 
                    return Ok(None)
                }
            }else {
                return Ok(None);
            }
        }else {
            return Ok(None);
        };

        let sample = Self::detect_sample(&kind, p).ok_or_else(|| {
            io::Error::new(
                io::ErrorKind::InvalidData,
                format!("Could not infer sample for path {}", p.display()),
            )
        })?;

        let experiment = Self::detect_experiment(scan_root, &kind, p);
        let path = match effective_path {
            Some(p) => p.to_string_lossy().to_string(),
            None => p.to_string_lossy().to_string()
        };

        let mut pf = ParsedFile {
            sample,
            experiment,
            kind,
            path,
            md5sum: None,
        };

        let _ = pf.ensure_md5sum()?; // files -> Some(md5), dirs -> None
        Ok(Some(pf))
    }

    fn tenx_triplet_dir_from_file(p: &Path) -> Option<PathBuf> {
        // We only trigger on matrix.mtx.gz (best single trigger)
        let fname = p.file_name()?.to_str()?;
        if fname != "matrix.mtx.gz" {
            return None;
        }
        p.parent().map(|d| d.to_path_buf())
    }


    // ---------- path helpers ----------

    pub fn geo_filename(&self) -> String {
        format!("{}_{}", self.experiment, self.basename() )
    }
    pub fn basename(&self) -> String {
        Path::new(&self.path)
            .file_name()
            .and_then(|s| s.to_str())
            .unwrap_or(&self.path)
            .to_string()
    }

    pub fn same_basename(&self, other: &ParsedFile) -> bool {
        self.basename() == other.basename()
    }

    pub fn is_file(&self) -> bool {
        fs::metadata(&self.path).map(|m| m.is_file()).unwrap_or(false)
    }

    pub fn is_dir(&self) -> bool {
        fs::metadata(&self.path).map(|m| m.is_dir()).unwrap_or(false)
    }

    // ---------- kind detection ----------

    fn looks_like_10x_triplet_dir(dir: &Path) -> io::Result<bool> {
        // matrix triplet signature
        let mtx = dir.join("matrix.mtx.gz");
        let bar = dir.join("barcodes.tsv.gz");
        let feat = dir.join("features.tsv.gz");
        let genes = dir.join("genes.tsv.gz");
        Ok(mtx.is_file() && bar.is_file() && (feat.is_file() || genes.is_file()))
    }

    fn parse_fastq_lane_role(p: &Path) -> io::Result<(String, String)> {
        let fname = p
            .file_name()
            .and_then(|s| s.to_str())
            .ok_or_else(|| io::Error::new(io::ErrorKind::InvalidData, "Non-utf8 filename"))?;

        let lower = fname.to_ascii_lowercase();

        let role = if Self::has_token(&lower, "r1") {
            "R1"
        } else if Self::has_token(&lower, "r2") {
            "R2"
        } else if Self::has_token(&lower, "i1") {
            "I1"
        } else if Self::has_token(&lower, "i2") {
            "I2"
        } else {
            return Err(io::Error::new(
                io::ErrorKind::InvalidData,
                format!("Could not determine read role (R1/R2/I1/I2) from FASTQ name: '{fname}'"),
            ));
        }
        .to_string();

        let lane = Self::find_lane_token(fname).unwrap_or_else(|| "1".to_string());
        Ok((lane, role))
    }

    fn has_token(lower: &str, tok: &str) -> bool {
        lower.contains(&format!("_{tok}")) || lower.contains(&format!("{tok}.")) || lower.contains(&format!("{tok}_"))
    }

    fn find_lane_token(fname: &str) -> Option<String> {
        // L001 style
        let b = fname.as_bytes();
        for i in 0..b.len().saturating_sub(4) {
            if b[i] == b'L'
                && b[i + 1].is_ascii_digit()
                && b[i + 2].is_ascii_digit()
                && b[i + 3].is_ascii_digit()
            {
                return Some(fname[i..i + 4].to_string());
            }
        }
        // numeric lane like "_1_" (example3_1_R1)
        // take first "_<digits>_" occurrence
        let parts: Vec<&str> = fname.split('_').collect();
        for part in parts {
            if !part.is_empty() && part.chars().all(|c| c.is_ascii_digit()) {
                return Some(part.to_string());
            }
        }
        None
    }

    // ---------- sample detection (keep your current rules, just moved here) ----------

    fn detect_sample(kind: &ParsedKind, p: &Path) -> Option<String> {
        match kind {
            ParsedKind::Fastq { .. } => Self::sample_from_fastq_name(p),
            ParsedKind::H5 | ParsedKind::TenX => Self::sample_from_parent_dir_or_stem(p),
        }
    }

    fn sample_from_parent_dir_or_stem(p: &Path) -> Option<String> {
        // your existing logic, but fixed/cleaned: "grandparent folder name else stem"
        if let Some(parent) = p
            .parent()
            .and_then(|pp| pp.parent())
            .and_then(|gp| gp.file_name())
            .and_then(|s| s.to_str())
        {
            if !parent.is_empty() {
                return Some(parent.to_string());
            }
        }
        p.file_stem().and_then(|s| s.to_str()).map(|s| s.to_string())
    }

    fn sample_from_fastq_name(p: &Path) -> Option<String> {
        // Default: cut at first marker token
        let fname = p.file_name()?.to_str()?;
        let cut = ["_S", "_L", "_R", "_I"]
            .iter()
            .filter_map(|tok| fname.find(tok))
            .min()
            .unwrap_or_else(|| fname.find('.').unwrap_or(fname.len()));
        Some(fname[..cut].to_string())
    }

    // ---------- experiment detection ----------

    fn detect_experiment(scan_root: &Path, kind: &ParsedKind, p: &Path) -> String {
        // You can make this more sophisticated per kind later (fastq marker, outs marker).
        // For now: best-effort marker search; fallback to first component under scan_root.
        let anchor = if p.extension().is_some() { p.parent().unwrap_or(p) } else { p };

        let exp = match kind {
            ParsedKind::Fastq { .. } => Self::folder_above_marker(anchor, "fastq")
                .or_else(|| Self::first_component_under_root(scan_root, anchor)),
            ParsedKind::TenX => Self::folder_above_marker(anchor, "outs")
                .or_else(|| Self::folder_above_leaf(anchor, "filtered_feature_bc_matrix"))
                .or_else(|| Self::first_component_under_root(scan_root, anchor)),
            ParsedKind::H5 => Self::folder_above_marker(anchor, "outs")
                .or_else(|| anchor.parent()
                    .and_then(|pp| pp.file_name())
                    .and_then(|s| s.to_str())
                    .map(|s| s.to_string()))
                .or_else(|| Self::first_component_under_root(scan_root, anchor)),
        };

        exp.unwrap_or_else(|| "exp1".to_string())
    }

    fn folder_above_marker(p: &Path, marker: &str) -> Option<String> {
        let comps: Vec<String> = p
            .components()
            .filter_map(|c| match c {
                Component::Normal(os) => Some(os.to_string_lossy().to_string()),
                _ => None,
            })
            .collect();
        let idx = comps.iter().rposition(|c| c.eq_ignore_ascii_case(marker))?;
        if idx >= 1 { Some(comps[idx - 1].clone()) } else { None }
    }

    fn folder_above_leaf(p: &Path, leaf: &str) -> Option<String> {
        let mut cur = Some(p);
        while let Some(pp) = cur {
            if pp.file_name()
                .and_then(|s| s.to_str())
                .map(|s| s.eq_ignore_ascii_case(leaf)) == Some(true)
            {
                return pp.parent()
                    .and_then(|par| par.file_name())
                    .and_then(|s| s.to_str())
                    .map(|s| s.to_string());
            }
            cur = pp.parent();
        }
        None
    }

    fn first_component_under_root(scan_root: &Path, p: &Path) -> Option<String> {
        let rel = p.strip_prefix(scan_root).ok().unwrap_or(p);
        rel.components().find_map(|c| match c {
            Component::Normal(os) => Some(os.to_string_lossy().to_string()),
            _ => None,
        })
    }

    // ---------- md5 (sidecar + compute) ----------

    fn md5_sidecar_path(&self) -> PathBuf {
        // robust: foo.fastq.gz -> foo.fastq.gz.md5sum
        PathBuf::from(format!("{}.md5sum", self.path))
    }

    pub fn ensure_md5sum(&mut self) -> io::Result<Option<&str>> {
        if self.md5sum.is_some() {
            return Ok(self.md5sum.as_deref());
        }

        let p = Path::new(&self.path);
        let md = fs::metadata(p)?;
        if md.is_dir() {
            return Ok(None);
        }

        let sidecar = self.md5_sidecar_path();
        if sidecar.exists() {
            if let Ok(file) = File::open(&sidecar) {
                let mut reader = BufReader::new(file);
                let mut line = String::new();
                if reader.read_line(&mut line).is_ok() {
                    let v = line.trim().to_string();
                    if !v.is_empty() {
                        self.md5sum = Some(v);
                        return Ok(self.md5sum.as_deref());
                    }
                }
            }
        }

        let md5 = Self::compute_file_md5_incremental(p)?;
        let _ = fs::write(&sidecar, format!("{md5}\n"));
        self.md5sum = Some(md5);
        Ok(self.md5sum.as_deref())
    }

    fn compute_file_md5_incremental(file_path: &Path) -> io::Result<String> {
        let mut f = File::open(file_path)?;
        let mut ctx = md5::Context::new();
        let mut buf = vec![0u8; 1024 * 1024];
        loop {
            let n = f.read(&mut buf)?;
            if n == 0 { break; }
            ctx.consume(&buf[..n]);
        }
        Ok(format!("{:x}", ctx.compute()))
    }
}
