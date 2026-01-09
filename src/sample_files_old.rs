use std::collections::HashMap;
use std::process::Command;
use std::fs::{self, File};
use std::io::{self, BufRead, BufReader, Write};
use std::path::Path;
use std::path::PathBuf;
use std::collections::HashSet;
use std::io::BufWriter;

use std::collections::BTreeMap;

#[derive(Debug, Default)]
pub struct LaneFastqs {
    pub reads: BTreeMap<String, String>,
}

impl LaneFastqs {
    pub fn add_read(&mut self, role: &str, path: String) {
        if let Some(existing) = self.reads.get(role) {
            panic!(
                "Duplicate read role '{}' for lane: already have '{}', tried to add '{}'",
                role, existing, path
            );
        }
        self.reads.insert(role.to_string(), path);
    }

    pub fn row_cells<F>(
        &self,
        roles: &[String],
        fmt: F,
    ) -> Vec<String>
    where
        F: Fn(&str) -> String,
    {
        roles
            .iter()
            .map(|role| {
                self.reads
                    .get(role)
                    .map(|p| fmt(p))
                    .unwrap_or_default()
            })
            .collect()
    }

}

#[derive(Debug, Default)]
pub struct SampleRecord {
    pub name: String,

    // 10x: you might have zero, one, or several bundles
    pub tenx: Option<String>,

    // h5: zero or more h5 files
    pub h5_files: Option<String>,

    // FASTQ lanes grouped by "technicality" WITHOUT R1/R2/I1
    // e.g. "S1_L001", "1", "L001"
    pub lanes: BTreeMap<String, LaneFastqs>,
}


impl SampleRecord{
    pub fn row_cells<F>(&self, roles: &[String], fmt: F )
    where
        F: Fn(&str) -> String,
    {
        let mut rows = Vec::new();

        // Precompute 10x + h5 cells (reused across lane rows)
        let tenx_cell = self.tenx.as_ref().map(|p| fmt(p)).unwrap_or_default();
        let h5_cell   = self.h5.as_ref().map(|p| fmt(p)).unwrap_or_default();
        rows.push( tenx_cell );
        rows.push( h5_cell );

         // --- one row per lane ---
        for (tech, lane) in &self.lanes {
            let mut cells: Vec<String> = Vec::new();
            rows.extend(lane.row_cells(fastq_roles, fmt));
        }

        rows
    }
    pub fn len(&self) -> usize{
        self.lanes.len()
    }
}


#[derive(Debug)]
pub struct SampleFiles {
    pub omit_md5: bool,
    pub samples: BTreeMap<String, SampleRecord>,
    processed_paths: HashSet<String>,   // canonical paths AFTER parsing
    types: HashSet<String>, // all roles seen: R1, R2, I1, MTX, ...
}

impl SampleFiles {
    pub fn new( omit_md5:bool ) -> Self {
        SampleFiles {
            omit_md5,
            samples: BTreeMap::new(),
            processed_paths: HashSet::new,
            types: HashSet::new(),
        }
    }

    pub fn len(&self) -> usize{
        self.filenames.len()
    }


    fn samples(&self) -> Vec<&String> {
        let mut sample_keys: Vec<&String> = self.filenames_by_sample.keys().collect();
        sample_keys.sort();
        sample_keys

    }


    /// Internal generic writer. `fmt` decides how to render the filename (full path, basename, etc.).
    fn write_sample_files_with<F>(&self, path: &str, fmt: F)
    where
        F: Fn(&str) -> String,
    {
        let file = File::create(path).unwrap_or_else(|e| {
            panic!("Could not create sample file {}:\n{}", path, e)
        });
        let mut w = BufWriter::new(file);

        // ------------------------------------------------------------
        // 1) Columns: from `self.types`
        // ------------------------------------------------------------
        let mut roles: Vec<String> = self.types.iter().cloned().collect();
        roles.sort(); // deterministic order

        let max_lanes = self
            .samples
            .values()
            .map(|s| s.len())
            .max()
            .unwrap_or(0);

        write!(w, "Sample_Lane").unwrap();
        write!(w, "{}TenX{}H5", sep as char, sep as char)?;

        for _block in 0..max_lanes {
            for role in &roles {
                write!(w, "{}{}", sep as char, role)?;
            }
        }
        writeln!(w)?;

        let expected = 1 + 2 + roles.len() * max_lanes;

        for sample in self.samples.values() {
            let mut row vec<String> = sample.row_cells( &roles, fmt );
            let missing = expected_cols.saturating_sub(row.len());
            row.extend(std::iter::repeat(String::new()).take(missing));
            w.writeln!(w, row.join('\t'))?;
        }
        Ok(())
    }

    /// Old behaviour: full paths
    pub fn write_sample_files(&self, path: &str) {
        self.write_sample_files_with(path, |p| p.to_string());
    }

    /// New variant: just basenames of the files
    pub fn write_sample_files_basename(&self, path: &str) {
        self.write_sample_files_with(path, |p| {
            self.extract_basename(p).unwrap_or_else(|| p.to_string())
        });
    }

    // Internal reusable writer: `fmt` decides how to render the file name.
    fn write_md5_files_with<F>(&self, path: &str, fmt: F) -> io::Result<()>
    where
        F: Fn(&str) -> String,
    {
        let file = File::create(path)?;
        let mut w = BufWriter::new(file);

        writeln!(w, "file_name\tmd5sum")?;

        // If you want deterministic order, sort by path:
        let mut entries: Vec<_> = self.filenames.iter().collect();
        entries.sort_by(|(p1, _, _), (p2, _, _)| p1.cmp(p2));

        for (file_path, _role, md5sum) in entries {
            let printed = fmt(file_path);
            writeln!(w, "{}\t{}", printed, md5sum)?;
        }

        Ok(())
    }

    /// Full paths + md5
    pub fn write_md5_files(&self, path: &str) -> io::Result<()> {
        self.write_md5_files_with(path, |p| p.to_string())
    }

    /// Basename only + md5
    pub fn write_md5_files_basename(&self, path: &str) -> io::Result<()> {
        self.write_md5_files_with(path, |p| {
            self.extract_basename(p).unwrap_or_else(|| p.to_string())
        })
    }

    /// Add a file into the SampleFiles structure.
    ///
    /// This function performs several steps in the correct order:
    /// 1. Reject unwanted files such as Undetermined/Unmapped FASTQs.
    /// 2. Parse the file name using `parse_filename_split()`, which:
    ///        - Detects and handles 10x CellRanger triplet folders.
    ///        - Possibly replaces the incoming path with a generated <sample>.10x.zip.
    ///        - Extracts sample name and technicalities for FASTQs.
    /// 3. Only AFTER parsing do we push the (possibly transformed) file into
    ///    the filenames list and compute md5.
    /// 4. Update the sample and sample+tech maps with the correct file index.
    ///
    /// This fixes the old broken logic:
    ///    - Previously md5 and filenames were created BEFORE parsing,
    ///      so 10x triplets stored invalid entries.
    ///    - Now the authoritative filename always comes from the parser first.
    pub fn add_file(&mut self, file_path: &str) {
        // Step 1: ignore junk FASTQs
        if let Some(basename) = self.extract_basename(file_path) {
            if basename.starts_with("Undetermined")
                || basename.starts_with("Unmapped")
                || basename.starts_with("Umapped")
            {
                return; // Not a useful file
            }
        }

        // Step 2: parse filename (may modify path e.g. for 10x → zip)
        let parsed = match self.parse_filename_split(file_path) {
            Some(v) => v,
            None => return, // Could not classify file → ignore gracefully
        };

        // parsed = (sample, technicalities, updated_path)
        let (sample, technicalities, updated_path) = parsed;

        // Canonicalize the final path
        let canonical = match std::fs::canonicalize(&file_path) {
            Ok(p) => p.to_string_lossy().to_string(),
            Err(_) => file_path.to_string(), // fallback
        };

        // NOW DO DEDUPLICATION INSIDE SampleFiles
        if !self.processed_paths.insert(canonical.clone()) {
            // Already seen — skip
            return;
        }

        // Step 3: Now and ONLY now compute md5 of the real file
        let md5sum = self.get_md5sum(&updated_path);

        // Step 4: Add the authoritative file path to main list
        let index = self.filenames.len();
        // keep track of all roles seen
        self.types.insert(technicalities.clone());
        eprintln!("I have found {} {} with {}", updated_path, technicalities, md5sum);

        self.filenames.push((updated_path.clone(), technicalities.clone(), md5sum));

        // Step 5: Update sample lookups
        self.filenames_by_sample
            .entry(sample.clone())
            .or_insert_with(Vec::new)
            .push(index);

        self.filenames_by_sample_tech
            .entry((sample, technicalities))
            .or_insert_with(Vec::new)
            .push(index);
    }

    fn parse_filename_split(&self, file_path: &str)
        -> Option<(String, String, String)>
    {
        // First: detect 10x – this may return a zip file
        if let Some(t) = self.parse_matrix_triplets(file_path) {
            return Some(t);
        }

        // FASTQs: same logic as before, but last field is the *original file path*
        let file_name = file_path.split('/').last()?;
        let file_name = file_name
            .strip_suffix(".fastq.gz")
            .or_else(|| file_name.strip_suffix(".fq.gz"))
            .unwrap_or(file_name);

        let parts: Vec<&str> = file_name.split('_').collect();
        if parts.is_empty() {
            return None;
        }

        let mut sample_parts = Vec::new();
        let mut tech_parts = Vec::new();

        for part in parts {
            if part.starts_with('S') && part[1..].chars().all(|c| c.is_digit(10)) {
                tech_parts.push(part.to_string());
            } else if part.starts_with('L') && part[1..].chars().all(|c| c.is_digit(10)) {
                tech_parts.push(part.to_string());
            } else if part.starts_with('R') || part == "1" || part == "2" || part.starts_with("I1") {
                // read info ignored
                tech_parts.push(part.to_string());
                break;
            } else {
                sample_parts.push(part.to_string());
            }
        }

        let sample = sample_parts.join("_");
        let technicalities = tech_parts.join("_");

        Some((sample, technicalities, file_path.to_string()))
    }



    /// When any of the 10x triplet files is encountered, check if the full
    /// triplet is present in that directory. If yes, zip them into
    /// <sample_name>.zip and return the zip path.
    ///
    /// sample_name is derived from the directory name that contains the triplet.
    fn handle_10x_triplet(&self, file_path: &str) -> Option<(String, PathBuf)> {
        use std::path::{Path, PathBuf};

        let path = Path::new(file_path);
        let fname = path.file_name()?.to_str()?;

        // Only trigger on the triplet files
        if !matches!(fname, "matrix.mtx.gz" | "features.tsv.gz" | "barcodes.tsv.gz") {
            return None;
        }

        // LEVEL 3 = directory containing the triplet files
        let level3 = path.parent()?;
        let level3_name = level3.file_name()?.to_string_lossy();

        // Check triplet completeness
        let has_triplet = level3.join("matrix.mtx.gz").exists()
            && level3.join("features.tsv.gz").exists()
            && level3.join("barcodes.tsv.gz").exists();

        if !has_triplet {
            return None;
        }

        // LEVEL 2 = parent directory
        let level2 = level3.parent()?;
        let level2_name = level2.file_name()?.to_string_lossy();

        // LEVEL 1 = parent of LEVEL 2
        let level1 = level2.parent()?;
        let level1_name = level1.file_name()?.to_string_lossy().to_string();

        // ---------------------------
        // Determine sample name (rules A + B)
        // ---------------------------
        let sample_name = if 
            // RULE A:
            // <sample>/outs/<filtered_feature_bc_matrix>/matrix.mtx.gz
            level3_name.ends_with("feature_matrix") || level3_name.ends_with("feature_bc_matrix")
                && level2_name == "outs"
        {
            level1_name
        }
        else if 
            // RULE B:
            // <sample>/outs/matrix.mtx.gz
            level2_name.starts_with("out")
        {
            level1_name
        }
        // ⭐ RULE C: GEO-style: sample/filtered_feature_bc_matrix (no outs at all)
        else if level3_name.ends_with("feature_matrix")
            || level3_name.ends_with("feature_bc_matrix")
            || level3_name.ends_with("filtered_feature_bc_matrix")
            || level3_name.ends_with("filtered_files")
        {
            // level2 is the sample
            level2_name.to_string()
        }
        else {
            eprintln!(
                "WARNING: Found 10x triplet at '{}' but directory structure \
                 does not match any rule → ignoring.",
                level3.display()
            );
            return None;   // --- EARLY EXIT ---
        };

        // ---------------------------
        // Zip creation or reuse (simple unified branch)
        // ---------------------------
        let zip_path = level2.join(format!("{sample_name}.10x.zip"));

        // If already exists, do not zip again
        if zip_path.exists() {
            return Some((sample_name, zip_path));
        }

        // Create zip from LEVEL 3 only
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

    /// fix the barcodes.tsv.gz, features.tsv.gz and matrix.mtx.gz files
    fn parse_matrix_triplets(&self, file_path: &str)
        -> Option<(String, String, String)>
    {
        if let Some((sample_name, zip_path)) = self.handle_10x_triplet(file_path) {
            return Some((
                sample_name,
                "10x".to_string(),
                zip_path.to_string_lossy().to_string(),
            ));
        }
        None
    }

    /// Retrieve filenames by sample name, sorted by filename lexicographically
    fn get_files(&self, ids: Option<&Vec<usize>> ) -> Vec<(String, String, String)> {
        let mut sorted_files = match ids {
            Some(id_s ) => id_s.into_iter().map(|id| self.filenames[*id].clone()).collect(),
            None => self.filenames.clone(),
        };
        // Sort by the filename (first element of the tuple)
        sorted_files.sort_by(|a, b| a.0.cmp(&b.0)); 
        sorted_files
    }

    // Retrieve filenames by sample name, sorted lexicographically
    fn get_files_by_sample(&self, sample: &str) -> Option<Vec<(String, String, String)>> {
        self.filenames_by_sample.get(sample).map(|indices| {
            // Sort the indices lexicographically
            self.get_files( Some(indices) )
        })
    }

    // Retrieve filenames by sample name + technicalities, sorted lexicographically
    fn get_files_by_sample_tech(&self, sample: &str, tech: &str) -> Option<Vec<(String,String,String)>> {
        self.filenames_by_sample_tech.get(&(sample.to_string(), tech.to_string())).map(|indices| {
            // Sort the indices lexicographically
            self.get_files( Some(indices) )
        })
    }

    // Helper function to extract the basename
    fn extract_basename(&self, file_path: &str ) -> Option<String> {
        Path::new(file_path).file_name() // Extract the file name
            .and_then(|name| name.to_str())               // Convert OsStr to &str
            .map(|s| s.to_string())                       // Convert &str to String
    }

    fn compute_file_md5_incremental( &self, file_path:&str ) -> io::Result<String> {
        if self.omit_md5 {
            return Ok("omitted".to_string());
        }
        // Run the md5sum command
        let output = Command::new("md5sum")
            .arg(file_path)
            .output()?;
        // Check if the command was successful
        if !output.status.success() {
            return Err(io::Error::new(io::ErrorKind::Other, "md5sum command failed"));
        }

        let hash = String::from_utf8_lossy(&output.stdout);
        Ok( format!("{}", hash.split_whitespace().next().unwrap() ) )
    }


    fn get_md5sum(&self, file_path: &str) -> String {
        let path = Path::new(file_path);
        let md5_file = path.with_extension("fastq.gz.md5sum");
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





}
