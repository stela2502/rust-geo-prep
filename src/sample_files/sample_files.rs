// src/sample_files/sample_files.rs
use std::collections::{BTreeMap, BTreeSet, HashMap, HashSet};
use std::io::{self, BufWriter, Write};
use std::path::{Path, PathBuf};

use std::fs::File;

use walkdir::WalkDir;

use crate::sample_files::lane_fastqs::LaneFastqs;
use crate::sample_files::sample_record::SampleRecord;
use crate::sample_files::parsed_file::{ParsedFile, ParsedKind};



#[derive(Debug, Clone, PartialEq, Eq, PartialOrd, Ord)]
pub struct SampleKey {
    pub experiment: String,
    pub sample: String,
}

#[derive(Debug, Default)]
pub struct SampleFiles {
    pub samples: BTreeMap<SampleKey, SampleRecord>,
    pub force_experiment_prefix_export: bool,

    // basename -> (md5 -> representative parsed file)
    seen: HashMap<String, HashMap<String, ParsedFile>>,
}

impl SampleFiles {
    pub fn new() -> Self {
        Self::default()
    }

    pub fn len(&self) -> usize {
        self.samples.len()
    }

    /// Walk a directory, parse relevant items into ParsedFile, dedup backups, and add into SampleRecords.
    pub fn ingest_dir<P: AsRef<Path>>(&mut self, scan_root: P) -> io::Result<()> {
        let scan_root = scan_root.as_ref();

        // loop protection for dirs + avoid silly duplicates by canonical path
        let mut visited_dirs: HashSet<(u64, u64)> = HashSet::new();
        let mut visited_paths: HashSet<PathBuf> = HashSet::new();

        for entry in WalkDir::new(scan_root).follow_links(true).into_iter().filter_map(Result::ok) {
            let p = entry.path();

            // directory loop protection
            if let Ok(md) = p.metadata() {
                if md.is_dir() {
                    #[cfg(unix)]
                    {
                        use std::os::unix::fs::MetadataExt;
                        let key = (md.dev(), md.ino());
                        if !visited_dirs.insert(key) {
                            continue;
                        }
                    }
                }
            }

            // avoid reprocessing same file path (canonicalized)
            let canon = std::fs::canonicalize(p).unwrap_or_else(|_| p.to_path_buf());
            if !visited_paths.insert(canon) {
                continue;
            }

            if let Some(mut parsed) = ParsedFile::from_path(scan_root, p)? {
                // make sure md5 for file artifacts is populated (dirs return None)
                let _ = parsed.ensure_md5sum()?;

                // global dedup / conflict logic (backup folders)
                if self.should_ignore_as_backup(&parsed) {
                    continue;
                }

                // If content differs for same basename across experiments, we will need exp-prefix export
                self.update_export_flags(&parsed);

                // now actually add it
                self.add_file(parsed);
            }
        }

        Ok(())
    }

    /// The central “add_file”: takes a ParsedFile and routes it into the correct SampleRecord.
    pub fn add_file(&mut self, parsed: ParsedFile) {
        let key = SampleKey {
            experiment: parsed.experiment.clone(),
            sample: parsed.sample.clone(),
        };

        let rec = self.samples.entry(key).or_insert_with(|| {
            let mut r = SampleRecord::default();
            r.name = parsed.sample.clone();
            r.experiment = parsed.experiment.clone(); // add this field to SampleRecord (recommended)
            r
        });

        match parsed.kind.clone() {
            ParsedKind::TenX => {
                // you can keep "one 10x per sample" rule
                if rec.tenx.is_some() {
                    eprintln!("Duplicate 10x bundle for {}:{} ignored: {}", rec.experiment, rec.name, parsed.path);
                } else {
                    rec.tenx = Some(parsed);
                }
            }
            ParsedKind::H5 => {
                if rec.h5_files.is_some() {
                    // if exact same path, ignore; otherwise warn
                    if rec.h5_files.as_ref().unwrap().path == parsed.path {
                        // ignore
                    } else {
                        eprintln!("Duplicate H5 for {}:{} ignored: {}", rec.experiment, rec.name, parsed.path);
                    }
                } else {
                    rec.h5_files = Some(parsed);
                }
            }
            ParsedKind::Fastq { lane, role } => {
                rec.lanes.entry(lane).or_default().add_read(&role, parsed);
            }
        }
    }

    // ---------- global policy ----------

    fn should_ignore_as_backup(&mut self, parsed: &ParsedFile) -> bool {
        let base = parsed.basename();

        // Only dedup file artifacts (need md5); directories can’t be deduped here
        let md5 = match parsed.md5sum.as_ref() {
            Some(m) => m,
            None => return false,
        };

        let by_md5 = self.seen.entry(base).or_default();

        // same basename + same md5 => backup duplicate (ignore)
        if by_md5.contains_key(md5) {
            // optional: log once
            // eprintln!("Backup duplicate ignored (same md5): {}", parsed.path);
            return true;
        }

        by_md5.insert(md5.to_string(), parsed.clone());
        false
    }

    fn update_export_flags(&mut self, parsed: &ParsedFile) {
        let base = parsed.basename();
        let md5 = match parsed.md5sum.as_ref() {
            Some(m) => m,
            None => return, // dirs
        };

        // Look for other variants with same basename but different md5
        if let Some(by_md5) = self.seen.get(&base) {
            if by_md5.len() >= 2 {
                // already a conflict; export must disambiguate
                self.force_experiment_prefix_export = true;
                return;
            }
            for (other_md5, other_pf) in by_md5 {
                if other_md5 != md5 {
                    // different content with same basename
                    if other_pf.experiment != parsed.experiment {
                        self.force_experiment_prefix_export = true;
                    } else {
                        // same experiment, same basename, different content => this is dangerous
                        self.force_experiment_prefix_export = true;
                        eprintln!(
                            "WARNING: same experiment '{}' has two different files with basename '{}' (md5 differs).",
                            parsed.experiment, base
                        );
                    }
                }
            }
        }
    }

    // ---------- naming helpers for writers ----------

    /// GEO upload filename to use for a source path.
    /// If force_experiment_prefix_export is true, prefixes with experiment id.
    pub fn geo_filename(&self, experiment: &str, src_path: &str) -> String {
        let base = Path::new(src_path)
            .file_name()
            .and_then(|s| s.to_str())
            .unwrap_or(src_path);

        if self.force_experiment_prefix_export {
            format!("{}_{}", experiment, base)
        } else {
            base.to_string()
        }
    }

    /// GEO sample name to use in tables (optional but recommended).
    pub fn geo_sample_name(&self, experiment: &str, sample: &str) -> String {
        if self.force_experiment_prefix_export {
            format!("{}_{}", experiment, sample)
        } else {
            sample.to_string()
        }
    }


    /// Iterate all ParsedFiles that are intended to be exported/copied.
    /// NOTE: If TenX is still stored as a directory, you probably want to zip first;
    /// this will still list it, but scripts will fail to copy dirs with cp/copy-item.
    fn iter_all_parsed_files(&self) -> Vec<&ParsedFile> {
        let mut out: Vec<&ParsedFile> = Vec::new();

        for (_key, rec) in &self.samples {
            if let Some(pf) = rec.tenx.as_ref() {
                out.push(pf);
            }
            if let Some(pf) = rec.h5_files.as_ref() {
                out.push(pf);
            }
            for lane in rec.lanes.values() {
                for pf in lane.reads.values() {
                    out.push(pf);
                }
            }
        }

        out
    }

    /// Like iter_all_parsed_files, but yields mutable refs (used for ensure_md5sum).
    fn iter_all_parsed_files_mut(&mut self) -> Vec<&mut ParsedFile> {
        let mut out: Vec<&mut ParsedFile> = Vec::new();

        for (_key, rec) in &mut self.samples {
            if let Some(pf) = rec.tenx.as_mut() {
                out.push(pf);
            }
            if let Some(pf) = rec.h5_files.as_mut() {
                out.push(pf);
            }
            for lane in rec.lanes.values_mut() {
                for pf in lane.reads.values_mut() {
                    out.push(pf);
                }
            }
        }

        out
    }

    /// Write md5 table using GEO filename (basename or exp-prefixed basename, depending on geo_filename()).
    pub fn write_md5_files_basename<P: AsRef<Path>>(&mut self, out_path: P) -> io::Result<()> {
        // Ensure md5 is computed for all file-path ParsedFiles that need it.
        for pf in self.iter_all_parsed_files_mut() {
            let _ = pf.ensure_md5sum()?; // dirs will return Ok(None)
        }

        // Collect rows: (geo_file_name, md5)
        let mut rows: Vec<(String, String)> = Vec::new();
        for pf in self.iter_all_parsed_files() {
            let geo_name = self.geo_filename(&pf.experiment, &pf.path);
            let md5 = pf.md5sum.clone().unwrap_or_else(|| "none".to_string());
            rows.push((geo_name, md5));
        }

        // Stable ordering
        rows.sort_by(|a, b| a.0.cmp(&b.0));

        let f = File::create(out_path)?;
        let mut w = BufWriter::new(f);

        writeln!(w, "file_name\tmd5sum")?;
        for (name, md5) in rows {
            writeln!(w, "{}\t{}", name, md5)?;
        }
        Ok(())
    }

    /// Generate bash script to copy all referenced files into DEST, using GEO filenames.
    pub fn write_collect_all_files_script_sh<P: AsRef<Path>>(
        &mut self,
        script_path: P,
        dest: &str,
    ) -> io::Result<()> {
        // Ensure md5 exists if you want scripts to be consistent with tables later
        // (optional, but cheap since you already computed earlier)
        for pf in self.iter_all_parsed_files_mut() {
            let _ = pf.ensure_md5sum()?;
        }

        // Collect pairs: (dst_name, src_path) for stable ordering
        let mut pairs: Vec<(String, String)> = Vec::new();
        for pf in self.iter_all_parsed_files() {
            let dst_name = self.geo_filename(&pf.experiment, &pf.path);
            pairs.push((dst_name, pf.path.clone()));
        }
        pairs.sort_by(|a, b| a.0.cmp(&b.0));

        let f = File::create(script_path)?;
        let mut w = BufWriter::new(f);

        writeln!(w, "#!/usr/bin/env bash")?;
        writeln!(w, "set -euo pipefail")?;
        writeln!(w, "DEST=\"{}\"", dest)?;
        writeln!(w, "mkdir -p \"$DEST\"")?;
        writeln!(w)?;
        writeln!(w, "COPY_CMD=(cp -f)")?;
        writeln!(w)?;

        for (dst_name, src) in pairs {
            writeln!(w, "\"${{COPY_CMD[@]}}\" \"{}\" \"$DEST/{}\"", src, dst_name)?;
        }

        Ok(())
    }

    /// Generate PowerShell script to copy all referenced files into DEST, using GEO filenames.
    pub fn write_collect_all_files_script_ps1<P: AsRef<Path>>(
        &mut self,
        script_path: P,
        dest: &str,
    ) -> io::Result<()> {
        for pf in self.iter_all_parsed_files_mut() {
            let _ = pf.ensure_md5sum()?;
        }

        let mut pairs: Vec<(String, String)> = Vec::new(); // (dst_name, src_path)
        for pf in self.iter_all_parsed_files() {
            let dst_name = self.geo_filename(&pf.experiment, &pf.path);
            pairs.push((dst_name, pf.path.clone()));
        }
        pairs.sort_by(|a, b| a.0.cmp(&b.0));

        let f = File::create(script_path)?;
        let mut w = BufWriter::new(f);

        writeln!(w, "Param()")?;
        writeln!(w, "$ErrorActionPreference = 'Stop'")?;
        writeln!(w, "$DEST = \"{}\"", dest)?;
        writeln!(w, "New-Item -ItemType Directory -Force -Path $DEST | Out-Null")?;
        writeln!(w)?;

        for (dst_name, src) in pairs {
            writeln!(
                w,
                "Copy-Item -LiteralPath \"{}\" -Destination (Join-Path $DEST \"{}\") -Force",
                src, dst_name
            )?;
        }

        Ok(())
    }

    /// Recreates your old sample table writer, now backed by ParsedFile.
    /// The table uses GEO upload filenames (geo_filename) for TenX/H5/FASTQ cells.
    pub fn write_sample_files_basename<P: AsRef<Path>>(&self, out_path: P) -> io::Result<()> {
        let mut f = BufWriter::new(File::create(out_path)?);

        // We need a stable global header: determine maximum #lanes and role order.
        // Approach: compute global max lanes and global role set.
        let mut global_roles: BTreeSet<String> = BTreeSet::new();
        let mut max_lanes: usize = 0;

        for (_key, rec) in &self.samples {
            let roles = Self::all_roles_sorted(rec);
            for r in roles {
                global_roles.insert(r);
            }
            let lane_count = rec.lanes.len();
            if lane_count > max_lanes {
                max_lanes = lane_count;
            }
        }

        // Prefer canonical ordering globally too
        let mut roles_vec: Vec<String> = {
            let mut tmp = Vec::new();
            for r in ["I1", "I2", "R1", "R2"] {
                if global_roles.remove(r) {
                    tmp.push(r.to_string());
                }
            }
            tmp.extend(global_roles.into_iter());
            tmp
        };

        if roles_vec.is_empty() {
            // still write a sane header if no fastqs found
            roles_vec = vec!["I1".into(), "R1".into(), "R2".into()];
        }

        // ---- header ----
        write!(f, "Source_Path(s)\tSample_Lane\tTenX\tH5")?;
        for _lane_idx in 0..max_lanes {
            for r in &roles_vec {
                write!(f, "\t{}", r)?;
            }
        }
        writeln!(f)?;

        // ---- rows ----
        // Sort output by (experiment, sample) to keep stable
        let mut keys: Vec<_> = self.samples.keys().cloned().collect();
        keys.sort();

        for key in keys {
            let rec = self.samples.get(&key).unwrap();

            let src_folders = Self::collect_source_folders_for_record(rec);
            let sample_name = rec.name.clone();

            // TenX/H5 cells: GEO upload name or empty
            let tenx_cell = rec
                .tenx
                .as_ref()
                .map(|pf| pf.geo_filename() )
                .unwrap_or_default();

            let h5_cell = rec
                .h5_files
                .as_ref()
                .map(|pf| pf.geo_filename() )
                .unwrap_or_default();

            write!(f, "{}\t{}\t{}\t{}", src_folders, sample_name, tenx_cell, h5_cell)?;

            // Render lanes in sorted lane-key order, but pad to max_lanes
            let lane_keys = Self::all_lane_keys_sorted(rec);

            for i in 0..max_lanes {
                if let Some(lk) = lane_keys.get(i) {
                    let lane = rec.lanes.get(lk).unwrap();
                    let fmt = |pf: &ParsedFile| pf.geo_filename();
                    let cells = lane.row_cells(&roles_vec, self.force_experiment_prefix_export, &fmt);
                    for c in cells {
                        write!(f, "\t{}", c)?;
                    }
                } else {
                    // pad missing lanes with empty cells
                    for _ in &roles_vec {
                        write!(f, "\t")?;
                    }
                }
            }

            writeln!(f)?;
        }

        Ok(())
    }
}
