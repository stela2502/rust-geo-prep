//sample_record.rs
use super::{LaneFastqs, ParsedFile};

use std::collections::{BTreeMap, BTreeSet};
use std::path::{Path, PathBuf};


#[derive(Debug, Default)]
pub struct SampleRecord {
    pub name: String,

    /// 10x bundle (zip), optional
    pub tenx: Option<ParsedFile>,

    /// keep a experiment hint in case of duplicate sample names!
    pub experiment: String,

    /// h5 file, optional
    pub h5_files: Option<ParsedFile>,

    /// FASTQ lanes grouped by lane key, each containing role→path (R1/R2/I1/...)
    pub lanes: BTreeMap<String, LaneFastqs>,
}

impl SampleRecord {
    
    pub fn fastq_source_folders(&self) -> String
    {
        let mut folders: BTreeSet<String> = BTreeSet::new();

        for (_lane_key, lane) in &self.lanes {
            for (_role, path) in &lane.reads {
                if let Some(parent) = Path::new(&path.path).parent() {
                    folders.insert(parent.to_string_lossy().to_string());
                }
            }
        }

        folders.into_iter().collect::<Vec<_>>().join(",")
    }

    /// Render a single flattened row for this sample: Sample + TenX + H5 + (lane blocks...)
    pub fn row_cells<F>(&self, roles: &[String], include_experiment:bool, fmt: &F, max_lanes: usize) -> Vec<String>
    where
        F: Fn(&ParsedFile) -> String,
    {
        let mut out = Vec::new();

        // first columns
        out.push(self.fastq_source_folders());
        out.push(self.name.clone());
        out.push(self.tenx.as_ref().map(|p| p.path.clone()).unwrap_or_default());
        out.push(self.h5_files.as_ref().map(|p| p.path.clone()).unwrap_or_default());

        // lane blocks (sorted by key)
        let mut lane_count = 0usize;
        for (_lane_key, lane) in &self.lanes {
            out.extend(lane.row_cells(roles, include_experiment, fmt));
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
    pub fn all_paths<'a>(&'a self) -> impl Iterator<Item = &'a ParsedFile> + 'a {
        let tenx = self.tenx.as_ref().into_iter();
        let h5   = self.h5_files.as_ref().into_iter();
        let fastqs = self
            .lanes
            .values()
            .flat_map(|lane| lane.reads.values())
            .map(|s| s);

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

    /// GEO sample name: prefix with experiment when conflicts exist.
    pub fn geo_sample_name(&self, force_experiment_prefix_export: bool ) -> String {
        if force_experiment_prefix_export {
            format!("{}_{}", self.experiment, self.name)
        } else {
            self.name.to_string()
        }
    }

    fn parent_dir_string(p: &str) -> Option<String> {
        Path::new(p)
            .parent()
            .map(|pp| pp.to_string_lossy().to_string())
    }


    /// Unique parent folders for all files referenced by this record, comma-separated.
    pub fn collect_source_folders_for_record(&self) -> String {
        let mut set: BTreeSet<String> = BTreeSet::new();

        if let Some(pf) = self.tenx.as_ref() {
            if let Some(par) = Self::parent_dir_string(&pf.path) {
                set.insert(par);
            }
        }
        if let Some(pf) = self.h5_files.as_ref() {
            if let Some(par) = Self::parent_dir_string(&pf.path) {
                set.insert(par);
            }
        }
        for lane in self.lanes.values() {
            for pf in lane.reads.values() {
                if let Some(par) = Self::parent_dir_string(&pf.path) {
                    set.insert(par);
                }
            }
        }

        set.into_iter().collect::<Vec<_>>().join(",")
    }

    /// Lane keys in stable order.
    pub fn lane_keys_sorted(&self) -> Vec<String> {
        let mut lanes: Vec<String> = self.lanes.keys().cloned().collect();
        lanes.sort();
        lanes
    }

    /// Role names in stable order (I1/I2/R1/R2 first, then the rest alphabetically).
    pub fn all_roles_sorted(&self) -> Vec<String> {
        let mut set: BTreeSet<String> = BTreeSet::new();
        for lane in self.lanes.values() {
            for role in lane.reads.keys() {
                set.insert(role.clone());
            }
        }

        let preferred = ["I1", "I2", "R1", "R2"];
        let mut out = Vec::new();
        for r in preferred {
            if set.remove(r) {
                out.push(r.to_string());
            }
        }
        out.extend(set.into_iter());
        out
    }
}
