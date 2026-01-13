//lane_fastqs.rs

use std::collections::{BTreeMap};
use crate::ParsedFile;

#[derive(Debug, Default)]
pub struct LaneFastqs {
    pub reads: BTreeMap<String, ParsedFile>,
}

impl LaneFastqs {
    /// Add a FASTQ for a lane under a specific role (R1/R2/I1/...)
    pub fn add_read(&mut self, role: &str, path: ParsedFile) {
        if let Some(existing) = self.reads.get(role) {
            eprintln!(
                "Duplicate read role '{}' for lane: already have '{}', tried to add '{}' - file is ignored!",
                role, existing.path, path.path
            );
        }else {
           self.reads.insert(role.to_string(), path); 
        }
        
    }

    /// Render FASTQ cells for this lane in the provided `roles` order.
    pub fn row_cells<F>(&self, roles: &[String], include_experiment: bool, fmt: &F) -> Vec<String>
    where
        F: Fn(&ParsedFile) -> String,
    {
        roles
            .iter()
            .map(|role| self.reads.get(role).map(|p| fmt(&p) ).unwrap_or_default())
            .collect()
    }
}