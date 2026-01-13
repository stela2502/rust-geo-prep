// src/sample_files/mod.rs
pub mod parsed_file;
pub mod sample_files;
pub mod lane_fastqs;
pub mod sample_record;

pub use parsed_file::{ParsedFile, ParsedKind};
pub use sample_files::{SampleFiles, SampleKey};
pub use lane_fastqs::LaneFastqs;
pub use sample_record::SampleRecord;