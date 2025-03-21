
use std::env;
use std::collections::HashMap;
use walkdir::WalkDir;

use rust_geo_prep::sample_files::SampleFiles;

fn main() {
    let args: Vec<String> = env::args().collect();
    
    let sample_file_path = "sample_collection_sample_lines.tsv";
    let files_file_path = "sample_collection_files_md5sum_lines.tsv";
    
    let sample_file_path_basename = "sample_collection_basename_sample_lines.tsv";
    let files_file_path_basename = "sample_collection_basename_files_md5sum_lines.tsv";

    
    let mut data = SampleFiles::new();
    
    // Parse files and group by sample name, technicalities, and read type
    for entry in WalkDir::new( "." ).into_iter().filter_map(Result::ok) {
        let file_path = entry.path();
        if let Some(file_name) = file_path.file_name().and_then(|n| n.to_str()) {
            if file_name.ends_with(".fastq.gz") {
                let path_str = file_path.to_string_lossy();
                data.add_file( &path_str );
            }
        }
    }

    data.write_sample_files(&sample_file_path);
    data.write_md5_files(&files_file_path);
    data.write_sample_files_basename(&sample_file_path_basename);
    data.write_md5_files_basename(&files_file_path_basename);

    println!("Data written to '{}', '{}', '{}' and '{}'", 
        &sample_file_path,
        &files_file_path, 
        &sample_file_path_basename,
        &files_file_path_basename
    );
}


