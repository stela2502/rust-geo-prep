
use std::env;
use rust_geo_prep::*;
use std::collections::HashMap;
use walkdir::WalkDir;

fn main() {
    let args: Vec<String> = env::args().collect();
    let path_to_search = args.get(1).cloned().unwrap_or_else(|| ".".to_string());
    let sampleid = args.get(2).cloned().unwrap_or_else(|| "".to_string());
    
    let path_sanitized = path_to_search.replace(|c: char| c.is_whitespace() || c == '/' || c == '\\', "_");
    
    let sample_file_path = format!("{}/sample_collection_sample_lines.tsv", path_sanitized);
    let files_file_path = format!("{}/sample_collection_files_md5sum_lines.tsv", path_sanitized);
    
    let sample_file_path_basename = format!("{}/sample_collection_basename_sample_lines.tsv", path_sanitized);
    let files_file_path_basename = format!("{}/sample_collection_basename_files_md5sum_lines.tsv", path_sanitized);

    let mut data: HashMap<String, HashMap<String, String>> = HashMap::new();
    
    //let re1 = Regex::new(r".*/(.*?)_(S.*L\d{3}).*(R[12]|I1).*\\.fastq\\.gz").unwrap();
    //let re2 = Regex::new(r".*/(.*?)_(S\d+).*(R[12]|I1).*\\.fastq\\.gz").unwrap();
    
    for entry in WalkDir::new(&path_to_search).into_iter().filter_map(Result::ok) {
        let file_path = entry.path();
        if let Some(file_name) = file_path.file_name().and_then(|n| n.to_str()) {
            if file_name.ends_with(".fastq.gz") {
                let path_str = file_path.to_string_lossy().to_string();
                //if let Some(captures) = parse_filename(&path_str, &re1, &re2) {
                if let Some( (sample, read_type) ) = parse_filename_split(&path_str) {
                    data.entry(format!("{}", sample))
                        .or_default()
                        .insert(read_type, path_str);
                }
            }
        }
    }
    let md5data = generate_md5_file_data( &data );
    write_sample_files(&sample_file_path, &data);
    write_md5_files(&files_file_path, &md5data);
    write_sample_files_basename(&sample_file_path_basename, &data);
    write_md5_files_basename(&files_file_path_basename, &md5data);
    println!("Data written to '{}', '{}', '{}' and '{}'", 
        &sample_file_path,
        &files_file_path, 
        &sample_file_path_basename,
        &files_file_path_basename
    );
}


