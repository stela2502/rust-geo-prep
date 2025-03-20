extern crate md5;

use std::collections::HashMap;
use std::env;
use std::fs::{self, File};
use std::io::{Read, BufRead, BufReader, Write};
use std::path::Path;
use md5::{Md5, Digest};
use walkdir::WalkDir;


fn main() {
    let args: Vec<String> = env::args().collect();
    let path_to_search = args.get(1).cloned().unwrap_or_else(|| ".".to_string());
    let sampleid = args.get(2).cloned().unwrap_or_else(|| "".to_string());
    
    let path_sanitized = path_to_search.replace(|c: char| c.is_whitespace() || c == '/' || c == '\\', "_");
    
    let sample_file_path = format!("{}/sample_collection_{}.sample_lines.tsv", path_sanitized, sampleid);
    let files_file_path = format!("{}/sample_collection_{}.files_md5sum_lines.tsv", path_sanitized, sampleid);
    
    let mut data: HashMap<String, HashMap<String, String>> = HashMap::new();
    
    //let re1 = Regex::new(r".*/(.*?)_(S.*L\d{3}).*(R[12]|I1).*\\.fastq\\.gz").unwrap();
    //let re2 = Regex::new(r".*/(.*?)_(S\d+).*(R[12]|I1).*\\.fastq\\.gz").unwrap();
    
    for entry in WalkDir::new(&path_to_search).into_iter().filter_map(Result::ok) {
        let file_path = entry.path();
        if let Some(file_name) = file_path.file_name().and_then(|n| n.to_str()) {
            if file_name.ends_with(".fastq.gz") {
                let path_str = file_path.to_string_lossy().to_string();
                //if let Some(captures) = parse_filename(&path_str, &re1, &re2) {
                if let Some(captures) = parse_filename_split(&path_str) {
                    let (sample, lane, read_type) = captures;
                    data.entry(format!("{}_{}", sample, lane))
                        .or_default()
                        .insert(read_type, path_str);
                }
            }
        }
    }
    
    write_sample_files(&sample_file_path, &data);
    write_md5_files(&files_file_path, &data);
    println!("Data written to '{}'", sample_file_path);
}


fn parse_filename_split(file_path: &str) -> Option<(String, String, String)> {
    let parts: Vec<&str> = file_path.split('/').last()?.split('_').collect();

    if parts.len() < 3 {
        return None;
    }

    // Read type is typically R1, R2, or I1, usually near the end
    let read_type = parts.last()?.to_string();

    // Find lane information, which often has 'S' and 'L' markers
    let lane_idx = parts.iter().position(|&s| s.starts_with('S'))?;
    let lane = parts.get(lane_idx)?.to_string();

    // The sample name is everything before the lane info
    let sample = parts[..lane_idx].join("_");

    Some((sample, lane, read_type))
}

/*fn parse_filename(file_path: &str, re1: &Regex, re2: &Regex) -> Option<(String, String, String)> {
    if let Some(caps) = re1.captures(file_path) {
        return Some((caps[1].to_string(), caps[2].to_string(), caps[3].to_string()));
    }
    if let Some(caps) = re2.captures(file_path) {
        return Some((caps[1].to_string(), caps[2].to_string(), caps[3].to_string()));
    }
    None
}*/


fn write_sample_files(path: &str, data: &HashMap<String, HashMap<String, String>>) {
    let mut file = File::create(path).expect("Could not create sample file");
    writeln!(file, "Sample_Lane\tR1\tR2\tI1").unwrap();
    for (sample_lane, reads) in data {
        writeln!(file, "{}\t{}\t{}\t{}", 
                 sample_lane,
                 reads.get("R1").unwrap_or(&"MISSING_R1".to_string()),
                 reads.get("R2").unwrap_or(&"MISSING_R2".to_string()),
                 reads.get("I1").unwrap_or(&"MISSING_I1".to_string())
        ).unwrap();
    }
}

fn write_md5_files(path: &str, data: &HashMap<String, HashMap<String, String>>) {
    let mut file = File::create(path).expect("Could not create md5 file");
    writeln!(file, "file_name\tmd5sum").unwrap();
    for reads in data.values() {
        for file_path in reads.values() {
            let md5sum = get_md5sum(file_path);
            writeln!(file, "{}\t{}", file_path, md5sum).unwrap();
        }
    }
}

fn get_md5sum(file_path: &str) -> String {
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
    if let Ok(mut file) = File::open(path) {
        let mut hasher = Md5::new();
        let mut buffer = vec![0; 8192];
        while let Ok(bytes) = file.read(&mut buffer) {
            if bytes == 0 { break; }
            hasher.update(&buffer[..bytes]);
        }
        let result = format!("{:x}", hasher.finalize());
        let _ = fs::write(&md5_file, &result);
        return result;
    }
    "none".to_string()
}

