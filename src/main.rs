
use std::env;
use std::collections::HashMap;
use walkdir::WalkDir;
use clap::Parser;

use rust_geo_prep::sample_files::SampleFiles;


#[derive(Parser)]
#[clap(version = "1.0.0", author = "Stefan L. <stefan.lang@med.lu.se>")]
struct Opts {
    /// the output prefix
    #[clap(short, long, default_value="sample_collection")]
    prefix: String,
    /// the column separator for the file
    #[clap(default_value= "\\t",short, long)]
    sep: String,
    /// the number of clusters
    #[clap(short, long)]
    clusters: usize,
    /// the starting temperature
    #[clap(default_value_t= 20.0,short, long)]
    temp: f32,
    /// the cooling factor
    #[clap(default_value_t= 0.9995,short, long)]
    cool: f32,
    ///max number of iterations
    #[clap(default_value_t= 1000*1000,short, long)]
    max_it: usize,
    /// the grouping outfile
    #[clap(short, long)]
    outfile: String,
}

fn main() {
    let opts: Opts = Opts::parse();
    
    let sample_file_path = format!("{}_sample_lines.tsv", opts.prefix);
    let files_file_path = format!("{}_files_md5sum_lines.tsv", opts.prefix);
    
    let sample_file_path_basename = format!("{}_basename_sample_lines.tsv", opts.prefix);
    let files_file_path_basename = format!("{}_basename_files_md5sum_lines.tsv", opts.prefix);

    
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


