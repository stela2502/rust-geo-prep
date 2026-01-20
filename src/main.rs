
use clap::Parser;
use std::path::{Path, PathBuf};

use rust_geo_prep::sample_files::SampleFiles;

/// Submitting data to GEO is complex. 
/// This tool helps by collecting the different fastq files and grouping them into samples groups.
/// It also calculates the md5sums and reports them for every fastq file.
#[derive(Parser)]
#[clap(version = "1.0.0", author = "Stefan L. <stefan.lang@med.lu.se>")]
struct Opts {
    /// the output prefix
    #[clap(short, long, default_value="sample_collection")]
    prefix: String,

    /// File suffixes treated as target files
    ///
    /// Can be specified multiple times:
    ///   --suffix .fastq.gz --suffix .fq.gz
    #[clap(
        short = 's',
        long = "suffix",
        multiple_occurrences = true,
        default_values = &[".fastq.gz", ".fq.gz"]
    )]
    suffixes: Vec<String>,

    /// path names to ignore
    ///
    /// Can be specified multiple times:
    ///   --exclude old_experiment --exclude geo_downloaded_data
    #[clap(
        short = 'e',
        long = "exclude",
        multiple_occurrences = true,
    )]
    exclude: Vec<String>,

    /// Root directory. Each direct subfolder is an experiment.
    #[clap(short, long )]
    input: Option<PathBuf>,

}


fn main(){
    let opts: Opts = Opts::parse();
    
    let sample_file_path = format!("{}.tsv", opts.prefix);
    let files_file_path = format!("{}_md5sum.tsv", opts.prefix);
    let pairs_file_path = format!("{}_pairs.tsv", opts.prefix);
    let collection_script_path = if cfg!(windows) {
        format!("{}_collection_script.ps1", opts.prefix)
    } else {
        format!("{}_collection_script.sh", opts.prefix)
    };
    let collection_dest = format!("{}_all_files_copied", opts.prefix);
    
    //let sample_file_path_basename = format!("{}_basename_sample_lines.tsv", opts.prefix);
    //let files_file_path_basename = format!("{}_basename_files_md5sum_lines.tsv", opts.prefix);

    println!("We are searching for files ending on either of these strings {:?}", opts.suffixes );

    let root = opts.input.as_deref().unwrap_or(Path::new("."));

    
    let mut data = SampleFiles::new();
    
    let (added, visited) = match data.ingest_dir(root, &opts.suffixes, &opts.exclude) {
        Err(e) => {
            eprintln!("\n❌ Failed while scanning input directories:");
            eprintln!("   {e}\n");
            std::process::exit(1);
        },
        Ok(i) => i,
    };

    let _ = data.write_sample_files_basename(&sample_file_path);
    let _ = data.write_md5_files_basename(&files_file_path);
    let _ = data.write_fastq_pairs_table(&pairs_file_path );
    let _ = if cfg!(windows) {
        data.write_collect_all_files_script_ps1(&collection_script_path, &collection_dest)
    } else {
        data.write_collect_all_files_script_sh(&collection_script_path, &collection_dest)
    };

    //let _ = data.write_sample_files_basename(&sample_file_path_basename);
    //let _ = data.write_md5_files_basename(&files_file_path_basename);

    let run_cmd = if cfg!(windows) {
        format!(".\\{}", collection_script_path)
    } else {
        format!("bash {}", collection_script_path)
    };



    println!(
        "\n{} target files detected; {} added and {} samples identified.\n\
         \nOutput files:\n\
         - Sample table      : {}\n\
         - MD5 checksum table: {}\n\
         - Pairs collection  : {}\n\
         - Collection script : {}\n\
         - Copy destination  : {}\n\
         \nNext steps:\n\
         1) Review the TSV files for correctness.\n\
         2) Run the collection script to gather all referenced files:\n\
            {}\n\
         3) Use the TSV files to fill in the official GEO submission spreadsheets.\n\
         \nNote: These files are intermediate manifests. Experimental metadata must be added manually.\n",
        visited, 
        added,
        data.len(),
        sample_file_path,
        files_file_path,
        pairs_file_path,
        collection_script_path,
        collection_dest,
        run_cmd
    );    
    if data.force_experiment_prefix_export{
        println!("Experiment names are part of the published file names as a sample id overlap was detected!")
    }
}


