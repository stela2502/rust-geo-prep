# rust-geo-prep

**A Rust CLI tool to prepare FASTQ and 10x data for GEO submission**

Submitting sequencing data to GEO is tedious: FASTQs and 10x outputs are
scattered across directories, samples must be grouped consistently, and
MD5 checksums must be reported correctly. `rust-geo-prep` automates this
process by scanning your project directory, grouping files into samples,
computing MD5 sums, and generating reproducible collection scripts.

It assumes that you have grouped your data into experiment folders and that 
both the FASTQ data and the count outputs are located in subfolders of this 
main experiment folder.
It further requires that files follow directory structures as created by 
the CellRanger tool:

- FASTQ files must contain the SampleID at the beginning of the filename and a read type identifier such as R1, R2, I1, or I2.
- matrix.mtx.gz files must be located in <sample_id>/outs/filtered_feature_bc_matrix/ and be accompanied by features.tsv.gz and barcodes.tsv.gz.
- filtered_feature_bc_matrix.h5 files must be located in <sample_id>/outs/.


The tool generates GEO-ready sample tables, MD5 checksum reports, FASTQ pairing tables, and a platform-specific collection script that safely gathers all referenced files into a single upload directory.

------------------------------------------------------------------------

## Features

-   Recursively scans experiment folders
-   Groups FASTQ files into GEO-ready sample groups
-   Supports multiple FASTQ suffixes (e.g. `.fastq.gz`, `.fq.gz`)
-   Supports 10x HDF5 / MTX triplets
-   Excludes arbitrary paths
-   Computes MD5 checksums for every file
-   Generates deterministic, reproducible outputs
-   Automatically resolves filename collisions during collection
-   Designed for large HPC / shared-storage projects

------------------------------------------------------------------------

## Installation

### From source

``` bash
cargo install --git https://github.com/stela2502/rust-geo-prep
```

Or build locally:

``` bash
cargo build --release
cp target/release/rust-geo-prep ~/bin/
```

------------------------------------------------------------------------

## Usage

``` bash
rust-geo-prep [OPTIONS]
```

### Options

  --------------------------------------------------------------------------
  Option                    Description
  ------------------------- ------------------------------------------------
  `-i, --input <DIR>`       Root directory. Each direct subfolder is treated
                            as one experiment

  `-e, --exclude <NAME>`    Path names to ignore (can be repeated)

  `-p, --prefix <PREFIX>`   Output file prefix (default:
                            `sample_collection`)

  `-s, --suffix <SUFFIX>`   File suffixes to include (can be repeated)

  `-h, --help`              Show help

  `-V, --version`           Show version
  --------------------------------------------------------------------------

Defaults:

``` text
suffixes: .fastq.gz .fq.gz
prefix:   sample_collection
```

------------------------------------------------------------------------

## Directory Layout Assumption

``` text
INPUT/
  experiment_1/
    sampleA_R1.fastq.gz
    sampleA_R2.fastq.gz
  experiment_2/
    ...
```

Each **direct subfolder of `--input` is treated as one experiment**.

------------------------------------------------------------------------

## FASTQ Example

``` bash
rust-geo-prep \
  --input /data/projects \
  --exclude old_runs \
  --exclude geo_downloaded_data \
  --suffix .fastq.gz \
  --suffix .fq.gz \
  --prefix geo_submission
```

------------------------------------------------------------------------

## 10x Matrix / HDF5 Example

To collect 10x CellRanger outputs:

``` bash
rust-geo-prep \
  --input /data/projects \
  --suffix filtered_feature_bc_matrix.h5 \
  --suffix matrix.mtx.gz \
  --prefix geo_10x
```
And collect the files using the also created copy script.

This allows you to prepare unique:

-   `<sample_id>_filtered_feature_bc_matrix.h5`
-   `<sample_id>.zip` (combining the 10x matrix triplets into one zip)

for GEO submission.

------------------------------------------------------------------------

## Generated Files

Typical outputs:

  File                        Purpose
  --------------------------- ------------------------------
  `*_sample_collection.tsv`   GEO sample table
  `*_files_md5sum.tsv`        MD5 checksum table
  `*_fastq_pairs.tsv`         FASTQ R1/R2 pairing table
  `*_collection_script.sh`    Bash collection script
  `*_collection_script.ps1`   PowerShell collection script

------------------------------------------------------------------------

## FASTQ Pair Table

The FASTQ pairs table contains one row per logical sample and groups:

-   R1
-   R2
-   I1 / I2 (if present)

into a single row. This makes it easy to inspect whether pairs are
complete and consistent before submission.

------------------------------------------------------------------------

## Collection Scripts

The generated scripts:

-   `*_collection_script.sh`
-   `*_collection_script.ps1`

copy all referenced files into a single destination directory.

### Automatic filename disambiguation

If two samples would result in identical filenames, the tool
automatically adds the experiment name to the **unique filenames** during collection while
keeping full traceability in the tables.

This guarantees:

-   No overwriting
-   GEO-safe flat upload directories
-   Stable reproducibility

You do not need to resolve collisions manually.

------------------------------------------------------------------------

## GEO Workflow

Recommended workflow:

1.  Run `rust-geo-prep`
2.  Inspect the generated TSV tables
3.  Run the collection script to gather all referenced files into one
    directory
4.  Upload the collected directory to GEO
5.  Use the MD5 table for GEO validation

------------------------------------------------------------------------

## Excluding Paths

``` bash
--exclude tmp --exclude backup --exclude geo_downloaded
```

Any path containing the excluded token is ignored.

------------------------------------------------------------------------

## Multiple Suffixes

``` bash
--suffix .fastq.gz --suffix .fq.gz
```

and to include the quantified data add

``` bash
--suffix filtered_feature_bc_matrix.h5 --suffix matrix.mtx.gz
```

All listed suffixes are treated as valid target files.

------------------------------------------------------------------------

## Platform Notes

-   Linux/macOS: use the generated `.sh` script
-   Windows: use the generated `.ps1` script
-   Paths are preserved exactly as discovered

------------------------------------------------------------------------

## Reproducibility

-   Files are sorted deterministically
-   Grouping is stable
-   MD5 sums are calculated on demand
-   Scripts are reproducible

------------------------------------------------------------------------

## Troubleshooting

### "Could not determine read role"

Your FASTQ filename does not follow standard R1/R2/I1/I2 naming
conventions.

### Duplicate filenames

Handled automatically by the collection script with unique renaming.

------------------------------------------------------------------------

## Philosophy

`rust-geo-prep` does **not** try to infer experimental biology.

It guarantees:

-   Correct grouping
-   Correct checksums
-   GEO-compatible file handling

Biological interpretation remains the responsibility of the researcher.

------------------------------------------------------------------------

## Author

Stefan Lang\
Division of Molecular Hematology, Lund University\
ORCID: 0000-0002-0854-2328

------------------------------------------------------------------------

## License

MIT License
