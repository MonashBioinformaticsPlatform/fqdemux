# fqdemux Nextflow Pipeline

## Introduction

This pipeline demultiplexes FASTQ files based on inline barcodes using `fqtk demux`. It takes multiplexed FASTQ files, along with definitions for barcodes, UMIs and sample names, and produces demultiplexed FASTQ files and a samplesheet compatible with the [nf-core/rnaseq](https://nf-co.re/rnaseq) pipeline.

## Usage

```bash
nextflow run <path/to/pipeline>/main.nf \
  --barcodes_samplesheet <path/to/barcodes_samplesheet.tsv> \
  --readstructure_samplesheet <path/to/readstructure_samplesheet.tsv> \
  --outdir <output_directory> \
  -profile apptainer  # e.g., local, docker, singularity
```

Replace placeholders like `<path/to/pipeline>`, `<path/to/barcodes_samplesheet.tsv>`, etc., with your actual paths.

## Input Parameters

*   `--barcodes_samplesheet` (Required): Path to a tab-separated values (TSV) file defining the barcodes and their corresponding sample names. Example format:
    ```tsv
    sample_id	barcode
    A01	AATGCGGCTTGACC
    A02	ATGGTAATCCAAGG
    A03	ACTGTGCAGCTGAG
    A04	AACCGCCGAGAAGG
    ```
*   `--readstructure_samplesheet` (Required): Path to a tab-separated values (TSV) file defining the input FASTQ files and their associated read structures. The pipeline expects columns:
    *   `filename`: Full path to the FASTQ file (e.g., `/path/to/your/data/example_R1.fastq.gz`).
    *   `read_structure`: The read structure string compatible with `fqtk` (e.g., `10B141T` for 10bp barcode followed by 141bp transcript).
    Example format:
    ```tsv
    filename	read_structure
    /path/to/BRBseq_1.fq.gz	14B14M122T
    /path/to/BRBseq_2.fq.gz	150T
    ```
    See the [`fqtk` documentation](https://github.com/fulcrumgenomics/fqtk) for more details on the read structure string format.
*   `--outdir`: Path to the directory where output files will be saved. Defaults to `./results`.

* `--strandedness` (Optional, default: 'auto') - set the strandedness column for the `samplesheet.csv` output ('auto', 'forward' or 'reverse'). 

## Output

The pipeline will create the specified `--outdir` and place the following files inside:

1.  **Demultiplexed FASTQ files**: FASTQ files split by sample based on the provided barcodes.
2.  **Samplesheet**: A `samplesheet.csv` CSV file compatible with `nf-core/rnaseq`.

