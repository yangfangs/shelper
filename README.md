# shelper: Sanger Sequencing Helper Tools

`shelper` is an R package designed to streamline Sanger sequencing analysis. It provides a comprehensive suite of tools for chromatogram visualization, variant validation, and primer design, all accessible through an intuitive Shiny application or R command-line functions.

## Features

- **Chromatogram Visualization**: View and analyze Sanger sequencing traces (`.ab1`, `.scf`) directly.
- **Variant Validation**: Automatically align Sanger reads to reference genomes (hg19/hg38) to identify and validate mutations (SNVs, Indels).
- **Primer Design**: Integrated interface for **Primer3** to design PCR primers for specific genomic targets.
- **Batch Processing**: Efficiently handle multiple sequencing files in parallel using `future` and `promises`.
- **Local Caching**: Optimized performance with intelligent caching of reference sequences.

## Installation

You can install the development version of `shelper` from GitHub:

```r
# install.packages("devtools")
devtools::install_github("yangfangs/shelper")
```

## System Requirements

To utilize the full functionality of `shelper`, the following external tools should be installed and accessible in your system path:

1.  **Samtools**: For efficient retrieval of reference genome sequences.
2.  **NCBI BLAST+ (`blastn`)**: For checking primer specificity.
3.  **Primer3 (`primer3_core`)**: For the primer design module.

## Configuration

Before launching the application, you can configure the paths to your local reference genomes and tools. This is especially important if these tools are not in your system PATH or if you use custom database locations.

```r
library(shelper)

# Configure paths (example)
set_shelper_config(
  hg19 = "/path/to/blastdb/ucsc.hg19.fasta",      # Path to hg19 fasta
  hg38 = "/path/to/hg38/hg38.fa",                 # Path to hg38 fasta
  blastn = "/usr/local/bin/blastn",               # Path to blastn executable
  blast_db = "/path/to/blastdb/hg19"              # Path to BLAST database
)
```

> **Note**: If you do not provide specific paths, the package will attempt to use default paths or look for tools in your system PATH.

## Usage

### Launching the Shiny App

The easiest way to use `shelper` is through its interactive web interface:

```r
library(shelper)

# Launch the app
run_shelper()
```

## License

This project is licensed under the MIT License.
