# KrakenClip

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Rust](https://img.shields.io/badge/rust-1.70%2B-blue.svg)](https://www.rust-lang.org/)
[![Docker Pulls](https://img.shields.io/docker/pulls/hecrp/krakenclip.svg)](https://hub.docker.com/r/hecrp/krakenclip)
[![Docker Image Size](https://img.shields.io/docker/image-size/hecrp/krakenclip.svg)](https://hub.docker.com/r/hecrp/krakenclip)
[![Build Status](https://img.shields.io/badge/build-passing-brightgreen)](https://github.com/hecrp/krakenclip)
[![Maintenance](https://img.shields.io/badge/Maintained%3F-yes-green.svg)](https://github.com/hecrp/krakenclip/graphs/commit-activity)

KrakenClip is a high-performance command-line utility written in [Rust](https://www.rust-lang.org/) for processing and analyzing [Kraken2](https://ccb.jhu.edu/software/kraken2/) bioinformatics software reports and log files. This toolkit focuses on fast and efficient processing of classifier outputs, making it an ideal choice for large datasets or bioinformatics pipelines as it works as a standalone binary with no external dependencies.

## Getting Started

### Prerequisites

- Rust programming language (latest stable version)
- Cargo (Rust's package manager)

### Installation

1. Clone the repository:
   ```
   git clone https://github.com/hecrp/krakenclip.git
   cd krakenclip
   ```

2. Build the project:
   ```
   cargo build --release
   ```

3. The executable will be available in `target/release/kraken2-parser`
    ```
    ./target/release/kraken2-parser --help
    ```

### Docker

You can also use the provided Dockerfile to build and run KrakenClip in a container:

1. Build the Docker image:
   ```
   docker build -t krakenclip .
   ```

2. Run the Docker container:
   ```
   docker run --rm krakenclip
   ```

3. To process files, mount volumes and run the toolkit:
   ```
   docker run --rm -v /path/to/local/data:/data krakenclip analyze /data/report.txt
   ```

## Usage

The basic syntax for using KrakenClip is shown below:

```
USAGE:
    krakenclip [SUBCOMMAND]

OPTIONS:
    -h, --help       Print help information
    -V, --version    Print version information

SUBCOMMANDS:
    analyze               Analyze Kraken2 report
    extract               Extract sequences based on Kraken2 results
    generate-test-data    Generate test data for performance testing
    help                  Print this message or the help of the given subcommand(s)
```

### Analyze Module

Used to analyze and extract information from Kraken2 reports:

```
USAGE:
    krakenclip analyze [OPTIONS] <REPORT>

ARGS:
    <REPORT>    Kraken2 report file

OPTIONS:
    -h, --help                 Print help information
        --json <JSON>          Generate JSON output
        --tax-id <TAXON_ID>    Search for a specific taxon by ID
```

### Extract Module

Used to extract sequences based on Kraken2 results:

```
USAGE:
    krakenclip extract [OPTIONS] --output <OUTPUT> --taxids <TAXIDS> <SEQUENCE> <LOG>

ARGS:
    <SEQUENCE>                Input FASTA/FASTQ file
    <LOG>                     Kraken2 log file

OPTIONS:
    -h, --help                Print help information
    -o, --output <OUTPUT>     Output file for extracted sequences
        --report <REPORT>     Kraken2 report file (required for hierarchy (--include) options)
        --taxids <TAXIDS>     Comma-separated list of taxids to extract
        --include-children    Include sequences from all descendant taxa
        --include-parents     Include sequences from all ancestor taxa
        --stats-output <FILE> Generate a statistics file with detailed extraction information
```

#### Hierarchical Extraction

The Extract module now supports hierarchical taxonomic extraction with two key options:

- **`--include-children`**: Extracts sequences from the specified taxids AND all their descendant taxa.

- **`--include-parents`**: Extracts sequences from the specified taxids AND all their ancestor taxa.

Both options require providing a Kraken2 report file with the `--report` option, as the taxonomic hierarchy information is extracted from there.

#### Statistics Report

The new `--stats-output` option generates a comprehensive markdown-formatted statistics file that includes:

- Total sequence counts (extracted vs. input)
- Breakdown of extracted sequences by taxid
- Percentage of sequences per taxid relative to total extracted and total input
- Distinction between original taxids and those added through hierarchical expansion (expanded)
- Summary statistics for original vs. expanded taxa

### Generate Test Data Module

Used to generate test data for benchmarking and testing:

```
USAGE:
    krakenclip generate-test-data [OPTIONS] --output <OUTPUT> --lines <LINES> --type <TYPE>

OPTIONS:
    -h, --help                Print help information
    -o, --output <OUTPUT>     Output file path
    -l, --lines <LINES>       Number of lines to generate
    -t, --type <TYPE>         Type of data to generate (wide, deep, fragments, dense, etc.)
```

## Performance

The KrakenClip is highly optimized for performance and memory efficiency. It implements several advanced techniques:

### Optimization Techniques
- **Block-based reading**: Uses optimized 512KB buffers for efficient file processing
- **Fast pattern searching**: Leverages `memchr` with SIMD instructions when available
- **Optimized numeric parsing**: Implements `fast-float` for zero-allocation number conversion
- **String caching**: Prevents memory duplication for repeated strings like taxonomic rank codes
- **Smart memory pre-allocation**: Minimizes reallocations during processing
- **Zero-copy processing**: Works directly with byte slices when possible

### Key Libraries
- **memchr**: Provides ultra-fast character searching using SIMD when available
- **fast-float**: High-performance numeric parsing
- **rayon**: Enables parallel processing for operations like sequence extraction
- **serde_json**: Efficient JSON serialization/deserialization

These optimizations result in 3-5x better performance compared to traditional parsers, with predictable memory usage even for extremely large files (millions of lines).

## License

This project is licensed under the MIT License - see the LICENSE file for details.
