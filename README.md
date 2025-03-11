# Kraken2 Toolkit

Kraken2 Toolkit is a command-line utility written in [Rust](https://www.rust-lang.org/) for processing and analyzing [Kraken2](https://ccb.jhu.edu/software/kraken2/) bioinformatics software reports and log files. This toolkit focuses on fast and efficient processing of classifier outputs, thus being a comprehensive option for large datasets or bioinformatics pipelines as it works as a binary file with no external dependencies.

## Getting Started

### Prerequisites

- Rust programming language (latest stable version)
- Cargo (Rust's package manager)

### Installation

1. Clone the repository:
   ```
   git clone https://github.com/yourusername/kraken2-toolkit.git
   cd kraken2-toolkit
   ```

2. Build the project:
   ```
   cargo build --release
   ```

3. The executable will be available in `target/release/kraken2-parser`
    ```
    ./target/release/kraken2-parser --help
    ```

## Usage

The basic syntax for using Kraken2 Toolkit is shown below:

```
USAGE:
    kraken2-parser [SUBCOMMAND]

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
    kraken2-parser analyze [OPTIONS] <REPORT>

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
    kraken2-parser extract [OPTIONS] --output <OUTPUT> --taxids <TAXIDS> <SEQUENCE> <LOG>

ARGS:
    <SEQUENCE>    Input FASTA/FASTQ file
    <LOG>         Kraken2 log file

OPTIONS:
    -h, --help                Print help information
    -o, --output <OUTPUT>     Output file for extracted sequences
        --report <REPORT>     Kraken2 report file
        --taxids <TAXIDS>     Comma-separated list of taxids to extract
```

### Generate Test Data Module

Used to generate test data for benchmarking and testing:

```
USAGE:
    kraken2-parser generate-test-data [OPTIONS] --output <OUTPUT> --lines <LINES> --type <TYPE>

OPTIONS:
    -h, --help                Print help information
    -o, --output <OUTPUT>     Output file path
    -l, --lines <LINES>       Number of lines to generate
    -t, --type <TYPE>         Type of data to generate (wide, deep, fragments, dense, etc.)
```

## Performance

The Kraken2 Toolkit is highly optimized for maximum performance and memory efficiency. The current implementation leverages advanced optimization techniques:

### Optimization Techniques

- **Block-based Reading**: Optimized 512KB buffer for reading files in optimal-sized chunks, reducing system calls
- **Accelerated Pattern Search**: Uses the `memchr` algorithm for ultra-fast character searching in byte arrays
- **Optimized Numeric Parsing**: Implements `fast_float` for text-to-number conversion without temporary allocations
- **String Caching**: Avoids memory duplication for repeated strings like taxonomic rank codes
- **Smart Memory Pre-allocation**: Intelligently reserves memory for vectors and structures, minimizing reallocations
- **Zero-copy Processing**: Works directly with byte slices in memory when possible
- **Optimized Hierarchical Algorithms**: Single-pass taxonomy construction

### Key Libraries

- **memchr**: Provides ultra-fast byte searching using SIMD instructions when available
- **fast-float**: High-performance numeric parsing library
- **rayon**: Enables parallel processing for operations like sequence extraction
- **serde_json**: Efficient JSON serialization/deserialization

These optimizations allow Kraken2 Toolkit to process massive taxonomic report files (millions of lines) in seconds, where simpler implementations might take minutes.

## License

This project is licensed under [LICENSE NAME] - see the LICENSE file for details.