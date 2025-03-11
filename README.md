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

### Docker

You can also use the provided Dockerfile to build and run Kraken2 Toolkit in a container:

1. Build the Docker image:
   ```
   docker build -t kraken2-toolkit .
   ```

2. Or for an even smaller image using Alpine Linux:
   ```
   docker build -f Dockerfile.alpine -t kraken2-toolkit:alpine .
   ```

3. Run the Docker container:
   ```
   docker run --rm kraken2-toolkit
   ```

4. To process files, mount volumes:
   ```
   docker run --rm -v /path/to/local/data:/data kraken2-toolkit analyze /data/report.txt
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
        --report <REPORT>     Kraken2 report file (optional but recommended)
        --taxids <TAXIDS>     Comma-separated list of taxids to extract
```

The Extract module requires both the Kraken log file and optionally the report file for different purposes:

1. **Kraken log file (required)**: Maps individual sequences to their assigned taxon IDs
2. **Kraken report file (optional)**: Provides the complete taxonomic hierarchy, enabling more sophisticated extractions like:
   - Extracting all descendants of a taxon (not just directly assigned sequences)
   - Extracting sequences by complete clades (e.g., entire family or genus)
   - Verifying hierarchical relationships between requested taxa

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

The Kraken2 Toolkit is highly optimized for performance and memory efficiency. It implements several advanced techniques:

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
