# Kraken2 Toolkit

Kraken2 Toolkit is a command-line utility written in [Rust](https://www.rust-lang.org/) for processing and analyzing [Kraken2](https://ccb.jhu.edu/software/kraken2/) bioinformatics software reports and log files. This toolkit focuses on fast and efficient processing of classifier outputs, thus being an comprehensive option for large datasets or bioinformatics pipelines as it works as a binary file with no external dependencies. This project is currently under development.

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

The basic syntax for using Kraken2 Toolkit is shown below ('kraken2-parser --help'):

   ```
   USAGE:
        kraken2-parser [SUBCOMMAND]

    OPTIONS:
        -h, --help       Print help information
        -V, --version    Print version information

    SUBCOMMANDS:
        analyze    Analyze Kraken2 report
        extract    Extract sequences based on Kraken2 results
        help       Print this message or the help of the given subcommand(s)
   ```
Analyze module:
   ```
   USAGE:
            kraken2-parser analyze [OPTIONS] <REPORT>

        ARGS:
            <REPORT>    Kraken2 report file

        OPTIONS:
                --alternative          Use alternative parser
                --compare <COMPARE>    Compare both parsers
                --filter <FILTER>      Filter results (e.g., by rank, percentage) (under development)
            -h, --help                 Print help information
                --info                 Display taxonomic or aggregated information (under development)
                --json <JSON>          Generate JSON output
                --search <SEARCH>      Search for a specific taxon by name (under development)
                --tax-id <TAXON_ID>    Search for a specific taxon by ID
   ```

Extract module:

   ```
    kraken2-parser-extract 
    Extract sequences based on Kraken2 results

    USAGE:
        kraken2-parser extract [OPTIONS] --output <OUTPUT> --taxids <TAXIDS> <SEQUENCE> <LOG>

    ARGS:
        <SEQUENCE>    Input FASTA/FASTQ file
        <LOG>         Kraken2 log file

    OPTIONS:
        -h, --help               Print help information
        -o, --output <OUTPUT>    Output file for extracted sequences
            --report <REPORT>    Kraken2 report file
            --taxids <TAXIDS>    Comma-separated list of taxids to extract
    ```