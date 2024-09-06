use clap::{App, Arg, SubCommand};
use crate::krk_parser::{self, KrakenReport, TaxonEntry};
use crate::alternative_krk_parser;
use std::time::Instant;
use memory_stats::memory_stats;
use std::fs::File;
use std::io::{BufReader, BufRead};
use crate::taxon_query::{find_taxon_info, print_taxon_info};
use std::collections::HashSet;
use crate::logkrk_parser;
use crate::sequence_processor;

pub fn run_cli() {
    let matches = App::new("Kraken2 Toolkit")
        .version("0.2.0")
        .author("Your Name")
        .about("A toolkit for processing Kraken2 reports, logs, and sequence files")
        .subcommand(SubCommand::with_name("analyze")
            .about("Analyze Kraken2 report")
            .arg(Arg::with_name("REPORT")
                .help("Kraken2 report file")
                .required(true)
                .index(1))
            .arg(Arg::with_name("JSON")
                .help("Generate JSON output")
                .long("json")
                .takes_value(true))
            .arg(Arg::with_name("TAXON_ID")
                .help("Search for a specific taxon by ID")
                .long("tax-id")
                .takes_value(true))
            .arg(Arg::with_name("SEARCH")
                .help("Search for a specific taxon by name")
                .long("search")
                .takes_value(true))
            .arg(Arg::with_name("FILTER")
                .help("Filter results (e.g., by rank, percentage)")
                .long("filter")
                .takes_value(true))
            .arg(Arg::with_name("INFO")
                .help("Display taxonomic or aggregated information")
                .long("info")
                .takes_value(false))
            .arg(Arg::with_name("ALTERNATIVE")
                .help("Use alternative parser")
                .long("alternative")
                .takes_value(false))
            .arg(Arg::with_name("COMPARE")
                .help("Compare both parsers")
                .long("compare")
                .takes_value(true)))
        .subcommand(SubCommand::with_name("extract")
            .about("Extract sequences based on Kraken2 results")
            .arg(Arg::with_name("SEQUENCE")
                .help("Input FASTA/FASTQ file")
                .required(true)
                .index(1))
            .arg(Arg::with_name("LOG")
                .help("Kraken2 log file")
                .required(true)
                .index(2))
            .arg(Arg::with_name("REPORT")
                .help("Kraken2 report file")
                .long("report")
                .takes_value(true))
            .arg(Arg::with_name("OUTPUT")
                .help("Output file for extracted sequences")
                .required(true)
                .short('o')
                .long("output")
                .takes_value(true))
            .arg(Arg::with_name("TAXIDS")
                .help("Comma-separated list of taxids to extract")
                .required(true)
                .long("taxids")
                .takes_value(true)))
        .get_matches();

    match matches.subcommand() {
        Some(("analyze", analyze_matches)) => {
            run_analyze(analyze_matches);
        }
        Some(("extract", extract_matches)) => {
            run_extract(extract_matches);
        }
        _ => {
            println!("Please specify a subcommand. Use --help for more information.");
        }
    }
}

fn run_analyze(matches: &clap::ArgMatches) {
    let input_file = matches.value_of("REPORT").unwrap();

    if matches.is_present("COMPARE") {
        let iterations = matches.value_of("COMPARE")
            .and_then(|v| v.parse().ok())
            .unwrap_or(1);
        compare_parsers(input_file, matches, iterations);
    } else {
        let (report, hierarchy_time) = if matches.is_present("ALTERNATIVE") {
            alternative_krk_parser::parse_kraken2_report(input_file)
        } else {
            krk_parser::parse_kraken2_report(input_file)
        };

        println!("Hierarchy build time: {:.6} seconds", hierarchy_time);
        process_report(matches, &report);

        if let Some(taxon_id) = matches.value_of("TAXON_ID").and_then(|id| id.parse().ok()) {
            if let Some(taxon_info) = find_taxon_info(&report.root, taxon_id) {
                print_taxon_info(&taxon_info);
            } else {
                println!("Taxon with ID {} not found", taxon_id);
            }
        }

        if let Some(json_output) = matches.value_of("JSON") {
            match krk_parser::write_json_report(&report, json_output) {
                Ok(_) => println!("JSON report written to {}", json_output),
                Err(e) => eprintln!("Error writing JSON report: {}", e),
            }
        }
    }
}

fn run_extract(matches: &clap::ArgMatches) {
    let sequence_file = matches.value_of("SEQUENCE").unwrap();
    let log_file = matches.value_of("LOG").unwrap();
    let output_file = matches.value_of("OUTPUT").unwrap();
    let taxids: HashSet<String> = matches.value_of("TAXIDS").unwrap()
        .split(',')
        .map(String::from)
        .collect();

    match logkrk_parser::parse_kraken_output(log_file, &taxids) {
        Ok(save_readids) => {
            match sequence_processor::process_sequence_files(&[sequence_file.to_string()], &save_readids, output_file) {
                Ok(_) => println!("Sequences extracted successfully to {}", output_file),
                Err(e) => eprintln!("Error processing sequence file: {}", e),
            }
        }
        Err(e) => eprintln!("Error parsing Kraken log: {}", e),
    }
}

fn compare_parsers(input_file: &str, matches: &clap::ArgMatches, iterations: u32) {
    println!("Comparing parsers ({} iterations)...", iterations);

    // Función para medir el uso de memoria
    let get_memory_usage = || {
        memory_stats().map(|usage: memory_stats::MemoryStats| usage.physical_mem).unwrap_or(0)
    };

    let mut original_stats = (0.0, 0.0, 0.0, 0.0, 0);
    let mut alternative_stats = (0.0, 0.0, 0.0, 0.0, 0);

    for _ in 0..iterations {
        // Original parser
        let start_memory = get_memory_usage();
        let start = Instant::now();
        
        let file = File::open(input_file).expect("Failed to open file");
        let reader = BufReader::new(file);
        let parse_start = Instant::now();
        let entries: Vec<TaxonEntry> = reader.lines()
            .filter_map(Result::ok)
            .filter_map(|line| krk_parser::parse_line(&line))
            .collect();
        let parse_time = parse_start.elapsed().as_secs_f64();
        
        let hierarchy_start = Instant::now();
        let root = krk_parser::build_hierarchy(entries);
        let hierarchy_time = hierarchy_start.elapsed().as_secs_f64();
        
        let original_report = KrakenReport {
            unclassified: None,
            root,
        };
        
        let original_total_time = start.elapsed().as_secs_f64();
        let end_memory = get_memory_usage();

        // Alternative parser
        let alt_start_memory = get_memory_usage();
        let start = Instant::now();
        
        let file = File::open(input_file).expect("Failed to open file");
        let reader = BufReader::with_capacity(1024 * 1024, file);
        let parse_start = Instant::now();
        let entries: Vec<TaxonEntry> = reader.split(b'\n')
            .filter_map(Result::ok)
            .filter_map(|line| alternative_krk_parser::parse_line(&line))
            .collect();
        let parse_time_alt = parse_start.elapsed().as_secs_f64();
        
        let hierarchy_start = Instant::now();
        let root = alternative_krk_parser::build_hierarchy(entries);
        let hierarchy_time_alt = hierarchy_start.elapsed().as_secs_f64();
        
        let alternative_report = KrakenReport {
            unclassified: None,
            root,
        };
        
        let alternative_total_time = start.elapsed().as_secs_f64();
        let alt_end_memory = get_memory_usage();

        // JSON writing time measurement
        let mut json_time_original = 0.0;
        let mut json_time_alternative = 0.0;

        if matches.is_present("JSON") {
            if let Some(output_file) = matches.value_of("OUTPUT") {
                let original_output = format!("{}_original.json", output_file);
                let alternative_output = format!("{}_alternative.json", output_file);

                let json_start = Instant::now();
                let _ = krk_parser::write_json_report(&original_report, &original_output);
                json_time_original = json_start.elapsed().as_secs_f64();

                let json_start = Instant::now();
                let _ = krk_parser::write_json_report(&alternative_report, &alternative_output);
                json_time_alternative = json_start.elapsed().as_secs_f64();
            }
        }

        original_stats.0 += original_total_time;
        original_stats.1 += parse_time;
        original_stats.2 += hierarchy_time;
        original_stats.3 += json_time_original;
        original_stats.4 += end_memory - start_memory;

        alternative_stats.0 += alternative_total_time;
        alternative_stats.1 += parse_time_alt;
        alternative_stats.2 += hierarchy_time_alt;
        alternative_stats.3 += json_time_alternative;
        alternative_stats.4 += alt_end_memory - alt_start_memory;
    }

    // Calculate averages
    let iterations_f64 = iterations as f64;
    original_stats = (
        original_stats.0 / iterations_f64,
        original_stats.1 / iterations_f64,
        original_stats.2 / iterations_f64,
        original_stats.3 / iterations_f64,
        original_stats.4 / iterations as usize,
    );

    alternative_stats = (
        alternative_stats.0 / iterations_f64,
        alternative_stats.1 / iterations_f64,
        alternative_stats.2 / iterations_f64,
        alternative_stats.3 / iterations_f64,
        alternative_stats.4 / iterations as usize,
    );

    // Print results
    println!("\nOriginal parser (average of {} iterations):", iterations);
    println!("  Total time: {:.6} seconds", original_stats.0);
    println!("  File parsing time: {:.6} seconds", original_stats.1);
    println!("  Hierarchy build time: {:.6} seconds", original_stats.2);
    println!("  JSON writing time: {:.6} seconds", original_stats.3);
    println!("  Memory usage: {} bytes", original_stats.4);

    println!("\nAlternative parser (average of {} iterations):", iterations);
    println!("  Total time: {:.6} seconds", alternative_stats.0);
    println!("  File parsing time: {:.6} seconds", alternative_stats.1);
    println!("  Hierarchy build time: {:.6} seconds", alternative_stats.2);
    println!("  JSON writing time: {:.6} seconds", alternative_stats.3);
    println!("  Memory usage: {} bytes", alternative_stats.4);

    // Compare results
    println!("\nComparison (Alternative vs Original):");
    println!("  Total time: {:.2}%", (alternative_stats.0 / original_stats.0) * 100.0);
    println!("  File parsing time: {:.2}%", (alternative_stats.1 / original_stats.1) * 100.0);
    println!("  Hierarchy build time: {:.2}%", (alternative_stats.2 / original_stats.2) * 100.0);
    println!("  JSON writing time: {:.2}%", (alternative_stats.3 / original_stats.3) * 100.0);
    println!("  Memory usage: {:.2}%", (alternative_stats.4 as f64 / original_stats.4 as f64) * 100.0);
}

fn process_report(matches: &clap::ArgMatches, _report: &krk_parser::KrakenReport) {
    if matches.is_present("INFO") {
        display_info(_report);
    }

    if let Some(search_term) = matches.value_of("SEARCH") {
        search_taxon(_report, search_term);
    }

    if let Some(filter) = matches.value_of("FILTER") {
        apply_filter(_report, filter);
    }
}

fn display_info(_report: &krk_parser::KrakenReport) {
    println!("Taxonomic Information:");
    // Add aggregated information as needed
}

fn search_taxon(_report: &krk_parser::KrakenReport, search_term: &str) {
    println!("Searching for taxon: {}", search_term);
    // Implement search functionality
}

fn apply_filter(_report: &krk_parser::KrakenReport, filter: &str) {
    println!("Applying filter: {}", filter);
    // Implement filter functionality
}
