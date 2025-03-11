use clap::{App, Arg, SubCommand};
use crate::krk_parser::{self, KrakenReport, TaxonEntry};
use std::time::Instant;
use memory_stats::memory_stats;
use std::fs::File;
use std::io::{BufReader, BufRead};
use crate::taxon_query::{find_taxon_info, print_taxon_info};
use std::collections::HashSet;
use crate::logkrk_parser;
use crate::sequence_processor;
use crate::generate_test_data;

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
                .takes_value(false)))
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
        .subcommand(SubCommand::with_name("generate-test-data")
            .about("Generate test data for performance testing")
            .arg(Arg::with_name("OUTPUT")
                .help("Output file for generated test data")
                .required(true)
                .short('o')
                .long("output")
                .takes_value(true))
            .arg(Arg::with_name("LINES")
                .help("Number of lines to generate")
                .long("lines")
                .default_value("100000")
                .takes_value(true))
            .arg(Arg::with_name("TYPE")
                .help("Type of test data to generate (wide, deep, fragments, dense, etc.)")
                .long("type")
                .default_value("random")
                .takes_value(true)))
        .get_matches();

    match matches.subcommand() {
        Some(("analyze", analyze_matches)) => {
            run_analyze(analyze_matches);
        }
        Some(("extract", extract_matches)) => {
            run_extract(extract_matches);
        }
        Some(("generate-test-data", generate_matches)) => {
            run_generate(generate_matches);
        }
        _ => {
            println!("Please specify a subcommand. Use --help for more information.");
        }
    }
}

fn run_analyze(matches: &clap::ArgMatches) {
    let input_file = matches.value_of("REPORT").unwrap();
    
    let start = Instant::now();
    let parse_start = Instant::now();
    
    // Obtener uso de memoria inicial
    let start_memory = memory_stats().map(|usage: memory_stats::MemoryStats| usage.physical_mem).unwrap_or(0);
    
    let (report, hierarchy_time) = krk_parser::parse_kraken2_report(input_file);
    
    let parse_time = parse_start.elapsed().as_secs_f64() - hierarchy_time;
    let total_time = start.elapsed().as_secs_f64();
    
    // Obtener uso de memoria final
    let end_memory = memory_stats().map(|usage: memory_stats::MemoryStats| usage.physical_mem).unwrap_or(0);
    let memory_used = end_memory - start_memory;

    println!("Hierarchy build time: {:.6} seconds", hierarchy_time);
    
    // Añadir métricas de rendimiento adicionales
    println!("Total time: {:.6} seconds", total_time);
    println!("File parsing time: {:.6} seconds", parse_time);
    println!("Memory usage: {} bytes", memory_used);
    
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

fn run_generate(matches: &clap::ArgMatches) {
    let output_file = matches.value_of("OUTPUT").unwrap();
    
    let num_lines = matches.value_of("LINES")
        .and_then(|v| v.parse().ok())
        .unwrap_or(100_000);
        
    let data_type = matches.value_of("TYPE").unwrap_or("random");
    
    match generate_test_data::generate_data(output_file, num_lines, data_type) {
        Ok(_) => println!("Test data generated successfully"),
        Err(e) => eprintln!("Error generating test data: {}", e),
    }
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

