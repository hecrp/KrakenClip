use clap::{App, Arg, SubCommand};
use crate::krk_parser::{self, KrakenReport, TaxonEntry};
use std::time::Instant;
use memory_stats::memory_stats;
use std::fs::File;
use std::io::{BufReader, BufRead, Write, BufWriter};
use crate::taxon_query::{find_taxon_info, print_taxon_info};
use std::collections::{HashSet, HashMap};
use crate::logkrk_parser;
use crate::sequence_processor;
use crate::generate_test_data;
use std::path::Path;
use chrono;

pub fn run_cli() {
    let matches = App::new("KrakenClip")
        .version("0.2.0")
        .author("Your Name")
        .about("A high-performance toolkit for processing Kraken2 reports, logs, and sequence files")
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
                .takes_value(true))
            .arg(Arg::with_name("INCLUDE_CHILDREN")
                .help("Include sequences from all descendant taxa")
                .long("include-children")
                .takes_value(false))
            .arg(Arg::with_name("INCLUDE_PARENTS")
                .help("Include sequences from all ancestor taxa")
                .long("include-parents")
                .takes_value(false))
            .arg(Arg::with_name("STATS_OUTPUT")
                .help("Generate a statistics file with detailed extraction information")
                .long("stats-output")
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
    
    // Get initial memory usage
    let start_memory = memory_stats().map(|usage: memory_stats::MemoryStats| usage.physical_mem).unwrap_or(0);
    
    let (report, hierarchy_time) = krk_parser::parse_kraken2_report(input_file);
    
    let parse_time = parse_start.elapsed().as_secs_f64() - hierarchy_time;
    let total_time = start.elapsed().as_secs_f64();
    
    // Get final memory usage
    let end_memory = memory_stats().map(|usage: memory_stats::MemoryStats| usage.physical_mem).unwrap_or(0);
    let memory_used = end_memory - start_memory;

    println!("Hierarchy build time: {:.6} seconds", hierarchy_time);
    
    // Add additional performance metrics
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
    
    // Calculate approximate number of taxids to pre-allocate memory
    let taxids_raw = matches.value_of("TAXIDS").unwrap();
    let taxids_count = taxids_raw.split(',').count();
    
    // Pre-allocate memory for the sets
    let mut taxids: HashSet<String> = HashSet::with_capacity(taxids_count * 2); // *2 for margin
    
    // Fill the taxids set
    for taxid in taxids_raw.split(',') {
        taxids.insert(String::from(taxid));
    }
    
    let include_children = matches.is_present("INCLUDE_CHILDREN");
    let include_parents = matches.is_present("INCLUDE_PARENTS");
    
    // Map to track which taxid each readid belongs to
    // We estimate an average of 1000 sequences per taxid
    let mut taxid_readid_map: HashMap<String, HashSet<String>> = HashMap::with_capacity(taxids_count * 2);
    
    // Save original taxids before expansion for tracking
    let original_taxids: HashSet<String> = taxids.clone();
    
    // If we need to include children or parents, we need the report file
    if (include_children || include_parents) && matches.value_of("REPORT").is_some() {
        let report_file = matches.value_of("REPORT").unwrap();
        println!("Reading taxonomy from report file: {}", report_file);
        
        // Parse the Kraken report to get the taxonomic tree
        let (report, _) = krk_parser::parse_kraken2_report(report_file);
        
        // Expand the taxids set based on options
        if include_children || include_parents {
            // Estimate initial size for expanded set
            // For large taxonomies, it could be 10x the original size
            let estimated_size = if taxids.len() < 10 {
                // For small sets of taxids, we might have many children
                taxids.len() * 100
            } else {
                // For larger sets, relative expansion is usually smaller
                taxids.len() * 10
            };
            
            let mut expanded_taxids = HashSet::with_capacity(estimated_size);
            
            for taxid in &taxids {
                // Add the original taxid
                expanded_taxids.insert(taxid.clone());
                
                if let Ok(taxid_num) = taxid.parse::<u64>() {
                    // If include_children is true, add all descendant taxids
                    if include_children {
                        collect_all_children(&report.root, taxid_num, &mut expanded_taxids);
                        println!("Including children for taxid: {}", taxid);
                    }
                    
                    // If include_parents is true, add all ancestor taxids
                    if include_parents {
                        collect_all_parents(&report.root, taxid_num, &mut expanded_taxids);
                        println!("Including parents for taxid: {}", taxid);
                    }
                }
            }
            
            // Replace original taxids with expanded set
            taxids = expanded_taxids;
        }
        
        let additional_count = taxids.len() - original_taxids.len();
        println!("Expanded to {} taxids (added {} through hierarchy)", taxids.len(), additional_count);
    } else if (include_children || include_parents) && matches.value_of("REPORT").is_none() {
        eprintln!("Warning: --include-children or --include-parents specified, but no report file provided. Hierarchy expansion disabled.");
    }

    match logkrk_parser::parse_kraken_output_with_taxids(log_file, &taxids, &mut taxid_readid_map) {
        Ok(save_readids) => {
            // Track total sequences in the input file
            let total_sequences = count_sequences_in_file(sequence_file);
            
            // Extract the sequences
            match sequence_processor::process_sequence_files(&[sequence_file.to_string()], &save_readids, output_file) {
                Ok(_) => {
                    println!("Sequences extracted successfully to {}", output_file);
                    println!("Extracted sequences matching {} taxids", taxids.len());
                    println!("Total extracted sequences: {}", save_readids.len());
                    
                    // Generate statistics file if requested
                    if let Some(stats_file) = matches.value_of("STATS_OUTPUT") {
                        match generate_statistics_file(
                            stats_file,
                            &taxid_readid_map,
                            &original_taxids,
                            &taxids,
                            save_readids.len(),
                            total_sequences
                        ) {
                            Ok(_) => println!("Statistics written to {}", stats_file),
                            Err(e) => eprintln!("Error writing statistics: {}", e),
                        }
                    }
                },
                Err(e) => eprintln!("Error processing sequence file: {}", e),
            }
        }
        Err(e) => eprintln!("Error parsing Kraken log: {}", e),
    }
}

// Count total sequences in a FASTA/FASTQ file
fn count_sequences_in_file(file_path: &str) -> usize {
    let file = match File::open(file_path) {
        Ok(file) => file,
        Err(_) => return 0,
    };
    
    // Use a large buffer (512 KB) for efficient reading
    let mut reader = BufReader::with_capacity(512 * 1024, file);
    let mut buffer = Vec::with_capacity(512 * 1024);
    
    // Detect file format by reading the first character
    let first_byte = match reader.fill_buf() {
        Ok(buf) if !buf.is_empty() => buf[0],
        _ => return 0, // Empty file or error
    };
    
    // Counters for statistics
    let mut count = 0;
    
    // Identify bytes for sequence identifiers
    let fasta_marker = b'>';
    let fastq_marker = b'@';
    
    match first_byte {
        b'>' => {
            // FASTA file
            // Count each line that starts with '>'
            loop {
                buffer.clear();
                let bytes_read = match reader.read_until(b'\n', &mut buffer) {
                    Ok(n) => n,
                    Err(_) => break,
                };
                
                if bytes_read == 0 {
                    break;
                }
                
                if !buffer.is_empty() && buffer[0] == fasta_marker {
                    count += 1;
                }
            }
        },
        b'@' => {
            // FASTQ file
            // In FASTQ each record has 4 lines, only count the first
            let mut line_num = 0;
            
            loop {
                buffer.clear();
                let bytes_read = match reader.read_until(b'\n', &mut buffer) {
                    Ok(n) => n,
                    Err(_) => break,
                };
                
                if bytes_read == 0 {
                    break;
                }
                
                if line_num % 4 == 0 && !buffer.is_empty() && buffer[0] == fastq_marker {
                    count += 1;
                }
                
                line_num += 1;
            }
        },
        _ => {
            // Unknown format
            eprintln!("Warning: Unrecognized sequence file format, assuming FASTA");
            // Try to process as FASTA anyway
            loop {
                buffer.clear();
                let bytes_read = match reader.read_until(b'\n', &mut buffer) {
                    Ok(n) => n,
                    Err(_) => break,
                };
                
                if bytes_read == 0 {
                    break;
                }
                
                if !buffer.is_empty() && buffer[0] == fasta_marker {
                    count += 1;
                }
            }
        }
    }
    
    count
}

// Generate comprehensive statistics file
fn generate_statistics_file(
    output_file: &str,
    taxid_readid_map: &HashMap<String, HashSet<String>>,
    original_taxids: &HashSet<String>,
    expanded_taxids: &HashSet<String>,
    total_extracted: usize,
    total_sequences: usize
) -> Result<(), Box<dyn std::error::Error>> {
    // Create a large buffer for writing
    let file = File::create(output_file)?;
    let mut writer = BufWriter::with_capacity(256 * 1024, file);
    
    // Calculate statistics once to avoid repeated calculations
    let extraction_percentage = (total_extracted as f64 / total_sequences as f64) * 100.0;
    let additional_taxids = expanded_taxids.len() - original_taxids.len();
    
    // Counters for summary statistics
    let mut original_count = 0;
    let mut expanded_count = 0;
    
    // Pre-calculate taxid counts and statistics at once
    // This avoids traversing the HashMap multiple times
    let mut taxid_stats: Vec<(String, bool, usize, f64, f64)> = Vec::with_capacity(taxid_readid_map.len());
    
    for (taxid, readids) in taxid_readid_map {
        let is_original = original_taxids.contains(taxid);
        let count = readids.len();
        let percent_of_extracted = (count as f64 / total_extracted as f64) * 100.0;
        let percent_of_total = (count as f64 / total_sequences as f64) * 100.0;
        
        // Update counters
        if is_original {
            original_count += count;
        } else {
            expanded_count += count;
        }
        
        taxid_stats.push((
            taxid.clone(),
            is_original,
            count,
            percent_of_extracted,
            percent_of_total
        ));
    }
    
    // Sort by number of sequences (descending)
    taxid_stats.sort_by(|a, b| b.2.cmp(&a.2));
    
    // Calculate percentages for the summary
    let original_percent = if total_extracted > 0 {
        (original_count as f64 / total_extracted as f64) * 100.0
    } else {
        0.0
    };
    
    let expanded_percent = if total_extracted > 0 {
        (expanded_count as f64 / total_extracted as f64) * 100.0
    } else {
        0.0
    };
    
    // Write metadata as comments (lines starting with #)
    // These will be recognized as comments by pandas, R, and other tools
    writeln!(writer, "# Extraction Statistics")?;
    writeln!(writer, "# Generated: {}", chrono::Local::now().format("%Y-%m-%d %H:%M:%S"))?;
    writeln!(writer, "# ")?;
    writeln!(writer, "# Summary:")?;
    writeln!(writer, "# Total sequences in input: {}", total_sequences)?;
    writeln!(writer, "# Total sequences extracted: {} ({:.2}%)", 
        total_extracted, 
        extraction_percentage
    )?;
    writeln!(writer, "# ")?;
    writeln!(writer, "# Taxonomic IDs summary:")?;
    writeln!(writer, "# Original taxids: {}", original_taxids.len())?;
    writeln!(writer, "# Expanded taxids: {}", expanded_taxids.len())?;
    writeln!(writer, "# Additional taxids from hierarchy: {}", additional_taxids)?;
    writeln!(writer, "# ")?;
    writeln!(writer, "# Original taxids yielded {} sequences ({:.2}% of extracted)", 
        original_count, original_percent)?;
    
    if expanded_taxids.len() > original_taxids.len() {
        writeln!(writer, "# Expanded taxids yielded {} sequences ({:.2}% of extracted)",
            expanded_count, expanded_percent)?;
    }
    
    // Write the CSV header - common CSV format with clear field names
    // Using descriptive column names compatible with data analysis software
    writeln!(writer, "taxid,origin,sequences,percent_of_extracted,percent_of_total")?;
    
    // Write detailed statistics as CSV rows
    for (taxid, is_original, count, percent_extracted, percent_total) in &taxid_stats {
        let origin = if *is_original { "Original" } else { "Expanded" };
        
        writeln!(writer, "{},{},{},{:.2},{:.2}",
            taxid, origin, count, percent_extracted, percent_total)?;
    }
    
    // Ensure all data is written
    writer.flush()?;
    
    Ok(())
}

// Recursively collect all child taxids
fn collect_all_children(node: &krk_parser::TaxonEntry, target_taxid: u64, taxids: &mut HashSet<String>) {
    // Pre-allocate space to avoid reallocations during collection
    // We estimate each node could have approximately 10 children with a tree depth of 10
    if taxids.capacity() < taxids.len() + 100 {
        taxids.reserve(100);
    }
    
    // Check if this node is the target
    if node.taxon_id == target_taxid {
        // If found, add all children recursively
        add_children_recursive(node, taxids);
        return;
    }
    
    // Otherwise, search in children
    for child in &node.children {
        collect_all_children(child, target_taxid, taxids);
    }
}

// Add all children of a node to the taxids set
fn add_children_recursive(node: &krk_parser::TaxonEntry, taxids: &mut HashSet<String>) {
    // Reserve space to avoid frequent reallocations
    let children_count = node.children.len();
    if taxids.capacity() < taxids.len() + children_count {
        taxids.reserve(children_count);
    }
    
    for child in &node.children {
        taxids.insert(child.taxon_id.to_string());
        add_children_recursive(child, taxids);
    }
}

// Recursively find all parent taxids
fn collect_all_parents(node: &krk_parser::TaxonEntry, target_taxid: u64, taxids: &mut HashSet<String>) {
    // Use a vector with pre-allocated capacity for the path
    // Most taxonomies don't exceed 30-40 levels of depth
    let mut path = Vec::with_capacity(40);
    find_parents_recursive(node, target_taxid, &mut path, taxids);
}

// Find parents using DFS with path tracking
fn find_parents_recursive(
    node: &krk_parser::TaxonEntry, 
    target_taxid: u64, 
    path: &mut Vec<u64>, 
    taxids: &mut HashSet<String>
) -> bool {
    // Add current node to the path
    path.push(node.taxon_id);
    
    // If this is the target node, add all parents from the path
    if node.taxon_id == target_taxid {
        // Reserve space for the number of parents to be added
        if taxids.capacity() < taxids.len() + path.len() {
            taxids.reserve(path.len());
        }
        
        // Add all nodes from the path except the last one (current node)
        for &parent_id in path.iter().take(path.len() - 1) {
            taxids.insert(parent_id.to_string());
        }
        
        path.pop(); // Restore path for backtracking
        return true;
    }
    
    // Search in children
    for child in &node.children {
        if find_parents_recursive(child, target_taxid, path, taxids) {
            path.pop(); // Restore path for backtracking
            return true;
        }
    }
    
    // Target node not found in this subtree
    path.pop(); // Restore path for backtracking
    false
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

