use std::fs::File;
use std::io::{BufRead, BufReader, Write};
use std::time::Instant;
use serde::Serialize;

#[derive(Debug, Clone, Serialize)]
pub struct TaxonEntry {
    pub level: usize,
    pub percentage: f32,
    pub clade_fragments: u64,
    pub direct_fragments: u64,
    pub rank_code: String,
    pub taxon_id: u64,
    pub name: String,
    pub children: Vec<TaxonEntry>,
}

#[derive(Debug, Serialize)]
pub struct KrakenReport {
    pub unclassified: Option<TaxonEntry>,
    pub root: Vec<TaxonEntry>,
}

pub fn parse_line(line: &str) -> Option<TaxonEntry> {
    let parts: Vec<&str> = line.split('\t').collect();
    if parts.len() < 6 {
        return None;
    }
    let level = parts[5].chars().take_while(|&c| c == ' ').count() / 2;
    let name = parts[5].trim_start().to_string();

    Some(TaxonEntry {
        level,
        percentage: parts[0].trim().parse().unwrap_or(0.0),
        clade_fragments: parts[1].trim().parse().unwrap_or(0),
        direct_fragments: parts[2].trim().parse().unwrap_or(0),
        rank_code: parts[3].trim().to_string(),
        taxon_id: parts[4].trim().parse().unwrap_or(0),
        name,
        children: Vec::new(),
    })
}

pub fn build_hierarchy(entries: Vec<TaxonEntry>) -> Vec<TaxonEntry> {
    let mut hierarchy: Vec<TaxonEntry> = Vec::new();
    let mut stack: Vec<*mut TaxonEntry> = Vec::new();

    for entry in entries {
        while let Some(&last) = stack.last() {
            if unsafe { (*last).level < entry.level } {
                break;
            }
            stack.pop();
        }

        let new_entry = if let Some(&parent) = stack.last() {
            unsafe {
                (*parent).children.push(entry);
                (*parent).children.last_mut().unwrap()
            }
        } else {
            hierarchy.push(entry);
            hierarchy.last_mut().unwrap()
        };

        stack.push(new_entry as *mut TaxonEntry);
    }

    hierarchy
}

pub fn parse_kraken2_report(file_path: &str) -> (KrakenReport, f64) {
    let start = Instant::now();
    let file = File::open(file_path).expect("Failed to open file");
    let reader = BufReader::new(file);
    let mut entries: Vec<TaxonEntry> = reader.lines()
        .filter_map(Result::ok)
        .filter_map(|line| parse_line(&line))
        .collect();

    let unclassified = if !entries.is_empty() && entries[0].level == 0 && entries[0].name == "unclassified" {
        Some(entries.remove(0))
    } else {
        None
    };

    let hierarchy_start = Instant::now();
    let root = build_hierarchy(entries);
    let hierarchy_time = hierarchy_start.elapsed().as_secs_f64();

    let total_time = start.elapsed().as_secs_f64();

    (KrakenReport {
        unclassified,
        root,
    }, hierarchy_time)
}

pub fn write_json_report(report: &KrakenReport, output_file: &str) -> std::io::Result<()> {
    let json = serde_json::to_string_pretty(report)?;
    let mut file = File::create(output_file)?;
    file.write_all(json.as_bytes())?;
    Ok(())
}

pub fn verify_hierarchy(entries: &[TaxonEntry]) -> bool {
    let root_sum: f32 = entries.iter()
        .map(|entry| entry.percentage)
        .sum();
    
    let result = (root_sum - 100.0).abs() < 0.01;
    
    result
}

fn print_hierarchy(entries: &[TaxonEntry], level: usize) {
    for entry in entries {
        println!("{}{}: {} ({}%)", "  ".repeat(level), entry.name, entry.clade_fragments, entry.percentage);
        print_hierarchy(&entry.children, level + 1);
    }
}