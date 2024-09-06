use crate::krk_parser::{KrakenReport, TaxonEntry};
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::time::Instant;
use memchr::memchr;
use fast_float;

// Elimina la redefinición de TaxonEntry

pub fn parse_line(line: &[u8]) -> Option<TaxonEntry> {
    let mut fields = [0; 6];
    let mut field_index = 0;
    let mut last_tab = 0;

    while field_index < 5 {
        match memchr(b'\t', &line[last_tab..]) {
            Some(pos) => {
                fields[field_index] = last_tab;
                field_index += 1;
                last_tab += pos + 1;
            }
            None => return None,
        }
    }

    fields[5] = last_tab;

    let level = line[fields[5]..].iter().take_while(|&&b| b == b' ').count() / 2;
    let name = std::str::from_utf8(&line[fields[5] + level * 2..]).unwrap_or("").trim();

    // Manejo especial para el campo de porcentaje
    let percentage_str = std::str::from_utf8(&line[fields[0]..fields[1] - 1]).unwrap_or("").trim();
    let percentage = percentage_str.parse::<f32>().unwrap_or(0.0);

    Some(TaxonEntry {
        level,
        percentage,
        clade_fragments: fast_float::parse::<f64, _>(std::str::from_utf8(&line[fields[1]..fields[2] - 1]).unwrap_or("")).unwrap_or(0.0) as u64,
        direct_fragments: fast_float::parse::<f64, _>(std::str::from_utf8(&line[fields[2]..fields[3] - 1]).unwrap_or("")).unwrap_or(0.0) as u64,
        rank_code: std::str::from_utf8(&line[fields[3]..fields[4] - 1]).unwrap_or("").trim().to_string(),
        taxon_id: fast_float::parse::<f64, _>(std::str::from_utf8(&line[fields[4]..fields[5] - 1]).unwrap_or("")).unwrap_or(0.0) as u64,
        name: name.to_string(),
        children: Vec::new(),
    })
}

pub fn build_hierarchy(entries: Vec<TaxonEntry>) -> Vec<TaxonEntry> {
    let mut stack: Vec<TaxonEntry> = Vec::new();
    let mut result: Vec<TaxonEntry> = Vec::new();

    for node in entries {
        while let Some(last) = stack.last() {
            if last.level >= node.level {
                let popped = stack.pop().unwrap();
                if let Some(parent) = stack.last_mut() {
                    parent.children.push(popped);
                } else {
                    result.push(popped);
                }
            } else {
                break;
            }
        }
        stack.push(node);
    }

    while let Some(node) = stack.pop() {
        if let Some(parent) = stack.last_mut() {
            parent.children.push(node);
        } else {
            result.push(node);
        }
    }

    result
}

pub fn build_hierarchy_unsafe(entries: Vec<TaxonEntry>) -> Vec<TaxonEntry> {
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
    let file = File::open(file_path).expect("Failed to open file");
    let reader = BufReader::with_capacity(1024 * 1024, file);
    let mut entries: Vec<TaxonEntry> = reader.split(b'\n')
        .filter_map(|line| line.ok())
        .filter_map(|line| parse_line(&line))
        .collect();

    let unclassified = if !entries.is_empty() && entries[0].level == 0 && entries[0].name == "unclassified" {
        Some(entries.remove(0))
    } else {
        None
    };

    let hierarchy_start = Instant::now();
    // let root = build_hierarchy(entries);
    let root = build_hierarchy_unsafe(entries);
    let hierarchy_time = hierarchy_start.elapsed().as_secs_f64();

    (KrakenReport {
        unclassified,
        root,
    }, hierarchy_time)
}

// Elimina la función write_json_report, ya que usaremos la de krk_parser