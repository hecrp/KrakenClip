use std::fs::File;
use std::io::Read;
use std::time::Instant;
use memchr::memchr_iter;
use fast_float;
use std::io::Write;
use std::path::Path;

// Optimized constants
const BUFFER_SIZE: usize = 512 * 1024; // 512KB is an optimal middle ground
const TAB_CHAR: u8 = b'\t';
const NL_CHAR: u8 = b'\n';
const SPACE_CHAR: u8 = b' ';
const NEWLINE_CHAR: u8 = b'\n';

// Structure representing a taxonomy node
#[derive(Debug, Clone)]
pub struct TaxonEntry {
    pub percentage: f32,   // Percentage of reads
    pub clade_reads: u64,  // Reads assigned to this clade
    pub direct_reads: u64, // Reads directly assigned
    pub rank: String,      // Taxonomic rank
    pub taxid: u32,        // Taxon ID
    pub name: String,      // Taxon name
    pub depth: usize,      // Depth in the tree
    pub children: Vec<TaxonEntry>, // Child nodes
    // Additional fields for compatibility
    pub level: usize,      // Hierarchy level (same as depth)
    pub taxon_id: u64,     // Taxon ID (same as taxid but as u64)
    pub clade_fragments: u64, // Alias for clade_reads
    pub direct_fragments: u64, // Alias for direct_reads
    pub rank_code: String, // Alias for rank
}

// Main structure for the Kraken report
#[derive(Debug)]
pub struct KrakenReport {
    pub root: TaxonEntry,
    pub taxon_map: std::collections::HashMap<u32, usize>,
    pub unclassified: Option<TaxonEntry>, // Additional field for unclassified data
}

impl TaxonEntry {
    // Create a new empty taxonomy node
    fn new(percentage: f32, clade_reads: u64, direct_reads: u64, rank: String, taxid: u32, name: String, depth: usize) -> Self {
        Self {
            percentage,
            clade_reads,
            direct_reads,
            rank: rank.clone(),
            taxid,
            name,
            depth,
            children: Vec::new(),
            // Initialize additional fields
            level: depth,
            taxon_id: taxid as u64,
            clade_fragments: clade_reads,
            direct_fragments: direct_reads,
            rank_code: rank,
        }
    }
    
    // Add a child to this node
    fn add_child(&mut self, child: TaxonEntry) {
        self.children.push(child);
    }
}

impl Default for KrakenReport {
    fn default() -> Self {
        let mut report = Self {
            root: TaxonEntry::new(0.0, 0, 0, "R".to_string(), 1, "root".to_string(), 0),
            taxon_map: std::collections::HashMap::new(),
            unclassified: None,
        };
        report.taxon_map.insert(1, 0); // Root always at position 0
        report
    }
}

// String cache to avoid memory duplication
struct StringCache {
    rank_codes: Vec<String>,
}

impl StringCache {
    fn new() -> Self {
        // Pre-populate with common rank codes
        let mut rank_codes = Vec::with_capacity(10);
        for code in ["D", "P", "C", "O", "F", "G", "S", "U"].iter() {
            rank_codes.push((*code).to_string());
        }
        Self { rank_codes }
    }

    fn get_rank_code(&mut self, code: &str) -> String {
        for existing in &self.rank_codes {
            if existing == code {
                return existing.clone();
            }
        }
        // If not in cache, add it
        let code_string = code.to_string();
        self.rank_codes.push(code_string.clone());
        code_string
    }
}

// Highly optimized line parser using pointers and avoiding unnecessary allocations
#[inline(always)]
pub fn parse_line(line: &[u8], string_cache: &mut StringCache) -> Option<TaxonEntry> {
    // Use fast tab search with memchr
    let tab_positions: Vec<usize> = memchr_iter(TAB_CHAR, line).take(5).collect();
    if tab_positions.len() < 5 {
        return None;
    }
    
    // Optimization: calculate field positions directly
    let field_starts = [0, tab_positions[0] + 1, tab_positions[1] + 1, tab_positions[2] + 1, tab_positions[3] + 1, tab_positions[4] + 1];
    let field_ends = [tab_positions[0], tab_positions[1], tab_positions[2], tab_positions[3], tab_positions[4], line.len()];
    
    // Fast calculation of indentation level
    let name_start = field_starts[5];
    let mut level = 0;
    let mut i = name_start;
    
    // Count spaces quickly in pairs
    while i + 1 < line.len() && line[i] == SPACE_CHAR && line[i + 1] == SPACE_CHAR {
        level += 1;
        i += 2;
    }
    
    // Extract name efficiently
    let name_bytes = &line[name_start + level * 2..];
    let name = std::str::from_utf8(name_bytes).unwrap_or("").trim();
    
    // Fast parsing of numeric values without unnecessary allocations
    let percentage_bytes = &line[field_starts[0]..field_ends[0]];
    let percentage = fast_float::parse::<f32, _>(std::str::from_utf8(percentage_bytes).unwrap_or("0.0")).unwrap_or(0.0);
    
    let clade_bytes = &line[field_starts[1]..field_ends[1]];
    let clade_fragments = std::str::from_utf8(clade_bytes).unwrap_or("0").parse::<u64>().unwrap_or(0);
    
    let direct_bytes = &line[field_starts[2]..field_ends[2]];
    let direct_fragments = std::str::from_utf8(direct_bytes).unwrap_or("0").parse::<u64>().unwrap_or(0);
    
    // Optimization for rank codes (avoid string duplication)
    let rank_bytes = &line[field_starts[3]..field_ends[3]];
    let rank_str = std::str::from_utf8(rank_bytes).unwrap_or("").trim();
    let rank_code = string_cache.get_rank_code(rank_str);
    
    let taxon_bytes = &line[field_starts[4]..field_ends[4]];
    let taxon_id = std::str::from_utf8(taxon_bytes).unwrap_or("0").parse::<u64>().unwrap_or(0);
    
    Some(TaxonEntry {
        level,
        percentage,
        clade_fragments,
        direct_fragments,
        rank_code: rank_code.clone(),
        taxid: taxon_id as u32,
        name: name.to_string(),
        depth: level,
        children: Vec::new(),
        taxon_id,
        clade_reads: clade_fragments,
        direct_reads: direct_fragments,
        rank: rank_code,
    })
}

// Building hierarchy with optimal memory reservation
pub fn build_hierarchy_optimized(entries: Vec<TaxonEntry>) -> Vec<TaxonEntry> {
    if entries.is_empty() {
        return Vec::new();
    }
    
    // Pre-allocate space for the final hierarchy
    let mut hierarchy = Vec::with_capacity(entries.len() / 8 + 1);
    
    // Use a stack with pre-allocated capacity (hierarchy rarely exceeds 20 levels)
    let mut stack: Vec<TaxonEntry> = Vec::with_capacity(20);
    
    for mut entry in entries {
        // Optimization to avoid chain reactions of push/pop
        while !stack.is_empty() && stack.last().unwrap().level >= entry.level {
            let popped = stack.pop().unwrap();
            if let Some(parent) = stack.last_mut() {
                // Reserve capacity if needed to avoid relocations
                if parent.children.is_empty() {
                    parent.children.reserve(4); // Most nodes have fewer than 4 children
                }
                parent.children.push(popped);
            } else {
                hierarchy.push(popped);
            }
        }
        stack.push(entry);
    }
    
    // Process remaining elements in the stack
    while let Some(entry) = stack.pop() {
        if let Some(parent) = stack.last_mut() {
            parent.children.push(entry);
        } else {
            hierarchy.push(entry);
        }
    }
    
    hierarchy
}

// Optimized buffer for block reading
struct OptimizedBuffer {
    buffer: Box<[u8]>,
    pos: usize,
    cap: usize,
    file: File,
}

impl OptimizedBuffer {
    fn new(file: File) -> Self {
        Self {
            buffer: vec![0; BUFFER_SIZE].into_boxed_slice(),
            pos: 0,
            cap: 0,
            file,
        }
    }
    
    fn fill_buffer(&mut self) -> std::io::Result<usize> {
        self.pos = 0;
        self.cap = self.file.read(&mut self.buffer)?;
        Ok(self.cap)
    }
    
    fn read_line(&mut self, line_buffer: &mut Vec<u8>) -> std::io::Result<bool> {
        line_buffer.clear();
        
        if self.pos >= self.cap && self.fill_buffer()? == 0 {
            return Ok(false);
        }
        
        loop {
            let mut i = self.pos;
            while i < self.cap {
                if self.buffer[i] == NEWLINE_CHAR {
                    line_buffer.extend_from_slice(&self.buffer[self.pos..i]);
                    self.pos = i + 1;
                    return Ok(true);
                }
                i += 1;
            }
            
            // No line break found, copy rest of buffer
            line_buffer.extend_from_slice(&self.buffer[self.pos..self.cap]);
            
            // Fill the buffer again
            if self.fill_buffer()? == 0 {
                return Ok(!line_buffer.is_empty());
            }
        }
    }
}

// Main analysis function with block processing
pub fn parse_kraken2_report(file_path: &str) -> (KrakenReport, f64) {
    let file = File::open(file_path).expect("Failed to open file");
    let file_size = file.metadata().map(|m| m.len() as usize).unwrap_or(0);
    
    let mut buffer = OptimizedBuffer::new(file);
    
    let start_time = Instant::now();
    let mut string_cache = StringCache::new();
    let mut line_buffer = Vec::with_capacity(1024);
    let mut entries = Vec::with_capacity(file_size / 50); // Approximate estimation
    
    while buffer.read_line(&mut line_buffer).unwrap_or(false) {
        if let Some(entry) = parse_line(&line_buffer, &mut string_cache) {
            entries.push(entry);
        }
    }
    
    // Optimized hierarchy construction
    let hierarchy = build_hierarchy_optimized(entries);
    
    // Determine if there's an "unclassified" node
    let unclassified = if !hierarchy.is_empty() && hierarchy[0].level == 0 && hierarchy[0].name == "unclassified" {
        Some(hierarchy[0].clone())
    } else {
        None
    };
    
    // Get root node
    let root = if hierarchy.len() > 1 || (hierarchy.len() == 1 && unclassified.is_some()) {
        let root_index = if unclassified.is_some() { 1 } else { 0 };
        if root_index < hierarchy.len() {
            hierarchy[root_index].clone()
        } else {
            TaxonEntry::new(0.0, 0, 0, "R".to_string(), 1, "root".to_string(), 0)
        }
    } else {
        TaxonEntry::new(0.0, 0, 0, "R".to_string(), 1, "root".to_string(), 0)
    };
    
    // Build taxon map
    let mut taxon_map = std::collections::HashMap::new();
    build_taxon_map(&root, &mut taxon_map, 0);
    
    let duration = start_time.elapsed().as_secs_f64();
    
    (KrakenReport {
        unclassified,
        root,
        taxon_map,
    }, duration)
}

// Build taxon map for quick access by ID
fn build_taxon_map(entry: &TaxonEntry, map: &mut std::collections::HashMap<u32, usize>, index: usize) -> usize {
    map.insert(entry.taxid, index);
    
    let mut next_index = index + 1;
    for child in &entry.children {
        next_index = build_taxon_map(child, map, next_index);
    }
    next_index
}

// Function to write the report in JSON format
pub fn write_json_report(report: &KrakenReport, output_path: &str) -> std::io::Result<()> {
    use serde_json;
    
    // Recursive function to convert a node to JSON format
    fn node_to_json(node: &TaxonEntry) -> serde_json::Value {
        let mut json = serde_json::json!({
            "name": node.name,
            "taxid": node.taxid,
            "rank": node.rank,
            "percentage": node.percentage,
            "clade_reads": node.clade_reads,
            "direct_reads": node.direct_reads,
            "level": node.level,
            "children": []
        });
        
        // Add children recursively
        let children = node.children.iter().map(node_to_json).collect::<Vec<_>>();
        json["children"] = serde_json::Value::Array(children);
        
        json
    }
    
    // Convert the entire report to JSON
    let mut json = serde_json::json!({
        "root": node_to_json(&report.root)
    });
    
    // Add unclassified node if it exists
    if let Some(ref unclassified) = report.unclassified {
        json["unclassified"] = node_to_json(unclassified);
    }
    
    // Write to file
    let file = std::fs::File::create(Path::new(output_path))?;
    let mut writer = std::io::BufWriter::new(file);
    
    serde_json::to_writer_pretty(&mut writer, &json)?;
    writer.flush()?;
    
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    
    #[test]
    fn test_parse_line() {
        let test_line = b"50.00\t1000\t500\tP\t123\t  Bacteria";
        let mut cache = StringCache::new();
        let result = parse_line(test_line, &mut cache);
        assert!(result.is_some());
        let entry = result.unwrap();
        assert_eq!(entry.level, 1);
        assert_eq!(entry.percentage, 50.0);
        assert_eq!(entry.clade_fragments, 1000);
        assert_eq!(entry.direct_fragments, 500);
        assert_eq!(entry.rank_code, "P");
        assert_eq!(entry.taxid, 123);
        assert_eq!(entry.name, "Bacteria");
    }
    
    #[test]
    fn test_build_hierarchy() {
        let entries = vec![
            TaxonEntry {
                level: 0,
                percentage: 100.0,
                clade_fragments: 1000,
                direct_fragments: 0,
                rank_code: "D".to_string(),
                taxid: 1,
                name: "Root".to_string(),
                depth: 0,
                children: Vec::new(),
            },
            TaxonEntry {
                level: 1,
                percentage: 80.0,
                clade_fragments: 800,
                direct_fragments: 200,
                rank_code: "P".to_string(),
                taxid: 2,
                name: "Bacteria".to_string(),
                depth: 1,
                children: Vec::new(),
            },
            TaxonEntry {
                level: 2,
                percentage: 60.0,
                clade_fragments: 600,
                direct_fragments: 100,
                rank_code: "C".to_string(),
                taxid: 3,
                name: "Proteobacteria".to_string(),
                depth: 2,
                children: Vec::new(),
            },
        ];
        
        let hierarchy = build_hierarchy_optimized(entries);
        assert_eq!(hierarchy.len(), 1);
        assert_eq!(hierarchy[0].name, "Root");
        assert_eq!(hierarchy[0].children.len(), 1);
        assert_eq!(hierarchy[0].children[0].name, "Bacteria");
        assert_eq!(hierarchy[0].children[0].children.len(), 1);
        assert_eq!(hierarchy[0].children[0].children[0].name, "Proteobacteria");
    }
} 