use std::collections::{HashSet, HashMap};
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::error::Error;
use memchr::memchr;

// Optimized buffer size (512 KB) for efficient I/O operations
// Larger buffers reduce system calls and improve throughput for sequential reads
const BUFFER_SIZE: usize = 512 * 1024;

// Constant characters for fast pattern matching operations
// Using raw byte literals for memory efficiency and faster comparisons
const TAB_CHAR: u8 = b'\t';
const LF_CHAR: u8 = b'\n';

/// Parse a Kraken output file and return a set of read IDs that match the specified tax IDs.
/// This version is kept for backward compatibility.
pub fn parse_kraken_output(kraken_output: &str, save_taxids: &HashSet<String>) -> Result<HashSet<String>, Box<dyn Error>> {
    let mut taxid_readid_map = HashMap::new();
    let result = parse_kraken_output_with_taxids(kraken_output, save_taxids, &mut taxid_readid_map)?;
    Ok(result)
}

/// Parse a Kraken output file and return a set of read IDs that match the specified tax IDs.
/// Additionally, track which tax ID each read ID corresponds to in the taxid_readid_map.
/// This version is optimized for performance using zero-copy operations when possible.
pub fn parse_kraken_output_with_taxids(
    kraken_output: &str, 
    save_taxids: &HashSet<String>,
    taxid_readid_map: &mut HashMap<String, HashSet<String>>
) -> Result<HashSet<String>, Box<dyn Error>> {
    // Open file and create a buffered reader with optimized buffer size
    let file = File::open(kraken_output)?;
    let mut reader = BufReader::with_capacity(BUFFER_SIZE, file);
    let mut save_readids = HashSet::new();
    
    // Pre-allocate space in HashSets to avoid reallocations during insertion
    // Capacity planning: assume each taxid may have ~1000 associated sequences
    // This significantly reduces hash table resizing operations which are costly
    save_readids.reserve(save_taxids.len() * 1000);
    
    // Reusable buffer to minimize heap allocations during file reading
    // This prevents allocating a new buffer for each line read
    let mut buffer = Vec::with_capacity(BUFFER_SIZE);
    
    loop {
        // Clear buffer before reading next line
        // This reuses allocated memory instead of creating a new buffer
        buffer.clear();
        
        // Read until newline character, returning number of bytes read
        // The `?` operator is Rust's syntax for error propagation
        let bytes_read = reader.read_until(LF_CHAR, &mut buffer)?;
        if bytes_read == 0 {
            // End of file reached
            break;
        }
        
        // Remove trailing newline if present for cleaner text processing
        if buffer.last() == Some(&LF_CHAR) {
            buffer.pop();
        }
        
        // Find tab positions using optimized memchr library
        // This is significantly faster than manual iteration or string splitting
        // The memchr crate uses SIMD instructions when available for vectorized searching
        let mut tab_positions = Vec::with_capacity(4); // Kraken logs typically have 4-5 fields
        let mut pos = 0;
        
        // Collect positions of tab characters for field extraction
        // This avoids string splits which would create multiple new allocations
        while let Some(offset) = memchr(TAB_CHAR, &buffer[pos..]) {
            let absolute_pos = pos + offset;
            tab_positions.push(absolute_pos);
            pos = absolute_pos + 1;
            
            // Early exit once we have enough tabs to find what we need
            // This avoids unnecessary work for very long lines
            if tab_positions.len() >= 3 {
                break;
            }
        }
        
        // We need at least 2 tabs to get readid and taxid
        // This guards against malformed input
        if tab_positions.len() < 2 {
            continue;
        }
        
        // Extract taxid (field 3)
        // Using byte slices and ranges for zero-copy operations
        let taxid_start = tab_positions[1] + 1;
        let taxid_end = if tab_positions.len() > 2 { 
            tab_positions[2] 
        } else { 
            buffer.len() 
        };
        
        // Convert byte slice to UTF-8 string without allocating a new String
        // This is a zero-copy operation until we need to store the value
        let taxid = std::str::from_utf8(&buffer[taxid_start..taxid_end])?;
        
        // Check if this taxid is one we're interested in
        if save_taxids.contains(taxid) {
            // Extract readid (field 2) using the same efficient byte slicing
            let readid_start = tab_positions[0] + 1;
            let readid_end = tab_positions[1];
            
            // Here we need to allocate a String because we'll store it in the HashSet
            // The to_string() method creates a new String from the &str slice
            let readid = std::str::from_utf8(&buffer[readid_start..readid_end])?.to_string();
            
            // Store the readid in our results set
            save_readids.insert(readid.clone());
            
            // Track which taxid this readid belongs to using a nested HashMap structure
            // This uses Rust's powerful entry API to conditionally insert values
            taxid_readid_map
                .entry(taxid.to_string())
                .or_insert_with(|| {
                    // This closure is called only when a new HashSet needs to be created
                    // It's a lazy initialization pattern in Rust for conditional allocation
                    let mut set = HashSet::new();
                    set.reserve(1000);  // Pre-allocate space for efficiency
                    set
                })
                .insert(readid);  // Insert into the existing or newly created HashSet
        }
    }

    Ok(save_readids)
}