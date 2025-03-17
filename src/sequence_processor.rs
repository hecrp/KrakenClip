use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::collections::HashSet;
use std::error::Error;
use rayon::prelude::*;
use memchr::memchr;
use std::sync::{Arc, Mutex};
use std::path::Path;

// Optimized buffer size constant for efficient I/O operations
const BUFFER_SIZE: usize = 1024 * 1024; // 1MB buffer

/// Process sequence files and extract those matching the specified read IDs
/// 
/// # Arguments
/// * `input_files` - List of input FASTA/FASTQ files to process
/// * `save_readids` - Set of read IDs to extract or exclude
/// * `output_file` - Path to the output file where matching sequences will be written
/// * `exclude` - If true, excludes the IDs in save_readids; if false, includes them
/// 
/// # Returns
/// * `Result<(), Box<dyn Error + Send + Sync>>` - Result of the operation
/// 
/// # Implementation Details
/// This function uses Rust's parallel processing capabilities through Rayon to process
/// multiple files concurrently. It employs a thread-safe shared buffer using Arc<Mutex<>>
/// to allow multiple threads to write to the same output file.
pub fn process_sequence_files(
    input_files: &[String], 
    save_readids: &HashSet<String>, 
    output_file: &str,
    exclude: bool
) -> Result<(), Box<dyn Error + Send + Sync>> {
    // Create a shared writer buffer that can be safely used across threads
    let output = File::create(output_file)?;
    let writer = Arc::new(Mutex::new(BufWriter::with_capacity(BUFFER_SIZE, output)));
    
    // Contador de secuencias extraídas
    let matched_count = Arc::new(Mutex::new(0usize));

    // Process files in parallel using Rayon's parallel iterator
    input_files.par_iter().try_for_each(|input_file| -> Result<(), Box<dyn Error + Send + Sync>> {
        let file = File::open(input_file)?;
        let mut reader = BufReader::with_capacity(BUFFER_SIZE, file);
        let mut buffer = Vec::with_capacity(BUFFER_SIZE);

        // Local buffer to minimize the number of mutex lock acquisitions
        // This improves performance by reducing thread contention
        let mut local_buffer = Vec::with_capacity(BUFFER_SIZE * 2);
        
        let mut read_count = 0;
        let mut local_matched_count = 0;
        
        loop {
            buffer.clear();
            let bytes_read = reader.read_until(b'\n', &mut buffer)?;
            if bytes_read == 0 {
                break;
            }

            if buffer[0] == b'>' || buffer[0] == b'@' {
                read_count += 1;
                let id = parse_id(&buffer);
                let should_write = if exclude {
                    !save_readids.contains(id)
                } else {
                    save_readids.contains(id)
                };

                if should_write {
                    local_matched_count += 1;
                    // Accumulate in the local buffer
                    local_buffer.extend_from_slice(&buffer);
                    
                    if buffer[0] == b'@' {
                        // FASTQ format (4 lines per record)
                        for _ in 0..3 {
                            buffer.clear();
                            reader.read_until(b'\n', &mut buffer)?;
                            local_buffer.extend_from_slice(&buffer);
                        }
                    } else {
                        // FASTA format (variable number of lines)
                        loop {
                            buffer.clear();
                            let bytes_read = reader.read_until(b'\n', &mut buffer)?;
                            if bytes_read == 0 || buffer[0] == b'>' {
                                break;
                            }
                            local_buffer.extend_from_slice(&buffer);
                        }
                    }
                    
                    // Si el buffer local está lleno, escribir a disco
                    if local_buffer.len() > BUFFER_SIZE {
                        let mut writer_guard = writer.lock().unwrap();
                        writer_guard.write_all(&local_buffer)?;
                        // Asegurarse de que el buffer se vacíe completamente
                        writer_guard.flush()?;
                        local_buffer.clear();
                    }
                }
            }
        }
        
        // Write any remaining data in the local buffer
        if !local_buffer.is_empty() {
            let mut writer_guard = writer.lock().unwrap();
            writer_guard.write_all(&local_buffer)?;
            writer_guard.flush()?;
        }
        
        // Actualizar el contador global
        {
            let mut count = matched_count.lock().unwrap();
            *count += local_matched_count;
        }
        
        Ok(())
    })?;

    // Asegurarse de que todos los datos se escriban en disco
    let arc_count = Arc::try_unwrap(matched_count).expect("Error al recuperar el contador").into_inner().unwrap();
    
    // Ensure all data is written to disk by unwrapping the Arc and Mutex
    // This is a Rust-specific pattern for safely reclaiming exclusive ownership
    // of a value that was previously shared between threads
    let mut final_writer = Arc::try_unwrap(writer)
        .expect("Failed to unwrap Arc - still has multiple owners")
        .into_inner()
        .expect("Failed to unwrap Mutex - still locked");
    
    final_writer.flush()?;
    
    Ok(())
}

/// Extract the sequence ID from a FASTA/FASTQ header line
/// 
/// # Arguments
/// * `line` - Header line that starts with '>' or '@'
/// 
/// # Returns
/// * `&str` - Extracted sequence ID
/// 
/// # Performance Note
/// This function uses the highly optimized memchr library for byte-level
/// searching, avoiding unnecessary UTF-8 validation until the final step.
fn parse_id(line: &[u8]) -> &str {
    if line.len() <= 1 {
        return "";
    }
    
    // Find the first space or tab after the initial character
    // Using memchr for optimized byte searching instead of iterating character by character
    let delimiter_pos = memchr(b' ', &line[1..])
        .or_else(|| memchr(b'\t', &line[1..]))
        .map(|pos| pos + 1) // Adjust for the offset from &line[1..]
        .unwrap_or(line.len() - 1);
    
    let id = std::str::from_utf8(&line[1..delimiter_pos])
        .unwrap_or("");
    
    id
}