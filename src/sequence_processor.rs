use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Write};
use std::collections::HashSet;
use std::error::Error;
use rayon::prelude::*;
use memchr::memchr;
use std::sync::{Arc, Mutex};

pub fn process_sequence_files(input_files: &[String], save_readids: &HashSet<String>, output_file: &str) -> Result<(), Box<dyn Error + Send + Sync>> {
    let output = File::create(output_file)?;
    let writer: Arc<Mutex<BufWriter<File>>> = Arc::new(Mutex::new(BufWriter::new(output)));

    input_files.par_iter().try_for_each(|input_file| -> Result<(), Box<dyn Error + Send + Sync>> {
        let writer: Arc<Mutex<BufWriter<File>>> = Arc::clone(&writer);
        let file = File::open(input_file)?;
        let mut reader: BufReader<File> = BufReader::with_capacity(1024 * 1024, file); // 1MB buffer
        let mut buffer: Vec<u8> = Vec::with_capacity(1024 * 1024); // 1MB buffer for sequences

        loop {
            buffer.clear();
            let bytes_read = reader.read_until(b'\n', &mut buffer)?;
            if bytes_read == 0 {
                break;
            }

            if buffer[0] == b'>' || buffer[0] == b'@' {
                let id = parse_id(&buffer);
                if save_readids.contains(id) {
                    writer.lock().unwrap().write_all(&buffer)?;
                    
                    if buffer[0] == b'@' {
                        // FASTQ format
                        for _ in 0..3 {
                            buffer.clear();
                            reader.read_until(b'\n', &mut buffer)?;
                            writer.lock().unwrap().write_all(&buffer)?;
                        }
                    } else {
                        // FASTA format
                        loop {
                            buffer.clear();
                            let bytes_read = reader.read_until(b'\n', &mut buffer)?;
                            if bytes_read == 0 || buffer[0] == b'>' {
                                break;
                            }
                            writer.lock().unwrap().write_all(&buffer)?;
                        }
                    }
                }
            }
        }
        Ok(())
    })?;

    Ok(())
}

fn parse_id(line: &[u8]) -> &str {
    let end = memchr(b' ', &line[1..]).unwrap_or(line.len() - 1);
    std::str::from_utf8(&line[1..=end]).unwrap_or("")
}