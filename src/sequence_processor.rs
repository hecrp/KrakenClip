use std::fs::File;
use std::io::{BufRead, BufReader, Write};
use std::collections::HashSet;
use std::error::Error;

pub fn process_sequence_files(input_files: &[String], save_readids: &HashSet<String>, output_file: &str) -> Result<(), Box<dyn Error>> {
    let mut output = File::create(output_file)?;

    for input_file in input_files {
        let file = File::open(input_file)?;
        let reader = BufReader::new(file);
        let mut lines = reader.lines();

        while let Some(header) = lines.next() {
            let header = header?;
            if header.starts_with('>') || header.starts_with('@') {
                let id = header.split_whitespace().next().unwrap().trim_start_matches(|c| c == '>' || c == '@');
                if save_readids.contains(id) {
                    writeln!(output, "{}", header)?;
                    if header.starts_with('@') {
                        // FASTQ format
                        for _ in 0..3 {
                            if let Some(line) = lines.next() {
                                writeln!(output, "{}", line?)?;
                            }
                        }
                    } else {
                        // FASTA format
                        while let Some(Ok(line)) = lines.next() {
                            if line.starts_with('>') {
                                break;
                            }
                            writeln!(output, "{}", line)?;
                        }
                    }
                }
            }
        }
    }

    Ok(())
}