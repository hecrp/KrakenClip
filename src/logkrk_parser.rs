use std::collections::{HashSet, HashMap};
use std::fs::File;
use std::io::{BufRead, BufReader};
use std::error::Error;

pub fn parse_kraken_output(kraken_output: &str, save_taxids: &HashSet<String>) -> Result<HashSet<String>, Box<dyn Error>> {
    let file = File::open(kraken_output)?;
    let reader = BufReader::new(file);
    let mut save_readids = HashSet::new();

    for line in reader.lines() {
        let line = line?;
        let fields: Vec<&str> = line.split('\t').collect();
        if fields.len() >= 3 {
            let taxid = fields[2];
            if save_taxids.contains(taxid) {
                save_readids.insert(fields[1].to_string());
            }
        }
    }

    Ok(save_readids)
}