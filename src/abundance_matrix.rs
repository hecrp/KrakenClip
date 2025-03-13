use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{BufWriter, Write};
use std::error::Error;
use crate::krk_parser::{KrakenReport, TaxonEntry};

// Constants for optimization
const BUFFER_SIZE: usize = 256 * 1024; // 256KB buffer for writing
const TAXON_LEVELS: &[(&str, &str)] = &[
    ("K", "kingdom"),
    ("P", "phylum"),
    ("C", "class"),
    ("O", "order"),
    ("F", "family"),
    ("G", "genus"),
    ("S", "species")
];

// Special name for unclassified entries in the matrix
const UNCLASSIFIED_NAME: &str = "Unclassified";

/// Structure to efficiently store the abundance matrix
#[derive(Default)]
pub struct AbundanceMatrix {
    // Map of taxa to their abundances by sample
    taxon_abundances: HashMap<String, HashMap<String, f64>>,
    // Set of all samples
    samples: HashSet<String>,
    // Current taxonomic level
    level: String,
    // Total reads by sample (for normalization)
    sample_totals: HashMap<String, f64>,
    // Specifically stores unclassified reads by sample
    unclassified_reads: HashMap<String, f64>,
    // Indicates whether to force inclusion of unclassified entries
    force_include_unclassified: bool,
}

impl AbundanceMatrix {
    /// Creates a new abundance matrix
    pub fn new(level: &str) -> Self {
        Self {
            taxon_abundances: HashMap::new(),
            samples: HashSet::new(),
            level: level.to_string(),
            sample_totals: HashMap::new(),
            unclassified_reads: HashMap::new(),
            force_include_unclassified: false,
        }
    }

    /// Sets whether unclassified reads should be forcibly included
    pub fn set_force_include_unclassified(&mut self, include: bool) {
        self.force_include_unclassified = include;
    }

    /// Adds a sample to the matrix
    pub fn add_sample(&mut self, report: &KrakenReport, sample_name: &str, min_abundance: f64, normalize: bool) {
        self.samples.insert(sample_name.to_string());
        
        // Calculate the total reads for this sample
        let total_reads = self.calculate_total_reads(report);
        self.sample_totals.insert(sample_name.to_string(), total_reads);
        
        // Store unclassified reads
        if let Some(ref unclassified) = report.unclassified {
            let unclassified_count = unclassified.clade_reads as f64;
            self.unclassified_reads.insert(sample_name.to_string(), unclassified_count);
            
            // If we're forcing the inclusion of unclassified reads, add them as a special taxon
            if self.force_include_unclassified {
                let abundance = if normalize {
                    (unclassified_count / total_reads) * 100.0
                } else {
                    unclassified_count
                };
                
                if abundance >= min_abundance {
                    let taxon_map = self.taxon_abundances
                        .entry(UNCLASSIFIED_NAME.to_string())
                        .or_default();
                    taxon_map.insert(sample_name.to_string(), abundance);
                }
            }
        }
        
        // Process the taxonomic tree recursively
        self.process_node(&report.root, sample_name, min_abundance, normalize);
    }

    /// Calculates the total reads in a sample
    fn calculate_total_reads(&self, report: &KrakenReport) -> f64 {
        let mut total = report.root.clade_reads as f64;
        if let Some(ref unclassified) = report.unclassified {
            total += unclassified.clade_reads as f64;
        }
        total
    }

    /// Processes a node in the taxonomic tree
    fn process_node(&mut self, node: &TaxonEntry, sample_name: &str, min_abundance: f64, normalize: bool) {
        // Check if the node is at the desired taxonomic level
        if node.rank == self.level {
            let abundance = if normalize {
                // Normalize by the total reads of the sample
                let total = self.sample_totals.get(sample_name).unwrap_or(&1.0);
                (node.clade_reads as f64 / total) * 100.0
            } else {
                node.clade_reads as f64
            };

            if abundance >= min_abundance {
                let taxon_map = self.taxon_abundances
                    .entry(node.name.clone())
                    .or_default();
                taxon_map.insert(sample_name.to_string(), abundance);
            }
        }

        // Process children recursively
        for child in &node.children {
            self.process_node(child, sample_name, min_abundance, normalize);
        }
    }

    /// Applies a transformation to convert absolute counts to proportions
    pub fn transform_to_proportions(&mut self) {
        // Only apply if we're not already normalized
        for (_, sample_abundances) in self.taxon_abundances.iter_mut() {
            for (sample, abundance) in sample_abundances.iter_mut() {
                let total = self.sample_totals.get(sample).unwrap_or(&1.0);
                *abundance = (*abundance / total) * 100.0;
            }
        }
    }

    /// Generates the abundance matrix in TSV format
    pub fn write_matrix(&self, output_file: &str) -> Result<(), Box<dyn Error>> {
        let file = File::create(output_file)?;
        let mut writer = BufWriter::with_capacity(BUFFER_SIZE, file);

        // Get sorted list of samples
        let mut samples: Vec<String> = self.samples.iter().cloned().collect();
        samples.sort();

        // Write header
        write!(writer, "Taxon")?;
        for sample in &samples {
            write!(writer, "\t{}", sample)?;
        }
        writeln!(writer)?;

        // Write data
        let mut taxons: Vec<String> = self.taxon_abundances.keys().cloned().collect();
        
        // Sort taxa but place "Unclassified" at the beginning if it exists
        taxons.sort_by(|a, b| {
            if a == UNCLASSIFIED_NAME {
                std::cmp::Ordering::Less
            } else if b == UNCLASSIFIED_NAME {
                std::cmp::Ordering::Greater
            } else {
                a.cmp(b)
            }
        });

        for taxon in taxons {
            write!(writer, "{}", taxon)?;
            let abundances = self.taxon_abundances.get(&taxon).unwrap();
            
            for sample in &samples {
                let abundance = abundances.get(sample).unwrap_or(&0.0);
                write!(writer, "\t{:.6}", abundance)?;
            }
            writeln!(writer)?;
        }

        writer.flush()?;
        Ok(())
    }
}

/// Validates the specified taxonomic level
pub fn validate_taxonomic_level(level: &str) -> bool {
    TAXON_LEVELS.iter().any(|(code, _)| *code == level)
}

/// Gets the full name of the taxonomic level
pub fn get_taxonomic_level_name(level: &str) -> Option<&str> {
    TAXON_LEVELS.iter()
        .find(|(code, _)| *code == level)
        .map(|(_, name)| *name)
} 