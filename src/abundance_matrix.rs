use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{BufWriter, Write};
use std::error::Error;
use std::path::Path;
use std::cmp::Ordering;
use crate::krk_parser::{KrakenReport, TaxonEntry};

/// Optimized buffer size for write operations
const BUFFER_SIZE: usize = 256 * 1024; // 256KB

/// Supported taxonomic levels
const TAXON_LEVELS: &[(&str, &str)] = &[
    ("K", "kingdom"),
    ("P", "phylum"),
    ("C", "class"),
    ("O", "order"),
    ("F", "family"),
    ("G", "genus"),
    ("S", "species")
];

/// Special name for unclassified entries in the matrix
const UNCLASSIFIED_NAME: &str = "Unclassified";

/// Specific errors for the abundance matrix module
#[derive(Debug)]
pub enum AbundanceMatrixError {
    IoError(std::io::Error),
    InvalidLevel(String),
    InvalidValue(String),
}

impl std::fmt::Display for AbundanceMatrixError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::IoError(e) => write!(f, "I/O error: {}", e),
            Self::InvalidLevel(s) => write!(f, "Invalid taxonomic level: {}", s),
            Self::InvalidValue(s) => write!(f, "Invalid value: {}", s),
        }
    }
}

impl Error for AbundanceMatrixError {}

impl From<std::io::Error> for AbundanceMatrixError {
    fn from(e: std::io::Error) -> Self {
        Self::IoError(e)
    }
}

/// Specialized result type for abundance matrix functions
pub type AbundanceResult<T> = Result<T, AbundanceMatrixError>;

/// Structure to efficiently store the abundance matrix
/// 
/// This implementation uses nested HashMaps to efficiently store
/// sparse abundance data, where many taxa may be present in some
/// samples but not others.
#[derive(Default)]
pub struct AbundanceMatrix {
    /// Map of taxa to their abundances by sample
    taxon_abundances: HashMap<String, HashMap<String, f64>>,
    /// Set of all samples
    samples: HashSet<String>,
    /// Current taxonomic level
    level: String,
    /// Total reads per sample (for normalization)
    sample_totals: HashMap<String, f64>,
    /// Specifically stores unclassified reads per sample
    unclassified_reads: HashMap<String, f64>,
    /// Indicates whether to force inclusion of unclassified entries
    force_include_unclassified: bool,
}

impl AbundanceMatrix {
    /// Creates a new abundance matrix for the specified taxonomic level
    /// 
    /// # Arguments
    /// * `level` - Taxonomic level (K, P, C, O, F, G, S)
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
    /// 
    /// # Arguments
    /// * `include` - If true, includes unclassified entries
    pub fn set_force_include_unclassified(&mut self, include: bool) {
        self.force_include_unclassified = include;
    }

    /// Adds a sample to the abundance matrix
    /// 
    /// # Arguments
    /// * `report` - Parsed Kraken2 report
    /// * `sample_name` - Name of the sample
    /// * `min_abundance` - Minimum abundance threshold to include a taxon
    /// * `normalize` - If true, normalizes the values during processing
    /// 
    /// # Implementation Details
    /// This method calculates the total reads for normalization and
    /// recursively processes the taxonomic tree to extract abundances
    /// at the specified taxonomic level.
    pub fn add_sample(&mut self, report: &KrakenReport, sample_name: &str, min_abundance: f64, normalize: bool) {
        self.samples.insert(sample_name.to_string());
        
        // Calculate the total reads for this sample
        let total_reads = self.calculate_total_reads(report);
        self.sample_totals.insert(sample_name.to_string(), total_reads);
        
        // Store unclassified reads
        if let Some(ref unclassified) = report.unclassified {
            let unclassified_count = unclassified.clade_reads as f64;
            self.unclassified_reads.insert(sample_name.to_string(), unclassified_count);
            
            // If we're forcing the inclusion of no classified, add them as a special taxon
            if self.force_include_unclassified {
                let abundance = if normalize {
                    (unclassified_count / total_reads) * 100.0
                } else {
                    unclassified_count
                };
                
                if abundance >= min_abundance {
                    self.taxon_abundances
                        .entry(UNCLASSIFIED_NAME.to_string())
                        .or_default()
                        .insert(sample_name.to_string(), abundance);
                }
            }
        }
        
        // Process the taxonomic tree recursively
        self.process_node(&report.root, sample_name, min_abundance, normalize);
    }

    /// Calculates the total reads in a sample
    /// 
    /// # Arguments
    /// * `report` - Parsed Kraken2 report
    /// 
    /// # Returns
    /// * `f64` - Total reads in the sample
    fn calculate_total_reads(&self, report: &KrakenReport) -> f64 {
        let mut total = report.root.clade_reads as f64;
        if let Some(ref unclassified) = report.unclassified {
            total += unclassified.clade_reads as f64;
        }
        total
    }

    /// Processes a node in the taxonomic tree
    /// 
    /// # Arguments
    /// * `node` - Current node to process
    /// * `sample_name` - Name of the sample
    /// * `min_abundance` - Minimum abundance threshold
    /// * `normalize` - If true, normalizes values during processing
    /// 
    /// # Implementation Details
    /// This recursive method traverses the taxonomic tree, extracting
    /// abundance data for nodes at the target taxonomic level.
    /// The recursive approach ensures we capture all taxa at the specified
    /// level, regardless of their position in the tree.
    fn process_node(&mut self, node: &TaxonEntry, sample_name: &str, min_abundance: f64, normalize: bool) {
        // Check if the node is at the desired taxonomic level
        if node.rank == self.level {
            let abundance = if normalize {
                // Normalize by the total reads in the sample
                let total = self.sample_totals.get(sample_name).unwrap_or(&1.0);
                (node.clade_reads as f64 / total) * 100.0
            } else {
                node.clade_reads as f64
            };

            if abundance >= min_abundance {
                self.taxon_abundances
                    .entry(node.name.clone())
                    .or_default()
                    .insert(sample_name.to_string(), abundance);
            }
        }

        // Process children recursively
        for child in &node.children {
            self.process_node(child, sample_name, min_abundance, normalize);
        }
    }

    /// Applies a transformation to convert absolute counts to proportions
    /// 
    /// # Implementation Details
    /// This method is used when the matrix was initially populated with
    /// absolute read counts but needs to be converted to proportions.
    /// It divides each taxon's count by the total reads in its sample
    /// and multiplies by 100 to get a percentage.
    pub fn transform_to_proportions(&mut self) {
        // Only apply if we're not already normalized
        for (_, sample_abundances) in self.taxon_abundances.iter_mut() {
            for (sample, abundance) in sample_abundances.iter_mut() {
                if let Some(total) = self.sample_totals.get(sample) {
                    if *total > 0.0 {
                        *abundance = (*abundance / total) * 100.0;
                    }
                }
            }
        }
    }

    /// Generates the abundance matrix in TSV format
    /// 
    /// # Arguments
    /// * `output_file` - Path to the output file
    /// 
    /// # Returns
    /// * `AbundanceResult<()>` - Result of the operation
    /// 
    /// # Implementation Details
    /// This method writes the abundance matrix to a TSV file with
    /// taxa as rows and samples as columns. It uses a BufWriter for
    /// efficient I/O operations and ensures "Unclassified" appears
    /// at the top of the matrix if present.
    pub fn write_matrix(&self, output_file: &str) -> AbundanceResult<()> {
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

        // Prepare the taxa
        let mut taxons: Vec<String> = self.taxon_abundances.keys().cloned().collect();
        
        // Sort taxa but place "Unclassified" at the beginning if it exists
        taxons.sort_by(|a, b| {
            if a == UNCLASSIFIED_NAME {
                Ordering::Less
            } else if b == UNCLASSIFIED_NAME {
                Ordering::Greater
            } else {
                a.cmp(b)
            }
        });

        // Write data
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
/// 
/// # Arguments
/// * `level` - Taxonomic level to validate (K, P, C, O, F, G, S)
/// 
/// # Returns
/// * `bool` - True if the level is valid, false otherwise
pub fn validate_taxonomic_level(level: &str) -> bool {
    TAXON_LEVELS.iter().any(|(code, _)| *code == level)
}

/// Gets the full name of the taxonomic level
/// 
/// # Arguments
/// * `level` - Taxonomic level code (K, P, C, O, F, G, S)
/// 
/// # Returns
/// * `Option<&str>` - Full name of the taxonomic level if valid
pub fn get_taxonomic_level_name(level: &str) -> Option<&str> {
    TAXON_LEVELS.iter()
        .find(|(code, _)| *code == level)
        .map(|(_, name)| *name)
} 