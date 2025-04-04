use std::collections::HashMap;
use std::fs::File;
use std::io::{BufWriter, Write};
use std::error::Error;
use serde_json::{json, Value};
use chrono;
use crate::krk_parser::{KrakenReport, TaxonEntry};

/// Optimized buffer size for write operations
const BUFFER_SIZE: usize = 256 * 1024; // 256KB

/// Specific errors for the BIOM format module
#[derive(Debug)]
pub enum BiomError {
    IoError(std::io::Error),
    JsonError(serde_json::Error),
    InvalidData(String),
}

impl std::fmt::Display for BiomError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::IoError(e) => write!(f, "I/O error: {}", e),
            Self::JsonError(e) => write!(f, "JSON error: {}", e),
            Self::InvalidData(s) => write!(f, "Invalid data: {}", s),
        }
    }
}

impl Error for BiomError {}

impl From<std::io::Error> for BiomError {
    fn from(e: std::io::Error) -> Self {
        Self::IoError(e)
    }
}

impl From<serde_json::Error> for BiomError {
    fn from(e: serde_json::Error) -> Self {
        Self::JsonError(e)
    }
}

/// Specialized result type for BIOM functions
pub type BiomResult<T> = Result<T, BiomError>;

/// Structure to efficiently store BIOM format data
#[derive(Default)]
pub struct BiomTable {
    /// Matrix of abundances (rows are taxa, columns are samples)
    pub data: Vec<Vec<f64>>,
    /// Row metadata (taxonomic information)
    pub rows: Vec<HashMap<String, String>>,
    /// Column metadata (sample information)
    pub columns: Vec<HashMap<String, String>>,
    /// Row IDs (taxon names)
    pub row_ids: Vec<String>,
    /// Column IDs (sample names)
    pub column_ids: Vec<String>,
    /// Matrix type (sparse or dense)
    pub matrix_type: String,
    /// Matrix element type
    pub element_type: String,
    /// Shape of the matrix
    pub shape: (usize, usize),
    /// Date of creation
    pub date: String,
    /// Generated by information
    pub generated_by: String,
    /// Table ID
    pub id: String,
    /// Table type
    pub table_type: String,
    /// Format URL
    pub format_url: String,
}

impl BiomTable {
    /// Creates a new BIOM table from a Kraken report
    /// 
    /// # Arguments
    /// * `report` - Parsed Kraken2 report
    /// * `sample_name` - Name of the sample
    /// * `normalize` - If true, normalizes abundances to percentages
    pub fn from_kraken_report(report: &KrakenReport, sample_name: &str, normalize: bool) -> Self {
        let mut table = Self::default();
        
        // Set basic BIOM metadata
        table.date = chrono::Local::now().to_rfc3339();
        table.generated_by = "KrakenClip".to_string();
        table.id = format!("krakenclip_{}", sample_name);
        table.table_type = "OTU table".to_string();
        table.format_url = "http://biom-format.org/documentation/format_versions/biom-1.0.html".to_string();
        table.matrix_type = "dense".to_string();
        table.element_type = "float".to_string();
        
        // Process the report
        let mut row_data = Vec::new();
        let mut row_metadata = Vec::new();
        let mut row_ids = Vec::new();
        
        // Process root and its children
        Self::process_node(&report.root, &mut row_data, &mut row_metadata, &mut row_ids, normalize);
        
        // Process unclassified if present
        if let Some(ref unclassified) = report.unclassified {
            Self::process_node(unclassified, &mut row_data, &mut row_metadata, &mut row_ids, normalize);
        }
        
        // Set the data
        table.data = vec![row_data];
        table.rows = row_metadata;
        table.row_ids = row_ids.clone();
        table.column_ids = vec![sample_name.to_string()];
        table.shape = (row_ids.len(), 1);
        
        table
    }
    
    /// Processes a taxonomic node and its children recursively
    fn process_node(
        node: &TaxonEntry,
        row_data: &mut Vec<f64>,
        row_metadata: &mut Vec<HashMap<String, String>>,
        row_ids: &mut Vec<String>,
        normalize: bool
    ) {
        // Add current node
        let abundance = if normalize {
            node.percentage as f64
        } else {
            node.clade_reads as f64
        };
        
        row_data.push(abundance);
        row_ids.push(node.name.clone());
        
        // Add metadata
        let mut metadata = HashMap::new();
        metadata.insert("taxid".to_string(), node.taxid.to_string());
        metadata.insert("rank".to_string(), node.rank.clone());
        metadata.insert("level".to_string(), node.level.to_string());
        row_metadata.push(metadata);
        
        // Process children
        for child in &node.children {
            Self::process_node(child, row_data, row_metadata, row_ids, normalize);
        }
    }
    
    /// Writes the BIOM table to a file in JSON format
    /// 
    /// # Arguments
    /// * `output_file` - Path to the output file
    /// 
    /// # Returns
    /// * `BiomResult<()>` - Result of the operation
    pub fn write_json(&self, output_file: &str) -> BiomResult<()> {
        let file = File::create(output_file)?;
        let mut writer = BufWriter::with_capacity(BUFFER_SIZE, file);
        
        // Create the BIOM JSON structure
        let biom_json = json!({
            "id": self.id,
            "format": "Biological Observation Matrix 1.0.0",
            "format_url": self.format_url,
            "type": self.table_type,
            "generated_by": self.generated_by,
            "date": self.date,
            "matrix_type": self.matrix_type,
            "matrix_element_type": self.element_type,
            "shape": self.shape,
            "data": self.data,
            "rows": self.rows.iter().map(|metadata| {
                json!({
                    "id": metadata.get("taxid").unwrap(),
                    "metadata": metadata
                })
            }).collect::<Vec<Value>>(),
            "columns": self.column_ids.iter().map(|id| {
                json!({
                    "id": id,
                    "metadata": {}
                })
            }).collect::<Vec<Value>>()
        });
        
        // Write to file
        serde_json::to_writer_pretty(&mut writer, &biom_json)?;
        writer.flush()?;
        
        Ok(())
    }
} 