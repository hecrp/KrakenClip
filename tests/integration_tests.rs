use std::process::Command;
use std::fs;
use std::path::Path;

// Test for Kraken2 report parser
#[test]
fn test_parse_kraken2_report() {
    let output = Command::new("./target/debug/krakenclip")
        .arg("analyze")
        .arg("data/kraken_report.txt")
        .output()
        .expect("Failed to execute analyze command");
    
    let stdout = String::from_utf8(output.stdout).unwrap();
    assert!(stdout.contains("Memory usage"), "The analysis should show information about memory usage");
    assert_eq!(output.status.code().unwrap(), 0, "The command should execute successfully");
}

// Test that verifies the help command
#[test]
fn test_cli_help_output() {
    let output = Command::new("./target/debug/krakenclip")
        .arg("--help")
        .output()
        .expect("Failed to execute command");
    
    let stdout = String::from_utf8(output.stdout).unwrap();
    assert!(stdout.contains("KrakenClip"), "The help should contain the program name");
    assert!(stdout.contains("analyze"), "The help should list the analyze command");
    assert!(stdout.contains("extract"), "The help should list the extract command");
    assert!(stdout.contains("abundance-matrix"), "The help should list the abundance-matrix command");
}

// Test for the analyze command
#[test]
fn test_analyze_command() {
    let output = Command::new("./target/debug/krakenclip")
        .arg("analyze")
        .arg("data/kraken_report.txt")
        .output()
        .expect("Failed to execute analyze command");
    
    assert_eq!(output.status.code().unwrap(), 0, "The command should execute successfully");
    let stdout = String::from_utf8(output.stdout).unwrap();
    assert!(stdout.contains("Hierarchy build time"), "The analysis should show information about hierarchy build time");
}

// Test for the extract command
#[test]
fn test_extract_command() {
    // Make sure the output directory exists
    let output_dir = Path::new("data/outputs/test_extract");
    if !output_dir.exists() {
        fs::create_dir_all(output_dir).expect("Could not create output directory");
    }
    
    let output_file = output_dir.join("extracted.fa");
    // Delete the output file if it exists
    if output_file.exists() {
        fs::remove_file(&output_file).expect("Could not delete output file");
    }
    
    let output = Command::new("./target/debug/krakenclip")
        .arg("extract")
        .arg("data/sample_reads.fastq")
        .arg("data/kraken_log.txt")
        .arg("--output")
        .arg(output_file.to_str().unwrap())
        .arg("--taxids")
        .arg("2")
        .output()
        .expect("Failed to execute extract command");
    
    assert_eq!(output.status.code().unwrap(), 0, "The command should execute successfully");
    assert!(output_file.exists(), "The output file should have been created");
    
    // Verify that the file is not empty
    let metadata = fs::metadata(output_file).expect("Could not read output file");
    assert!(metadata.len() > 0, "The output file should not be empty");
}

// Test for the abundance-matrix command
#[test]
fn test_abundance_matrix_command() {
    // Make sure the output directory exists
    let output_dir = Path::new("data/outputs/test_matrix");
    if !output_dir.exists() {
        fs::create_dir_all(output_dir).expect("Could not create output directory");
    }
    
    let output_file = output_dir.join("abundance.tsv");
    // Delete the output file if it exists
    if output_file.exists() {
        fs::remove_file(&output_file).expect("Could not delete output file");
    }
    
    let output = Command::new("./target/debug/krakenclip")
        .arg("abundance-matrix")
        .arg("data/kraken_report.txt")
        .arg("-o")
        .arg(output_file.to_str().unwrap())
        .output()
        .expect("Failed to execute abundance-matrix command");
    
    assert_eq!(output.status.code().unwrap(), 0, "The command should execute successfully");
    assert!(output_file.exists(), "The output file should have been created");
    
    // Verify that the file is not empty
    let metadata = fs::metadata(output_file).expect("Could not read output file");
    assert!(metadata.len() > 0, "The output file should not be empty");
}

// Test for error handling - non-existent file
#[test]
fn test_error_nonexistent_file() {
    let output = Command::new("./target/debug/krakenclip")
        .arg("analyze")
        .arg("data/nonexistent_file.txt")
        .output()
        .expect("Failed to execute command");
    
    assert_ne!(output.status.code().unwrap(), 0, "The command should fail with a non-existent file");
    let stderr = String::from_utf8(output.stderr).unwrap();
    assert!(stderr.contains("Error"), "It should display an error message");
}

// Test for error handling - incorrect parameters
#[test]
fn test_error_invalid_parameters() {
    let output = Command::new("./target/debug/krakenclip")
        .arg("abundance-matrix")
        .arg("data/kraken_report.txt")
        .arg("--level")
        .arg("Z") // Invalid taxonomic level
        .output()
        .expect("Failed to execute command");
    
    let stderr = String::from_utf8(output.stderr).unwrap();
    assert!(stderr.contains("error:"), "It should display an error message");
}

// Test for abundance matrix normalization
#[test]
fn test_abundance_matrix_normalization() {
    // Make sure the output directory exists
    let output_dir = Path::new("data/outputs/test_matrix");
    if !output_dir.exists() {
        fs::create_dir_all(output_dir).expect("Could not create output directory");
    }
    
    let output_file = output_dir.join("abundance_normalized.tsv");
    // Delete the output file if it exists
    if output_file.exists() {
        fs::remove_file(&output_file).expect("Could not delete output file");
    }
    
    let output = Command::new("./target/debug/krakenclip")
        .arg("abundance-matrix")
        .arg("data/kraken_report.txt")
        .arg("-o")
        .arg(output_file.to_str().unwrap())
        .arg("--normalize")
        .output()
        .expect("Failed to execute abundance-matrix command with normalization");
    
    assert_eq!(output.status.code().unwrap(), 0, "The command should execute successfully");
    assert!(output_file.exists(), "The output file should have been created");
    
    // Verify that the file is not empty
    let metadata = fs::metadata(output_file).expect("Could not read output file");
    assert!(metadata.len() > 0, "The output file should not be empty");
}

// Test for extraction with unclassified sequences inclusion
#[test]
fn test_abundance_matrix_include_unclassified() {
    // Make sure the output directory exists
    let output_dir = Path::new("data/outputs/test_matrix");
    if !output_dir.exists() {
        fs::create_dir_all(output_dir).expect("Could not create output directory");
    }
    
    let output_file = output_dir.join("abundance_unclassified.tsv");
    // Delete the output file if it exists
    if output_file.exists() {
        fs::remove_file(&output_file).expect("Could not delete output file");
    }
    
    let output = Command::new("./target/debug/krakenclip")
        .arg("abundance-matrix")
        .arg("data/kraken_report.txt")
        .arg("-o")
        .arg(output_file.to_str().unwrap())
        .arg("--include-unclassified")
        .output()
        .expect("Failed to execute abundance-matrix command with include-unclassified");
    
    assert_eq!(output.status.code().unwrap(), 0, "The command should execute successfully");
    assert!(output_file.exists(), "The output file should have been created");
    
    // Verify that the file is not empty
    let metadata = fs::metadata(output_file).expect("Could not read output file");
    assert!(metadata.len() > 0, "The output file should not be empty");
}