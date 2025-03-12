use kraken2_parser::krk_parser::parse_kraken2_report;
use std::process::Command;

#[test]
fn test_parse_kraken2_report() {
    let entries = parse_kraken2_report("data/example_report.txt");
    assert!(!entries.0.root.is_empty());
    // Add more assertions here
}

#[test]
fn test_cli_help_output() {
    // Verificar que la ayuda contiene el nombre correcto
    let output = Command::new("./target/debug/krakenclip")
        .arg("--help")
        .output()
        .expect("Failed to execute command");
    
    let stdout = String::from_utf8(output.stdout).unwrap();
    assert!(stdout.contains("KrakenClip"));
}