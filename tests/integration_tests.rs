use kraken2_parser::krk_parser::parse_kraken2_report;

#[test]
fn test_parse_kraken2_report() {
    let entries = parse_kraken2_report("data/example_report.txt");
    assert!(!entries.0.root.is_empty());
    // Add more assertions here
}