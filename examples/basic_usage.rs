use kraken2_parser::krk_parser::parse_kraken2_report;

fn main() {
    let entries = parse_kraken2_report("examples/sample_report.txt");
    println!("Parsed {} entries", entries.0.root.len());
    // Add more example usage here
}