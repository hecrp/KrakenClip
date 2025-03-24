mod krk_parser;
mod taxon_query;
mod logkrk_parser;
mod sequence_processor;
mod generate_test_data;
mod abundance_matrix;
mod cli;
mod biom;

fn main() {
    println!("KrakenClip - High-performance Kraken2 processing toolkit");
    cli::run_cli();
}