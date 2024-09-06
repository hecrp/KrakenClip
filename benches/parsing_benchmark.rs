use criterion::{black_box, criterion_group, criterion_main, Criterion};
use kraken2_parser::krk_parser::parse_kraken2_report;

fn benchmark_parsing(c: &mut Criterion) {
    c.bench_function("parse kraken2 report", |b| {
        b.iter(|| parse_kraken2_report(black_box("benches/sample_report.txt")))
    });
}

criterion_group!(benches, benchmark_parsing);
criterion_main!(benches);