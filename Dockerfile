# First stage: Build the application
FROM rust:slim as builder

WORKDIR /app

# Copy the Cargo files for dependency caching
COPY Cargo.toml Cargo.lock ./

# Create a dummy main.rs to build dependencies
RUN mkdir -p src && \
    echo "fn main() {}" > src/main.rs && \
    cargo build --release && \
    rm -rf src

# Copy the actual source code
COPY src/ src/

# Build the application
RUN cargo build --release

# Second stage: Create the minimal runtime image
FROM debian:slim

WORKDIR /app

# Install necessary runtime libraries (add any runtime dependencies here)
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
    ca-certificates \
    && rm -rf /var/lib/apt/lists/*

# Copy the binary from the builder stage
COPY --from=builder /app/target/release/kraken2-parser /app/kraken2-parser

# Set the binary as the entrypoint
ENTRYPOINT ["/app/kraken2-parser"]

# Default command (can be overridden)
CMD ["--help"] 