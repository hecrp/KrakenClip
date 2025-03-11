# First stage: Build the application
FROM rust:latest as builder

WORKDIR /app

# Copy all source code first
COPY . .

# Build the application
RUN cargo build --release

# Second stage: Create the minimal runtime image
FROM debian:bookworm-slim

WORKDIR /app

# Install necessary runtime libraries
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