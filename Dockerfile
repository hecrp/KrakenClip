# First stage: Build the application
FROM rust:alpine as builder

WORKDIR /app

# Install build dependencies
RUN apk add --no-cache musl-dev libc-dev

# Copy all source code
COPY . .

# Build the application with optimization for size
RUN cargo build --release

# Second stage: Create the minimal runtime image
FROM alpine:latest

WORKDIR /app

# Install necessary runtime libraries (minimal)
RUN apk add --no-cache ca-certificates

# Copy the binary from the builder stage
COPY --from=builder /app/target/release/kraken2-parser /app/kraken2-parser

# Set the binary as the entrypoint
ENTRYPOINT ["/app/kraken2-parser"]

# Default command (can be overridden)
CMD ["--help"] 