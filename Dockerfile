# First stage: Build the application with maximum optimization for size
FROM rust:alpine AS builder

WORKDIR /app

# Install build dependencies
RUN apk add --no-cache musl-dev libc-dev

# Copy all source code
COPY . .

# Build the optimized application
RUN cargo build --release && \
    strip /app/target/release/krakenclip

# Second stage: Create the minimal runtime image
FROM alpine:latest

WORKDIR /app

# No need for extra runtime libraries as we're using a statically linked binary
# Copy only the binary from the builder stage
COPY --from=builder /app/target/release/krakenclip /app/krakenclip

# Set the binary as the entrypoint
ENTRYPOINT ["/app/krakenclip"]

# Default command (can be overridden)
CMD ["--help"] 