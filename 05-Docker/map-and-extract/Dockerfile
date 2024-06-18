# Use the latest official Ubuntu image as a base
FROM ubuntu:latest

# Install dependencies
RUN apt-get update && apt-get install -y \
    minimap2 \
    samtools \
    procps \
    wget \
    openjdk-11-jre-headless \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Download and install BBTools
RUN wget -O bbmap.tar.gz https://sourceforge.net/projects/bbmap/files/BBMap_38.90.tar.gz/download && \
    mkdir -p /opt/bbmap && \
    tar -xzvf bbmap.tar.gz -C /opt/bbmap --strip-components=1 && \
    ln -s /opt/bbmap/*.sh /usr/local/bin/ && \
    rm bbmap.tar.gz

# Set up environment variables
ENV PATH="/usr/local/bin:/opt/bbmap:${PATH}"

# Default command
CMD ["bash"]
