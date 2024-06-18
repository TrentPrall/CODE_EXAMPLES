# Use an ARM64-friendly base image
FROM ubuntu:latest

# Set environment variables
ENV DEBIAN_FRONTEND=noninteractive

# Update the package list and install necessary dependencies
RUN apt-get update -y && apt-get install -y \
    build-essential \
    wget \
    curl \
    bzip2 \
    libz-dev \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

# Install Miniforge for ARM64
RUN wget --quiet https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-aarch64.sh -O /tmp/miniforge.sh && \
    /bin/bash /tmp/miniforge.sh -b -p /opt/conda && \
    rm /tmp/miniforge.sh && \
    /opt/conda/bin/conda clean -a -y

# Set environment variables for conda
ENV PATH="/opt/conda/bin:$PATH"

# Add bioconda and conda-forge channels
RUN conda config --add channels conda-forge && \
    conda config --add channels bioconda && \
    conda config --set channel_priority strict

# Create conda environment and install hifiasm
RUN conda install -y hifiasm

# Set the working directory
WORKDIR /workdir

# Define the default command to run when the container starts
CMD ["/bin/bash"]
