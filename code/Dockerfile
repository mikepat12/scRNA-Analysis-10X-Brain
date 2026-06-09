# Base image with R 4.5.0
FROM rocker/r-ver:4.5.0

# System dependencies
RUN apt-get update && apt-get install -y \
  libcurl4-openssl-dev \
  libssl-dev \
  libxml2-dev \
  libhdf5-dev \
  libgit2-dev \
  libharfbuzz-dev \
  libfribidi-dev \
  libfreetype6-dev \
  libpng-dev \
  libtiff5-dev \
  libjpeg-dev \
  libglpk-dev \
  libgsl-dev \
  libcairo2-dev \
  libxt-dev \
  libfontconfig1-dev \
  libudunits2-dev \
  pandoc \
  git \
  && rm -rf /var/lib/apt/lists/*

# Install BiocManager and dependency packages
RUN R -e "install.packages('BiocManager', repos='https://cloud.r-project.org')" && \
  R -e "BiocManager::install(version = '3.21', ask = FALSE, update = FALSE)" && \
  R -e "options(repos = BiocManager::repositories())" && \
  R -e "BiocManager::install(c( \
    'batchelor', \
    'BiocGenerics', \
    'DelayedArray', \
    'DelayedMatrixStats', \
    'limma', \
    'lme4', \
    'S4Vectors', \
    'SingleCellExperiment', \
    'SummarizedExperiment', \
    'HDF5Array', \
    'ggrastr' \
  ), ask = FALSE, update = FALSE)"

# Install CRAN packages
RUN R -e "install.packages('devtools')"
RUN R -e "remotes::install_version('tidyverse', '2.0.0')"
RUN R -e "remotes::install_version('Seurat', version='4.3.0.1')"

# Install Monocle3 and BPCells from gitbub
RUN R -e "remotes::install_github('bnprks/BPCells/r')"
RUN R -e "remotes::install_github('cole-trapnell-lab/monocle3')"

# Set default command
CMD ["R"]