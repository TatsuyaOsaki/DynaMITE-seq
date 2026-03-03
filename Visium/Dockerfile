#FROM rocker/r-ver:latest
#FROM rocker/tidyverse:latest
# Use rocker/geospatial which includes sf, terra, and spatial dependencies
FROM rocker/geospatial:latest

LABEL maintainer="iddryg@gmail.com"

# Install system dependencies for R packages
RUN apt-get update && apt-get install -y --no-install-recommends \
    libmagick++-dev \
    libxml2-dev \
    libssl-dev \
    libcurl4-openssl-dev \
    libhdf5-dev \
    libfftw3-dev \
    libgsl-dev \
    libpng-dev \
    libjpeg-dev \
    libtiff5-dev \
    libgit2-dev \
    libglpk-dev \
    libudunits2-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libfontconfig1-dev \
    libfreetype6-dev \
    pkg-config \
    pandoc \
    python3-pip \
    python3-venv \
    && rm -rf /var/lib/apt/lists/*

# Configure writable R library path
ENV R_LIBS_USER=/usr/local/lib/R/site-library
RUN mkdir -p ${R_LIBS_USER} && chmod 777 ${R_LIBS_USER}
RUN mkdir -p /usr/local/lib/R/library && chmod 777 /usr/local/lib/R/library
RUN echo "options(repos = c(CRAN = 'https://cloud.r-project.org'), download.file.method = 'libcurl')" >> /usr/local/lib/R/etc/Rprofile.site && \
    echo ".libPaths(c('${R_LIBS_USER}', .libPaths()))" >> /usr/local/lib/R/etc/Rprofile.site

# Install Jupyter Lab with PEP 668 workaround
RUN pip3 install jupyterlab --break-system-packages

# Install R kernel
RUN Rscript -e "install.packages('IRkernel', repos='https://cloud.r-project.org')" && \
    Rscript -e "IRkernel::installspec(user = FALSE, sys_prefix=TRUE)"

# Install concaveman (sf is already installed in rocker/geospatial)
# Install this separately first to catch any errors
RUN Rscript -e "install.packages('concaveman', repos='https://cloud.r-project.org', verbose=TRUE)" || exit 1

# Install CRAN packages with error checking
# Tidyverse, including ggplot2, stringr, etc are already installed with FROM rocker/tidyverse:latest
# RUN Rscript -e "install.packages('tidyverse', repos='https://tidyverse.r-universe.dev')"
RUN Rscript -e "install.packages(c( \
    'Seurat', 'Matrix', 'doParallel', 'ggplot2', \
    'msigdbr', 'pheatmap', 'readxl', 'corrplot', 'ggbeeswarm', \
    'openxlsx', 'RColorBrewer', 'rlang', 'GGally', 'data.table', \
    'scales', 'ggrepel', 'devtools', 'patchwork', \
    'shiny', 'ashr', 'magick', 'paletteer', 'scCustomize', 'umap', \
    'hdf5r', 'spatstat', 'arrow', 'alphahull', 'jsonlite', \
    'future', 'furrr'), repos='https://cloud.r-project.org', verbose=TRUE)"


# Bioconductor is installed with FROM rocker/tidyverse:latest
# Install Bioconductor
#RUN Rscript -e "install.packages('BiocManager', repos='https://cloud.r-project.org')" && \
#    Rscript -e "BiocManager::install(version = '3.21')"

# Install Bioconductor packages
RUN Rscript -e "BiocManager::install(c( \
    'SpatialExperiment', 'GSVA', 'clusterProfiler', 'org.Hs.eg.db', 'dittoSeq', \
    'BaseSet', 'fgsea', 'enrichplot', 'iSEE', 'limma', \
    'tximport', 'DESeq2', 'PCAtools', 'EnhancedVolcano', 'standR', \
    'escape', 'UCell' \
    ), ask=FALSE, update=TRUE, checkBuilt=TRUE)"

# Create necessary directories and set permissions
RUN mkdir /.local && chmod 777 /.local && mkdir /.jupyter && chmod 777 /.jupyter \
    && chmod 777 ${R_LIBS_USER} && chmod 777 /usr/local/lib/R/library \
    && mkdir -p /.cache/R && chmod -R 777 /.cache \
    && mkdir -p /tmp/.cache && chmod -R 777 /tmp/.cache \
    && mkdir -p /home/.cache && chmod -R 777 /home/.cache

# Set environment variables for R cache
ENV R_USER_CACHE_DIR=/tmp/.cache
ENV XDG_CACHE_HOME=/tmp/.cache

# Install msigdbf to get full dataset
RUN Rscript -e "install.packages('msigdbdf', repos = 'https://igordot.r-universe.dev')"

WORKDIR /home

CMD ["jupyter", "lab", "--ip=0.0.0.0", "--port=8888", "--allow-root"]
