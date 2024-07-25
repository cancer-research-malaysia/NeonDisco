FROM --platform=linux/x86_64 mambaorg/micromamba

LABEL maintainer="Suffian Azizan"
LABEL version="1.0"
LABEL description="container image of Arriba program v2.3.0 for CRM"

# change to root user
USER root
# update Debian OS packages and install additional Linux system utilities with procps; also install R, then finally remove cached package lists
RUN apt-get update && DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends procps curl r-base \
&& rm -rf /var/lib/apt/lists/*

# change user
USER $MAMBA_USER

# Configure Micromamba to use a single thread for package extraction
RUN micromamba config set extract_threads 1

# copy the env file into the container 
COPY --chown=$MAMBA_USER:$MAMBA_USER arriba/context/base_env.yaml /tmp/base_env.yaml

# Create a new base environment based on the YAML file
RUN micromamba install -y -f /tmp/base_env.yaml && \
micromamba clean --all --yes

# activate the environment during container startup
ARG MAMBA_DOCKERFILE_ACTIVATE=1

# add conda bins to PATH
ENV PATH="/opt/conda/bin:/opt/conda/condabin:/opt/conda/var/lib/arriba:$PATH"


# R package installation
USER root

RUN R -e 'install.packages("BiocManager", repos="http://cran.us.r-project.org")'
RUN R -e 'BiocManager::install("argparse")'
RUN R -e 'BiocManager::install(c("GenomicRanges", "GenomicAlignments"))'
RUN R -e 'install.packages("circlize", repos="https://cran.csiro.au/")'

USER $MAMBA_USER
WORKDIR /home

