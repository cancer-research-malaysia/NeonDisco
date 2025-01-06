FROM mambaorg/micromamba:git-911a014-bookworm-slim
USER root
LABEL maintainer="Suffian Azizan"
LABEL version="1.5"
LABEL description="container image of tools for NGS reads preprocessing (SAMTools, Picard, Yara)"

# change to root user
USER root
# update Debian OS packages and install additional Linux system utilities, then finally remove cached package lists
RUN apt-get update && DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends \
build-essential tar wget curl pigz gzip zip unzip gcc g++ bzip2 procps git cmake locales \
&& rm -rf /var/lib/apt/lists/* \
&& echo "en_US.UTF-8 UTF-8" > /etc/locale.gen && locale-gen

# change user
USER $MAMBA_USER

# Configure Micromamba to use a single thread for package extraction
RUN micromamba config set extract_threads 1

# copy the env file into the container 
COPY --chown=$MAMBA_USER:$MAMBA_USER ngs-preproc/context/base_env.yaml /tmp/base_env.yaml

# Create a new base environment based on the YAML file
RUN micromamba install -y -f /tmp/base_env.yaml && \
micromamba clean --all --yes

# activate the environment during container startup
ARG MAMBA_DOCKERFILE_ACTIVATE=1

# add conda bins to PATH
ENV PATH="/opt/conda/bin:/opt/conda/condabin:$PATH"

# change user to root
USER root

# install yara
RUN git clone https://github.com/seqan/seqan.git
RUN mkdir yara-build && cd yara-build && cmake ../seqan -DSEQAN_BUILD_SYSTEM=APP:yara -DCMAKE_CXX_COMPILER=/usr/bin/g++-12 && make all && cp bin/yara* /usr/local/bin

# download hla reference files from OptiType repo for yara
RUN curl -o /tmp/hla_ref_for_yara_DNA.fa ftp://ftp.ebi.ac.uk/pub/databases/ipd/imgt/hla/hla_gen.fasta && curl -o /tmp/hla_ref_for_yara_RNA.fa ftp://ftp.ebi.ac.uk/pub/databases/ipd/imgt/hla/hla_nuc.fasta

# Docker suffers from absolutely atrocious way of consolidating the paradigm of restricting privileges when running containers (rootless mode) with writing outputs to bound host volumes without using Docker volumes or other convoluted workarounds.

# Fortunately there is this tool that removes this altogether and helps matches the UID and GID of whoever is running the container image on a host machine

# Install MatchHostFsOwner. Using version 1.0.1
# See https://github.com/FooBarWidget/matchhostfsowner/releases
ADD https://github.com/FooBarWidget/matchhostfsowner/releases/download/v1.0.1/matchhostfsowner-1.0.1-x86_64-linux.gz /sbin/matchhostfsowner.gz
RUN gunzip /sbin/matchhostfsowner.gz && \
  chown root: /sbin/matchhostfsowner && \
  chmod +x /sbin/matchhostfsowner
RUN addgroup --gid 9999 app && \
  adduser --uid 9999 --gid 9999 --disabled-password --gecos App app

# set workdir
WORKDIR /home/app

# build yara index
RUN mkdir -p /home/app/refs/HLA-yara_index/dnaseq/ /home/app/refs/HLA-yara_index/rnaseq/ 
RUN yara_indexer -o /home/app/refs/HLA-yara_index/dnaseq/hla_ref_for_yara_DNA /tmp/hla_ref_for_yara_DNA.fa && yara_indexer -o /home/app/refs/HLA-yara_index/rnaseq/hla_ref_for_yara_RNA /tmp/hla_ref_for_yara_RNA.fa

# transfer download script to download reference data for STAR and Arriba
COPY ngs-preproc/src/download_refs_for_star.sh /tmp/download_refs_for_star.sh

ENTRYPOINT ["/usr/local/bin/_entrypoint.sh", "/sbin/matchhostfsowner"]
