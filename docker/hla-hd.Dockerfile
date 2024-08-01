FROM --platform=linux/x86_64 mambaorg/micromamba
USER root
LABEL maintainer="Suffian Azizan"
LABEL version="1.0"
LABEL description="container image of HLA-HD program v1.7.0"

# change to root user
USER root
# update Debian OS packages and install additional Linux system utilities, then finally remove cached package lists
RUN apt-get update && DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends \
tar wget curl pigz gzip zip unzip gcc g++ bzip2 procps \
&& rm -rf /var/lib/apt/lists/*

# change user
USER $MAMBA_USER

# Configure Micromamba to use a single thread for package extraction
RUN micromamba config set extract_threads 1

# copy the env file into the container 
COPY --chown=$MAMBA_USER:$MAMBA_USER hla-hd/context/base_env.yaml /tmp/base_env.yaml

# Create a new base environment based on the YAML file
RUN micromamba install -y -f /tmp/base_env.yaml && \
micromamba clean --all --yes

# activate the environment during container startup
ARG MAMBA_DOCKERFILE_ACTIVATE=1

# add conda bins to PATH
ENV PATH="/opt/conda/bin:/opt/conda/condabin:$PATH"

USER root
# copy the HLA-HD source code to the container image
COPY hla-hd/src/hlahd.1.7.0.tar.gz /tmp/hlahd.1.7.0.tar.gz

# unpack tar
RUN tar -zvxf /tmp/hlahd.1.7.0.tar.gz && rm /tmp/hlahd.1.7.0.tar.gz && cd /tmp/hlahd.1.7.0 && sh install.sh

# export to PATH
ENV PATH="$PATH:/tmp/hlahd.1.7.0/bin"

# Set user and group
ARG user=appuser
ARG group=appuser
ARG uid=1000
ARG gid=1000
RUN groupadd -g ${gid} ${group}
RUN useradd -u ${uid} -g ${group} -s /bin/sh -m ${user} 
# the '-m' create a user home directory

# Switch to user
USER ${uid}

