FROM --platform=linux/x86_64 mambaorg/micromamba
USER root
LABEL maintainer="Suffian Azizan"
LABEL version="1.0"
LABEL description="container image of FusionCatcher program v1.33"

# change to root user
USER root
# update Debian OS packages and install additional Linux system utilities, then finally remove cached package lists

# RUN apt-get update && DEBIAN_FRONTEND=noninteractive apt-get install -y build-essential libncurses5-dev default-jdk gawk gcc g++ bzip2 \
# make cmake automake gzip zip unzip zlib1g-dev zlib1g wget curl \
# pigz tar parallel libtbb-dev libtbbmalloc2 \

RUN apt-get update && DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends \
tar wget curl pigz gzip zip unzip gcc g++ bzip2 \
&& rm -rf /var/lib/apt/lists/*

# change user
USER $MAMBA_USER

# Configure Micromamba to use a single thread for package extraction
RUN micromamba config set extract_threads 1

# copy the env file into the container 
COPY --chown=$MAMBA_USER:$MAMBA_USER fusioncatcher/context/base_env.yaml /tmp/base_env.yaml

# Create a new base environment based on the YAML file
RUN micromamba install -y -f /tmp/base_env.yaml && \
micromamba clean --all --yes

# activate the environment during container startup
ARG MAMBA_DOCKERFILE_ACTIVATE=1

# install pip packages
RUN pip install xlrd

# add conda bins to PATH
ENV PATH="/opt/conda/bin:/opt/conda/condabin:$PATH"

USER $MAMBA_USER

# # copy database download script
# COPY --chown=$MAMBA_USER:$MAMBA_USER --chmod=0755 fusioncatcher/src/download-human-db-patched.sh /tmp/dl-human-db-patched.sh
# # change dir
# WORKDIR /tmp
# RUN ./dl-human-db-patched.sh

# change start dir
WORKDIR /home/$MAMBA_USER
