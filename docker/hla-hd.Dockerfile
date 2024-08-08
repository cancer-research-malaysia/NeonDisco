FROM --platform=linux/x86_64 mambaorg/micromamba:git-911a014-bookworm-slim
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


# change user to root
USER root
# copy the HLA-HD source code to the container image
COPY hla-hd/src/hlahd.1.7.0.tar.gz /tmp/hlahd.1.7.0.tar.gz
# also copy the specific dat file of HLA dictionary for this version and a custom update.dict.sh
COPY hla-hd/src/IMGTHLA-3.50.0/ /tmp/IMGTHLA-3.50.0/
COPY hla-hd/src/scripts/update.dictionary.custom.sh /tmp/update.dictionary.custom.sh

# unpack tar
RUN tar -zvxf /tmp/hlahd.1.7.0.tar.gz && rm /tmp/hlahd.1.7.0.tar.gz && cd /tmp/hlahd.1.7.0 && sh install.sh

# now remove the original update.dictionary.sh script
RUN rm /tmp/hlahd.1.7.0/update.dictionary.sh && mv /tmp/update.dictionary.custom.sh /tmp/hlahd.1.7.0/

# rename the custom script and move hla dict dat file to this directory
RUN cd /tmp/hlahd.1.7.0/ && mv update.dictionary.custom.sh update.dictionary.sh && mv /tmp/IMGTHLA-3.50.0/hla.dat .

# export to PATH
ENV PATH="$PATH:/tmp/hlahd.1.7.0/bin"

# now install HLA dictionary; use the custom update.dictionary.sh 
RUN cd /tmp/hlahd.1.7.0/ && sh update.dictionary.sh

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
WORKDIR /work

ENTRYPOINT ["/usr/local/bin/_entrypoint.sh", "/sbin/matchhostfsowner"]
