FROM mambaorg/micromamba:git-911a014-bookworm-slim
USER root
LABEL maintainer="Suffian Azizan"
LABEL version="2.0"
LABEL description="container image of Spladder program v3.1.0"

# change to root user
USER root
# update Debian OS packages and install additional Linux system utilities, then finally remove cached package lists
RUN apt-get update && DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends \
git tar wget curl pigz gzip zip unzip gcc g++ bzip2 procps coreutils gawk grep sed less nano \
&& rm -rf /var/lib/apt/lists/*

# change user
USER $MAMBA_USER

# Configure Micromamba to use a single thread for package extraction
RUN micromamba config set extract_threads 1

# copy the env file into the container 
COPY --chown=$MAMBA_USER:$MAMBA_USER spladder/context/base_env.yaml /tmp/base_env.yaml

# Create a new base environment based on the YAML file
RUN micromamba install -y -f /tmp/base_env.yaml && \
micromamba clean --all --yes

# activate the environment during container startup
ARG MAMBA_DOCKERFILE_ACTIVATE=1

# clone and install spladder from source
RUN pip install spladder==3.1.0

# add conda bins to PATH
ENV PATH="/opt/conda/bin:/opt/conda/condabin:$PATH"

# Docker suffers from absolutely atrocious way of consolidating the paradigm of restricting privileges when running containers (rootless mode) with writing outputs to bound host volumes without using Docker volumes or other convoluted workarounds.

# Fortunately there is this tool that removes this altogether and helps matches the UID and GID of whoever is running the container image on a host machine

USER root
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

ENTRYPOINT ["/usr/local/bin/_entrypoint.sh", "/sbin/matchhostfsowner"]
