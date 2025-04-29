FROM mambaorg/micromamba:git-911a014-bookworm-slim
USER root
LABEL maintainer="Suffian Azizan"
LABEL version="2.0"
LABEL description="container image of AGFusion program v1.4.1"

# change to root user
USER root
# update Debian OS packages and install additional Linux system utilities, then finally remove cached package lists
RUN apt-get update && DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends \
tar wget curl pigz gzip zip unzip gcc g++ bzip2 procps coreutils gawk grep sed \
&& rm -rf /var/lib/apt/lists/*

# change user
USER $MAMBA_USER

# Configure Micromamba to use a single thread for package extraction
RUN micromamba config set extract_threads 1

# copy the env file into the container 
COPY --chown=$MAMBA_USER:$MAMBA_USER agfusion/context/base_env.yaml /tmp/base_env.yaml

# Create a new base environment based on the YAML file
RUN micromamba install -y -f /tmp/base_env.yaml && \
micromamba clean --all --yes

# activate the environment during container startup
ARG MAMBA_DOCKERFILE_ACTIVATE=1

# install pip packages
RUN pip install agfusion
RUN pip install mysqlclient

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


# change user
USER app

# change user back to root
USER root

# add pfam file
COPY agfusion/src/Pfam-A.clans.tsv /tmp/Pfam-A.clans.tsv

# replace the cli.py of agfusion with my edited one that handles running agfusion cli with no arguments
COPY agfusion/src/cli-v2.py /tmp/cli-v2.py
RUN mv /tmp/cli-v2.py /tmp/cli.py && rm -rf /opt/conda/lib/python3.12/site-packages/agfusion/cli.py && mv /tmp/cli.py /opt/conda/lib/python3.12/site-packages/agfusion/


# install pyensembl
RUN pyensembl install --release 113 --species homo_sapiens
# build agfusion db
RUN agfusion build -d . -s homo_sapiens -r 113 --pfam Pfam-A.clans.tsv

# set workdir
WORKDIR /home/app

ENTRYPOINT ["/usr/local/bin/_entrypoint.sh", "/sbin/matchhostfsowner"]
