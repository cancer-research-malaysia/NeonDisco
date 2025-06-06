FROM mambaorg/micromamba
USER root
LABEL maintainer="Suffian Azizan"
LABEL version="2.0"
LABEL description="container image of FusionCatcher program v1.33"

# change to root user
USER root
# update Debian OS packages and install additional Linux system utilities, then finally remove cached package lists

RUN apt-get update && DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends \
tar wget curl pigz gzip zip unzip gcc g++ bzip2 procps coreutils gawk grep sed locales nano less \
&& rm -rf /var/lib/apt/lists/* \
&& echo "en_US.UTF-8 UTF-8" > /etc/locale.gen && locale-gen

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

# also store the bash script to download references as root
COPY --chown=mambauser:mambauser fusioncatcher/src/scripts/download-human-db-patched.sh /home/app/download-human-db-patched.sh
RUN chmod +x /home/app/download-human-db-patched.sh

# replace get_pcawg.py with the one called v2 (removed codes that access inactive URLs)
COPY --chown=mambauser:mambauser fusioncatcher/src/scripts/get_pcawg-v2.py /tmp/get_pcawg-v2.py
RUN chmod +x /tmp/get_pcawg-v2.py
RUN rm /opt/conda/bin/get_pcawg.py && mv /tmp/get_pcawg-v2.py /opt/conda/bin/get_pcawg.py

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
