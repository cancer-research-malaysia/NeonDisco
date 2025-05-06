FROM trinityctat/starfusion
USER root
LABEL maintainer="Suffian Azizan"
LABEL version="1.0"
LABEL description="container image of STAR Fusion from TrinityCTAT"

# change to root user
USER root
# update Debian OS packages and install additional Linux system utilities, then finally remove cached package lists
RUN apt-get update && apt-get install -y --no-install-recommends \
tar wget curl pigz gzip zip unzip gcc g++ bzip2 procps coreutils gawk grep sed less nano \
&& rm -rf /var/lib/apt/lists/*

## copy files
COPY starfusion/src/fusion_lib.Mar2021.dat.gz /tmp/fusion_lib.Mar2021.dat.gz
COPY starfusion/src/AnnotFilterRule.pm /tmp/AnnotFilterRule.pm

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

# Set workdir
WORKDIR /home/app

# Create the entrypoint directory if it doesn't exist
RUN mkdir -p /usr/local/bin

# Create the entrypoint script
RUN echo '#!/bin/sh' > /usr/local/bin/_entrypoint.sh && \
    echo 'exec "$@"' >> /usr/local/bin/_entrypoint.sh && \
    chmod +x /usr/local/bin/_entrypoint.sh

# Verify the entrypoint script exists and is executable
RUN ls -la /usr/local/bin/_entrypoint.sh && \
    cat /usr/local/bin/_entrypoint.sh

ENTRYPOINT ["/usr/local/bin/_entrypoint.sh", "/sbin/matchhostfsowner"]
