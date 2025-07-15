FROM mambaorg/micromamba:git-911a014-bookworm-slim

LABEL maintainer="Suffian Azizan"
LABEL version="4.0"
LABEL description="minimal container image of Python env for CRMY - revision removed MatchHostFsOwner"

# change to root user
USER root

# update Debian OS packages and install additional Linux system utilities with procps; also install R, then finally remove cached package lists
RUN apt-get update && DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends build-essential procps curl wget tar pigz gzip zip unzip gcc g++ bzip2 procps coreutils gawk grep sed nano less \
&& rm -rf /var/lib/apt/lists/*


ARG NEW_MAMBA_USER=ec2-user
ARG NEW_MAMBA_USER_ID=1000
ARG NEW_MAMBA_USER_GID=1000

RUN if grep -q '^ID=alpine$' /etc/os-release; then \
      # alpine does not have usermod/groupmod
      apk add --no-cache --virtual temp-packages shadow; \
    fi && \
    usermod "--login=${NEW_MAMBA_USER}" "--home=/home/${NEW_MAMBA_USER}" \
        --move-home "-u ${NEW_MAMBA_USER_ID}" "${MAMBA_USER}" && \
    groupmod "--new-name=${NEW_MAMBA_USER}" \
        "-g ${NEW_MAMBA_USER_GID}" "${MAMBA_USER}" && \
    if grep -q '^ID=alpine$' /etc/os-release; then \
      # remove the packages that were only needed for usermod/groupmod
      apk del temp-packages; \
    fi && \
    # Update the expected value of MAMBA_USER for the
    # _entrypoint.sh consistency check.
    echo "${NEW_MAMBA_USER}" > "/etc/arg_mamba_user" && \
    :
ENV MAMBA_USER=$NEW_MAMBA_USER

# change user
USER $MAMBA_USER

# Configure Micromamba to use a single thread for package extraction
RUN micromamba config set extract_threads 1

# copy the env file into the container 
COPY --chown=$MAMBA_USER:$MAMBA_USER pyenv/context/base_env.yaml /tmp/base_env.yaml

# Create a new base environment based on the YAML file
RUN micromamba install -y -f /tmp/base_env.yaml && \
micromamba clean --all --yes

# activate the environment during container startup
ARG MAMBA_DOCKERFILE_ACTIVATE=1

# add conda bins to PATH
ENV PATH="/opt/conda/bin:/opt/conda/condabin:$PATH"

# Docker suffers from absolutely atrocious way of consolidating the paradigm of restricting privileges when running containers (rootless mode) with writing outputs to bound host volumes without using Docker volumes or other convoluted workarounds.

# Fortunately there is this tool that removes this altogether and helps matches the UID and GID of whoever is running the container image on a host machine

# Install MatchHostFsOwner. Using version 1.0.1
# See https://github.com/FooBarWidget/matchhostfsowner/releases
# ADD https://github.com/FooBarWidget/matchhostfsowner/releases/download/v1.0.1/matchhostfsowner-1.0.1-x86_64-linux.gz /sbin/matchhostfsowner.gz
# RUN gunzip /sbin/matchhostfsowner.gz && \
#   chown root: /sbin/matchhostfsowner && \
#   chmod +x /sbin/matchhostfsowner
# RUN addgroup --gid 9999 app && \
#   adduser --uid 9999 --gid 9999 --disabled-password --gecos App app

# set workdir
WORKDIR /home/ec2-user

ENTRYPOINT ["/usr/local/bin/_entrypoint.sh"]
