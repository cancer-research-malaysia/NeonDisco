FROM mambaorg/micromamba:git-911a014-bookworm-slim
USER root
LABEL maintainer="Suffian Azizan"
LABEL version="2.0"
LABEL description="container image of tools for HLA typer ArcasHLA v0.6.0 - revision removed MatchHostFsOwner"

# change to root user
USER root

# update Debian OS packages and install additional Linux system utilities, then finally remove cached package lists
RUN apt-get update && DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends \
build-essential tar wget curl pigz gzip zip unzip gcc g++ bzip2 procps git cmake locales coreutils gawk grep sed nano less git-lfs \
&& rm -rf /var/lib/apt/lists/* \
&& echo "en_US.UTF-8 UTF-8" > /etc/locale.gen && locale-gen

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
COPY --chown=$MAMBA_USER:$MAMBA_USER arcashla/context/base_env.yaml /tmp/base_env.yaml

# Create a new base environment based on the YAML file
RUN micromamba install -y -f /tmp/base_env.yaml && \
micromamba clean --all --yes

# activate the environment during container startup
ARG MAMBA_DOCKERFILE_ACTIVATE=1

# add conda bins to PATH
ENV PATH="/opt/conda/bin:/opt/conda/condabin:$PATH"

# install dependencies for arcasHLA
RUN pip install --no-cache-dir biopython==1.83 numpy==1.26.3 scipy==1.12.0

### delete and then copy reference static script into the container
###### UPDATE: the 0.6.0-2 version of arcasHLA has fixed the reference command so there is no need to monkey fix it below
# RUN rm -f /opt/conda/share/arcas-hla-0.6.0-1/scripts/reference.py
# COPY --chown=$MAMBA_USER:$MAMBA_USER arcashla/src/reference-static.py /opt/conda/share/arcas-hla-0.6.0-1/scripts/reference.py

# download the IMGTHLA reference file you want to use
RUN wget -O IMGTHLA-3.61.0-alpha.tar.gz https://github.com/ANHIG/IMGTHLA/archive/refs/tags/v3.61.0-alpha.tar.gz && \
    tar -xzf IMGTHLA-3.61.0-alpha.tar.gz && \
    mv IMGTHLA-3.61.0-alpha /opt/conda/share/arcas-hla-0.6.0-2/dat/IMGTHLA && \
    rm IMGTHLA-3.61.0-alpha.tar.gz
# unzip the hla.dat.zip file
RUN cd /opt/conda/share/arcas-hla-0.6.0-2/dat/IMGTHLA && unzip hla.dat.zip
# update arcasHLA reference (version 3.61.0)
RUN arcasHLA reference --update

# change user to root
# USER root

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
