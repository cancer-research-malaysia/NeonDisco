FROM mambaorg/micromamba
LABEL maintainer="Suffian Azizan"
LABEL version="3.0"
LABEL description="container image of FusionCatcher program v1.33 - revision removed MatchHostoffsOwner"

# change to root user
USER root

# update Debian OS packages and install additional Linux system utilities, then finally remove cached package lists
RUN apt-get update && DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends \
build-essential tar wget curl pigz gzip zip unzip gcc g++ bzip2 procps coreutils gawk grep sed locales nano less \
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
COPY --chown=$MAMBA_USER:$MAMBA_USER fusioncatcher/src/scripts/download-human-db-patched.sh /home/ec2-user/download-human-db-patched.sh
RUN chmod +x /home/ec2-user/download-human-db-patched.sh

# replace get_pcawg.py with the one called v2 (removed codes that access inactive URLs)
COPY --chown=$MAMBA_USER:$MAMBA_USER fusioncatcher/src/scripts/get_pcawg-v2.py /tmp/get_pcawg-v2.py
RUN chmod +x /tmp/get_pcawg-v2.py
RUN rm /opt/conda/bin/get_pcawg.py && mv /tmp/get_pcawg-v2.py /opt/conda/bin/get_pcawg.py

# set workdir
WORKDIR /home/ec2-user

ENTRYPOINT ["/usr/local/bin/_entrypoint.sh"]
