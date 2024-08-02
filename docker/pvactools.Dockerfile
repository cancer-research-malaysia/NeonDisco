FROM --platform=linux/x86_64 mambaorg/micromamba
USER root
LABEL maintainer="Suffian Azizan"
LABEL version="1.0"
LABEL description="container image of pVacTools program v4.0.4"

# change to root user
USER root

# update Debian OS packages and install additional Linux system utilities, then finally remove cached package lists
RUN apt-get update && DEBIAN_FRONTEND=noninteractive apt-get install -y --no-install-recommends \
tar wget curl pigz gzip zip unzip gcc g++ bzip2 procps tcsh gawk git \
&& rm -rf /var/lib/apt/lists/*

# change user
USER $MAMBA_USER

# Configure Micromamba to use a single thread for package extraction
RUN micromamba config set extract_threads 1

# copy the env file into the container 
COPY --chown=$MAMBA_USER:$MAMBA_USER pvactools/context/base_env.yaml /tmp/base_env.yaml

# Create a new base environment based on the YAML file
RUN micromamba install -y -f /tmp/base_env.yaml && \
micromamba clean --all --yes

# activate the environment during container startup
ARG MAMBA_DOCKERFILE_ACTIVATE=1

# install pip packages
RUN pip install pvactools==4.0.4
RUN pip install git+https://github.com/griffithlab/bigmhc.git#egg=bigmhc
RUN pip install git+https://github.com/griffithlab/deepimmuno.git#egg=deepimmuno

# add conda bins to PATH
ENV PATH="/opt/conda/bin:/opt/conda/condabin:$PATH"

#### Now install IEDB binding prediction class I and class I archives
RUN wget https://downloads.iedb.org/tools/mhci/3.1.5/IEDB_MHC_I-3.1.5.tar.gz && \
tar -zxvf IEDB_MHC_I-3.1.5.tar.gz && \
cd mhc_i && \
./configure

RUN wget https://downloads.iedb.org/tools/mhcii/3.1.11/IEDB_MHC_II-3.1.11.tar.gz && \
tar -zxvf IEDB_MHC_II-3.1.11.tar.gz && cd mhc_ii
RUN ./configure.py

# download mhcflurry datasets and trained models
RUN mhcflurry-downloads fetch

