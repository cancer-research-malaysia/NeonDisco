FROM --platform=linux/x86_64 mambaorg/micromamba
USER root
LABEL maintainer="Suffian Azizan"
LABEL version="1.0"
LABEL description="container image of FusionInspector program v2.8.0"

# update Linux OS packages and install additional Linux system utilities with procps and also add parallel and finally remove cached package lists
RUN apt-get update && DEBIAN_FRONTEND=noninteractive apt-get install -y sudo && apt-get install -y gcc g++ perl automake make parallel procps wget curl libdb-dev bzip2 zlib1g zlib1g-dev unzip libbz2-dev liblzma-dev gfortran libreadline-dev libcurl4-openssl-dev libx11-dev \
libxt-dev x11-common libcairo2-dev libpng-dev libjpeg-dev pkg-config \
libxml2-dev libssl-dev libcurl4-openssl-dev pbzip2 git \
libharfbuzz-dev libfribidi-dev libfreetype6-dev libpng-dev libtiff5-dev libjpeg-dev \
&& apt-get install -y r-base r-base-dev && \
rm -rf /var/lib/apt/lists/*

## perl lib installations

RUN curl -L https://cpanmin.us | perl - App::cpanminus

RUN cpanm install PerlIO::gzip
RUN cpanm install Set::IntervalTree
RUN cpanm install DB_File
RUN cpanm install URI::Escape
RUN cpanm install Carp::Assert
RUN cpanm install JSON::XS.pm

# change user
USER $MAMBA_USER

# Configure Micromamba to use a single thread for package extraction
RUN micromamba config set extract_threads 1

# copy the env files into the container 
COPY --chown=$MAMBA_USER:$MAMBA_USER fusioninspector/context/base_env.yaml /tmp/base_env.yaml

# Create a new base environment based on the YAML file
RUN micromamba install -y -f /tmp/base_env.yaml && \
micromamba clean --all --yes

# activate the environment during container startup
ARG MAMBA_DOCKERFILE_ACTIVATE=1

# add conda bins to PATH
ENV PATH="/opt/conda/bin:/opt/conda/condabin:$PATH"

# install pip packages
# RUN pip install --no-cache-dir requests igv-reports==1.8.0

USER root
## set up tool config and deployment area
ENV SRC /usr/local/src
ENV BIN /usr/local/bin
ENV DATA /usr/local/data
RUN mkdir $DATA

## R installation:
# WORKDIR $SRC
# ENV R_VERSION=R-4.4.0

# RUN curl https://cran.r-project.org/src/base/R-4/$R_VERSION.tar.gz -o $R_VERSION.tar.gz && tar xvf $R_VERSION.tar.gz && cd $R_VERSION && ./configure && make && make install
    
RUN R --error -e 'install.packages("BiocManager", repos="http://cran.us.r-project.org")'
RUN R --error -e 'BiocManager::install("argparse")'
RUN R --error -e 'install.packages("tidyverse", repos="https://cran.csiro.au/")'
RUN R --error -e 'BiocManager::install("cowplot")'
RUN R --error -e 'BiocManager::install("ranger")'


######################
## Tool installations: #Only uncomment if they need to be custom installed and are not installed via micromamba
######################

## Bowtie2
# WORKDIR $SRC
# RUN wget https://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.3.3.1/bowtie2-2.3.3.1-linux-x86_64.zip/download -O bowtie2-2.3.3.1-linux-x86_64.zip && \
# unzip bowtie2-2.3.3.1-linux-x86_64.zip && \
# mv bowtie2-2.3.3.1-linux-x86_64/bowtie2* $BIN && \
# rm *.zip && \
# rm -r bowtie2-2.3.3.1-linux-x86_64

## Jellyfish
WORKDIR $SRC
RUN wget https://github.com/gmarcais/Jellyfish/releases/download/v2.2.7/jellyfish-2.2.7.tar.gz && \
tar xvf jellyfish-2.2.7.tar.gz && \
cd jellyfish-2.2.7/ && \
./configure && make && make install

## Picard tools
# WORKDIR $SRC
# RUN wget https://github.com/broadinstitute/picard/releases/download/2.20.3/picard.jar
# ENV PICARD_HOME $SRC

########
## Minimap2 (original author used v2.26)

WORKDIR $SRC
RUN curl -L https://github.com/lh3/minimap2/releases/download/v2.28/minimap2-2.28_x64-linux.tar.bz2 | tar -jxvf - && \
mv ./minimap2-2.28_x64-linux/minimap2 $BIN/

# K8 (original author used 0.2.4)
RUN curl -L https://github.com/attractivechaos/k8/releases/download/v1.2/k8-1.2.tar.bz2 | tar -jxf - && \
cp k8-1.2/k8-x86_64-`uname -s` $BIN/k8

# ## Salmon
# WORKDIR $SRC
# ENV SALMON_VERSION=1.10.3
# # author originally used 1.5.2
# RUN wget https://github.com/COMBINE-lab/salmon/releases/download/v${SALMON_VERSION}/Salmon-${SALMON_VERSION}_linux_x86_64.tar.gz && \
# tar xvf Salmon-${SALMON_VERSION}_linux_x86_64.tar.gz && \
# ln -sf $SRC/salmon-${SALMON_VERSION}_linux_x86_64/bin/salmon $BIN/.

# finally, copy a custom script made by FusionInspector authors to the image bin
COPY fusioninspector/src/sam_readname_cleaner.py $BIN/

# FusionInspector
# ENV FI_VERSION=2.9.0
# ENV FI_HASH=a43480df8dac6cfae0c01c2b636fd11de0d7bb98

ENV FI_VERSION=2.8.0
ENV FI_HASH=f798c9d9b51ddfdbe44e24094ad7dfb7f42b598c

RUN git clone --recursive https://github.com/FusionInspector/FusionInspector.git && \
cd FusionInspector/ && \
git checkout ${FI_HASH} && \
git submodule init && git submodule update && \
make && \
mv * $BIN

USER $USER
