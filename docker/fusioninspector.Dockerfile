FROM --platform=linux/x86_64 mambaorg/micromamba:git-911a014-bookworm-slim
USER root
LABEL maintainer="Suffian Azizan"
LABEL version="1.0"
LABEL description="container image of FusionInspector program v2.8.0"

# update Linux OS packages and install additional Linux system utilities with procps and also add parallel and finally remove cached package lists
RUN apt-get update && DEBIAN_FRONTEND=noninteractive apt-get install -y sudo && apt-get install -y \
automake \
build-essential \
bzip2 \
cmake \
curl \
default-jre \
fort77 \
ftp \
gcc \
g++ \
gfortran \
git \
libblas-dev \
libbz2-dev \
libcairo2-dev \
libcurl4-openssl-dev \
libdb-dev \
libfreetype6-dev \
libfribidi-dev \
libghc-zlib-dev \
libharfbuzz-dev \
libjpeg-dev \
liblzma-dev \
libncurses-dev \
libncurses5-dev \
libpcre3-dev \
libpng-dev \
libreadline-dev \
libssl-dev \
libtbb-dev \
libtiff5-dev \
libx11-dev \
libxml2-dev \
libxt-dev \
libzmq3-dev \
make \
nano \
pbzip2 \
perl \
pkg-config \
procps \
r-base \
r-base-dev \
rsync \
texlive \
texlive-latex-base \
texlive-latex-extra \
tzdata \
unzip \
wget \
x11-common \
zlib1g-dev \
&& \
rm -rf /var/lib/apt/lists/* /var/log/dpkg.log

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

################ R installation:
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

## autoconf 2.69 needed for trinity installed version of htslib
WORKDIR $SRC
RUN wget http://ftp.gnu.org/gnu/autoconf/autoconf-2.69.tar.gz && \
      tar xf autoconf* && \
      cd autoconf-2.69 && \
      sh configure --prefix /usr/local && \
      make install

#### INSTALL TRINITY
WORKDIR $SRC
##### Docker build note: The original trinityrnaseq github repo hardcoded the destination dir for the install in the python script trinity_installer.py, which is /usr/local/bin, and this actually clutters the bin dir, and prevents FusionInspector from being installed at the same location due to the similarity of several package directory names, so mv command will return an error. I forked the repo to my personal github, and modified the script so that the package is installed in one directory under $BIN. Note that this version of Trinity is the last maintained version as of 2024 so there is no possibility of upgrading and thus breaking this code.
ENV TRINITY_VERSION="2.15.2-alpha.1"
ENV TRINITY_CO=b63de51cfd600e69161fcd1c163114d79355415c

RUN git clone --recursive https://github.com/sufyazi/trinityrnaseq.git && \
    cd trinityrnaseq && \
    git checkout ${TRINITY_CO} && \
    git submodule init && git submodule update && \
    git submodule foreach --recursive git submodule init && \
    git submodule foreach --recursive git submodule update && \
    rm -rf ./trinity_ext_sample_data && \
    make && make plugins && \
    make install && \
    cd ../ && rm -r trinityrnaseq && \
    mv /usr/local/bin/Trinity-pkg /usr/local/bin/Trinity-${TRINITY_VERSION}-pkg

# set Trinity executable in PATH
ENV TRINITY_HOME /usr/local/bin/Trinity-${TRINITY_VERSION}-pkg
ENV PATH=${TRINITY_HOME}:${PATH}

## GMAP
ENV GSNAP_VER 2021-07-23
WORKDIR $SRC
RUN GMAP_URL="http://research-pub.gene.com/gmap/src/gmap-gsnap-$GSNAP_VER.tar.gz" && \
wget $GMAP_URL && \
tar xvf gmap-gsnap-$GSNAP_VER.tar.gz && \
cd gmap-$GSNAP_VER && ./configure --prefix=`pwd` && make && make install && \
cp bin/* $BIN/

## Jellyfish
WORKDIR $SRC
RUN wget https://github.com/gmarcais/Jellyfish/releases/download/v2.2.7/jellyfish-2.2.7.tar.gz && \
tar xvf jellyfish-2.2.7.tar.gz && \
cd jellyfish-2.2.7/ && \
./configure && make && make install

###############################################
## Minimap2 (original author used v2.26)
WORKDIR $SRC
RUN curl -L https://github.com/lh3/minimap2/releases/download/v2.28/minimap2-2.28_x64-linux.tar.bz2 | tar -jxvf - && \
mv ./minimap2-2.28_x64-linux/minimap2 $BIN/

# K8 (original author used 0.2.4)
RUN curl -L https://github.com/attractivechaos/k8/releases/download/v1.2/k8-1.2.tar.bz2 | tar -jxf - && \
cp k8-1.2/k8-x86_64-`uname -s` $BIN/k8

####### NOW INSTALL FUSIONINSPECTOR ################
ENV FI_VERSION=2.9.0
ENV FI_HASH=a43480df8dac6cfae0c01c2b636fd11de0d7bb98
# ENV FI_VERSION=2.8.0
# ENV FI_HASH=f798c9d9b51ddfdbe44e24094ad7dfb7f42b598c

RUN git clone --recursive https://github.com/FusionInspector/FusionInspector.git && \
cd FusionInspector/ && \
git checkout ${FI_HASH} && \
git submodule init && git submodule update && \
make && \
mkdir ${BIN}/FusionInspector_${FI_VERSION}-pkg && \
mv * ${BIN}/FusionInspector_${FI_VERSION}-pkg/

# copy a custom script made by FusionInspector authors to the image bin
COPY fusioninspector/src/scripts/sam_readname_cleaner.py ${BIN}/FusionInspector_${FI_VERSION}-pkg/

# cleanup tar.gz
WORKDIR ${SRC}
RUN rm *.tar.gz

###########################################################################################################
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
