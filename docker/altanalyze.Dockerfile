FROM python:2
LABEL maintainer="Suffian Azizan"
LABEL version="1.0"
LABEL description="container image of AltAnalyze and SNAF splicing neoantigen pipeline by Frank Li (frankligy123)"

# change to root user
USER root

# copy script
COPY altanalyze/src/AltAnalyze.sh /home/app/AltAnalyze.sh

# install AltAnalyze
WORKDIR /home/app
RUN apt-get clean \
    && apt-get update \
    && apt-get install -y parallel nano less build-essential tar wget curl pigz gzip zip unzip gcc g++ bzip2 procps git cmake locales coreutils gawk grep sed
RUN cd /home/app \
    && curl -o requirements.txt https://raw.githubusercontent.com/frankligy/SNAF/main/AltAnalyze/requirements_slim_test.txt \
    && apt-get install -y libigraph0-dev \
    && pip install python-igraph==0.7.1.post6 \
    && pip install lxml \
    && pip install patsy \
    && pip install fastcluster \
    && pip install --no-cache-dir -r requirements.txt \
    && git clone https://github.com/nsalomonis/altanalyze.git

RUN cd altanalyze \
    && python AltAnalyze.py --species Hs --update Official --version EnsMart91 --additional all \
    && curl -o Hs.bed https://raw.githubusercontent.com/frankligy/SNAF/main/AltAnalyze/Hs.bed \
    && mv Hs.bed AltDatabase/EnsMart91/ensembl/Hs \
    && curl -o ../prune.py https://raw.githubusercontent.com/frankligy/SNAF/main/AltAnalyze/prune.py \
    && chmod 777 /home/app/prune.py \
    && chmod 777 /home/app/AltAnalyze.sh \
    && chmod -R 777 /home/app/altanalyze 

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
WORKDIR /home/app

ENTRYPOINT ["/sbin/matchhostfsowner", "/home/app/AltAnalyze.sh"]
