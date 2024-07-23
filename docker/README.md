# Building a Minimal Docker Container Image for the Pipeline
## Preface

As most of the bioinformatics tools used in this neoantigen identification workflow are CLI programs, we would like to package them all into one container so the whole pipeline can be shipped around without much installation fuss. As AWS instances support Docker, we decided to use Docker for containerization. Singularity images would still work with the Nextflow pipelines in general, and Singularity is able to convert Docker images on the go so starting off with a Docker image is the best course of action. 

## The Steps:

Follow the steps below to replicate the containerization processes. Note that some binaries/programs need to be downloaded first prior to these steps.

1. Write a Dockerfile

### Strategy
a) Pull a minimal Linux image with micromamba installed from [micromamba.org](https://micromamba-docker.readthedocs.io), and then set up a base environment for FusionInspector and Arriba. 

b) Then switch to another environment specially made for `fusioncatcher` program and finish downloading its requisite database. 
