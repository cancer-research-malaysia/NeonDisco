# Building a Software Stack for the Pipeline
## Preface

As most of the bioinformatics tools used in this neoantigen identification workflow are CLI programs, it would be more convenient and allows a more replicable pipeline if the software components needed to run the pipeline are bundled together as a set of Docker container images. I decided to custom-build the Docker images as several software components that make up the stack described in the workflow are now a few official releases ahead from the versions that were used in the original fusion transcript detection workflow. This does not incur additional time cost at all as many of the tools used in the workflow are trivially installable via Conda. I decided to use the slim Debian Linux image that is shipped with Micromamba environment/package manager as base image for every Docker container image I built in the software stack. The ease of Docker builder also allows us to rebuild each Docker image if we decide to change or upgrade the version of the tool used as part of the software stack for the workflow. Docker is chosen as the image building platform due to its ease of use (with Dockerfiles), and the fact that Apptainer/Singularity can make use of Docker images on the fly, would also enable direct compatibility of this software stack for running on high performance computing (HPC) clusters. 

### Software Stack (currently used versions as of DATE_PLACEHOLDER)

|Purpose| |
|:----|:----|
|Tool| |
|Input Data| |
|Some Dependencies(refer documentations for the full list)| |
|Installation typeDetermination of HLA-alleles| |
|HLA-HD (v1.7.0)| |
|.bam files of MHC region reads and unmapped reads extracted from WES or RNA-Seq data| |
|HLA-HD database| |
|Downloaded pre-compiled binary; added to PATHDetection of fusion transcripts| |
|Arriba (v2.3.0)| |
|.fastq file of RNA-Seq data (raw file, or converted from bam but make sure that unmapped reads are not discarded)| |
|STARCTAT_Genome_Lib & Arriba database| |
|Downloaded pre-compiled binary; added to PATHFusionCatcher (v1.33)| |
|Python 2.7.x (>2.6, <3.0)FusionCatcher database| |
|Downloaded pre-compiled binary; added to PATH (Python based – executed within conda env)In silico validation of fusion calls| |
|FusionInspector (v2.8.0)| |
|.fastq or RNA-Seq data (raw file);list of fusions (e.g. A--B, C--D)| |
|CTAT_Genome_Lib| |
|Singularity imagePrediction of fusion transcript coding potential| |
|AGFusion (v1.4.1)| |
|Raw output file from Arriba / FusionCatcher| |
|Python >=3.7AGFusion database| |
|Downloaded pre-compiled binary; added to PATH (Python based – executed within conda env)Prediction of immunogenic fusion neoantigen| |
|pVacFuse (v4.0.4)| |
|Output folder from AGFusion| |
|--| |
|Docker image |



## The Steps:

Follow the steps below to replicate the containerization processes. Note that some binaries/programs need to be downloaded first prior to these steps.

1. Write a Dockerfile

### Strategy
a) Pull a minimal Linux image with micromamba installed from [micromamba.org](https://micromamba-docker.readthedocs.io), and then set up a base environment for FusionInspector and Arriba. 

b) Then switch to another environment specially made for `fusioncatcher` program and finish downloading its requisite database. 
