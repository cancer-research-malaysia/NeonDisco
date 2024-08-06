# Building a Software Stack for the Pipeline
## Preface

As most of the bioinformatics tools used in this neoantigen identification workflow are CLI programs, it would be more convenient and allows a more replicable pipeline if the software components needed to run the pipeline are bundled together as a set of Docker container images. I decided to custom-build the Docker images as several software components that make up the stack described in the workflow are now a few official releases ahead from the versions that were used in the original fusion transcript detection workflow. This does not incur additional time cost at all as many of the tools used in the workflow are trivially installable via Conda. I decided to use the slim Debian Linux image that is shipped with Micromamba environment/package manager as base image for every Docker container image I built in the software stack. The ease of Docker builder also allows us to rebuild each Docker image if we decide to change or upgrade the version of the tool used as part of the software stack for the workflow. Docker is chosen as the image building platform due to its ease of use (with Dockerfiles), and the fact that Apptainer/Singularity can make use of Docker images on the fly, would also enable direct compatibility of this software stack for running on high performance computing (HPC) clusters. 

### Software Stack (currently used versions as of DATE_PLACEHOLDER)

<div class="divTable">
    <div class="row">
        <div class="cell">Purpose </div>
        <div class="cell">Tool </div>
        <div class="cell">Databases  </div>
        <div class="cell">Input Data </div>
    </div>
    <div class="row">
        <div class="cell">Determination of HLA-alleles </div>
        <div class="cell">HLA-HD (v1.7.0) </div>
        <div class="cell">HLA-HD database </div>
        <div class="cell">.bam files of MHC region reads and unmapped reads extracted from WES or RNA-Seq data </div>
    </div>
    <div class="row">
        <div class="cell">Detection of fusion transcripts </div>
        <div class="cell">Arriba (v2.3.0) </div>
        <div class="cell">STAR Index Directory, CTAT Genome Library &amp; Arriba database  </div>
        <div class="cell">.fastq file of RNA-Seq data (raw file, or converted from bam but make sure that unmapped reads are not discarded) </div>
    </div>
    <div class="row">
        <div class="cell"></div>
        <div class="cell">FusionCatcher (v1.33) </div>
        <div class="cell">FusionCatcher database</div>
        <div class="cell"></div>
    </div>
    <div class="row">
        <div class="cell">In silico validation of fusion calls </div>
        <div class="cell">FusionInspector (v2.8.0) </div>
        <div class="cell">CTAT Genome Library</div>
        <div class="cell">.fastq or RNA-Seq data (raw file);  </div>
    </div>
    <div class="row">
        <div class="cell">Prediction of fusion transcript coding potential </div>
        <div class="cell">AGFusion (v1.4.1) </div>
        <div class="cell">AGFusion database </div>
        <div class="cell">Raw output file from Arriba / FusionCatcher </div>
    </div>
    <div class="row">
        <div class="cell">Prediction of immunogenic fusion neoantigen </div>
        <div class="cell">pVacFuse (v4.0.4) </div>
        <div class="cell">Several prediction tool databases (installed by default into the Docker image)</div>
        <div class="cell">Output folder from AGFusion</div>
    </div>
</div>



## The Steps:

Follow the steps below to replicate the containerization processes. Note that some binaries/programs need to be downloaded first prior to these steps.

1. Write a Dockerfile

### Strategy
a) Pull a minimal Linux image with micromamba installed from [micromamba.org](https://micromamba-docker.readthedocs.io), and then set up a base environment for FusionInspector and Arriba. 

b) Then switch to another environment specially made for `fusioncatcher` program and finish downloading its requisite database. 
