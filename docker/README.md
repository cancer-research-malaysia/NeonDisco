# Building a Software Stack for the Pipeline
## Preface

As most of the bioinformatics tools used in this neoantigen identification workflow are CLI programs, it would be more convenient and allows a more replicable pipeline if the software components needed to run the pipeline are bundled together as a set of Docker container images. I decided to custom-build the Docker images as several software components that make up the stack described in the workflow are now a few official releases ahead from the versions that were used in the original fusion transcript detection workflow. This does not incur additional time cost at all as many of the tools used in the workflow are trivially installable via Conda. I decided to use the slim Debian Linux image that is shipped with Micromamba environment/package manager as base image for every Docker container image I built in the software stack. The ease of Docker builder also allows us to rebuild each Docker image if we decide to change or upgrade the version of the tool used as part of the software stack for the workflow. Docker is chosen as the image building platform due to its ease of use (with Dockerfiles), and the fact that Apptainer/Singularity can make use of Docker images on the fly, would also enable direct compatibility of this software stack for running on high performance computing (HPC) clusters. 

### Software Stack (currently used versions as of DATE_PLACEHOLDER)

<table>
    <tr>
        <td>Purpose </td>
        <td>Tool </td>
        <td>Databases  </td>
        <td>Input Data </td>
    </tr>
    <tr>
        <td>Determination of HLA-alleles </td>
        <td>HLA-HD (v1.7.0) </td>
        <td>HLA-HD database </td>
        <td>.bam files of MHC region reads and unmapped reads extracted from WES or RNA-Seq data </td>
    </tr>
    <tr>
        <td>Detection of fusion transcripts </td>
        <td>Arriba (v2.3.0) </td>
        <td>STAR Index Directory, CTAT Genome Library &amp; Arriba database  </td>
        <td>.fastq file of RNA-Seq data (raw file, or converted from bam but make sure that unmapped reads are not discarded) </td>
    </tr>
    <tr>
        <td></td>
        <td>FusionCatcher (v1.33) </td>
        <td>FusionCatcher database</td>
        <td></td>
    </tr>
    <tr>
        <td>In silico validation of fusion calls </td>
        <td>FusionInspector (v2.8.0) </td>
        <td>CTAT Genome Library</td>
        <td>.fastq or RNA-Seq data (raw file);  </td>
    </tr>
    <tr>
        <td>Prediction of fusion transcript coding potential </td>
        <td>AGFusion (v1.4.1) </td>
        <td>AGFusion database </td>
        <td>Raw output file from Arriba / FusionCatcher </td>
    </tr>
    <tr>
        <td>Prediction of immunogenic fusion neoantigen </td>
        <td>pVacFuse (v4.0.4) </td>
        <td>Several prediction tool databases (installed by default into the Docker image)</td>
        <td>Output folder from AGFusion</td>
    </tr>
</table>



## The Steps:

Follow the steps below to replicate the containerization processes. Note that some binaries/programs need to be downloaded first prior to these steps.

1. Write a Dockerfile

### Strategy
a) Pull a minimal Linux image with micromamba installed from [micromamba.org](https://micromamba-docker.readthedocs.io), and then set up a base environment for FusionInspector and Arriba. 

b) Then switch to another environment specially made for `fusioncatcher` program and finish downloading its requisite database. 
