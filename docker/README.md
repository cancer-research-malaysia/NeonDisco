# Building a Software Stack for the Pipeline
## Preface

As most of the bioinformatics tools used in this neoantigen identification workflow are CLI programs, it would be more convenient and allows a more replicable pipeline if the software components needed to run the pipeline are bundled together as a set of Docker container images. I decided to custom-build the Docker images as several software components that make up the stack described in the workflow are now a few official releases ahead from the versions that were used in the original fusion transcript detection workflow. This does not incur additional time cost at all as many of the tools used in the workflow are trivially installable via Conda. I decided to use the slim Debian Linux image that is shipped with Micromamba environment/package manager as base image for every Docker container image I built in the software stack. The ease of Docker builder also allows us to rebuild each Docker image if we decide to change or upgrade the version of the tool used as part of the software stack for the workflow. Docker is chosen as the image building platform due to its ease of use (with Dockerfiles), and the fact that Apptainer/Singularity can make use of Docker images on the fly, would also enable direct compatibility of this software stack for running on high performance computing (HPC) clusters. 

### Software Stack (currently used versions as of 2024-08-06)

<table>
    <tr>
        <th><b>Purpose</b></th>
        <th><b>Tool</b></th>
        <th><b>Databases</b></th>
        <th><b>Input Data</b></th>
    </tr>
    <tr>
        <td>Determination of HLA-alleles </td>
        <td>HLA-HD (v1.7.0) </td>
        <td>HLA-HD database </td>
        <td>.bam files of MHC region reads and unmapped reads extracted from WES or RNA-Seq data </td>
    </tr>
    <tr>
        <td rowspan="2">Detection of fusion transcripts </td>
        <td>Arriba (v2.3.0) </td>
        <td>STAR Index Directory, CTAT Genome Library &amp; Arriba database  </td>
        <td rowspan="2">.fastq file of RNA-Seq data (raw file, or converted from bam but make sure that unmapped reads are not discarded) </td>
    </tr>
    <tr>
        <td>FusionCatcher (v1.33) </td>
        <td>FusionCatcher database</td>
    </tr>
    <tr>
        <td>In silico validation of fusion calls </td>
        <td>FusionInspector (v2.8.0) </td>
        <td>CTAT Genome Library</td>
        <td>.fastq of RNA-Seq data (raw file);  </td>
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

## Software Stack Building Strategy

Follow the steps below to replicate the containerization processes. Note that some binaries/programs need to be downloaded first as packages locally prior to container image building. This section serves primarily to document what I have carried out to build the software stack, as ideally, an end-user of the Nextflow pipeline would only need to pull the built container images off my personal Dockerhub. 

1. Write App-specific Dockerfiles under One Build Context

Docker `build` allows direct reference of the Dockerfile from which a particular image is to be built with the `-f` flag so we make make use of this by writing several Dockerfiles even under the same build context (i.e. one docker directory). 

All of the Docker container images has a minimal Linux image (`debian-slim`) with `micromamba` installed from [micromamba.org](https://micromamba-docker.readthedocs.io). The currently used version is `mambaorg/micromamba:git-911a014-bookworm-slim`. 

Each Dockerfile follows roughly the same layering structure, with minimal Linux package updating and installation followed by micromamba installation and setup. Under each app-specific subcontext lies an app-specific `base_env.yml` that is used for the Micromamba base environment setup. Several apps require additional database downloads (*pvactools*, *HLA-HD*) so this process is run following Micromamba environment setup. 

Finally, each Dockerfile's framework incorporate a tiny program called **Matchhostfsowner** to enable instant matching of user IDs on the host machine where each containerized app is run, with the internal container's file permissions. This will be expounded later. 

2. Build Application Images for the Stack

Using Docker's build engine, these Dockerfiles can act as the blueprint for the building of the containerized app images. Simply run:

```bash
docker build --no-cache -f APP.Dockerfile  -t DOCKERHUB-USER/APP-IMAGE .
```

> Note: The `.` implies '*current directory*' so run this command inside the `docker` directory of this repo. Additionally, the option `--no-cache` can be omitted when rebuilding images after editing Dockerfiles, which would only build changed layers and circumvent redownloading large resources again.

3. Check Application Images and Push to Dockerhub

With `docker run --rm -it sufyazi/fusioninspector-crm` for example, we can check the built image locally and see if all the required tools have been installed within the image. If all is fine then the image can then be pushed to Dockerhub by just doing `docker push sufyazi/fusioninspector-crm` for instance. Make sure that you have previously logged onto Dockerhub with your credentials and you are pushing to your own hub.


### Bonus Guide: How to run a command in a containerized instance

These dockerized images will be used to run commands in a containerized environment. This means we need to first spin up a container instance using `docker run`, and then running `docker exec` to send commands you want to run into the container. 

As we need to run as our own UID (Docker runs as root by default) we need to instantiate a restricted container to make use of the *Matchhostofowner* application we installed into each app image. To do so we need to run `docker run` with a modified host UID and GID and then run `docker exec` with the `-u` flag that sets up an instance-internal ID called `app`. 




