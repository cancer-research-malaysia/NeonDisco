# Building a Minimal Docker Container Image for the Pipeline
## Preface

As most of the bioinformatics tools used in this neoantigen identification workflow are CLI programs, we would like to package them all into one container so the whole pipeline can be shipped around without much installation fuss. As AWS instances support Docker, we decided to use Docker for containerization. 

## The Steps:

Follow the steps below to replicate the containerization processes. Note that some binaries/programs need to be downloaded first prior to these steps.

1. Write a Dockerfile
