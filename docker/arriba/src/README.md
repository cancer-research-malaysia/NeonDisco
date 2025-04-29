## NOTES FOR REFERENCE GENOME PREPROCESSING FOR STAR AND ARRIBA

Before running the pipeline, you need to prepare the reference genome for STAR and ARIBA **INSIDE** a `arriba-crmy` container, saving the output in the local host. This is done by running the following command:

``` docker run ... -c bash /tmp/download_refs_for_star.sh [genome+annotation] [-optional argument: path/to/local_ref_genome.gz] [-optional argument: path/to/local_annotation.gz]```

This will download the reference genome and annotation of interest and create the necessary indexes for STAR and ARIBA. A list of available genome+annotation key names can be seen by running the script without any argument. The reference genome will be downloaded from the Ensembl FTP server and the viral sequences will be downloaded from the NCBI FTP server. The STAR index will be saved in the `/home/app/ref/STAR_index_GRCh38viral+ENSEMBL104` directory, so make sure to mount a volume to `/home/app/ref` when running the container.

Note that you should ensure that `RefSeq_viral_genomes_v2.3.0.fa.gz` is present in the script directory. This file is required for the `download_refs_for_star.sh` script to work.

The STAR index should ideally be generated in a `arriba-crmy` container instance with at least 64GB of RAM prior to running the whole Nextflow pipeline. The index generation process can take several hours to complete and would generate a directory of STAR index. 