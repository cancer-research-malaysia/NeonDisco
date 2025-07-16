## NOTES FOR REFERENCE GENOME PREPROCESSING FOR STAR AND ARRIBA

Before running the pipeline, you need to prepare the reference genome for STAR and ARIBA **INSIDE** a `arriba-crmy` container, saving the output in the local host. This is done by running the following command:

```
docker run -e "MHF_HOST_UID=$(id -u)" -e "MHF_HOST_GID=$(id -g)" -v /home/ec2-user/repos/NeonDisco/docker/arriba/src:/scripts -v /home/ec2-user/refs/star-db/GRCh38viral_ENSEMBL113:/home/app/refs sufyazi/arriba-crmy bash -c "/scripts/download_refs_for_star.sh GRCh38viral+ENSEMBL113 /home/app/refs/local_sources/Homo_sapiens.GRCh38.113.dna.primary_assembly.fa.gz /home/app/refs/local_sources/Homo_sapiens.GRCh38.113.chr.gtf.gz"

```
**NOTE** The `-v` parameters specify binding of two locations: where the `download_refs_for_star.sh` is located (together with the local RefSeq viral genomes fa.gz) and where you want the output of this script is saved to. Make sure the virtual path for the second binding be `/home/app/refs` to ensure the script can run.

This will download the reference genome and annotation of interest and create the necessary indexes for STAR and ARIBA. A list of available genome+annotation key names can be seen by running the script without any argument. The reference genome will be downloaded from the Ensembl FTP server and the viral sequences will be downloaded from the NCBI FTP server. The STAR index will be saved in the `/home/app/refs/STAR_index_XXXXXXX` directory, so make sure to mount a volume to `/home/app/refs` when running the container.

Note that you should ensure that `RefSeq_viral_genomes_v2.3.0.fa.gz` is present in the script directory. This file is required for the `download_refs_for_star.sh` script to work.

The STAR index should ideally be generated in a `arriba-crmy` container instance with at least 64GB of RAM prior to running the whole Nextflow pipeline. The index generation process can take several hours to complete and would generate a directory of STAR index. 