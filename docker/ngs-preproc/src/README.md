## NOTES FOR REFERENCE GENOME PREPROCESSING FOR STAR AND ARRIBA

Before running the pipeline, you need to preprocess the reference genome for STAR and ARIBA INSIDE a `ngs-preproc-crm` container. This is done by running the following commands:

``` bash /tmp/download_refs_for_star.sh GRCh38viral+ENSEMBL104 ```

This will download the reference genome and create the necessary indexes for STAR and ARIBA. The reference genome will be downloaded from the Ensembl FTP server and the viral sequences will be downloaded from the NCBI FTP server. The STAR index will be saved in the `/home/app/refs/STAR_index_GRCh38viral+ENSEMBL104` directory.