## NeonDisco Version 6.0 alpha

### Implementations


- ✅ Add option to either run NeonDisco with `--inputDir` (for locally stored input files) or `--manifestPath` (for s3 remote files or locally stored files)
- ✅ Streamline sample/individual-centric processing up until aggregate fusion calling postprocessing
- ✅ Replace `HLAHD` with `arcasHLA` for HLA typing; this removes the need to use WES or DNA-seq data as extra inputs
- ☑️ Incorporate "panel of normals" filtering step on the raw aggregate fusion transcript list
- ☑️ Incorporate "detectedInCCLE" column on the FILTERED fusion transcript list
- ☑️ Add a "detectedBy" "BOTH" column on the FILTERED FT list
- ☑️ Implement parser and executor script for FusionInspector
- ☑️ Upgrade FusionCatcher database used (originally ensembl v102) to v113. v113 is chosen because Fusioncatcher is only able to fetch the latest release so to harmonize references and annotations used with Arriba, we have to use v113 (also Fusioncatcher is using Gencode v47 GTF annotation; full details is in `version.txt` in `fuscat-db` dir in `refs`). Unfortunately AGFusion uses annotation from `pyensembl`, which has not caught up with v113 (v111 is the latest one); so AGFusion is now using v111 only
- ☑️ Implement read filtering step to reduce computational cost during fusion calling. This idea is adopted from Easyfuse pipeline recently published. This is implemented in a v2 workflow (currently named `main.nf`)
- ☑️ Implement STAR-Fusion (using custom build ENSEMBL v113 annotation) to test aggregate fusion calling
- ☑️ 
