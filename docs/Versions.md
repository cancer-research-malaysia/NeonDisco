## NeonDisco Version 6.0 alpha

### Implementations


- ✅ Add option to either run NeonDisco with `--inputDir` (for locally stored input files) or `--manifestPath` (for s3 remote files or locally stored files)
- ✅ Streamline sample/individual-centric processing up until combined Arriba+FusionCatcher output fusion transcripts (UNFILTERED)
- ☑️ Replace `HLAHD` with `arcasHLA` for HLA typing; this removes the need to use WES or DNA-seq data as extra inputs
- ☑️ Incorporate "panel of normals" filtering step on the UNFILTERED Arriba+FusionCatcher output fusion transcript list
- ☑️ Implement parser and executor script for FusionInspector
- ☑️ Upgrade AGFusion database used (originally Ensembl 95) to 103 to match the one used with Arriba (Ensembl 103)