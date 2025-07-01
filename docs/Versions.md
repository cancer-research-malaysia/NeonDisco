## NeonDisco Version 8.0 alpha

This massive update primarily consists of optimization of the gene fusion discovery module of the base NeonDisco.

### Current Implementations
#### Execution-related
- ✅ Add option to either run NeonDisco with `--inputDir` (for locally stored input files) or `--manifestPath` (for s3 remote files or locally stored files)
- ✅ Add staged file (when using S3-stored objects as inputs) and intermediate file options to minimize storage space consumption during runtime
- ✅ Redesign sample/individual-centric processing up until neopeptide prediction step

#### Logic-related
- ✅ Replace `HLAHD` with `arcasHLA` for HLA typing; this removes the need to use WES or DNA-seq data as inputs, further reducing data type prerequisites for running NeonDisco. As we are interested mostly in population-wide HLA allele set, we can tolerate some degree of HLA typing inaccuracy at sample level. HLAHD might be implemented as optional alternative
- ✅ Implement read filtering step to reduce computational cost during fusion calling. This idea is adopted from Easyfuse pipeline recently published. This is implemented in a v2 workflow (currently named `main.nf`)
- ✅ Implement FusionInspector in silico validation step in NeonDisco
- ✅ Implement STAR-Fusion (using custom build ENSEMBL v113 annotation) to improve on aggregate fusion calling approach
- ✅ Implement pVACfuse neopeptide prediction step in NeonDisco just on FI-validated fusion gene pairs
- ✅ Improve the "detectedBy" column's processing on the normal-filtered fusion transcript list

#### Database-related
- ✅ Incorporate "panel of normals" (PoN) filtering step on the raw aggregate fusion transcript list. Currently our PoN consists of adjacent normal TCGA-BRCA dataset 
- ✅ Incorporate "detectedInCCLE" column on the FILTERED fusion transcript list based on a limited CCLE-based fusion transcripts detected in-house using the same analysis workflow


### Future Improvement
#### Execution–Related
- ✅ Add parallel execution capability of *sample-specific* and *cohort-wide* neopeptide prediction using both sample-specific HLA types and cohort-determined HLA types based on a population frequency cutoff
- ✅ Add parsing capability on manifest files to process inputs based on `sampleType` (Tumor or Normal)
- ✅ Add `S3-manifest-file-generator.sh` utility script (in `bin/utils/`)
- ✅ Implement a subworkflow called `generate-fusion-pons.nf` that can be run prior to main NeonDisco run, to process Normal datasets and add fusions detected in Normals for filtering in the main pipeline
- ☑️ Wrap the `Nextflow run` command into a CLI-based program written in Golang

#### Program Logic–Related 
- ☑️ Implement HLAHD as optional alternative HLA typing module that is selectable by a running parameter
- ☑️ Expand PoN to include Babiceanu et al., Gao et al., and Klijn et al., ~~**GTEx-based**~~ chimeric transcript list (grant proposed – to sequence matched normal of MyBrCa cohort to be added to PoN)
- ~~☑️ Improve false-positive filtering module to also use annotation information in Arriba and FusionCatcher output files (STARFusion already does such filtering internally)~~
