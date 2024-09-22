# Fusion Transcript Neoantigen Identification Pipeline (CRM)

## Background

Cancer vaccines are an emerging therapeutic option for tumor diseases with a potential to convert
immune-“cold” tumours to immune-“hot” tumours in combination with immune checkpoint inhibitor immunotherapy. They are also
of particular interest in the quest for the generation of
universal cancer vaccines that can be used “off-the-shelf", in contrast to other highly
personalized approaches such as adoptive cell therapies that are likely to be too resource-
intensive and expensive to be scaled up for general use outside first-world countries.

The hunt for a universal cancer vaccine continues to be a major focus area for cancer
immunotherapy. Various biotarget candidates have been proposed, including cancer neoantigens, tumour-associated antigens, tumour-associated viral
antigens, and more recently, “dark” transient epigenetically driven antigens. One of the main challenges towards developing universal cancer vaccines is that most
neoantigen targets appear to be highly individualized, with little overlap between different
patients. Another main challenge is that acquired resistance to individual targets develops
quickly due to tumour evolution. Thus, a multi-target, highly combinatorial
approach is likely to work best towards finding a generalized set of neoantigens for cancer vaccine development. 

[Nextflow](https://www.nextflow.io/) is a free and open source software project which makes it easier to run a computational workflow consisting of a series of interconnected steps.

This repository documents the porting process of an in-house cancer neoantigen identification workflow at CRM using multiple bioinformatics tools into a Nextflow pipeline to enable portability and reproducibility. The main Nextflow script is written in the [Nextflow DSL2](https://www.nextflow.io/docs/latest/dsl2.html) syntax. The file structure of this repository was soft-cloned from [here](https://github.com/FredHutch/workflow-template-nextflow) but the current file content and scripts have been completely customized and rewritten to fit the goal of this pipeline.

### Repository Structure

The essential components of the workflow repository are as follows:
- `main.nf`: Contains the primary workflow code which pulls in all additional code from the repository
- `modules/`: Contains all of the sub-workflows which are used to organize large chunks of analysis
- `bin/`: Contains all of the code which is executed in each individual step of the workflow

### User Input of Parameters

There are three ways by which users can easily provide their own inputs to a workflow; (1) with command-line flags, (2) with a params file and (3) by setting default values in `nextflow.config`.

On the command line, parameters are provided using two dashes before the parameter name, e.g. `--param_name value`. One limitation of this approach is that the provided value will be interpreted as a string. The best example of this is the edge case of the the negative boolean (`false`), which will be interpreted by Nextflow as a string (`'false'`). The second limitation is that the command line string starts to become rather long. Another consideration of providing parameters on the command line is that they may be interpreted by the shell before execution. For example, in the context of a BASH script `--param_name *.fastq.gz` will first be expanded into a list of files which match that pattern (e.g., `--param_name 1.fastq.gz 2.fastq.gz 3.fastq.gz`), which may not be the intention. This behavior can be prevented explicitly with single-quotes in BASH, with `--param_name '*.fastq.gz'` being unaltered by the shell before execution.

By using a params file, the user is able to more explicitly define the set of parameters which will be provided. The params file can be formatted as JSON or YAML, with the example below shown in JSON.

```
{
    "param_name": "*.fastq.gz",
    "second_param": false,
    "third_param": 5
}
```

The params file is provided by the user with the `-params-file` flag. While this approach requires the user to create an additional file, it also provides a method for defining variables without worrying about the nuances of the shell interpreter. If both methods are used for providing parameters, the command line flags will take precedence over the params file ([docs](https://www.nextflow.io/docs/latest/config.html)).

The third way to pass parameters to the workflow is to set up the default values in `nextflow.config` in the `params` scope (e.g. `params{param_name = 'default_value'}`). If a user passes in a value on the command line, then the configured default `params.param_name` will be overridden. The really useful thing about `params` is that they are inherited by every sub-workflow and process that is invoked. In other words, without having to do _anything_ else, you can use `${params.param_name}` in one of the script files in `bin/`, and you know that it will contain the value that was initially set.

### Templates

One of the options for defining the code that is run inside a Nextflow process is to use their [template syntax](https://www.nextflow.io/docs/latest/process.html#template). The advantage of this approach is that the code can be defined in a separate file with the appropriate file extension which can be recognized by your favorite IDE and linter. Any variables from Nextflow will be interpolated using an easy `${var_name}` syntax, and all other code will be native to the desired language. 

One important caveat to note while writing the template structure is that backslashes are used to escape Nextflow interpolation (e.g. internal BASH variables can be specified with `\$INTERNAL_VAR_NAME`), thus any use of backslashes for special characters must have two backslashes. Put simply, if you want to strip the newline character in Python, you would need to write `str.strip('\\n')` instead of `str.strip('\n')`.

### Software Containers

Each individual step in a workflow should be run inside a container (using either Docker or Singularity) which has the required dependencies. For this project, I have manually custom-built Docker container images of the complete software stack needed to run the prediction pipeline. This adds modularity to the stack so users can choose to run only the necessary tools in the stack if needs be. 

Software containers should be defined as parameters in `main.nf`, which allows the value to propagate automatically to all imported sub-workflows, while also being able to be overridden easily by the user if necessary. Practically speaking, this means that every process should have a `container` declared which follows the pattern `container "${params.container__toolname}"`, and which was set in `nextflow.config` with `params{container__toolname = "quay.io/org/image:tag"}`. It is crucial that the parameter be set _before_ the subworkflows are imported, as shown in this example workflow.

------
### Workflow Style Guide

While a workflow could be made in almost any way imaginable, there are some tips and tricks which make debugging and development easier. This is a highly opinionated list, and should not be taken as a hard-and-fast rule.

- Never use file names to encode metadata (like specimen name, `.trimmed`, etc.)
- Always publish files with `mode: 'copy', overwrite: true`
- Use `.toSortedList()` instead of `.collect()` for reproducible ordering
- Add `set -Eeuo pipefail` to the header of any BASH script
- Every process uses a `container`, which is defined as a `param.container__toolname` in `main.nf`
- Never use `.baseName` to remove file extension, instead use (e.g.) `.name.replaceAll('.fastq.gz', '')`

--------------------

## Overview of the Neoantigen Identification Workflow

The flowchart belows describes the workflow underlying the Nextflow pipeline. Each step has been individually tested with local cohort data.

```mermaid
flowchart TD
    subgraph "`*Fusion transcript calling*`"
    A[RNA-seq raw reads] -->|Arriba| B(Fusion transcripts)
    A[RNA-seq raw reads] -->|FusionCatcher| B(Candidate fusion transcripts)
    end
    subgraph "`*HLA genotyping*`"
    A2[WES aligned reads] --> |HLA-HD| B2(HLA allele haplotypes)
    end
    B2 --> |shared allele filtering| C(Top 5% HLA alleles in MyBRCA)
    B --> |AGFusion| C2(Annotated fusion transcripts)

    subgraph "`*neoantigen prediction*`"
    C ---> D(Common HLA alleles in MyBRCA + cancer-specific coding fusion transcripts)
    C2 ---> D
    D ---> |pVacFuse| E(Predicted immunogenic fusion neoantigens)
    end
    subgraph "`*neoepitope shortlisting*`"
    E ---> F(Fusion transcripts)
    F ---> |median IC50 <500nM| G(Fusion neoepitopes with coding potential)
    G ---> H1(Sample occurrence = 1)
    G ---> H2(Sample occurrence >= 2)
    H2 ---> I1(Reference proteome match)
    H2 --->I2(Novel peptide)
    end

```
