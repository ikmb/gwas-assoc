# IKMB GWAS Association Testing Pipeline

## Prerequisites
- Nextflow: https://www.nextflow.io/
- Java 8 or higher
- Singularity 3.4 or higher
- A dataset in chromosome-wise .vcf.gz format. The dataset *must have* a DS tag for dosage data or a GT tag for genotyped data. If both are present, DS is chosen. 
- Annotations and FAM files as generated by [The QC Pipeline](https://github.com/ikmb/gwas-qc/) (i.e. `dataset_individuals_annotation.txt` and `dataset.fam`)
- Optionally, a whitespace-separated file with additional covariates for association testing (see below)

Please ensure that you have 16 GB RAM installed on the computer where you intend to run the pipeline (i.e. your local computer or your HPC compute nodes).

## Quick Start

1. Run the [Quality Control Pipeline on the example](https://github.com/ikmb/gwas-qc/blob/master/Readme.md#quick-start) first
    - all files necessary for the association testing pipeline are automatically generated
2. Run the wrapper script (included in example package): `bash run-assoc.sh`


## How to Start

→ If you are working on the UKSH medcluster, please see [Advanced Configuration](https://github.com/ikmb/gwas-qc/#uksh-medcluster-configuration) first.

Fortunately, the association testing pipeline itself requires very little configuration. If you intend to run the pipeline on an HPC cluster, please review the advanced configuration items first.

You will need:
- A set of `.vcf.gz` files with the following specifics:
    - at most one chromosome per `.vcf.gz` file. Multiple chromosomes per files are not supported. If a file happens to have multiple chromosomes, only the first will be analyzed.
    - any chromosome codes are supported (i.e. `chrX`, `X`, `23`, `chr23`, `chromosomeX` are just fine)
    - your VCF files require dosage data (DS tag) for imputed data or genotype calls (GT tag) for genotyped data. If both tags are found, DS is chosen.
    - the `INFO` column in the VCF files should contain an imputation score. This is used to filter the input variants to create a good null model for SAIGE. For topmed imputations, we found `R2>0.8` to yield good results. If the given tag is not found in the VCF file, the default value 1.0 is assumed, making all variants pass the filter. This filter is not used for Plink-based association testing. 
- A FAM file to update sex and phenotype from. Only those files in the FAM file will be used from the VCFs. Note that *all* samples from the FAM file must be present in the VCF files. **Note: if the FAM family ID is 0, the sample name in the VCF should be the individual ID. If the FID is not 0, the sample name should be the family ID and the individual ID with an underscore (i.e. FID 1234 IID 9876 should be 1234_9876 in the VCF file)**
- The genome build of the input data, which will result in SAIGE handling the sex chromosomes according to the coordinates of the pseudoautosomal regions. Possible values are 37 and 38.
- Optionally, you can specify additional covariates to be used in association testing. By default, 10 principal components are automatically generated and used. If you want additional covariates, have a whitespace-separated file with a header at hand. The first column should be the sample ID in `FID_IID` format, any futher column is treated as a covariate. Specify the covars file with `--more_covars $FILE` and the columns to be used with `--more_covars_cols AGE,SEX`, where `AGE` and `SEX` are the respective column headers from the covar file that you wish to be included.

The following wrapper script is also included in [the example](https://github.com/ikmb/gwas-qc/blob/master/Readme.md#quick-start)
```bash
# For additional clarity, we use variables here.

# a shell-like glob expression for specifying VCF.GZ file sets.
# Note the additional quoting, we do not want the shell to expand the '*'
INPUT="$HOME/example/output/Example/QCed/"'*.noATCG.vcf.gz'

# Filter by these FAM files. In our example, they are the same as in the VCF
FAM="$HOME/example/output/Example/QCed/1000G_QCed.fam"

# Individuals annotation as generated by the QC Pipeline
ANNO="$HOME/example/output/Example/QCed/1000G_QCed.annotation.txt"

# Name prefix of the output files
NAME="1000G"

# Target folder
OUTPUT="output/Example/Assoc"

# Actual call
# The Assoc.config defines computation resources and can be fine-tuned if necessary. You can
# use the one from the example package or the generic one from 
# https://raw.githubusercontent.com/ikmb/gwas-assoc/master/Assoc.config

if [ ! -f Assoc.config ]; then
    wget https://raw.githubusercontent.com/ikmb/gwas-assoc/master/Assoc.config
fi

nextflow run -c Assoc.config ikmb/gwas-assoc \
    --input_imputed_glob "$INPUT" \
    --fam "$FAM" \
    --collection_name "$NAME" \
    --output "$OUTPUT" \
    --build 37
```

The pipeline output and reports will be written to the ```Assoc_output``` directory.

### Parameters

The following list covers all parameters that may be specified for the Association Testing Pipeline:

```
--input_imputed_glob [glob]     A glob expression to specify the .vcf.gz files that should be used
                                for association analysis
--fam [file.fam]                 A Plink-style FAM file that will be used to select a
                                subset of samples from the provided VCFs
--collection_name [name]        Output filename prefix
--output [directory]            Output directory
--more_covars [covars.txt]      Optional, whitespace-separated list of covariates. See above.
--more_covars_cols [COL1,COL2]  Optional, comma-separated list of covar column header names
--null_filter [filter]          Optional, bcftools-style formatted INFO filter for generation of
                                the SAIGE null model. Default: "R2>0.8"
--build [build]                 Genome build, 37 or 38
--trait [type]                  Trait type to analyze. May be 'binary' (default) or 'quantitative')
-resume                         Restart where the pipeline was cancelled or aborted. May or may
                                not work, depending on your filesystem specifics
```


## Advanced Configuration
Please refer to the [Advanced Configuration section of the QC Pipeline](https://github.com/ikmb/gwas-qc/#advanced-configuration) The same principles apply.
