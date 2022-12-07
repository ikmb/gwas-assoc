# Usage information

## Quick Start
1. Run the [Quality Control Pipeline](https://github.com/ikmb/gwas-qc/blob/master/Readme.md#quick-start) on the example first.
    - All files necessary for the association testing pipeline are automatically generated.
2. Run the gwas-assoc pipeline with like so:
   ```
   nextflow run ikmb/gwas-assoc -r DSL2 \
    --input_imputed_glob "gwas-qc/example/output/Example/QCed/"'*.noATCG.vcf.gz'" \
    --fam "gwas-qc/example/output/Example/SNPQCII/Example_QCed.fam" \
    --collection_name "EXAMPLE" \
    --output "output/Example/Assoc" \
    --build 37 
    ```

## How to Start

You will need to:
- A set of `.vcf.gz` files with the following specifics:
    - at most one chromosome per `.vcf.gz` file. Multiple chromosomes per files are not supported. If a file happens to have multiple chromosomes, only the first will be analyzed.
    - any chromosome codes are supported (i.e. `chrX`, `X`, `23`, `chr23`, `chromosomeX` are just fine)
    - **your VCF files require dosage data (DS tag) for imputed data or genotype calls (GT tag) for genotyped data. If both tags are found, DS is chosen.**
    - the `INFO` column in the VCF files should contain an imputation score. This is used to filter the input variants to create a good null model for SAIGE. For topmed imputations, we found `R2>0.8` to yield good results. If the given tag is not found in the VCF file, the default value 1.0 is assumed, making all variants pass the filter. This filter is not used for Plink-based association testing. 
- A FAM file to update sex and phenotype from. Only those files in the FAM file will be used from the VCFs. Note that *all* samples from the FAM file must be present in the VCF files. **Note: if the FAM family ID is 0, the sample name in the VCF should be the individual ID. If the FID is not 0, the sample name should be the family ID and the individual ID with an underscore (i.e. FID 1234 IID 9876 should be 1234_9876 in the VCF file)**
- The genome build of the input data, which will result in SAIGE handling the sex chromosomes according to the coordinates of the pseudoautosomal regions. Possible values are 37 and 38.
- Optionally, you can specify additional covariates to be used in association testing. By default, 10 principal components are automatically generated and used. If you want additional covariates, have a whitespace-separated file with a header at hand. The first column should be the sample ID in `FID_IID` format, any futher column is treated as a covariate. Specify the covars file with `--more_covars $FILE` and the columns to be used with `--more_covars_cols AGE,SEX`, where `AGE` and `SEX` are the respective column headers from the covar file that you wish to be included.
- The pipeline output and reports will be written to the `--output` directory.

### Parameters

The following list covers all parameters that may be specified for the Association Testing Pipeline:

```
--input_imputed_glob [glob]     [REQUIRED] A glob expression to specify the .vcf.gz files that should be used
                                    for association analysis
--fam [file.fam]                [REQUIRED] A Plink-style FAM file that will be used to select a
                                    subset of samples from the provided VCFs
--build [build]                 [REQUIRED] Genome build, 37 (default) or 38
--trait [type]                  [ADVISED] Trait type to analyze. May be 'binary' (default) or 'quantitative')
--collection_name [name]        [ADVISED] Output filename prefix
--output [directory]            [ADVISED] Output directory
--more_covars [covars.txt]      [OPTIONAL] whitespace-separated list of covariates. See above.
--more_covars_cols [COL1,COL2]  [OPTIONAL] comma-separated list of covar column header names
--null_filter [filter]          [OPTIONAL] bcftools-style formatted INFO filter for generation of
                                    the SAIGE null model. Default: "R2>0.8"
-resume                         [OPTIONAL] Restart where the pipeline was cancelled or aborted. May or may
                                    not work, depending on your filesystem specifics

--phenofile                     [OPTIONAL] Phenotype file for multiple phenotype/traits-testing with regenie. 
                                    Tab separated file with columnsheader "FID IID Phenotype1 Phenotype2" Entries must be "0" for FID, "FID_IID" for IID and all phenotypes must be either binary or quantitaive, don't mix! Missing Samples will be ignored. Binary traits should be specified as control=1,case=2,missing=NA.
--additional_regenie_parameter  [OPTIONAL] Add additional parameters to step2 of regenie e.g. annotation and mask parameters 
                                    for gene-based testing.
```
