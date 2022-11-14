

process regenie_step1 {
    tag "${params.collection_name}"
    publishDir params.output, mode: 'copy'
	//scratch params.scratch
    scratch false
    label 'regenie'

    input:
    tuple path(bed), path(bim), path(fam), path(logfile)// from for_saige_null
    tuple path(covars), path(covars_cols)// from for_saige_null_covars
    
    output:
    //tuple path("${params.collection_name}.rda"), path("${params.collection_name}.varianceRatio.txt") //into for_saige_assoc, nulldump
    tuple path('fit_bin_out_*.loco*'), path('fit_bin_out_pred.list')
   	
shell:
    '''
sed 's/^chr//' !{bim.baseName}.bim >tmp.bim

ln -s !{bed} tmp.bed
ln -s !{fam} tmp.fam

export R_LIBS_USER=/dev/null
if [ "!{params.trait}" == "binary" ]; then
    TRAIT_ARGS="--bt --cc12"
elif [ "!{params.trait}" == "quantitative" ]; then
    TRAIT_ARGS="--qt --apply-rint"# --tauInit=!{params.tauInit}
else
    echo "Unsupported trait type. Only 'binary' and 'quantitative' traits are supported." >/dev/stderr
    exit 1
fi

awk -F '\t' 'NR==1  {print "FID\tIID\tPhenotype"}{gsub(/"/,""); print $1,$2,$6}' !{fam} > phenotype.txt

regenie \
  --step 1 \
  --bed tmp \
  --covarFile !{covars} \
  --covarCol PC{1:!{params.pca_dims}} \
  --phenoFile phenotype.txt \
  --phenoCol "Phenotype" \
  --bsize 100 \
  $TRAIT_ARGS \
   --lowmem \
  --lowmem-prefix tmp_rg \
  --out fit_bin_out \
  --gz
    '''
}


process regenie_step2 {
    //scratch params.scratch
    scratch false
    tag "${params.collection_name}"
    label 'regenie'

    input:
    // <-- single VCF and chromosome number
//    tuple file(vcf), file(tbi), val(chrom), val(filetype),file(chunk) from for_saige_assoc_imp.transpose() // turn [chr1, [chunk1, chunk2, chunk3]] into [chr1, chunk1], [chr1, chunk2], [chr1, chunk3]
//    tuple file(modelfile), file(varratio) from for_saige_assoc
    
    
    //tuple file(vcf), file(tbi), file(field), val(chrom), val(filetype), file(chunk), file(modelfile), file(varratio)// from saige_jobs
    tuple path(bed), path(bim), path(fam), path(logfile)// from for_saige_null
    
    //each path(inc_fam)

    tuple path(covars), path(covars_cols)// from for_saige_null_covars
    tuple path(locofiles), path(predlist)
    output:
    //file("${chrom}.${chunk.name}.SAIGE.stats")// into for_merge_sumstats
    
shell:
'''
if [ "!{params.trait}" == "binary" ]; then
    TRAIT_ARGS="--bt --cc12"
elif [ "!{params.trait}" == "quantitative" ]; then
    TRAIT_ARGS="--qt --apply-rint"# --tauInit=!{params.tauInit}
else
    echo "Unsupported trait type. Only 'binary' and 'quantitative' traits are supported." >/dev/stderr
    exit 1
fi


sed 's/^chr//' !{bim.baseName}.bim >tmp.bim

ln -s !{bed} tmp.bed
ln -s !{fam} tmp.fam

awk -F '\t' 'NR==1  {print "FID\tIID\tPhenotype"}{gsub(/"/,""); print $1,$2,$6}' !{fam} > phenotype.txt

regenie \
  --step 2 \
  --bed tmp \
  --covarFile !{covars} \
  --covarCol PC{1:!{params.pca_dims}} \
  --phenoFile phenotype.txt \
  --bsize 200 \
  $TRAIT_ARGS \
  --firth --approx \
  --pThresh 0.01 \
  --pred !{predlist} \
  --out test_bin_out_firth \
  --gz
'''
}
