process phenofile_from_fam {
    tag "${params.collection_name}"
	scratch params.scratch
    label 'base'

    input:
        path(assocfam)
    output:
        path('phenotype.txt')
   	
    shell:
    '''
    #gawk  'NR==1  {print "FID\tIID\tPhenotype"}{print "0\t"$1"_"$2"\t"$6}' !{assocfam} > phenotype.txt #-v 'OFS= '  -F '\t'
    gawk  'NR==1  {print "FID\tIID\tPhenotype"}{print $1"\t"$2"\t"$6}' !{assocfam} > phenotype.txt #-v 'OFS= '  -F '\t'
    '''
}

process regenie_step1 {
    tag "${params.collection_name}"
	scratch params.scratch
    label 'regenie'

    input:
        tuple path(bed), path(bim), path(fam), path(logfile)
        tuple path(covars), path(covars_cols)
        path(phenofile)
    output:
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

regenie \
  --step 1 \
  --bed tmp \
  --covarFile !{covars} \
  --covarCol PC{1:!{params.pca_dims}} \
  --phenoFile !{phenofile} \
  --use-relative-path \
  --bsize 100 \
  $TRAIT_ARGS \
  --lowmem \
  --loocv	\
  --lowmem-prefix tmp_rg \
  !{params.additional_regenie_parameter} \
  --out fit_bin_out \
  --gz
    '''
}
//--phenoCol "Phenotype" \

process regenie_step2 {
    scratch params.scratch
    tag "${params.collection_name}"
    label 'regenie'
    publishDir params.output, mode: 'copy'


    input:
        tuple path(bed), path(bim), path(fam), path(logfile)
        tuple path(covars), path(covars_cols)
        tuple path(locofiles), path(predlist)
        path(phenofile)
    output:
        path("${params.collection_name}_regenie_firth*")
    
    shell:
        outprefix = params.collection_name + '_regenie_firth'
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

regenie \
  --step 2 \
  --bed tmp \
  --covarFile !{covars} \
  --covarCol PC{1:!{params.pca_dims}} \
  --phenoFile !{phenofile} \
  --bsize 200 \
  $TRAIT_ARGS \
  --firth --approx \
  --pThresh 0.01 \
  --loocv	\
  --pred !{predlist} \
  --out !{outprefix} \
  !{params.additional_regenie_parameter} \
  --gz
'''
}
