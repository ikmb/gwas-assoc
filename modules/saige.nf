
/* Create a SAIGE null model from b38 Plink build and PCs */
process saige_null {
    tag "${params.collection_name}"
    publishDir params.output, mode: 'copy'
	scratch params.scratch
    label 'saige'

    input:
    tuple path(bed), path(bim), path(fam), path(logfile)// from for_saige_null
    tuple path(pheno), path(covar_list)// from for_saige_null_covars
    
    output:
    tuple path("${params.collection_name}.rda"), path("${params.collection_name}.varianceRatio.txt") //into for_saige_assoc, nulldump
   	
shell:
    '''
sed 's/^chr//' !{bim.baseName}.bim >tmp.bim
ln -s !{bed} tmp.bed
ln -s !{fam} tmp.fam
export R_LIBS_USER=/dev/null
if [ "!{params.trait}" == "binary" ]; then
    TRAIT_ARGS="--traitType=binary"
elif [ "!{params.trait}" == "quantitative" ]; then
    TRAIT_ARGS="--traitType=quantitative --invNormalize=TRUE --tauInit=!{params.tauInit}"
else
    echo "Unsupported trait type. Only 'binary' and 'quantitative' traits are supported." >/dev/stderr
    exit 1
fi
step1_fitNULLGLMM.R \
    --plinkFile=tmp \
    --phenoFile=!{pheno} \
    --phenoCol=Pheno \
    --covarColList=$(cat !{covar_list}) \
    --sampleIDColinphenoFile=IID \
    $TRAIT_ARGS        \
    --outputPrefix=!{params.collection_name} \
    --nThreads=!{task.cpus} \
    --LOCO=FALSE
    '''
}

/* Perform SAIGE Assoc test on imputed VCFs */
process saige_assoc {
    scratch params.scratch
    tag "${params.collection_name}.${chrom}.${chunk}"
    label 'saige'

    input:
    // <-- single VCF and chromosome number
//    tuple file(vcf), file(tbi), val(chrom), val(filetype),file(chunk) from for_saige_assoc_imp.transpose() // turn [chr1, [chunk1, chunk2, chunk3]] into [chr1, chunk1], [chr1, chunk2], [chr1, chunk3]
//    tuple file(modelfile), file(varratio) from for_saige_assoc
    tuple file(vcf), file(tbi), file(field), val(chrom), val(filetype), file(chunk), file(modelfile), file(varratio)// from saige_jobs
    each path(inc_fam)

    output:
    file("${chrom}.${chunk.name}.SAIGE.stats")// into for_merge_sumstats
    
shell:
'''
FIRSTPOS=$(head -n1 <!{chunk})
LASTPOS=$(tail -n1 <!{chunk})
#if [ "!{params.chrchr}" = "0" ]; then
    CHR=!{chrom}
#else
#    CHR=chr!{chrom}
#fi
FIELD="DS"
AVAILABLE=$(cat !{field})
# Choose between DS and GT, prefer DS over GT
if [[ ! $AVAILABLE =~ "DS" ]]; then
    if [[ ! $AVAILABLE =~ "GT" ]]; then
        echo "Error: Neither dosage nor genotypes are present in the VCF files." >/dev/stderr
        exit 1
    else
        echo "Warning: Could not find dosage field in VCFs, will use genotype field." >/dev/stderr
        FIELD="GT"
    fi
fi
# Choose PAR region coordinates based on genome build
case !{params.build} in
    19 | 37 | hg19 | hg37 | GRCh37)
        XPAR="60001-2699520,154931044-155260560"
        ;;
    38 | hg38 | GRCh38)
        XPAR="10001-2781479,155701383-156030895"
        ;;
    *)
        XPAR=""
        ;;
esac
# Check whether we need to handle males in the X chromosome
case $CHR in
    23 | X | chr23 | chrX)
        # set up special X handling for SAIGE
        EXTRA_ARGS="--is_rewrite_XnonPAR_forMales=TRUE --X_PARregion=$XPAR --sampleFile_male=males.txt"
        # Make list of males with double-ID
        awk '{if($1!="0") {$2=$1"_"$2; $1="0";}  if($5=="1") {print $2}}' !{inc_fam} >males.txt
        if [ ! -s males.txt ]; then
	   EXTRA_ARGS=""
	fi
        
;;
    *)
        EXTRA_ARGS=""
esac
export R_LIBS_USER=/dev/null
step2_SPAtests.R \
    --vcfFile="!{vcf}" \
    --vcfFileIndex="!{tbi}" \
    --vcfField=$FIELD \
    --chrom=$CHR \
    --minMAF=0.0001 \
    --minMAC=1 \
    --GMMATmodelFile="!{modelfile}" \
    --varianceRatioFile="!{varratio}" \
    --SAIGEOutputFile=temp.stats \
    --numLinesOutput=2 \
    --start=$FIRSTPOS \
    --end=$LASTPOS \
    --IsOutputAFinCaseCtrl=TRUE \
    --IsOutputNinCaseCtrl=TRUE \
    --IsOutputHetHomCountsinCaseCtrl=TRUE \
    $EXTRA_ARGS
# add odds ratio
<temp.stats awk 'NR==1{print $0 " OR";next} {print $0 " " exp($10)}' \
    >"!{chrom}.!{chunk.name}.SAIGE.stats"
'''
}
