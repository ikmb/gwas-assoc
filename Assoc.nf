
params.tauInit="1,0"
params.collection_name = "collection_name"
params.output = "output"
params.more_covars = "."
params.more_covars_cols = ""
params.liftover = 0
params.trait = "binary"
params.chrchr = 0 // should not be used anymore

params.null_filter = "R2>0.8"

params.ucsc_liftover = "" // points to base path where chain files and liftOver binary are found

inc_fam = file(params.fam)

def get_file_details(filename) {
    def m = filename =~ /\/([^\/]+)(\.vcf\.gz|\.bgen)$/
    return [ m[0][1], m[0][2] ]
}

def get_chromosome_code(filename) {
    def m = filename =~ /\/([^\/]+).filtered.vcf.gz$/
    return m[0][1]
}

for_saige_imp = Channel.fromPath(params.input_imputed_glob, checkIfExists: true).map { it ->
    def match = get_file_details(it)
    [it, match[1]] }

for_saige_null_covars = Channel.create()
for_merge_sumstats_covars = Channel.create()


for_split_vcf = Channel.create()
for_extract_dosage = Channel.create()
for_gen_r2 = Channel.create()
for_plink_imp = Channel.create()


process option_check {
    errorStrategy 'terminate'
    input:
    output:

shell:
'''
#!/bin/bash
ERROR=0

# FAM file
if [ ! -f "!{params.fam}" ]; then
    echo "Cannot open FAM file. Please fix the filename given to --fam. You also might want to check the documentation about mounting external paths:" >/dev/stderr
    echo "   https://github.com/ikmb/gwas-qc/#mounting-paths-into-the-singularity-container" >/dev/stderr
    ERROR=1
fi

# Covars
if [ "!{params.more_covars}" != "." ] && [ ! -f "!{params.more_covars}" ]; then
    echo "Cannot open covariate file.  Please fix the filename given to --fam. You also might want to check the documentation about mounting external paths:" >/dev/stderr
    echo "   https://github.com/ikmb/gwas-qc/#mounting-paths-into-the-singularity-container" >/dev/stderr
    ERROR=1
fi

if [ "!{params.more_covars}" != "." ] && [ -f "!{params.more_covars}" ] && [ "!{params.more_covars_cols}" == "" ]; then
    echo "A covariate file has been specified with --more_covars but no covar columns were given using --more_covars_cols." >/dev/stderr
    ERROR=1
fi


# Genome build
case "!{params.build}" in
    37) ;;
    38) ;;
    *) echo "Unsupported genome build. Please specify --build 37 or --build 38." >/dev/stderr
       ERROR=1
       ;;
esac

# Trait type
case "!{params.trait}" in
    binary) ;;
    quantitative) ;;
    *) echo "Unsupported trait type. Please specify 'binary' or 'quantitative' for --trait. " >/dev/stderr
       ERROR=1
       ;;
esac


exit $ERROR
'''
}

// prefilter
// ---------
// Prepares the imputed files for futher analysis:
// - generates a tabix index
// - filters samples according to the FAM file given via --fam.
// Creates new VCF/GZ sets.
process prefilter {
    tag "${params.collection_name}"
    label "big_mem"
    label "long_running"
    cpus 4
    input:
    tuple file(vcf), val(filetype) from for_saige_imp
    file inc_fam

    output:
    tuple file("*.filtered.vcf.gz"), file("*.filtered.vcf.gz.tbi"), val(filetype) into for_extract_chromosome_code //for_split_vcf, for_extract_dosage, for_gen_r2, for_plink_imp

shell:
'''
tabix -p vcf !{vcf}

# Gather lists of samples with proper double-ID handling
bcftools query -l !{vcf} | sort >vcf-samples
gawk '{if($1!="0") {$2=$1"_"$2; $1="0";} print $2}' !{inc_fam} | sort >fam-samples

# Keep only those samples that are in VCF and in FAM file
comm -12 vcf-samples fam-samples >keep-samples

# Count those that we ought to remove
REMOVE_COUNT=$(comm -23 vcf-samples fam-samples | wc -l)

bcftools query -f '%CHROM\\n' !{vcf} | head -n1 >chrom
CHROM=$(cat chrom)

# Ok, we need to remove some samples. tabix afterwards
# If there's nothing to remove, skip this time-consuming step.
if [ "$REMOVE_COUNT" != "0" ]; then
    bcftools view --threads 4 -S keep-samples !{vcf} -o $CHROM.filtered.vcf.gz -O z
    tabix $CHROM.filtered.vcf.gz
else
    ln -s !{vcf} $CHROM.filtered.vcf.gz
    ln -s !{vcf}.tbi $CHROM.filtered.vcf.gz.tbi
fi

'''
}

for_extract_chromosome_code.map { it -> [it[0], it[1], get_chromosome_code(it[0]), it[2]] }.into { for_split_vcf; for_extract_dosage; for_gen_r2; for_plink_imp }

/* Prepare a list of r2-filtered variants for SAIGE step1 model generation */
process gen_r2_list {
    tag "${params.collection_name}.${chrom}"
    input:
    tuple file(vcf), file(tbi), val(chrom), val(filetype) from for_gen_r2

    output:
    file("r2-include.${chrom}") into for_prune_r2

shell:
'''
set +e
bcftools query -i'!{params.null_filter}' -f '%CHROM:%POS:%REF:%ALT\\n' !{vcf} >r2-include.!{chrom}

if [ $? -ne 0 ]; then
    echo Filter not found or genotyped-only data available.
    bcftools query -f '%CHROM:%POS:%REF:%ALT\\n' !{vcf} >r2-include.!{chrom}
fi

exit 0

'''
}

process merge_r2 {
    tag "${params.collection_name}"
    input:
    file r2 from for_prune_r2.collect()

    output:
    file("r2-include.sorted") into for_prune_merged

shell:
'''
cat r2-include.* >r2-include
export TMPDIR=.
gawk '$0 !~ /^chr/ {$1="chr"$1} {print}' r2-include | sort >r2-include.sorted
'''
}

/* Convert chromosome-wise imputed VCFs to Plink bim/bed/fam sets, also fills
   phenotype and sex columns from QC FAM file. */
process make_plink {
//    publishDir params.output, mode: 'copy'
    label "big_mem"
    tag "${params.collection_name}.$chrom"

    input:
    tuple file(vcf), file(tbi), val(chrom), val(filetype) from for_plink_imp
    file fam from inc_fam

    output:
    tuple file("${params.collection_name}.${chrom}.bed"), file("${params.collection_name}.${chrom}.bim"), file("${params.collection_name}.${chrom}.fam"), file("${params.collection_name}.${chrom}.log"),val(chrom) into for_liftover, for_merge_b38

shell:
'''
# Generate double-id FAM
MEM=!{task.memory.toMega()-1000}
gawk '{if($1!="0") {$2=$1"_"$2; $1="0";} print $0}' !{fam} >new-fam
/opt/plink2 --vcf !{vcf} --const-fid --memory $MEM --allow-no-sex --pheno new-fam --mpheno 4 --update-sex new-fam 3 --output-chr chrM --make-bed --keep-allele-order --out !{params.collection_name}.!{chrom}

mv !{params.collection_name}.!{chrom}.bim old_bim
gawk '$1 !~ /^chr/ {$1="chr"$1} {$2=$1":"$4":"$6":"$5; print}' <old_bim >!{params.collection_name}.!{chrom}.bim

# Might need some "chr" prefixing here
'''
}

process merge_plink {
    tag "${params.collection_name}"
    publishDir params.output, mode: 'copy'
    label 'big_mem'
    cpus 4
    label 'long_running'
    input:
    file(filelist) from for_merge_b38.collect()

    output:
    tuple file("${params.collection_name}.bed"), file("${params.collection_name}.bim"), file("${params.collection_name}.fam"), file("${params.collection_name}.log") into for_prune //, for_liftover


    shell:
'''
ls *.bed | xargs -i -- basename {} .bed | tail -n +2 >merge-list
MEM=!{task.memory.toMega()-1000}
FIRSTNAME=$(ls *.bed | xargs -i -- basename {} .bed | head -n1)
module load Plink/1.9
plink --bfile $FIRSTNAME --threads 4 --memory $MEM --merge-list merge-list --make-bed --allow-no-sex --indiv-sort none --keep-allele-order --out !{params.collection_name}
#mv !{params.collection_name}.bim tmp
#mawk '{$2=$1":"$4":"$6":"$5; print $0}' tmp >!{params.collection_name}.bim
'''
}


process prune {
    tag "${params.collection_name}"
    label "big_mem"
    publishDir params.output, mode: 'copy'
    cpus 4
    input:
    tuple file(bed), file(bim), file(fam), file(logfile) from for_prune
    file r2 from for_prune_merged

    output:
    tuple file("${params.collection_name}.pruned.bed"), file("${params.collection_name}.pruned.bim"), file("${params.collection_name}.pruned.fam"), file("${params.collection_name}.pruned.log") into for_saige_null, for_liftover_pruned, for_generate_pcs

shell:
'''
module load Plink/1.9
echo Generating PCA SNP List file for variant selection

R2=!{r2}
MEM=!{task.memory.toMega()-1000}

if [ ! -s $R2 ]; then
    # No entries in include list. This means that either we have genotyped-only
    # data where no R2 score is present or the filter is malformed. We're assuming
    # genotyped-only data here. In this case, populate the r2 list.
    cut -f2 -d" " !{bim} | sort >new-r2
    R2=new-r2
fi


/opt/plink2 --memory $MEM --threads 4 --bed !{bed} --bim !{bim} --fam !{fam} --extract $R2 --indep-pairwise 50 5 0.2 --out _prune
/opt/plink2 --memory $MEM --threads 4 --bed !{bed} --bim !{bim} --fam !{fam} --extract _prune.prune.in --maf 0.05 --make-bed --out intermediate

python -c 'from SampleQCI_helpers import *; write_snps_autosomes_noLDRegions_noATandGC_noIndels("!{bim}", "include_variants")'
sort include_variants >include_variants.sorted

comm -12 include_variants.sorted $R2 >include-r2-variants
/opt/plink2 --memory $MEM --threads 4 --bfile intermediate --chr 1-22 --extract include-r2-variants --output-chr chrM --make-bed --out !{params.collection_name}.pruned
rm intermediate*

'''
}

process generate_pcs {
    cpus 16
    tag "${params.collection_name}"
    label "big_mem"

input:
    tuple file(bed), file(bim), file(fam), file(log) from for_generate_pcs

output:
    file "pcs.txt" into for_saige_covars

shell:
'''
module load FlashPCA
MEM=!{task.memory.toMega()-1000}

flashpca2 -d 10 --bfile !{bim.baseName} --memory $MEM --numthreads 16 --outload loadings.txt --outmeansd meansd.txt

'''
}


/* Create a SAIGE null model from b38 Plink build and PCs */
process saige_null {
    tag "${params.collection_name}"
    publishDir params.output, mode: 'copy'
    container "docker://wzhou88/saige:0.42.1"
    cpus 16
    label 'long_running'

    input:
    tuple file(bed), file(bim), file(fam), file(logfile) from for_saige_null
    tuple file(pheno), file(covar_list) from for_saige_null_covars

    output:
    tuple file("${params.collection_name}.rda"), file("${params.collection_name}.varianceRatio.txt") into for_saige_assoc, nulldump

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
    --nThreads=16 \
    --LOCO=FALSE
    '''
}


//params.saige_chunk_size = 50
params.saige_chunk_size = 20000

process split_vcf_ranges {
    tag "${params.collection_name}.${chrom}"
    input:
    tuple file(vcf), file(tbi), val(chrom), val(filetype) from for_split_vcf
    output:
    tuple file(vcf), file(tbi), file(field), val(chrom), val(filetype), file("0*") into for_saige_assoc_imp, testdump

    shell:
    '''
    bcftools query -f "%POS\\n" !{vcf} >positions
    bcftools view -H !{vcf} | head -n1 | cut -f9 >field
    split -d --suffix-length=8 --lines=!{params.saige_chunk_size} positions '0'
    '''
}



testdump.transpose().combine(nulldump.toList()).map{ it -> [ it[0], it[1], it[2], it[3], it[4], it[5], it[6][0], it [6][1] ] }.set{ saige_jobs }

/* Perform SAIGE Assoc test on imputed VCFs */
process saige_assoc {
    tag "${params.collection_name}.${chrom}.${chunk}"
    container "docker://wzhou88/saige:0.42.1"
    label 'long_running'
    input:
    // <-- single VCF and chromosome number
//    tuple file(vcf), file(tbi), val(chrom), val(filetype),file(chunk) from for_saige_assoc_imp.transpose() // turn [chr1, [chunk1, chunk2, chunk3]] into [chr1, chunk1], [chr1, chunk2], [chr1, chunk3]
//    tuple file(modelfile), file(varratio) from for_saige_assoc
    tuple file(vcf), file(tbi), file(field), val(chrom), val(filetype), file(chunk), file(modelfile), file(varratio) from saige_jobs
    file inc_fam

    output:
    file("${chrom}.${chunk.name}.SAIGE.stats") into for_merge_sumstats

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


//for_merge_sumstats.collect().dump()

process merge_saige_results {
    tag "${params.collection_name}"
    publishDir params.output, mode: 'copy'
    input:
    file(sumstats) from for_merge_sumstats.collect()
    tuple file(pheno), file(cols) from for_merge_sumstats_covars
//    file in_covars


    output:
    file("${params.collection_name}.SAIGE.stats") into for_lift_sumstats_saige
shell:
'''

ls -1 *.*.stats | sort -n >allfiles


head -n1 !{sumstats[0]} | tr -s '\t ' ' '>tmp
while read -r line; do
    tail -n +2 $line >>tmp
done <allfiles


<tmp gawk 'NR==1{print} NR>1{if(substr($1,1,3)!="chr"){$1="chr"$1} $3=$1":"$2":"$4":"$5; print}' >>!{params.collection_name}.SAIGE.stats

'''
}


process extract_dosage {
    tag "${params.collection_name}.${chrom}"
    label 'long_running'

    input:
    tuple file (gz), file(tbi), val(chrom), val(filetype) from for_extract_dosage

    output:
    tuple file("${chrom}.PLINKdosage.map"), file("${chrom}.PLINKdosage.gz"), val(chrom) into for_plink

shell:
'''

VCF=!{gz}
TARGET=!{chrom}.PLINKdosage

# Theoratically, gawk or awk would also work but are much slower
AWK=gawk
BGZIP=bgzip

FIFO=$TARGET.fifo
rm -f $FIFO
mkfifo $FIFO

# Generate awk script
cat <<'AwkProg2' >awkprog2.awk
BEGIN{
    printf "CHR BP SNP A1 A2 INFO"
    gpfield = 0
    gtfield = 0
}
{
    # skip header lines
    if (($1 ~ /^##/)) {
    } else {
        # Copy header line
        if (($1 ~ /^#CHROM/)) {
            # Process one FID_IID-part into "0 FID_IID", so Plink thinks there
            # are really two values
            for (i=10; i<=NF; i++)
                printf " "0" "$i
            printf "\\n"
        } else {
            # Find field index of GP and GT
            n = split($9, fields, ":")
            for(i=1; i in fields; i++) {
                if(fields[i]=="GP") {
                    gpfield = i
                }
                if(fields[i]=="GT") {
                    gtfield = i
                }
            }

            if(gtfield == 0 && gpfield == 0) {
                print "Error: Cannot find a GT or GP field in input VCF." >"/dev/stderr"
                exit 1
            }

            # Extract r^2 score from INFO field, set INFO to 1.0 if none found
            {
                n = split($8, fields, /(;|=)/)
                info = "R2=1.0" # default to 1 if no usable info score is found
                for(i=1; i in fields; i++) {
                    if(fields[i] == "INFO" || fields[i] == "DR2" || fields[i] == "R2") {
                        info = fields[i+1]
                        break
                    }
                }
            }

            # keep a dictionary of variants to skip duplicates
            if(!schonmalgesehen[$1":"$2":"$4":"$5]) {
                printf $1" "$2" "$1":"$2":"$4":"$5" "$4" "$5" "info

                # split value field
                for (i=10; i<=NF; i++) {
                    n=split($i,array,":")
                    # Choose GP over GT
                    if(gpfield == 0) {
                        n=split($i,array,/[\\/|\\|]/)
                        if(n==2) {
                            # "convert" hard calls to probabilities
                            if(array[1]+array[2]==0) { printf " 1 0 0" }
                            else if(array[1]+array[2]==1) { printf " 0 1 0" }
                            else { printf " 0 0 1" }
                        
                        } else if(n==1) {
                            # male samples on chrX, only one genotype
                            printf " "$i " 0 "$i
                        } else {
                            print "VCF Error: GP or GT field not present or malformatted" >"/dev/stderr"
                            exit 1
                        }
                    } else {

                        # Missing fields are specified as ".", not splittable by ","
                        if(array[gpfield] == ".") {
                            printf " . . ."
                        } else {
                            n=split(array[gpfield],array_GP,",")
                            if(n==3) {
                                printf " "array_GP[1]" "array_GP[2]" "array_GP[3]
                            } else {
                                # males on chrX
                                printf " "array_GP[1]" 0 "array_GP[2]
                            }

                        }
                    }
                }

                printf "\\n"

                schonmalgesehen[$1":"$2":"$4":"$5] = 1
            }
        }

    }
}
AwkProg2

(<$FIFO tail -n +2 | $AWK '{print $1, $3, 0, $2}' | uniq) >$TARGET.map &

# unpack vcf.gz only once
$BGZIP -c -d $VCF \
    | $AWK -f awkprog2.awk \
    | tee $FIFO \
    | gzip >$TARGET.gz

wait

rm -f $FIFO
'''
}



/* Create SAIGE-compatible covars from PCA eigenvectors */
process make_saige_covars {

    publishDir params.output, mode: 'copy'

    input:
    // file in_covars
    file (evec) from for_saige_covars

    file inc_fam

    output:
    tuple file("${params.collection_name}.covars"), file("${params.collection_name}.covar_cols") into for_saige_null_covars, for_merge_sumstats_covars, for_plink_covars
    file "${params.collection_name}.double-id.fam" into for_plink_fam

    shell:
'''
# Re-format sample ID and family ID to VCF rules
gawk '{if($1!="0") {$2=$1"_"$2; $1="0";} print $0}' !{inc_fam}  | sort >new-fam

# Also re-format, but keep evec header
# mawk 'NR==1{print $0;next} {if($1!="0") {$2=$1"_"$2; $1="0";} print $0}' {in_covars}  | sort >evec

# Take evec file as a whole, SAIGE does not care about additional columns.
# Filter evec file according to samples contained in FAM (if not already done)

gawk 'FNR==NR{samples[$2];next} {if($2 in samples) { for(i=1;i<=(12);i++) {printf "%s%s", $i, (i<12?OFS:ORS)}}}' new-fam !{evec} >filtered-evec

if [ -f "!{params.more_covars}" ]; then
    gawk 'FNR==NR{samples[$2];next} {if($2 in samples) { for(i=1;i<=(NF);i++) {printf "%s%s", $i, (i<NF?OFS:ORS)}}}' new-fam !{params.more_covars} | sort >filtered-covars
    cut -f3- !{params.more_covars} >covars-column
    < filtered-covars cut -f3- -d" " >>covars-column
fi

EVEC_LINES=$(wc -l <filtered-evec)
FAM_LINES=$(wc -l <new-fam)
if [ $EVEC_LINES -ne $FAM_LINES ]; then
    echo I could not find covariates in !{evec} for every sample specified in !{inc_fam}. Please check.
    echo Aborting.
    exit 1
fi

echo "FID IID PC1 PC2 PC3 PC4 PC5 PC6 PC7 PC8 PC9 PC10" >evec.double-id.withheader
cat filtered-evec >>evec.double-id.withheader

# Take phenotype info from FAM, translate to SAIGE-encoding
echo "Pheno" >pheno-column
<new-fam  gawk '{if($6=="2") {$6="1";} else if($6=="1") {$6="0";} print $6}' >>pheno-column

mv new-fam !{params.collection_name}.double-id.fam

# Merge both, replace space runs with single tabs for SAIGE
touch covars-column

if [ ! -d "!{params.more_covars}" ]; then
    echo PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10,!{params.more_covars_cols} >!{params.collection_name}.covar_cols
else
    echo PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10 >!{params.collection_name}.covar_cols
fi

paste -d" " evec.double-id.withheader pheno-column covars-column | tr -s ' ' \\\\t >!{params.collection_name}.covars

'''
}


process plink_assoc {
    tag "${params.collection_name}.${chrom}"
    label 'long_running'
    label 'big_mem'

    input:
    tuple file(dosagemap), file(dosage), val(chrom) from for_plink
    tuple file(pheno), file(cols) from for_plink_covars
    file fam from for_plink_fam

    output:
    file "${chrom}.assoc.dosage" into for_plink_results

shell:
'''
module load Plink/1.9
MEM=!{task.memory.toMega()-1000}

plink --fam !{fam} --map !{dosagemap} --dosage !{dosage} skip0=2 skip1=0 skip2=1 format=3 case-control-freqs \\
    --covar !{pheno} --covar-name $(cat !{cols}) \\
    --allow-no-sex --ci 0.95 \\
    --out !{chrom} --memory $MEM

'''
}

process merge_plink_results {
    tag "${params.collection_name}"

    publishDir params.output, mode: 'copy'
    input:
    file(stats) from for_plink_results.collect()

    output:
    file("${params.collection_name}.Plink.stats") into for_lift_sumstats_plink

shell:
'''
# extract first line, convert tabs to space
head -n1 !{stats[0]} | tr -s '\t ' ' ' | xargs >!{params.collection_name}.Plink.stats

ls !{stats} | sort -n | xargs -n1 tail -n +2 | gawk '{if(substr($1,1,3)!="chr"){$1="chr"$1} $2=$1":"$3":"$4":"$5; print}'>>!{params.collection_name}.Plink.stats
'''
}

// ----------------------------------------------------------------------------

process liftover {
    tag "${params.collection_name}.${chrom}"
    cpus 2
    label 'big_mem'
//    publishDir params.output, mode: 'copy'
    input:
    tuple file(bed), file(bim), file(fam), file(logfile),val(chrom) from for_liftover
    output:
    tuple file("${bed.baseName}_lifted.bed"), file("${bim.baseName}_lifted.bim"), file("${fam.baseName}_lifted.fam"), file("${logfile.baseName}_lifted.log") optional true into for_merge_lifted
    file ("postlift.${chrom}") optional true into for_merge_tables
//    tuple file ("postlift.${chrom}"), file("unmapped-variants"), file("duplicates") into for_merge_lift
//    tuple file("postlift.${chrom}"), file("unmapped-variants") into for_gen_liftdb
shell:
'''
MEM=!{task.memory.toMega()-1000}

if [ "!{params.build}" == "37" ]; then
    SUFFIX="hg38to37"
    TEXT="hg38 to hg37"
    NEWBUILD="38"
else
    SUFFIX="hg37to38"
    TEXT="hg37 to hg38"
    NEWBUILD="37"
fi

LIFTOVER=!{params.ucsc_liftover}/liftOver
BASENAME=!{bed.getBaseName()}_lifted

if [ "!{params.build}" == "37" ]; then
    CHAIN=!{params.ucsc_liftover}/hg19ToHg38.over.chain.gz
elif [ "!{params.build}" == "38" ]; then
    CHAIN=!{params.ucsc_liftover}/hg38ToHg19.over.chain.gz
else
    echo "Error: Input must be build 37 or build 38 but is specified as !{params.build}" >&2
    exit 0
fi

if [ ! -x "$LIFTOVER" ]; then
    echo "Error: Cannot execute ${LIFTOVER}. File does not exist or is not executable." >&2
    exit 0
fi



## Translate to chr:pos-pos and chr23 -> chrX
#mawk '{print "chr"$1":"$4"-"$4}' !{bim} \\
#    | sed 's/chr23/chrX/g' \\
#    | sed 's/chr24/chrY/g' \\
#    | sed 's/chr25/chrX/g' \\
#    | sed 's/chr26/chrMT/g' \\
#    > prelift.pos

# Create UCSC BED file from BIM
# Columns: chrX pos-1 pos rsID
<!{bim} gawk '{ print $1, $4-1, $4, $4, $2 }' \\
    | sed 's/^chr23/chrX/' \\
    | sed 's/^chr24/chrY/' \\
    | sed 's/^chr25/chrX/' \\
    | sed 's/^chr26/chrMT/' >prelift.bed

# lift-over
$LIFTOVER prelift.bed $CHAIN postlift.bed unmapped.bed

# Generate exclude list for unmapped variants
<unmapped.bed gawk '$0~/^#/{next} {print $5}' >unmapped-variants &

# Generate update-pos list for Plink, exclude double variants
<postlift.bed gawk '{print $5,$3}' >new-pos &

wait

<postlift.bed gawk '{if($1 != "chr!{chrom}" && $1 != "!{chrom}") {print $5}}' >chromosome-switchers &
<new-pos sort | uniq -d | cut -f1 -d" ">duplicates &

wait

cat chromosome-switchers >>duplicates

/opt/plink2 --bed !{bed} --bim !{bim} --fam !{fam} --memory $MEM --exclude duplicates --make-pgen --out no-duplicates
/opt/plink2 --pfile no-duplicates --memory $MEM --exclude unmapped-variants --update-map new-pos --make-pgen --sort-vars  --out ${BASENAME}
/opt/plink2 --pfile ${BASENAME} --memory $MEM --make-bed --out ${BASENAME}.tmp

plink --bfile ${BASENAME}.tmp --memory $MEM --merge-x no-fail --make-bed --out ${BASENAME}
mv ${BASENAME}.bim tmp
gawk '{$2="chr"$1":"$4":"$6":"$5; print $0}' tmp >${BASENAME}.bim
mv postlift.bed postlift.!{chrom}
'''
}

// Merge b37 Plink chromosome-wise sets into a single Plink set
process merge_plink_lifted {
    tag "${params.collection_name}"
    publishDir params.output, mode: 'copy'
    cpus 4
    label 'big_mem'
    label 'long_running'

    input:
    file(filelist) from for_merge_lifted.collect()

    output:
    tuple file("${params.collection_name}_lifted*.bed"), file("${params.collection_name}_lifted*.bim"), file("${params.collection_name}_lifted*.fam"), file("${params.collection_name}_lifted*.log")

    shell:
'''
MEM=!{task.memory.toMega()-1000}


if [ "!{params.build}" == "37" ]; then
    SUFFIX="hg38to37"
    TEXT="hg38 to hg37"
    NEWBUILD="38"
else
    SUFFIX="hg37to38"
    TEXT="hg37 to hg38"
    NEWBUILD="37"
fi

ls *.bed | xargs -i -- basename {} .bed | tail -n +2 >merge-list
FIRSTNAME=$(ls *.bed | xargs -i -- basename {} .bed | head -n1)
module load Plink/1.9
plink --bfile $FIRSTNAME --memory $MEM --threads 4 --merge-list merge-list --make-bed --allow-no-sex --indiv-sort none --keep-allele-order --out !{params.collection_name}_lifted_b$NEWBUILD
mv !{params.collection_name}_lifted_b${NEWBUILD}.bim tmp
gawk '{$2="chr"$1":"$4":"$6":"$5; print $0}' tmp >!{params.collection_name}_lifted_b${NEWBUILD}.bim
'''
}


process liftover_pruned {
    tag "${params.collection_name}"
    cpus 2
    publishDir params.output, mode: 'copy'
    input:
    tuple file(bed), file(bim), file(fam), file(logfile) from for_liftover_pruned
    output:
    tuple file("${params.collection_name}.pruned_lifted*.bed"), file("${params.collection_name}.pruned_lifted*.bim"), file("${params.collection_name}.pruned_lifted*.fam"), file("${params.collection_name}.pruned_lifted*.log") optional true
    //file "postlift.${chrom}" into for_merge_lift
shell:
'''
module load Plink/1.9
MEM=!{task.memory.toMega()-1000}

LIFTOVER="!{params.ucsc_liftover}/liftOver"

if [ "!{params.build}" == "37" ]; then
    SUFFIX="hg38to37"
    TEXT="hg38 to hg37"
    NEWBUILD="38"
else
    SUFFIX="hg37to38"
    TEXT="hg37 to hg38"
    NEWBUILD="37"
fi

BASENAME="!{bed.getBaseName()}_lifted_b$NEWBUILD"

if [ "!{params.build}" == "37" ]; then
    CHAIN=!{params.ucsc_liftover}/hg19ToHg38.over.chain.gz
elif [ "!{params.build}" == "38" ]; then
    CHAIN=!{params.ucsc_liftover}/hg38ToHg19.over.chain.gz
else
    echo "Error: Input must be build 37 or build 38 but is specified as !{params.build}" >&2
    exit 0
fi

if [ ! -x "$LIFTOVER" ]; then
    echo "Error: Cannot execute ${LIFTOVER}. File does not exist or is not executable." >&2
    exit 0
fi


## Translate to chr:pos-pos and chr23 -> chrX
#mawk '{print "chr"$1":"$4"-"$4}' !{bim} \\
#    | sed 's/chr23/chrX/g' \\
#    | sed 's/chr24/chrY/g' \\
#    | sed 's/chr25/chrX/g' \\
#    | sed 's/chr26/chrMT/g' \\
#    > prelift.pos

# Create UCSC BED file from BIM
# Columns: chrX pos-1 pos rsID
<!{bim} mawk '{ print $1, $4-1, $4, $4, $2 }' \\
    | sed 's/^chr23/chrX/' \\
    | sed 's/^chr24/chrY/' \\
    | sed 's/^chr25/chrX/' \\
    | sed 's/^chr26/chrMT/' >prelift.bed

# lift-over
$LIFTOVER prelift.bed $CHAIN postlift.bed unmapped.bed

# Generate exclude list for unmapped variants
<unmapped.bed mawk '$0~/^#/{next} {print $5}' >unmapped-variants &

# Generate update-pos list for Plink, exclude double variants
<postlift.bed mawk '{print $5,$3}' >new-pos &
wait


# Find chromosome switchers
<postlift.bed mawk '{newchr=substr($1,4); split($5,oldparts,":");oldchr=substr(oldparts[1],4); if(newchr!=oldchr){print $5}}' >chromosome-switchers &
# for chromosome-wise only: <postlift.bed mawk '{if($1 != "chr{chrom}") {print $5}}' >chromosome-switchers
<new-pos sort | uniq -d | cut -f1 -d" ">duplicates &
wait


cat chromosome-switchers >>duplicates

/opt/plink2 --bed !{bed} --bim !{bim} --fam !{fam} --memory $MEM --exclude duplicates --make-pgen --out no-duplicates
/opt/plink2 --pfile no-duplicates --exclude unmapped-variants --memory $MEM --update-map new-pos --make-pgen --sort-vars  --out ${BASENAME}
/opt/plink2 --pfile ${BASENAME} --make-bed --output-chr chrM --memory $MEM --out ${BASENAME}.tmp

plink --bfile ${BASENAME}.tmp --merge-x no-fail --output-chr chrM --memory $MEM --make-bed --out ${BASENAME}
mv ${BASENAME}.bim tmp
mawk '{$2="chr"$1":"$4":"$6":"$5; print $0}' tmp >${BASENAME}.bim
rm postlift.bed ${BASENAME}.tmp*

'''
}


process merge_liftover_tables {
    tag "${params.collection_name}"
    publishDir params.output, mode: 'copy'

    input:
    file(postlift) from for_merge_tables.collect()

    output:
    file ("${params.collection_name}.liftover.table.*") into (for_lift_sumstats_table_plink, for_lift_sumstats_table_saige)

    shell:
'''

if [ "!{params.build}" == "37" ]; then
    SUFFIX="hg37to38"
else
    SUFFIX="hg38to37"
fi

cat postlift.* >"!{params.collection_name}.liftover.table.${SUFFIX}"

'''
}


process lift_plink_sumstats {
    tag "${params.collection_name}"
    publishDir params.output, mode: 'copy'
    label 'big_mem'
    label 'long_running'

    input:
    file(sumstats) from for_lift_sumstats_plink
    file(lifttable) from for_lift_sumstats_table_plink

    output:
    file("${sumstats}.lifted.*")

shell:
'''

if [ "!{params.build}" == "37" ]; then
    SUFFIX="hg38to37"
    TEXT="hg38 to hg37"
    NEWBUILD="38"
else
    SUFFIX="hg37to38"
    TEXT="hg37 to hg38"
    NEWBUILD="37"
fi


echo "Lifted from ${TEXT}. Liftover protocol: " >new-meta
stats2hg19_plink.pl !{lifttable} !{sumstats} >>new-meta
mv !{sumstats}.b37 !{sumstats}.lifted.${NEWBUILD}
mv new-meta !{sumstats}.lifted.${NEWBUILD}.txt
'''
}


process lift_saige_sumstats {
    tag "${params.collection_name}"
    publishDir params.output, mode: 'copy'
    label 'big_mem'
    label 'long_running'

    input:
    file(sumstats) from for_lift_sumstats_saige
    file(lifttable) from for_lift_sumstats_table_saige

    output:
    file("${sumstats}.lifted.*")

shell:
'''

if [ "!{params.build}" == "37" ]; then
    SUFFIX="hg38to37"
    TEXT="hg38 to hg37"
    NEWBUILD="38"
else
    SUFFIX="hg37to38"
    TEXT="hg37 to hg38"
    NEWBUILD="37"
fi


echo "Lifted from ${TEXT}. Liftover protocol: " >new-meta
stats2hg19_saige.pl !{lifttable} !{sumstats} >>new-meta
mv !{sumstats}.b37 !{sumstats}.lifted.${NEWBUILD}
mv new-meta !{sumstats}.lifted.${NEWBUILD}.txt
'''
}

