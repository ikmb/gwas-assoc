
// returns a file handle or bails if the file does not exist
fileExists = { fn ->
              if (fn.exists())
                  return fn;
              else
                  error("File not found: $fn")
}

/* expected parameters:
 - params.input_imputed_glob='/path/to/[1-9][0-9].vcf.gz' (both vcf.gz and bgen are acceptable)
   - chromosome number must be the very first part of the file name.
 - params.input_qced_glob='/path/to/files.{bim,bed,fam}'
 - params.fam=/path/to/qced.fam
 - params.covars=/path/to/pca.evec
*/

// inputs
for_saige_imp = Channel.fromPath(params.input_imputed_glob, checkIfExists: true).map { it -> 
    def match = get_file_details(it)
    [it, match[0], match[1]] }

//inc_qced = Channel.fromPath(params.input_qced_glob, checkIfExists: true)
inc_fam = fileExists(file(params.fam))

anno = fileExists(file(params.anno))

// defaults
params.collection_name = "dataset"
// params.nxfdir = "/home/sukmb388/nextflow-staging"
params.nxfdir = "${workflow.projectDir}"
params.output = "."
params.more_covars = ""
params.liftover = 0

//more_covars = file(params.more_covars)

// Use "chr" as chromosome prefix
params.chrchr = 1
// params.null_filter = 'INFO/INFO>0.9'
params.null_filter='R2>0.80'

def get_file_details(filename) {
    def m = filename =~ /\/([^\/]+)(\.vcf\.gz|\.bgen)$/
    println filename
    return [ m[0][1], m[0][2] ]
}

// prefilter
// ---------
// Prepares the imputed files for futher analysis:
// - generates a tabix index
// - filters samples according to the FAM file given via --fam.
// Creates new VCF/GZ sets.
process prefilter {
    tag "${params.collection_name}.${chrom}"
//    publishDir params.output, mode: 'copy'
    time {48.h * task.attempt}
    cpus 4
    input:
    tuple file(vcf), val(chrom), val(filetype) from for_saige_imp
    file inc_fam

    output:
    tuple file("${chrom}.filtered.vcf.gz"), file("${chrom}.filtered.vcf.gz.tbi"), val(chrom), val(filetype) into for_split_vcf, for_gen_r2, for_plink_imp, for_extract_dosage

shell:
'''
tabix -p vcf !{vcf}

# Gather lists of samples with proper double-ID handling
bcftools query -l !{vcf} | sort >vcf-samples
mawk '{if($1!="0") {$2=$1"_"$2; $1="0";} print $2}' !{inc_fam} | sort >fam-samples

# Keep only those samples that are in VCF and in FAM file
comm -12 vcf-samples fam-samples >keep-samples

# Count those that we ought to remove
REMOVE_COUNT=$(comm -23 vcf-samples fam-samples | wc -l)

# Ok, we need to remove some samples. tabix afterwards
# If there's nothing to remove, skip this time-consuming step.
if [ "$REMOVE_COUNT" != "0" ]; then
    bcftools view --threads 4 -S keep-samples !{vcf} -o !{chrom}.filtered.vcf.gz -O z
    tabix !{chrom}.filtered.vcf.gz
else
    ln -s !{vcf} !{chrom}.filtered.vcf.gz
    ln -s !{vcf}.tbi !{chrom}.filtered.vcf.gz.tbi
fi

'''
}

process extract_dosage {
    tag "${params.collection_name}.${chrom}"

    input:
    tuple file(vcf), file(tbi), val(chrom), val(filetype) from for_extract_dosage

    output:
    tuple file("${chrom}.PLINKdosage.gz"), file("${chrom}.PLINKdosage.map"), val(chrom) into for_plink_assoc

shell:
'''
VCF=!{vcf}
TARGET=!{chrom}.PLINKdosage

# Theoratically, gawk or awk would also work but are much slower
AWK=mawk
BGZIP=bgzip
FIFO=$TARGET.fifo
rm -f $FIFO
mkfifo $FIFO

# Generate awk script
cat <<'AwkProg2' >awkprog2.awk
BEGIN{
    printf "CHR BP SNP A1 A2 INFO"
}
{
    # skip header lines
    if (($1 ~ /^##/)) {
    } else {
        # Copy header line
        if (($1 ~ /^#CHROM/)) {
            for (i=10; i<=NF; i++)
                printf " "$i" "$i
            printf "\n"
        } else {
            # Extract r^2 score from INFO field
            {
                split($8, fields, /(;|=)/)
                for(i=1; i in fields; i++) {
                    if(fields[i] == "INFO" || fields[i] == "DR2" || fields[i] == "R2") {
                            info = fields[i+1]
                            break
                    }
                }
            }

            if(!schonmalgesehen[$1":"$2":"$4":"$5]) {
                printf $1" "$2" "$1":"$2":"$4":"$5" "$4" "$5" "info
                for (i=10; i<=NF; i++) {
                    split($i,array,":")
                    # Missing fields are specified as ".", not splittable by ","
                    if(array[4] == ".") {
                        printf " . . ."
                    } else {
                        n=split(array[4],array_GP,",")
                        if(n==3) {
                            printf " "array_GP[1]" "array_GP[2]" "array_GP[3]
                        } else {
                            printf " "array_GP[1]" 0 "array_GP[2]
                    }
                    }
                }
                printf "\n"

                schonmalgesehen[$1":"$2":"$4":"$5] = 1
            }
        }

    }
}
AwkProg2


(<$FIFO tail -n +2 | $AWK '{print $1, $3, 0, $2}' | uniq) >$TARGET.map &
#cat $FIFO  &


$BGZIP -c -d $VCF \
    | $AWK -f awkprog2.awk \
    | tee $FIFO \
    | gzip >$TARGET.gz

wait
#rm -f awkprog2.awk
rm -f $FIFO
'''
}


/* Prepare a list of r2-filtered variants for SAIGE step1 model generation */
process gen_r2_list {
    tag "${params.collection_name}.${chrom}"
    input:
    tuple file(vcf), file(tbi), val(chrom), val(filetype) from for_gen_r2

    output:
    file("r2-include.${chrom}") into for_prune_r2

shell:
'''
bcftools query -i'!{params.null_filter}' -f '%CHROM:%POS:%REF:%ALT\\n' !{vcf} >r2-include.!{chrom}
'''
}

/* Convert chromosome-wise imputed VCFs to Plink bim/bed/fam sets, also fills
   phenotype and sex columns from QC FAM file. */
process make_plink {
//    publishDir params.output, mode: 'copy'
    memory 50.GB
    tag "${params.collection_name}.$chrom"

    input:
    tuple file(vcf), file(tbi), val(chrom), val(filetype) from for_plink_imp
    file fam from inc_fam

    output:
    tuple file("${params.collection_name}.${chrom}.bed"), file("${params.collection_name}.${chrom}.bim"), file("${params.collection_name}.${chrom}.fam"), file("${params.collection_name}.${chrom}.log"),val(chrom) into for_liftover, for_merge_b38

shell:
'''
# Generate double-id FAM
mawk '{if($1!="0") {$2=$1"_"$2; $1="0";} print $0}' !{fam} >new-fam
/opt/plink2 --vcf !{vcf} --const-fid --memory 46000 --allow-no-sex --pheno new-fam --mpheno 4 --update-sex new-fam 3 --make-bed --keep-allele-order --out !{params.collection_name}.!{chrom}
mv !{params.collection_name}.!{chrom}.bim bim
mawk '{$2=$1":"$4":"$6":"$5; print $0}' <bim >!{params.collection_name}.!{chrom}.bim
'''
}

/* Use UCSC liftOver to convert Plink b38 to Plink b37 */
process liftover {
    tag "${params.collection_name}.${chrom}"

    when:
    params.liftover == 1

    input:
    tuple file(bed), file(bim), file(fam), file(logfile),val(chrom) from for_liftover

    output:
    tuple file("${bed.baseName}_b37.bed"), file("${bim.baseName}_b37.bim"), file("${fam.baseName}_b37.fam"), file("${logfile.baseName}_b37.log") into for_merge_b37
    file ("postlift.${chrom}") into for_merge_lift

//    tuple file ("postlift.${chrom}"), file("unmapped-variants"), file("duplicates") into for_merge_lift
//    tuple file("postlift.${chrom}"), file("unmapped-variants") into for_gen_liftdb
shell:
'''
module load Plink/1.9

LIFTOVER=!{params.nxfdir}/liftover
BASENAME=!{bed.getBaseName()}_b37
CHAIN=$LIFTOVER/assets/liftover/hg38ToHg19.over.chain.gz

## Translate to chr:pos-pos and chr23 -> chrX
#mawk '{print "chr"$1":"$4"-"$4}' !{bim} \\
#    | sed 's/chr23/chrX/g' \\
#    | sed 's/chr24/chrY/g' \\
#    | sed 's/chr25/chrX/g' \\
#    | sed 's/chr26/chrMT/g' \\
#    > prelift.pos

# Create UCSC BED file from BIM
# Columns: chrX pos-1 pos rsID
<!{bim} mawk '{ print "chr"$1, $4-1, $4, $4, $2 }' \\
    | sed 's/^chr23/chrX/' \\
    | sed 's/^chr24/chrY/' \\
    | sed 's/^chr25/chrX/' \\
    | sed 's/^chr26/chrMT/' >prelift.bed

# lift-over
$LIFTOVER/bin/liftOver prelift.bed $CHAIN postlift.bed unmapped.bed

# Generate exclude list for unmapped variants
<unmapped.bed mawk '$0~/^#/{next} {print $5}' >unmapped-variants

# Generate update-pos list for Plink, exclude double variants
<postlift.bed mawk '{print $5,$3}' >new-pos
<postlift.bed mawk '{if($1 != "chr!{chrom}") {print $5}}' >chromosome-switchers
<new-pos sort | uniq -d | cut -f1 -d" ">duplicates

cat chromosome-switchers >>duplicates

/opt/plink2 --bed !{bed} --bim !{bim} --fam !{fam} --exclude duplicates --make-pgen --out no-duplicates
/opt/plink2 --pfile no-duplicates --exclude unmapped-variants --update-map new-pos --make-pgen --sort-vars  --out ${BASENAME}
/opt/plink2 --pfile ${BASENAME} --make-bed --out ${BASENAME}.tmp

plink --bfile ${BASENAME}.tmp --merge-x no-fail --make-bed --out ${BASENAME}
mv ${BASENAME}.bim tmp
mawk '{$2="chr"$1":"$4":"$6":"$5; print $0}' tmp >${BASENAME}.bim
mv postlift.bed postlift.!{chrom}

'''
}

/* Merge b37 Plink chromosome-wise sets into a single Plink set */
process merge_plink_b37 {
    tag "${params.collection_name}"
    publishDir params.output, mode: 'copy'
    memory 64000
    cpus 4
    time {8.h * task.attempt}

    when:
    params.liftover == 1

    input:
    file(filelist) from for_merge_b37.collect()

    output:
    tuple file("${params.collection_name}_b37.bed"), file("${params.collection_name}_b37.bim"), file("${params.collection_name}_b37.fam"), file("${params.collection_name}_b37.log")

    shell:
    '''
ls *.bed | xargs -i -- basename {} .bed | tail -n +2 >merge-list
FIRSTNAME=$(ls *.bed | xargs -i -- basename {} .bed | head -n1)
module load Plink/1.9
plink --bfile $FIRSTNAME --memory 60000 --threads 4 --merge-list merge-list --make-bed --allow-no-sex --indiv-sort none --keep-allele-order --out !{params.collection_name}_b37

    '''

}

/* Merge b38 Plink chromosome-wise sets into a single Plink set */
process merge_plink_b38 {
    tag "${params.collection_name}"
    publishDir params.output, mode: 'copy'
    memory 64000
    cpus 4
    time {8.h * task.attempt}
    input:
    file(filelist) from for_merge_b38.collect()

    output:
    tuple file("${params.collection_name}.bed"), file("${params.collection_name}.bim"), file("${params.collection_name}.fam"), file("${params.collection_name}.log") into for_prune


    shell:
    '''
ls *.bed | xargs -i -- basename {} .bed | tail -n +2 >merge-list
FIRSTNAME=$(ls *.bed | xargs -i -- basename {} .bed | head -n1)
module load Plink/1.9
plink --bfile $FIRSTNAME --threads 4 --memory 60000 --merge-list merge-list --make-bed --allow-no-sex --indiv-sort none --keep-allele-order --out !{params.collection_name}
mv !{params.collection_name}.bim tmp
mawk '{$2="chr"$1":"$4":"$6":"$5; print $0}' tmp >!{params.collection_name}.bim
'''
}

/* Prune b38 for SAIGE null model */
process prune {
    tag "${params.collection_name}"
    publishDir params.output, mode: 'copy'
    memory 60000
    cpus 4
    input:
    tuple file(bed), file(bim), file(fam), file(logfile) from for_prune
    file "r2-include.*" from for_prune_r2.collect()

    output:
    tuple file("${params.collection_name}.pruned.bed"), file("${params.collection_name}.pruned.bim"), file("${params.collection_name}.pruned.fam"), file("${params.collection_name}.pruned.log") into for_saige_null, for_liftover_pruned, for_generate_pcs

shell:
'''
module load Plink/1.9
echo Generating PCA SNP List file for variant selection

cat r2-include.* >r2-include

# Check if "chr" prefix is used in VCFs. If not, prefix the r2 list now
if [ "!{params.chrchr}" = "0" ]; then
    mawk '{print "chr"$0}' r2-include | sort >r2-include.sorted
else
    sort r2-include >r2-include.sorted
fi
/opt/plink2 --memory 60000 --threads 4 --bed !{bed} --bim !{bim} --fam !{fam} --extract r2-include.sorted --indep-pairwise 50 5 0.2 --out _prune
/opt/plink2 --memory 60000 --threads 4 --bed !{bed} --bim !{bim} --fam !{fam} --extract _prune.prune.in --maf 0.05 --make-bed --out intermediate

python -c 'from SampleQCI_helpers import *; write_snps_autosomes_noLDRegions_noATandGC_noIndels("!{bim}", "include_variants")'
sort include_variants >include_variants.sorted

comm -12 include_variants.sorted r2-include.sorted >include-r2-variants
/opt/plink2 --memory 60000 --threads 4 --bfile intermediate --chr 1-22 --extract include-r2-variants --output-chr chrM --make-bed --out !{params.collection_name}.pruned
rm intermediate*

'''


}

process generate_pcs {
    cpus 16
    tag "${params.collection_name}"
    memory 60000

    publishDir params.output, mode: 'copy'
input:
    tuple file(bed), file(bim), file(fam), file(log) from for_generate_pcs

output:
    file "pcs.txt" into for_draw_pcs, for_saige_covars, for_plink_covars

shell:
'''
module load FlashPCA
flashpca2 -d 10 --bfile !{bim.baseName} --numthreads 16 --outload loadings.txt --outmeansd meansd.txt
'''
}

process fix_annotations {
input:
file(anno)

output:
file("fixed-annotations.txt") into for_draw_pcs_anno

shell:
'''
mawk 'FNR==1{print $0;next} {if($1!="0") {$2=$1"_"$2; $1="0";} print $0}' !{anno} >fixed-annotations.txt
'''
}

process draw_pcs {

    publishDir params.output, mode: 'copy'
    validExitStatus 0,1
input:
    file(pcs) from for_draw_pcs
    file(ann) from for_draw_pcs_anno
output:
    file "pcs.txt.batch.png" optional true
    file "pcs.txt.diagnosis.png" optional true
shell:
'''

(
cat <<'EOF_R'
library(ggplot2)
library(dplyr)
library(data.table)
library(ggrepel)

name <- "!{params.collection_name}"
file.pca.evec <- "!{pcs}"
file.individuals_annotation <- "!{ann}"
pca.evec <- fread(file.pca.evec, sep = "\\t", header= TRUE)
individuals_annotation <- fread(file.individuals_annotation, sep = "\\t", header= TRUE)

data <- pca.evec %>%
  left_join(individuals_annotation, ., by=c("individualID"="IID"))

ggplot(data, aes(x=PC1, y=PC2)) +

    geom_point(aes(color=batch)) +

    theme_bw() +
    theme(
    legend.position="top",
    legend.title = element_text(size = 3),
    legend.text = element_text(size = 3)
    ) +

    ggtitle(name) + ylab("PC2") + xlab("PC1")

ggsave(paste(file.pca.evec, ".batch.png", sep=""), device = "png", width = 5, height = 5)

ggplot(data, aes(x=PC1, y=PC2)) +

    geom_point(aes(color=diagnosis)) +

    theme_bw() +
    theme(
    legend.position="top",
    legend.title = element_text(size = 3),
    legend.text = element_text(size = 3)
    ) +

    ggtitle(name) + ylab("PC2") + xlab("PC1")

ggsave(paste(file.pca.evec, ".diagnosis.png", sep=""), device = "png", width = 5, height = 5)
EOF_R
) >das_r_script.R

R --no-save <das_r_script.R
'''
}

process liftover_pruned {
    tag "${params.collection_name}"
    publishDir params.output, mode: 'copy'

    when:
    params.liftover == 1

    input:
    tuple file(bed), file(bim), file(fam), file(logfile) from for_liftover_pruned

    output:
    tuple file("${params.collection_name}.pruned_b37.bed"), file("${params.collection_name}.pruned_b37.bim"), file("${params.collection_name}.pruned_b37.fam"), file("${params.collection_name}.pruned_b37.log")
    //file "postlift.${chrom}" into for_merge_lift
shell:
'''
module load Plink/1.9

LIFTOVER=!{params.nxfdir}/liftover
BASENAME=!{params.collection_name}.pruned_b37
CHAIN=$LIFTOVER/assets/liftover/hg38ToHg19.over.chain.gz

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
$LIFTOVER/bin/liftOver prelift.bed $CHAIN postlift.bed unmapped.bed

# Generate exclude list for unmapped variants
<unmapped.bed mawk '$0~/^#/{next} {print $5}' >unmapped-variants

# Generate update-pos list for Plink, exclude double variants
<postlift.bed mawk '{print $5,$3}' >new-pos

# Find chromosome switchers
<postlift.bed mawk '{newchr=substr($1,4); split($5,oldparts,":");oldchr=substr(oldparts[1],4); if(newchr!=oldchr){print $5}}' >chromosome-switchers
# for chromosome-wise only: <postlift.bed mawk '{if($1 != "chr{chrom}") {print $5}}' >chromosome-switchers
<new-pos sort | uniq -d | cut -f1 -d" ">duplicates

cat chromosome-switchers >>duplicates

/opt/plink2 --bed !{bed} --bim !{bim} --fam !{fam} --exclude duplicates --make-pgen --out no-duplicates
/opt/plink2 --pfile no-duplicates --exclude unmapped-variants --update-map new-pos --make-pgen --sort-vars  --out ${BASENAME}
/opt/plink2 --pfile ${BASENAME} --make-bed --out ${BASENAME}.tmp

plink --bfile ${BASENAME}.tmp --merge-x no-fail --make-bed --out ${BASENAME}
mv ${BASENAME}.bim tmp
mawk '{$2="chr"$1":"$4":"$6":"$5; print $0}' tmp >${BASENAME}.bim
rm postlift.bed

'''
}

/* Create SAIGE-compatible covars from PCA eigenvectors */
process make_saige_covars {
    input:
    file (evec) from for_saige_covars
    file inc_fam

    output:
    tuple file("pheno"), file("covar_cols.txt") into for_saige_null_covars, for_merge_sumstats_covars

    shell:
'''
# Re-format sample ID and family ID to VCF rules
mawk '{if($1!="0") {$2=$1"_"$2; $1="0";} print $0}' !{inc_fam}  | sort >new-fam

# Also re-format, but keep evec header
# mawk 'NR==1{print $0;next} {if($1!="0") {$2=$1"_"$2; $1="0";} print $0}' {in_covars}  | sort >evec

# Take evec file as a whole, SAIGE does not care about additional columns.
# Filter evec file according to samples contained in FAM (if not already done)

mawk 'FNR==NR{samples[$2];next} {if($2 in samples) { for(i=1;i<=(12);i++) {printf "%s%s", $i, (i<12?OFS:ORS)}}}' new-fam !{evec} >filtered-evec

if [ -e "!{params.more_covars}" ]; then
    mawk 'FNR==NR{samples[$2];next} {if($2 in samples) { for(i=1;i<=(6);i++) {printf "%s%s", $i, (i<6?OFS:ORS)}}}' new-fam !{params.more_covars} | sort >filtered-covars
    cut -f3-6 !{params.more_covars} >covars-column
    < filtered-covars cut -f3-6 -d" " >>covars-column
fi

EVEC_LINES=$(wc -l <filtered-evec)
FAM_LINES=$(wc -l <new-fam)
if [ $EVEC_LINES -ne $FAM_LINES ]; then
    echo I could not find covariates in !{evec} for every sample specified in !{inc_fam}. Please check.
    echo Aborting.
    exit 1
fi

echo "FID IID PC1 PC2 PC3 PC4 PC5 PC6 PC7 PC8 PC9 PC10" >evec.double-id.withheader
#echo "FID IID PC1 PC2 PC3 PC4 PC5 PC6 PC7 PC8 PC9 PC10 PC11 PC12 PC13 PC14 PC15 PC16 PC17 PC18 PC19 PC20" >evec.double-id.withheader
cat filtered-evec >>evec.double-id.withheader

# Take phenotype info from FAM, translate to SAIGE-encoding
echo "Pheno" >pheno-column
<new-fam  mawk '{if($6=="2") {$6="1";} else if($6=="1") {$6="0";} print $6}' >>pheno-column

# Merge both, replace space runs with single tabs for SAIGE
touch covars-column
if [ -e "!{params.more_covars}" ]; then
    echo PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC6,PC7,PC8,PC9,PC10,!{params.more_covars_cols} >covar_cols.txt
else
    echo PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC6,PC7,PC8,PC9,PC10 >covar_cols.txt
fi

paste -d" " evec.double-id.withheader pheno-column covars-column | tr -s ' ' \\\\t >pheno

'''
}

/* Create a SAIGE null model from b38 Plink build and PCs */
process saige_null {
    tag "${params.collection_name}"
    publishDir params.output, mode: 'copy'
    container "docker://wzhou88/saige:0.38"
    cpus 16
    time {2.d * task.attempt}

    input:
    tuple file(bed), file(bim), file(fam), file(logfile) from for_saige_null
    tuple file(pheno), file(covar_list) from for_saige_null_covars

    output:
    tuple file("${params.collection_name}.rda"), file("${params.collection_name}.varianceRatio.txt") into for_saige_assoc

shell:
    '''

sed 's/^chr//' !{bim.baseName}.bim >tmp.bim
ln -s !{bed} tmp.bed
ln -s !{fam} tmp.fam

step1_fitNULLGLMM.R \
    --plinkFile=tmp \
    --phenoFile=!{pheno} \
    --phenoCol=Pheno \
    --covarColList=$(cat !{covar_list}) \
    --sampleIDColinphenoFile=IID \
    --traitType=binary        \
    --outputPrefix=!{params.collection_name} \
    --nThreads=16 \
    --LOCO=FALSE
    '''
}


params.saige_chunk_size = 20000

process split_vcf_ranges {
    tag "${params.collection_name}.${chrom}"

    input:
    tuple file(vcf), file(tbi), val(chrom), val(filetype) from for_split_vcf
    output:
    tuple file(vcf), file(tbi), val(chrom), val(filetype), file("0*") into for_saige_assoc_imp

    shell:
    '''
    bcftools query -f "%POS\\n" !{vcf} >positions
    split -d --suffix-length=8 --lines=!{params.saige_chunk_size} positions '0'
    '''
}

/* Perform SAIGE Assoc test on imputed VCFs */
process saige_assoc {
    tag "${params.collection_name}.${chrom}.${chunk}"
    container "docker://wzhou88/saige:0.38"
    time {2.h * task.attempt}

    input:
    // <-- single VCF and chromosome number
    tuple file(vcf), file(tbi),  val(chrom), val(filetype),file(chunk) from for_saige_assoc_imp.transpose() // turn [chr1, [chunk1, chunk2, chunk3]] into [chr1, chunk1], [chr1, chunk2], [chr1, chunk3]
    tuple file(modelfile), file(varratio) from for_saige_assoc

    output:
    file("${params.collection_name}.${chrom}.${chunk.name}.SAIGE.stats") into for_merge_sumstats
shell:
'''
FIRSTPOS=$(head -n1 <!{chunk})
LASTPOS=$(tail -n1 <!{chunk})

if [ "!{params.chrchr}" = "0" ]; then
    CHR=!{chrom}
else
    CHR=chr!{chrom}
fi

step2_SPAtests.R \
    --vcfFile="!{vcf}" \
    --vcfFileIndex="!{tbi}" \
    --vcfField=DS \
    --chrom=$CHR \
    --minMAF=0.001 \
    --minMAC=1 \
    --GMMATmodelFile="!{modelfile}" \
    --varianceRatioFile="!{varratio}" \
    --SAIGEOutputFile=temp.stats \
    --numLinesOutput=2 \
    --start=$FIRSTPOS \
    --end=$LASTPOS \
    --IsOutputAFinCaseCtrl=TRUE \
    --IsOutputNinCaseCtrl=TRUE \
    --IsOutputHetHomCountsinCaseCtrl=TRUE

# add odds ratio
<temp.stats mawk 'NR==1{print $0 " OR";next} {print $0 " " exp($10)}' \
    >"!{params.collection_name}.!{chrom}.!{chunk.name}.SAIGE.stats"
'''
}

process merge_saige_results {
    tag "${params.collection_name}"
    publishDir params.output, mode: 'copy'

    input:
    file(sumstats) from for_merge_sumstats.collect()
    tuple file(pheno), file(cols) from for_merge_sumstats_covars

    output:
    tuple file("${params.collection_name}.SAIGE.*.stats"), file("${params.collection_name}.SAIGE.*.stats.txt") into for_lift_sumstats_sumstats

shell:
'''

head -n1 !{params.collection_name}.1.000000000.SAIGE.stats >!{params.collection_name}.SAIGE.stats

# for i in $(ls...) would not work anymore with thousands of files... too long command line.
# List all files into file, read it line by line.
#ls -1 !{params.collection_name}.*.*.stats >allfiles
#sort -t. -k2 allfiles >allfiles.sorted

for i in $(ls -1 !{params.collection_name}.*.*.stats | sort -t. -k2)
do 
    tail -n +2 $i >>!{params.collection_name}.SAIGE.stats
done

#while read -r line; do
#    tail -n +2 $line >>!{params.collection_name}.SAIGE.stats
#done <allfiles.sorted



#REAL_COVARS=$(readlink -f {in_covars})
#COVAR_TIME=$(stat -c %y $REAL_COVARS)

echo Collection name: !{params.collection_name}
#echo Covariates source: $REAL_COVARS >meta
#echo Last modified: $COVAR_TIME >>meta
echo Covars used: $(cat !{cols}) >>meta
echo FAM file used: !{inc_fam} >>meta
if [ -e "!{params.more_covars}" ]; then
    echo Additional covar file: !{params.more_covars} >>meta
else
    echo No additional covars. >>meta
fi
echo Null model calibration filter: !{params.null_filter} >>meta

mv !{params.collection_name}.SAIGE.stats !{params.collection_name}.SAIGE.stats
mv meta !{params.collection_name}.SAIGE.stats.txt
'''
}

process merge_liftover_tables {
    tag "${params.collection_name}"
    publishDir params.output, mode: 'copy'

    when:
    params.liftover == 1

    input:
    file(postlift) from for_merge_lift.collect()

    output:
    file "${params.collection_name}.hg38toHg19.table" into for_lift_sumstats_table

    shell:
    '''
cat postlift.* >"!{params.collection_name}.hg38toHg19.table"
    '''
}

process lift_sumstats {
    tag "${params.collection_name}"
    publishDir params.output, mode: 'copy'
    memory 40.GB
    time { 12.h * task.attempt }

    input:
    tuple file(sumstats), file(meta) from for_lift_sumstats_sumstats
    file(lifttable) from for_lift_sumstats_table

    output:
    file("${sumstats}.b37")
    file("${sumstats}.b37.txt")

shell:
'''
cp !{meta} new-meta
echo "Lifted from hg38 to hg19. Liftover protocol: " >>new-meta
stats2hg19.pl !{lifttable} !{sumstats} >>new-meta
mv new-meta !{sumstats}.b37.txt
'''
}


process plink_assoc {

input:
tuple file(gz), file(map), val(chrom) from for_plink_assoc
file inc_fam
file covars from for_plink_covars

output:
tuple file("${params.collection_name}.${chrom}.assoc.dosage"), file("${params.collection_name}.${chrom}.log") for_merge_dosage

shell:
'''
module load Plink/1.9

mawk '{ $1 = $2; print $0 }' !{inc_fam} >pcs.txt
head -n 1 !{covars} >pcs.txt
tail -n +2 | mawk '{print $2, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11, $12}' >>pcs.txt

plink --fam !{inc_fam} --map !{map} \\
      --dosage !{gz} skip0=2 skip1=0 skip2=1 format=3 case-control-freqs \\
      --covar pcs.txt --covar-name PC1-PC10 \\
      --allow-no-sex --ci 0.95 \\
      --out !{params.collection_name}.!{chrom}
'''
}

process merge_plink_stats {
input:
tuple file(dosage), file(log) from for_merge_dosage.collect()

output:
file("${params.collection_name}.plink.stats")

shell:
'''
cat !{params.collection_name}.!{chrom}.assoc.dosage >!{params.collection_name}.dosage
for (( i=2; i<=25; i++ ))
do
    if [ -e !{params.collection_name}.$i.assoc.dosage ]; then
        tail -n +2 !{params.collection_name}.$i.assoc.dosage >>!{params.collection_name}.dosage
    fi
done

mawk '{ if($8 >= 0.6) print }' !{params.collection_name}.dosage >!{params.collection_name}.dosage.RSQ0.6

'''

}
