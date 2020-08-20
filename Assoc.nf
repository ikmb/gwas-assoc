

params.collection_name = "collection_name"
params.output = "."
params.more_covars = "."
params.liftover = 0

params.chrchr = 1
params.null_filter = "R2>0.8"

def get_file_details(filename) {
    def m = filename =~ /\/([^\/]+)(\.vcf\.gz|\.bgen)$/
    println filename
    return [ m[0][1], m[0][2] ]
}

for_saige_imp = Channel.fromPath(params.input_imputed_glob, checkIfExists: true).map { it -> 
    def match = get_file_details(it)
    [it, match[0], match[1]] }


// does not really do anything besides pulling a docker container into the
// local cache that will be needed for saige.
// The saige pipeline master will be run on a compute node and might not
// have internet access for pulling from the docker registry.
process get_containers {
    container "docker://wzhou88/saige:0.38"
output:
    val("ok") into for_saige_gate
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
    tuple file("${chrom}.filtered.vcf.gz"), file("${chrom}.filtered.vcf.gz.tbi"), val(chrom), val(filetype) into for_saige, for_extract_dosage, for_gen_r2, for_plink_imp

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

process merge_r2 {
    tag "${params.collection_name}.${chrom}"
    input:
    file r2 from for_prune_r2.collect()

    output:
    file("r2-include.sorted") into for_prune_merged

shell:
'''
cat r2-include.* >r2-include

# Check if "chr" prefix is used in VCFs. If not, prefix the r2 list now
if [ "!{params.chrchr}" = "0" ]; then
    mawk '{print "chr"$0}' r2-include | sort >r2-include.sorted
else
    sort r2-include >r2-include.sorted
fi
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


/* Prune b38 for SAIGE null model */
process prune {
    tag "${params.collection_name}"
    publishDir params.output, mode: 'copy'
    memory 60000
    cpus 4
    input:
    tuple file(bed), file(bim), file(fam), file(logfile), val (chrom) from for_prune
    file r2 from for_prune_merged

    output:
    tuple file("${params.collection_name}.pruned.${chrom}.bed"), file("${params.collection_name}.pruned.${chrom}.bim"), file("${params.collection_name}.pruned.${chrom}.fam"), file("${params.collection_name}.pruned.${chrom}.log") into for_saige_null, for_liftover_pruned, for_generate_pcs

shell:
'''
module load Plink/1.9
echo Generating PCA SNP List file for variant selection

/opt/plink2 --memory 60000 --threads 4 --bed !{bed} --bim !{bim} --fam !{fam} --extract !{r2} --indep-pairwise 50 5 0.2 --out _prune
/opt/plink2 --memory 60000 --threads 4 --bed !{bed} --bim !{bim} --fam !{fam} --extract _prune.prune.in --maf 0.05 --make-bed --out intermediate

python -c 'from SampleQCI_helpers import *; write_snps_autosomes_noLDRegions_noATandGC_noIndels("!{bim}", "include_variants")'
sort include_variants >include_variants.sorted

comm -12 include_variants.sorted !{r2} >include-r2-variants
/opt/plink2 --memory 60000 --threads 4 --bfile intermediate --chr 1-22 --extract include-r2-variants --output-chr chrM --make-bed --out !{params.collection_name}.pruned
rm intermediate*

'''


}

process saige {
    tag "${params.collection_name}"

    input:
    tuple file(gz), file(tbi), val(chrom), val(filetype) from for_saige.collect()

    // wait for process to acquire container needed for SAIGE, just in case
    // the compute nodes, where process saige could be spawned, does not have
    // access to the docker registry
    val container_ok from for_saige_gate
    
shell:
'''

'''
}

process extract_dosage {
    tag "${params.collection_name}.${chrom}"

    input:
    tuple file (gz), file(tbi), val(chrom), val(filetype)

    output:
    tuple file("${chrom}.PLINKdosage.map"), file("${chrom}.PLINKdosage.gz"), val(chrom) into for_plink

shell:
'''

VCF=!{gz}
TARGET=$(basename $VCF .vcf.gz).PLINKdosage

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


$BGZIP -c -d $VCF \
    | $AWK -f awkprog2.awk \
    | tee $FIFO \
    | gzip >$TARGET.gz

wait

rm -f $FIFO
'''
}

process plink {
    tag "${params.collection_name}"

    input:
    tuple file(dosagemap), file(dosage), val(chrom) from for_plink
    file inc_fam

shell:
'''

'''
}
