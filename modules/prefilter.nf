// prefilter
// ---------
// Prepares the imputed files for futher analysis:
// - generates a tabix index
// - filters samples according to the FAM file given via --fam.
// Creates new VCF/GZ sets.
process prefilter {
    tag "${params.collection_name}"
    label "prefilter"
	scratch params.scratch

    input:
    tuple path(vcf), val(filetype)
    each path(inc_fam)

    output:
    tuple path("*.ap_prf.vcf.gz"), path("*.ap_prf.vcf.gz.tbi"), val(filetype) //into for_extract_chromosome_code //for_split_vcf, for_extract_dosage, for_gen_r2, for_plink_imp

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
bcftools query -f '%CHROM\\n' !{vcf} >chromlong

head  chromlong -n 1 >chrom

CHROM=$(cat chrom)
# Ok, we need to remove some samples. tabix afterwards
# If there's nothing to remove, skip this time-consuming step.
if [ "$REMOVE_COUNT" != "0" ]; then
    bcftools view --threads !{task.cpus} -S keep-samples !{vcf} -o $CHROM.ap_prf.vcf.gz -O z
    tabix $CHROM.ap_prf.vcf.gz
else
    ln -s !{vcf} $CHROM.ap_prf.vcf.gz
    ln -s !{vcf}.tbi $CHROM.ap_prf.vcf.gz.tbi
fi

'''
}
//TODO: include maf filter for prefiltering step?
//Like so: bcftools view -q 0.01:minor
//-i 'MAF>0.01'