
process merge_plink {
    tag "${params.collection_name}"
	 scratch params.scratch
    publishDir params.output, mode: 'copy'
    label 'plink'

    input:
    file(filelist)

    output:
    tuple file("${params.collection_name}.bed"), file("${params.collection_name}.bim"), file("${params.collection_name}.fam"), file("${params.collection_name}.log")// into for_prune //, for_liftover


    shell:
	'''
		ls *.bed | xargs -i -- basename {} .bed | tail -n +2 >merge-list
		MEM=!{task.memory.toMega()-1000}
		FIRSTNAME=$(ls *.bed | xargs -i -- basename {} .bed | tail -n +1 | head -n 1)
		plink --bfile $FIRSTNAME --threads !{task.cpus} --memory $MEM --merge-list merge-list --make-bed --allow-no-sex --indiv-sort none --keep-allele-order --out !{params.collection_name}
	'''
}

process plink_assoc {
    tag "${params.collection_name}.${chrom}"
    label 'plink'
    scratch params.scratch
    input:
    tuple path(dosagemap), path(dosage), val(chrom), path(pheno), path(cols), path(fam)

    output:
    path("${chrom}.assoc.dosage")// into for_plink_results

shell:
'''
MEM=!{task.memory.toMega()-1000}
plink --fam !{fam} --map !{dosagemap} --dosage !{dosage} skip0=2 skip1=0 skip2=1 format=3 case-control-freqs \\
    --covar !{pheno} --covar-name $(cat !{cols}) \\
    --allow-no-sex --ci 0.95 \\
    --out !{chrom} --memory $MEM
'''
}

// Merge b37 Plink chromosome-wise sets into a single Plink set
process merge_plink_lifted {
    tag "${params.collection_name}"
    publishDir params.output, mode: 'copy'
	scratch params.scratch
    label 'plink'

    input:
    path(filelist)// from for_merge_lifted.collect()

    output:
    tuple path("${params.collection_name}_lifted*.bed"), path("${params.collection_name}_lifted*.bim"), path("${params.collection_name}_lifted*.fam"), path("${params.collection_name}_lifted*.log")

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
FIRSTNAME=$(ls *.bed | xargs -i -- basename {} .bed | tail -n +1 | head -n1)
plink --bfile $FIRSTNAME --memory $MEM --threads !{task.cpus} --merge-list merge-list --make-bed --allow-no-sex --indiv-sort none --keep-allele-order --out !{params.collection_name}_lifted_b$NEWBUILD
mv !{params.collection_name}_lifted_b${NEWBUILD}.bim tmp
#this was gawk originally:
awk '{$2="chr"$1":"$4":"$6":"$5; print $0}' tmp >!{params.collection_name}_lifted_b${NEWBUILD}.bim
'''
}
