process merge_saige_results {
    tag "${params.collection_name}"
    publishDir params.output, mode: 'copy'
	scratch params.scratch
	label 'base'

    input:
    path(sumstats) //from for_merge_sumstats.collect()
    tuple path(pheno), path(cols) //from for_merge_sumstats_covars
//    file in_covars


    output:
    path("${params.collection_name}.SAIGE.stats")// into for_lift_sumstats_saige
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

process merge_liftover_tables {
    tag "${params.collection_name}"
    publishDir params.output, mode: 'copy'
	scratch params.scratch

    input:
    file(postlift)// from for_merge_tables.collect()

    output:
    file ("${params.collection_name}.liftover.table.*")// into (for_lift_sumstats_table_plink, for_lift_sumstats_table_saige)

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

process merge_plink_results {
    tag "${params.collection_name}"
	label 'base'
	scratch params.scratch
    publishDir params.output, mode: 'copy'

    input:
    path(stats)// from for_plink_results.collect()

    output:
    path("${params.collection_name}.Plink.stats") //into for_lift_sumstats_plink

shell:
'''
# extract first line, convert tabs to space
head -n1 !{stats[0]} | tr -s '\t ' ' ' | xargs >!{params.collection_name}.Plink.stats
ls !{stats} | sort -n | xargs -n1 tail -n +2 | gawk '{if(substr($1,1,3)!="chr"){$1="chr"$1} $2=$1":"$3":"$4":"$5; print}'>>!{params.collection_name}.Plink.stats
'''
}

process merge_r2 {
    tag "${params.collection_name}"
	scratch params.scratch
	label 'base'

    input:
    path(r2)

    output:
    path("r2-include.sorted")

	shell:
	'''
		cat r2-include.* >r2-include
		export TMPDIR=.
		gawk '$0 !~ /^chr/ {$1="chr"$1} {print}' r2-include | sort >r2-include.sorted
	'''
}
