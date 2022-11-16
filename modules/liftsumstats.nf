process lift_plink_sumstats {
    tag "${params.collection_name}"
    publishDir params.output, mode: 'copy'
    label 'perl'

    input:
    file(sumstats)// from for_lift_sumstats_plink
    file(lifttable)// from for_lift_sumstats_table_plink

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
    label 'perl'

    input:
    file(sumstats)// from for_lift_sumstats_saige
    file(lifttable)// from for_lift_sumstats_table_saige

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
