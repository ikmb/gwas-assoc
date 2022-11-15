//needs splitting to individual processes
process liftover_pruned {
    tag "${params.collection_name}"
    label 'liftover'
    publishDir params.output, mode: 'copy'
	scratch params.scratch

    input:
        tuple path(bed), path(bim), path(fam), path(logfile) 
    output:
        tuple path("${params.collection_name}.pruned_lifted*.bed"), path("${params.collection_name}.pruned_lifted*.bim"), path("${params.collection_name}.pruned_lifted*.fam"), path("${params.collection_name}.pruned_lifted*.log"), optional: true //emit: outputcollection,
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

process liftover {
    tag "${params.collection_name}.${chrom}"
    label 'liftover'
	scratch params.scratch
//    publishDir params.output, mode: 'copy'
    input:
    tuple path(bed), path(bim), path(fam), path(logfile),val(chrom)// from for_liftover
    output:
    tuple file("${bed.baseName}_lifted.bed"), file("${bim.baseName}_lifted.bim"), file("${fam.baseName}_lifted.fam"), file("${logfile.baseName}_lifted.log"), emit: mergelifted, optional: true // into for_merge_lifted
    path("postlift.${chrom}"), emit: for_merge_tables, optional: true //into for_merge_tables
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

