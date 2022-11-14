process prune_python_helper {
	tag "${params.collection_name}"
    label 'python2'
    scratch params.scratch

    input:
    tuple path(bed), path(bim), path(fam), path(logfile)
	path(r2)
    output:
    path(outputfile)
	shell:
	outputfile = "include-r2-variants"
	'''
	R2=!{r2}
	python -c 'from SampleQCI_helpers import *; write_snps_autosomes_noLDRegions_noATandGC_noIndels("!{bim}", "include_variants")'
	sort include_variants >include_variants.sorted
	comm -12 include_variants.sorted $R2 > !{outputfile}
	'''
}