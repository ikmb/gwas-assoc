process generate_pcs {
    tag "${params.collection_name}"
    label 'flashpca2'
	scratch params.scratch

input:
    tuple file(bed), file(bim), file(fam), file(log)

output:
    file "pcs.txt"

shell:
'''
MEM=!{task.memory.toMega()-1000}
#flashpca2 -d 10 --bfile !{bim.baseName} --memory $MEM --numthreads !{task.cpus} --outload loadings.txt --outmeansd meansd.txt
/home/flashpca-user/flashpca/flashpca -d !{params.pca_dims} --bfile !{bim.baseName} --memory $MEM --numthreads !{task.cpus} --outload loadings.txt --outmeansd meansd.txt

'''
}
