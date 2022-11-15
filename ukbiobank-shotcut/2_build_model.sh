#!/bin/bash

source env

<$COVARFILE gawk 'FNR==1{print $0;next} {$2--; print $0}' >covars.txt
<$COVARFILE gawk 'BEGIN{OFS=" "; print "AGExAGE AGExSEX"} NR>1{print $4*$4,$4*$3}' >more-covars.txt
paste -d" " covars.txt more-covars.txt | tr -s ' ' >all-covars.txt

rm -f ukb_covid19.saige.varianceRatio.txt
rm -f ukb_covid19.saige.rda
$SINGULARITY exec -B /work_ifs $CONTAINER step1_fitNULLGLMM.R \
	--plinkFile=ukb_covid19_pruned \
	--phenoFile=all-covars.txt \
	--phenoCol=PHENO \
	--sampleIDColinphenoFile=IID \
	--covarColList=SEX,AGE,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10 \
	--traitType=binary \
	--outputPrefix=ukb_covid19.saige \
	--nThreads=32 \
	--LOCO=FALSE
