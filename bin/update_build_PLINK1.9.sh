#!/bin/sh

#A script for updating a binary ped file using one of Will's strand files
#NRR 17th Jan 2012

#V2 13th Feb 2012. Added code to retain only SNPs in the strand file
#V3 16th Aug 2018 dellinghaus. Update to Plink 1.9
#V4 17th Aug 2018 jkaessens. Place temporaries in SLURM scratch dir (if available)

#Required parameters:
#1. The original bed stem (not including file extension suffix)
#2. The strand file to apply
#3. The new stem for output
#Result: A new bed file (etc) using the new stem

#Unpack the parameters into labelled variables
stem=$1
strand_file=$2
outstem=$3
echo Input stem is $stem
echo Strand file is $strand_file
echo Output stem is $outstem

TEMP_DIR="."
if [ -d /scratch/SlurmTMP/`whoami`.$SLURM_JOB_ID ]; then
    echo Using /scratch/SlurmTMP/`whoami`.$SLURM_JOB_ID/ for temporaries
fi

tmp_base=$TEMP_DIR/$(basename $strand_file)

#Cut the strand file into a series of Plink slices
chr_file=$tmp_base.chr
pos_file=$tmp_base.pos
flip_file=$tmp_base.flip
cat $strand_file | cut -f 1,2 > $chr_file
cat $strand_file | cut -f 1,3 > $pos_file
cat $strand_file | awk '{if ($5=="-") print $0}' | cut -f 1 > $flip_file

#Because Plink only allows you to update one attribute at a time, we need lots of temp
#Plink files
temp_prefix=$TEMP_DIR/TEMP_FILE_XX72262628_
temp1="$temp_prefix-after-update-chr"
temp2="$temp_prefix-after-update-pos"
temp3="$temp_prefix-after-flip"

#1. Apply the chr
#plink --noweb --allow-no-sex --bfile $stem --update-map $chr_file --update-chr --make-bed --out $temp1
echo Trying to update chromosomes
plink --allow-no-sex --bfile $stem --update-chr $chr_file --make-bed --out $temp1
#2. Apply the pos
echo Trying to update positions
#plink --noweb --allow-no-sex --bfile $temp1 --update-map $pos_file --make-bed --out $temp2
plink --allow-no-sex --bfile $temp1 --update-map $pos_file --make-bed --out $temp2
#3. Apply the flip
echo Not flipping.
#plink --noweb --allow-no-sex --bfile $temp2 --flip $flip_file --make-bed --out $temp3
#plink --allow-no-sex --bfile $temp2 --flip $flip_file --make-bed --out $temp3
#4. Extract the SNPs in the pos file, we don't want SNPs that aren't in the strand file
echo Trying to filter out those SNPs without strand information
#plink --noweb --allow-no-sex --bfile $temp3 --extract $pos_file --make-bed --out $outstem
plink --allow-no-sex --bfile $temp2 --extract $pos_file --make-bed --out $outstem

#Now delete any temporary artefacts produced
# rm -f $temp_prefix*

