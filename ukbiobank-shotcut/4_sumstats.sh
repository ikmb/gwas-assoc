#!/bin/bash
#SBATCH --time=1-0
#SBATCH --output=4_sumstats.log

# Parameters

source env

echo $(date)

rm -f ukb_covid19.stats.files

echo "Generating merge file list..."
for chr in $CHROMOSOMES
do
    ls tmp_$chr/*.stats -1 | sort -n -t. -k3 >>ukb_covid19.stats.files
done

echo "Merging files..."
head -n1 $(head -n1 ukb_covid19.stats.files) >ukb_covid19.stats
while read line; do
    tail -n +2 $line >>ukb_covid19.stats
done <ukb_covid19.stats.files
echo "Done"

