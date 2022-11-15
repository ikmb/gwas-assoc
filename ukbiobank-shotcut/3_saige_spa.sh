#!/bin/bash

source env


echo $(date)
BASEDIR="$IMPUTED_DATA"

for chr in $CHROMOSOMES
do
	FILENAME="$BASEDIR/ukb_imp_chr${chr}_v3.bgen"
	TARGET=$(basename $FILENAME .bgen)
	
	echo "Preparing for chromosome $chr..."
	mkdir -p "tmp_$chr"
	cd tmp_$chr
	rm *.log *.stats
	CHUNKSIZE=10000
    $SQLITE3 ${FILENAME}.bgi "SELECT DISTINCT rsid FROM Variant" | tee ../tmp_${chr}_positions | split -a9 --numeric-suffixes=1 -l $CHUNKSIZE - chunk
	cd ..
	CHUNKS=$(ls tmp_$chr -1 | wc -l)
	
	echo "Scheduling $CHUNKS jobs..."
	
	cat >tmp.sbatch <<EOT
#!/bin/bash
#SBATCH -c1
#SBATCH --job-name=$TARGET.$chr
#SBATCH --array=1-$CHUNKS%200
#SBATCH --output=tmp_${chr}/%a.log
#SBATCH --time=6:00:00
#SBATCH --mem=4G

    NAME=\$(printf "chunk%09d" \$SLURM_ARRAY_TASK_ID)
	date
	echo "This is task \$SLURM_ARRAY_TASK_ID. Processing tmp_$chr/\$NAME"
	
    $SINGULARITY exec -B/work_ifs $CONTAINER	step2_SPAtests.R \
		--bgenFile=$FILENAME \
		--bgenFileIndex=${FILENAME}.bgi \
		--sampleFile=samples.1col \
		--minMAF=0.0001 \
		--minMAC=1 \
		--GMMATmodelFile=ukb_covid19.saige.rda \
		--varianceRatioFile=ukb_covid19.saige.varianceRatio.txt \
		--SAIGEOutputFile=tmp_$chr/${TARGET}.${chr}.\${SLURM_ARRAY_TASK_ID}.stats \
		--numLinesOutput=2 \
		--IsOutputAFinCaseCtrl=TRUE \
		--IsOutputNinCaseCtrl=TRUE \
		--IsOutputHetHomCountsinCaseCtrl=TRUE \
		--IsDropMissingDosages=TRUE \
		--idstoIncludeFile=tmp_$chr/\$NAME 
		
	date
EOT
    echo $(date)
	sbatch --wait --export=ALL tmp.sbatch
done

echo $(date)
echo Finished.


	
