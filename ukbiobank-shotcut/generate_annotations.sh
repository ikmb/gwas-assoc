#!/bin/bash

# Extracts annotations from UKB phenotype db and Covid19 table dumps to
# generate covariates for association testing.

# Jan Kässens <j.kaessens@ikmb.uni-kiel.de>

if [ "$#" -ne 4 ]; then
    echo "Syntax: $0 <ukb.fam> <blacklist.csv> <covid19_results.txt> <ukbXXXX.tab>"
    echo "  where:"
    echo "    <ukb.fam>       Current PLINK FAM in the UKB genotype calls folder"
    echo "    <blacklist.csv> Current patient consent withdrawal list"
    echo "    <covid19_res..> UKB table dump for covid19 results"
    echo "    <ukbXXXX.tab>   Path to the tab-separated UKB phenotype database"
    echo "Output files will be generated:"
    echo "   our.covid.covars.txt with id, age, sex, phenotype and PCs in Plink encoding"
    echo "   our.covid.fam with sex and phenotype for further processing with plink"
    echo "   our.covid.pca.evec with PCs in flashpca2-compatible output"
    echo "Plink-compatible encoding in output: Sex (1=Male, 2=Female), Pheno: (1=Control,2=Case)"
    exit 1
fi

FAM=$1
BLACKLIST=$2
COVID=$3
TSV=$4

# Year of birth, plain integer, e.g. 1980
COL_BYEAR=6

# 0=Female, 1=Male
# Self-reported is column 4, genetic sex is column 605
COL_SEX=606

# Column of first PC, 40 PCs in total
COL_PCFIRST=614
NUM_PCS=10

COL_BATCH=605

COL_KINSHIP=656

echo "Parameters:"
echo "  UKB FAM database:  $FAM"
echo "  COVID-19 database: $COVID"
echo "  Phenotypes TSV:    $TSV"
echo "  Number of PCs:     $NUM_PCS"


# First, create an intersection of our sample database and the list of Covid-19 results
echo Scanning samples...
tail -n +2 "$COVID" | cut -f1 | sort | uniq >covid.samples
cut -f1 -d" " <"$FAM" | sort | uniq >our.samples
comm -12 covid.samples our.samples >our.covid
echo "  Found" $(wc -l <our.covid) of $(wc -l <covid.samples) samples in our database.
#rm -f our.samples covid.samples

# Collect Case/Control status from Covid results table
#   Samples might appear multiple times with different status. If it was at least once 1 (positive),
#   the negative status will be overriden.
echo Collecting case/control status...
gawk 'NR==FNR{ids[$1];next} {if($1 in ids && $6=="1"){print $1}}' our.covid "$COVID" | sort | uniq >our.covid.with-status.1
echo "   Found" $(wc -l <our.covid.with-status.1) cases

echo Processing patient consent withdrawal list...
gawk 'NR==FNR{blacklist[$1];next} {if( !($1 in blacklist)) print $0}' $BLACKLIST $TSV >filtered.tab
echo "   UKB samples before: " $(tail -n +2 <$TSV | wc -l) UKB samples after: $(tail -n +2 <filtered.tab | wc -l)
echo Extracting age, sex, batch and $NUM_PCS principal components...

echo -n "IID PHENO SEX AGE BATCH KINSHIP " >our.covid.covars.txt.tmp
for (( i=1; i<=$NUM_PCS; i++ ))
do 
    echo -n "PC$i " >>our.covid.covars.txt.tmp
done

echo >>our.covid.covars.txt.tmp

gawk -v BYEAR=$COL_BYEAR -v SEX=$COL_SEX -v FIRSTPC=$COL_PCFIRST -v NUM_PC=$NUM_PCS -v BATCH=$COL_BATCH -v KINSHIP=$COL_KINSHIP \
    'NR==FNR{status[$1];next}
    FNR>1 {
        covid=1
        if($1 in status) {
            covid=2
        }
        pcs=""
        for(i=FIRSTPC;i<FIRSTPC+NUM_PC;i++) {
            pcs=pcs" "$i
        }
        if(!($1 in blacklist) && $SEX != "NA") {
            # Change to Plink-compatible encoding for sex
            if($SEX=="0")
                $SEX="2"
            print $1,covid,$SEX,(2020-$BYEAR),$BATCH,$KINSHIP,pcs
        }
    }' \
    our.covid.with-status.1 \
    filtered.tab \
    >>our.covid.covars.txt.tmp

gawk 'NR==FNR{codes[$1]=$2;next} FNR==1{print $0;next} {$5=codes[$5]; print $0}' Coding22000.txt our.covid.covars.txt.tmp >our.covid.covars.txt

echo "   Done. Covars for "$(tail -n +2 <our.covid.covars.txt | wc -l)" samples written to our.covid.covars.txt"
echo "Cases left:   " $(cut -f2 -d" " <our.covid.covars.txt | grep 2 | wc -l)
echo "Controls left:" $(cut -f2 -d" " <our.covid.covars.txt | grep 1 | wc -l)
#rm filtered.tab
#rm our.covid.with-status.1 our.covid

echo "Generating our.covid.fam"
gawk 'NR>1{print $1,$1,0,0,$3,$2}' our.covid.covars.txt >our.covid.fam
echo "Generating flashpca2-compatible PC output our.covid.pca.evec"
gawk ''

echo "FID IID " >fidiid
gawk 'NR>1{print $1,$1}' our.covid.covars.txt >>fidiid
cut -f5- -d" " our.covid.covars.txt >pcs
paste fidiid pcs | tr -s ' ' "\t">our.covid.pca.evec

