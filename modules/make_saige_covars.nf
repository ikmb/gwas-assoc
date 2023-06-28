process make_saige_covars {
	scratch params.scratch
    publishDir params.output, mode: 'copy'
	label 'base'
    input:
    // file in_covars
    path(evec)
    path(inc_fam)

    output:
    tuple path("${params.collection_name}.covars"), path("${params.collection_name}.covar_cols"), emit: covars/// into for_saige_null_covars, for_merge_sumstats_covars, for_plink_covars
    path("${params.collection_name}.double-id.fam"), emit: plink_fam// into for_plink_fam
	tuple path("${params.collection_name}.covars"), path("${params.collection_name}.covar_cols"), path("${params.collection_name}.double-id.fam"), emit: for_plink

    shell:
'''
# Re-format sample ID and family ID to VCF rules
gawk '{if($1!="0") {$2=$1"_"$2; $1="0";} print $0}' !{inc_fam}  | sort >new-fam
# Also re-format, but keep evec header
# mawk 'NR==1{print $0;next} {if($1!="0") {$2=$1"_"$2; $1="0";} print $0}' {in_covars}  | sort >evec
# Take evec file as a whole, SAIGE does not care about additional columns.
# Filter evec file according to samples contained in FAM (if not already done)
gawk 'FNR==NR{samples[$2];next} {if($2 in samples) { for(i=1;i<=(12);i++) {printf "%s%s", $i, (i<12?OFS:ORS)}}}' new-fam !{evec} >filtered-evec
if [ -f "!{params.more_covars}" ]; then
    gawk 'FNR==NR{samples[$2];next} {if($2 in samples) { for(i=1;i<=(NF);i++) {printf "%s%s", $i, (i<NF?OFS:ORS)}}}' new-fam !{params.more_covars} | sort >filtered-covars
    
    # identify column indices of selected covars
    # NOTE: We checked above that all selected columns are available, so the list won't be empty
    COLSTR=$(head -n1 "!{params.more_covars}"  | gawk 'BEGIN { split("'"!{params.more_covars_cols}"'",cols,","); colstr="" } { for (c in cols) { for (i=3; i<=NF; i++) { if ($i == cols[c]) {colstr=colstr "," i} }}} END { print substr(colstr,2) }')
    
    # extract header
    cut -f$COLSTR -d" " !{params.more_covars} | awk 'NR<=1' >covars-column
    # fill values
    cut -f$COLSTR -d" " filtered-covars >>covars-column
fi
EVEC_LINES=$(wc -l <filtered-evec)
FAM_LINES=$(wc -l <new-fam)
if [ $EVEC_LINES -ne $FAM_LINES ]; then
    echo I could not find covariates in !{evec} for every sample specified in !{inc_fam}. Please check.
    echo Aborting.
    exit 1
fi
#echo "FID IID PC1 PC2 PC3 PC4 PC5 PC6 PC7 PC8 PC9 PC10" >evec.double-id.withheader
echo "FID IID" PC{1..!{params.pca_dims}} >evec.double-id.withheader

cat filtered-evec >>evec.double-id.withheader
# Take phenotype info from FAM, translate to SAIGE-encoding
echo "Pheno" >pheno-column
<new-fam  gawk '{if($6=="2") {$6="1";} else if($6=="1") {$6="0";} print $6}' >>pheno-column
mv new-fam !{params.collection_name}.double-id.fam
# Merge both, replace space runs with single tabs for SAIGE
touch covars-column
if [ ! -d "!{params.more_covars}" ]; then
    echo PC{1..10}, !{params.more_covars_cols} | sed 's/\\ ,/\\,/g' | tr ' ' ,  | sed 's/\\,\\,/\\,/g' >!{params.collection_name}.covar_cols
else
    echo PC{1..!{params.pca_dims}} | tr ' ' , >!{params.collection_name}.covar_cols
fi
paste -d" " evec.double-id.withheader pheno-column covars-column | tr -s ' ' \\\\t >!{params.collection_name}.covars
'''
}

