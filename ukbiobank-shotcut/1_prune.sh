#!/bin/bash

source env

$PLINK --bfile $SOURCE --memory 25000 --threads 8 --indep-pairwise 50 5 0.2 --out prunedata
$PLINK --bfile $SOURCE --memory 25000 --threads 8 --extract prunedata.prune.in --make-bed --out ukb_covid19_pruned
rm prunedata*

