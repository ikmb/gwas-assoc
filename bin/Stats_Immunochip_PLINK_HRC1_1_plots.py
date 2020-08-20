#!/usr/bin/python
import os
from os.path import *
import string
import re
import gzip
import math
import decimal
import datetime
import subprocess
from os import listdir

import sys

sys.path.append(join(sys.path[0], "../lib/all_scripts"))
from all_common import *

sys.path.append(join(sys.path[0], "../lib/lib"))
from plink_classes import *




######################################
#
#  functions
#
######################################

def aliasing(sample_with_params):
    """ find shorter variable names for parsed parameter list """

    # parameters
    global batch_mode
    global path_source_code_env
    global path_project_snpchip_data_env
    global path_project_collections_env

    # parameters
    global collection_name
    global disease_data_set_prefix_release
    global disease_data_set_suffix_release_imputed
    global disease_data_set_prefix_release_statistics
    global PLINK_dir
    global Imputation_orig_dir

    global QQplot_null_variants
    global QQplot_SNPexcludeList
    global QQplot_SNPexcludeList_exists
    global Multiallelic_SNPexcludeList
    global Multiallelic_SNPexcludeList_exists

    global maf
    global clumpp1
    global clumpp2
    global clumpr2
    global clumpkb
    global ldwindowkb
    global ldwindow
    global ldwindowr2

    global numofPCs
    global gwascatalog
    global numoftophits
    global build

    global timestamp

    # True or False
    batch_mode = False
    if sample_with_params['batch_mode'] == "True":
        batch_mode = True

    timestamp = "tmp" + "".join(str(x) for x in datetime.datetime.now().timetuple())

    # paths from bash env variables
    path_source_code_env = os.getenv(sample_with_params['path_source_code_env'])
    path_project_snpchip_data_env = os.getenv(sample_with_params['path_project_snpchip_data_env'])
    path_project_collections_env  = os.getenv(sample_with_params['path_project_collections_env'])

    # name of the new collection of batches
    collection_name              = sample_with_params['collection_name']

    # QCed release
    disease_data_set_prefix_release = join(path_project_collections_env,\
        collection_name,\
        sample_with_params['disease_data_set_dir_release'],\
        sample_with_params['disease_data_set_prefix_release'])
    if not os.path.isfile(disease_data_set_prefix_release + ".fam"):
        print >> sys.stderr, "abort: file \"" + disease_data_set_prefix_release + ".fam" + "\" does not exist!"
        sys.exit(1)

    disease_data_set_prefix_release_statistics = join(path_project_collections_env,\
        collection_name,\
        sample_with_params['disease_data_set_dir_release'],\
        sample_with_params['disease_data_set_prefix_release_statistics'])
    if not os.path.isfile(disease_data_set_prefix_release_statistics + ".ped"):
        print >> sys.stderr, "abort: file \"" + disease_data_set_prefix_release_statistics + ".ped" + "\" does not exist!"
        sys.exit(1)

    ## Output original dir for HRC imputation dosage files from Imputation server
    #Imputation_orig_dir = join(path_project_collections_env,\
    #    collection_name,\
    #    sample_with_params['disease_data_set_dir_release'],\
    #    "HRC1.1_imputation")

    # TODO aendere pfad zum imputation dir
    # Output original dir for HRC imputation dosage files from Imputation server
    Imputation_orig_dir = join(path_project_collections_env,\
        collection_name,\
        "HRC1.1_imputation")

    disease_data_set_suffix_release_imputed = sample_with_params['disease_data_set_suffix_release_imputed']
    if not os.path.isfile(join(Imputation_orig_dir, "1."+ disease_data_set_suffix_release_imputed +".PLINKdosage.gz")):
        print >> sys.stderr, "abort: file \"" + join(Imputation_orig_dir, "1."+ disease_data_set_suffix_release_imputed +".PLINKdosage.gz")+ "\" does not exist!"
        sys.exit(1)

    # Output dir for PLINK
    PLINK_dir = join(path_project_collections_env,\
        collection_name,\
        "PLINK_HRC1.1_" + basename(disease_data_set_prefix_release_statistics))
    os.system("mkdir -p %s" %(PLINK_dir))

    QQplot_null_variants = join(path_source_code_env, "annotation", "QQplots", sample_with_params['QQplot_null_variants'])
    if not os.path.isfile(QQplot_null_variants):
        print >> sys.stderr, "abort: file \"" +QQplot_null_variants+ "\" does not exist!"
        sys.exit(1)

    # Optional: remove variants when calculating null SNPs
    if os.path.isfile(join(path_project_snpchip_data_env, str(sample_with_params['QQplot_SNPexcludeList']))):
        QQplot_SNPexcludeList = join(path_project_snpchip_data_env, str(sample_with_params['QQplot_SNPexcludeList']))
        QQplot_SNPexcludeList_exists = True
    else:
        QQplot_SNPexcludeList_exists = False

    # Optional: remove variants from imputation
    if os.path.isfile(join(Imputation_orig_dir, str(sample_with_params['Multiallelic_SNPexcludeList']))):
        Multiallelic_SNPexcludeList = join(Imputation_orig_dir, str(sample_with_params['Multiallelic_SNPexcludeList']))
        Multiallelic_SNPexcludeList_exists = True
    else:
        Multiallelic_SNPexcludeList_exists = False

    maf = float(sample_with_params['maf'])
    clumpp1 = float(sample_with_params['clumpp1'])
    clumpp2 = float(sample_with_params['clumpp2'])
    clumpr2 = float(sample_with_params['clumpr2'])
    clumpkb = float(sample_with_params['clumpkb'])

    # LD parameters
    ldwindowkb = str(sample_with_params['ldwindowkb'])
    ldwindow   = str(sample_with_params['ldwindow'])
    ldwindowr2 = str(sample_with_params['ldwindowr2'])

    numofPCs = int(sample_with_params['numofPCs'])

    # gwascatalog with known hits
    gwascatalog                  = str(sample_with_params['gwascatalog'])
    # number of hits to be plotted
    numoftophits = str(sample_with_params['numoftophits'])
    # NCBI build
    build = str(sample_with_params['build'])




def determine_numof_cases_controls(prefix):
    """ determine number of cases and controls from assoc log file """

    try:
        fh_r = file(prefix + ".log", "r")
    except IOError, e:
        print e
        sys.exit(1)

    numof_cases    = 1000
    numof_controls = 1000

    special_line_pattern = re.compile("^Among remaining phenotypes, .* are cases and .* are controls.*$")
    line = fh_r.readline()
    while line:

        if(special_line_pattern.search(line)):
            list = re.split("\s+",line)
            numof_cases    = int(list[3])
            numof_controls = int(list[7])
            break

        line = fh_r.readline()

    fh_r.close()
    return (numof_cases, numof_controls)




def extractQQplot_null_variants(assoc_input, assoc_output):
    """ extract null variants from assoc file """

    try:
        fh_r = file(QQplot_null_variants, "r")
    except IOError, e:
      print e
      sys.exit(1)

    null = {}
    line = fh_r.readline().rstrip('\n')
    while line:
        list = re.split("\s+",line)
        rs   = list[0]
        null[rs] = True
        line = fh_r.readline().rstrip('\n')

    fh_r.close()

    try:
        fh_r = file(assoc_input, "r")
        fh_w = file(assoc_output, "w")
    except IOError, e:
      print e
      sys.exit(1)

    line = fh_r.readline().rstrip('\n')
    fh_w.writelines(line + "\n")
    line = fh_r.readline().rstrip('\n')
    while line:
        list = re.split("\s+",line)
        if list[0] == "":
            del list[0]
        rs = list[1]
        if null.has_key(rs):
            fh_w.writelines(line + "\n")
        line = fh_r.readline().rstrip('\n')

    fh_r.close()
    fh_w.close()




def excludeQQplot_variants(assoc_input, assoc_output):
    """ exclude specific regions from assoc null file """

    try:
        fh_r = file(QQplot_SNPexcludeList, "r")
    except IOError, e:
      print e
      sys.exit(1)

    exclude_regions = []
    line = fh_r.readline().rstrip('\n')
    while line:
        list = re.split("\s+",line)
        chr    = int(list[0])
        left   = int(list[1])
        right  = int(list[2])
        exclude_regions.append((chr, left, right))
        line = fh_r.readline().rstrip('\n')

    fh_r.close()

    try:
        fh_r = file(assoc_input, "r")
        fh_w = file(assoc_output, "w")
    except IOError, e:
      print e
      sys.exit(1)

    line = fh_r.readline().rstrip('\n')
    fh_w.writelines(line + "\n")
    line = fh_r.readline().rstrip('\n')
    while line:
        list = re.split("\s+",line)
        if list[0] == "":
            del list[0]
        chr = int(list[0])
        pos = int(list[2])
        exclude_this_region = False
        for region in exclude_regions:
            if chr == region[0] and pos >= region[1] and pos <= region[2]:
                exclude_this_region = True
                break
        if exclude_this_region:
            pass
        else:
            fh_w.writelines(line + "\n")
        line = fh_r.readline().rstrip('\n')

    fh_r.close()
    fh_w.close()




def qqplot_null(file_PLINK_null, numof_cases, numof_controls):
    """ extract null variants """

    # use normal qqplot scripts instead of null qqplot scripts because usage of
    # ld pruned variants as null SNPs
    print "\n    qq-plotpval for null SNVs from association results ...\n\n"
    #os.system("/home/sukmb113/progs/R-2.12.1/bin/R --slave --args %s %s %s %s < %s" \
    os.system("/home/sukmb113/progs/R-2.12.1/bin/R --no-save --args %s %s %s %s < %s" \
        %(file_PLINK_null,\
        str(numof_cases),\
        str(numof_controls),\
        collection_name,\
        join(path_source_code_env, "sources", "qqplotP.r")))
        #join(path_source_code_env, "sources", "qqplotP_null-snps.r")))

    print "\n    qq-plotpval for null SNVs from association results ...\n\n"
    #os.system("/home/sukmb113/progs/R-2.12.1/bin/R --slave --args %s %s %s %s < %s" \
    os.system("/home/sukmb113/progs/R-2.12.1/bin/R --no-save --args %s %s %s %s < %s" \
        %(file_PLINK_null,\
        str(numof_cases),\
        str(numof_controls),\
        collection_name,\
        join(path_source_code_env, "sources", "qqplotP_v2.r")))

    print "\n    qq-plot for null SNVs from association results ...\n\n"
    #os.system("/home/sukmb113/progs/R-2.12.1/bin/R --slave --args %s %s %s %s < %s" \
    os.system("/home/sukmb113/progs/R-2.12.1/bin/R --no-save --args %s %s %s %s < %s" \
        %(file_PLINK_null,\
        str(numof_cases),\
        str(numof_controls),\
        collection_name,\
        join(path_source_code_env, "sources", "qqplot_assoc_pval2CHISQ.r")))

    if QQplot_SNPexcludeList_exists:

        excludeQQplot_variants(assoc_input=file_PLINK_null, \
                               assoc_output=file_PLINK_null + "_excludeRegions")

        print "\n    qq-plotpval for null SNVs with specific regions excluded from association results ...\n\n"
        #os.system("/home/sukmb113/progs/R-2.12.1/bin/R --slave --args %s %s %s %s < %s" \
        os.system("/home/sukmb113/progs/R-2.12.1/bin/R --no-save --args %s %s %s %s < %s" \
            %(file_PLINK_null + "_excludeRegions",\
            str(numof_cases),\
            str(numof_controls),\
            collection_name,\
            join(path_source_code_env, "sources", "qqplotP.r")))
            #join(path_source_code_env, "sources", "qqplotP_null-snps.r")))

        os.system("/home/sukmb113/progs/R-2.12.1/bin/R --no-save --args %s %s %s %s < %s" \
            %(file_PLINK_null + "_excludeRegions",\
            str(numof_cases),\
            str(numof_controls),\
            collection_name,\
            join(path_source_code_env, "sources", "qqplotP_v2.r")))




def qqplots_manhattanplots(file_PLINK_geno_imp_rsq0_4, file_PLINK_geno_imp_rsq0_8):
    """ generate qq plots and manhattan plots """

    file_PLINK_geno_imp_rsq0_4_null = file_PLINK_geno_imp_rsq0_4 +".null"
    file_PLINK_geno_imp_rsq0_8_null = file_PLINK_geno_imp_rsq0_8 +".null"

    file_PLINK_geno_imp_rsq0_4_noxMHC = file_PLINK_geno_imp_rsq0_4 +".noxMHC"
    file_PLINK_geno_imp_rsq0_8_noxMHC = file_PLINK_geno_imp_rsq0_8 +".noxMHC"

    numof_cases, numof_controls = \
        determine_numof_cases_controls(prefix=join(PLINK_dir, basename(disease_data_set_prefix_release_statistics)))
    print(numof_cases)
    print(numof_controls)

    # ------------------- #
    # -- null variants -- #
    # ------------------- #
    # Too few null variants available for QQplots because only high density
    # regions and MHC/KIR were extracted here

    extractQQplot_null_variants(assoc_input=file_PLINK_geno_imp_rsq0_4, \
                                assoc_output=file_PLINK_geno_imp_rsq0_4_null)
    qqplot_null(file_PLINK_null=file_PLINK_geno_imp_rsq0_4_null,\
                numof_cases=numof_cases,\
                numof_controls=numof_controls)

    extractQQplot_null_variants(assoc_input=file_PLINK_geno_imp_rsq0_8, \
                                assoc_output=file_PLINK_geno_imp_rsq0_8_null)
    qqplot_null(file_PLINK_null=file_PLINK_geno_imp_rsq0_8_null,\
                numof_cases=numof_cases,\
                numof_controls=numof_controls)

    # ---------------------------- #
    # -- filter for non-xMHC    -- #
    # -- xMHC: chr6, [25,34[ Mb -- #
    # ---------------------------- #
    os.system("gawk '{ if (!($1 == 6 && ($3 >= 25000000 && $3 < 34000000))) print }' %s > %s"\
              %(file_PLINK_geno_imp_rsq0_4, file_PLINK_geno_imp_rsq0_4_noxMHC))
    qqplot_xMHC_noxMHC(file_PLINK=file_PLINK_geno_imp_rsq0_4,\
                       file_PLINK_noxMHC=file_PLINK_geno_imp_rsq0_4_noxMHC,\
                       numof_cases=numof_cases,\
                       numof_controls=numof_controls)
    os.system("rm -f %s" %(file_PLINK_geno_imp_rsq0_4_noxMHC))

    os.system("gawk '{ if (!($1 == 6 && ($3 >= 25000000 && $3 < 34000000))) print }' %s > %s"\
              %(file_PLINK_geno_imp_rsq0_8, file_PLINK_geno_imp_rsq0_8_noxMHC))
    qqplot_xMHC_noxMHC(file_PLINK=file_PLINK_geno_imp_rsq0_8,\
                       file_PLINK_noxMHC=file_PLINK_geno_imp_rsq0_8_noxMHC,\
                       numof_cases=numof_cases,\
                       numof_controls=numof_controls)
    os.system("rm -f %s" %(file_PLINK_geno_imp_rsq0_8_noxMHC))

    # -------------------------------- #
    # -- manhattan plot for chr1-22 -- #
    # -------------------------------- #

    print "\n    manhattan-plot from association results ...\n\n"

    os.system("gawk ' $1 !~ /23/ { print $0 } ' %s > %s_tmp1"\
              %(file_PLINK_geno_imp_rsq0_4, file_PLINK_geno_imp_rsq0_4))
    os.system("gawk ' $1 !~ /24/ { print $0 } ' %s_tmp1 > %s_tmp2"\
              %(file_PLINK_geno_imp_rsq0_4, file_PLINK_geno_imp_rsq0_4))
    os.system("gawk ' $1 !~ /25/ { print $0 } ' %s_tmp2 > %s_tmp3"\
              %(file_PLINK_geno_imp_rsq0_4, file_PLINK_geno_imp_rsq0_4))
    os.system("gawk ' $1 !~ /26/ { print $0 } ' %s_tmp3 > %s_noChr23"\
              %(file_PLINK_geno_imp_rsq0_4, file_PLINK_geno_imp_rsq0_4))
    os.system("rm -f %s_tmp*" %(file_PLINK_geno_imp_rsq0_4))

    os.system("R --slave --args %s_noChr23 %s < %s" \
        %(file_PLINK_geno_imp_rsq0_4,\
        collection_name,\
        join(path_source_code_env, "sources", "manhattan.r")))
    os.system("R --slave --args %s_noChr23 %s < %s" \
        %(file_PLINK_geno_imp_rsq0_4,\
        collection_name,\
        join(path_source_code_env, "sources", "manhattan2.r")))

    os.system("gawk ' $1 !~ /23/ { print $0 } ' %s > %s_tmp1"\
              %(file_PLINK_geno_imp_rsq0_8, file_PLINK_geno_imp_rsq0_8))
    os.system("gawk ' $1 !~ /24/ { print $0 } ' %s_tmp1 > %s_tmp2"\
              %(file_PLINK_geno_imp_rsq0_8, file_PLINK_geno_imp_rsq0_8))
    os.system("gawk ' $1 !~ /25/ { print $0 } ' %s_tmp2 > %s_tmp3"\
              %(file_PLINK_geno_imp_rsq0_8, file_PLINK_geno_imp_rsq0_8))
    os.system("gawk ' $1 !~ /26/ { print $0 } ' %s_tmp3 > %s_noChr23"\
              %(file_PLINK_geno_imp_rsq0_8, file_PLINK_geno_imp_rsq0_8))
    os.system("rm -f %s_tmp*" %(file_PLINK_geno_imp_rsq0_8))

    os.system("R --slave --args %s_noChr23 %s < %s" \
        %(file_PLINK_geno_imp_rsq0_8,\
        collection_name,\
        join(path_source_code_env, "sources", "manhattan.r")))
    os.system("R --slave --args %s_noChr23 %s < %s" \
        %(file_PLINK_geno_imp_rsq0_8,\
        collection_name,\
        join(path_source_code_env, "sources", "manhattan2.r")))

    # -- excluderegions -- #
    os.system("gawk ' $1 !~ /23/ { print $0 } ' %s > %s_tmp1"\
              %(file_PLINK_geno_imp_rsq0_4 +"_excludeRegions", file_PLINK_geno_imp_rsq0_4 +"_excludeRegions"))
    os.system("gawk ' $1 !~ /24/ { print $0 } ' %s_tmp1 > %s_tmp2"\
              %(file_PLINK_geno_imp_rsq0_4 +"_excludeRegions", file_PLINK_geno_imp_rsq0_4 +"_excludeRegions"))
    os.system("gawk ' $1 !~ /25/ { print $0 } ' %s_tmp2 > %s_tmp3"\
              %(file_PLINK_geno_imp_rsq0_4 +"_excludeRegions", file_PLINK_geno_imp_rsq0_4 +"_excludeRegions"))
    os.system("gawk ' $1 !~ /26/ { print $0 } ' %s_tmp3 > %s_noChr23"\
              %(file_PLINK_geno_imp_rsq0_4 +"_excludeRegions", file_PLINK_geno_imp_rsq0_4 +"_excludeRegions"))
    os.system("rm -f %s_tmp*" %(file_PLINK_geno_imp_rsq0_4 +"_excludeRegions"))

    os.system("R --slave --args %s_noChr23 %s < %s" \
        %(file_PLINK_geno_imp_rsq0_4 +"_excludeRegions",\
        collection_name,\
        join(path_source_code_env, "sources", "manhattan.r")))
    os.system("R --slave --args %s_noChr23 %s < %s" \
        %(file_PLINK_geno_imp_rsq0_4 +"_excludeRegions",\
        collection_name,\
        join(path_source_code_env, "sources", "manhattan2.r")))

    os.system("gawk ' $1 !~ /23/ { print $0 } ' %s > %s_tmp1"\
              %(file_PLINK_geno_imp_rsq0_8 +"_excludeRegions", file_PLINK_geno_imp_rsq0_8 +"_excludeRegions"))
    os.system("gawk ' $1 !~ /24/ { print $0 } ' %s_tmp1 > %s_tmp2"\
              %(file_PLINK_geno_imp_rsq0_8 +"_excludeRegions", file_PLINK_geno_imp_rsq0_8 +"_excludeRegions"))
    os.system("gawk ' $1 !~ /25/ { print $0 } ' %s_tmp2 > %s_tmp3"\
              %(file_PLINK_geno_imp_rsq0_8 +"_excludeRegions", file_PLINK_geno_imp_rsq0_8 +"_excludeRegions"))
    os.system("gawk ' $1 !~ /26/ { print $0 } ' %s_tmp3 > %s_noChr23"\
              %(file_PLINK_geno_imp_rsq0_8 +"_excludeRegions", file_PLINK_geno_imp_rsq0_8 +"_excludeRegions"))
    os.system("rm -f %s_tmp*" %(file_PLINK_geno_imp_rsq0_8 +"_excludeRegions"))

    os.system("R --slave --args %s_noChr23 %s < %s" \
        %(file_PLINK_geno_imp_rsq0_8 +"_excludeRegions",\
        collection_name,\
        join(path_source_code_env, "sources", "manhattan.r")))
    os.system("R --slave --args %s_noChr23 %s < %s" \
        %(file_PLINK_geno_imp_rsq0_8 +"_excludeRegions",\
        collection_name,\
        join(path_source_code_env, "sources", "manhattan2.r")))




def extract_Rsq_variants(assoc_logistic_input, assoc_dosage_input, assoc_merge_output_rsq0_4, assoc_merge_output_rsq0_8):
    """ extract imputed variants based on Rsq value """

    try:
        fh_w1 = file(assoc_merge_output_rsq0_8, "w")
        fh_w2 = file(assoc_merge_output_rsq0_4, "w")
    except IOError, e:
        print e
        sys.exit(1)

    store_lines_rsq0_4 = []
    store_lines_rsq0_8 = []

    snp2line = {}

    # -- assoc_logistic_input - all genotyped marker -- #
    try:
        fh_r = file(assoc_logistic_input, "r")
    except IOError, e:
        print e
        sys.exit(1)
    # header
    line = fh_r.readline().rstrip('\n')
    line = fh_r.readline().rstrip('\n')
    while line:
        list = re.split("\s+",line)
        if list[0] == "":
            del list[0]
        if list[-1] == "":
            del list[-1]

        # rs numbers instead of chr:pos in column SNP
        snp2line[list[0]+":"+list[2]] = line
        chr  = decimal.Decimal(list[0])
        pos  = decimal.Decimal(list[2])
        store_lines_rsq0_4.append( (chr, pos, line) )
        store_lines_rsq0_8.append( (chr, pos, line) )
        line = fh_r.readline().rstrip('\n')
    fh_r.close()

    # -- assoc_dosage_input - all imputed marker -- #
    try:
        fh_r = file(assoc_dosage_input, "r")
    except IOError, e:
        print e
        sys.exit(1)

    # header
    header = fh_r.readline().rstrip('\n')
    fh_w1.writelines(header+"\n")
    fh_w2.writelines(header+"\n")
    line = fh_r.readline().rstrip('\n')
    while line:
        list = re.split("\s+",line)
        if list[0] == "":
            del list[0]
        if list[-1] == "":
            del list[-1]
        # if not already genotyped --> add
        if not snp2line.has_key(list[0]+":"+list[2]):
            if float(list[7]) >= 0.4:
                chr  = decimal.Decimal(list[0])
                pos  = decimal.Decimal(list[2])
                store_lines_rsq0_4.append( (chr, pos, line) )
                if float(list[7]) >= 0.8:
                    store_lines_rsq0_8.append( (chr, pos, line) )
        line = fh_r.readline().rstrip('\n')
    fh_r.close()

    # sort merged files by chr, pos
    s_rsq0_4 = sorted(store_lines_rsq0_4, key=lambda tupel: tupel[1])
    t_rsq0_4 = sorted(s_rsq0_4, key=lambda tupel: tupel[0], reverse=False)
    s_rsq0_8 = sorted(store_lines_rsq0_8, key=lambda tupel: tupel[1])
    t_rsq0_8 = sorted(s_rsq0_8, key=lambda tupel: tupel[0], reverse=False)

    for tupel in t_rsq0_4:
        line = tupel[2]
        fh_w2.writelines(line+"\n")

    for tupel in t_rsq0_8:
        line = tupel[2]
        fh_w1.writelines(line+"\n")

    fh_w1.close()
    fh_w2.close()

    # save another version with chr:pos as SNPids for Locuszoom
    try:
        fh_w1 = file(assoc_merge_output_rsq0_8+".locuszoom", "w")
        fh_w2 = file(assoc_merge_output_rsq0_4+".locuszoom", "w")
    except IOError, e:
        print e
        sys.exit(1)

    fh_w1.writelines(header+"\n")
    fh_w2.writelines(header+"\n")

    check_duplicates_rsq0_4 = {}
    check_duplicates_rsq0_8 = {}

    for tupel in t_rsq0_4:
        line = tupel[2]
        list = re.split("\s+",line)
        if not check_duplicates_rsq0_4.has_key(list[0]+":"+list[2]):
            fh_w2.writelines(list[0]+"\t"+\
                         list[0]+":"+list[2] +"\t"+\
                         list[2]+"\t"+\
                         list[3]+"\t"+\
                         list[4]+"\t"+\
                         list[5]+"\t"+\
                         list[6]+"\t"+\
                         list[7]+"\t"+\
                         list[8]+"\t"+\
                         list[9]+"\t"+\
                         list[10]+"\n")
            check_duplicates_rsq0_4[list[0]+":"+list[2]] = True

    for tupel in t_rsq0_8:
        line = tupel[2]
        list = re.split("\s+",line)
        if not check_duplicates_rsq0_8.has_key(list[0]+":"+list[2]):
            fh_w1.writelines(list[0]+"\t"+\
                         list[0]+":"+list[2] +"\t"+\
                         list[2]+"\t"+\
                         list[3]+"\t"+\
                         list[4]+"\t"+\
                         list[5]+"\t"+\
                         list[6]+"\t"+\
                         list[7]+"\t"+\
                         list[8]+"\t"+\
                         list[9]+"\t"+\
                         list[10]+"\n")
            check_duplicates_rsq0_8[list[0]+":"+list[2]] = True

    fh_w1.close()
    fh_w2.close()




def qqplot_xMHC_noxMHC(file_PLINK, file_PLINK_noxMHC, numof_cases, numof_controls):
    """ extract variants based on Rsq value """

    # ---------------------------- #
    # -- filtered for non-xMHC  -- #
    # -- xMHC: chr6, [25,34[ Mb -- #
    # ---------------------------- #

    print "\n    qq-plot from association results ...\n\n"
    os.system("/home/sukmb113/progs/R-2.12.1/bin/R --slave --args %s %s %s %s < %s" \
        %(file_PLINK_noxMHC,\
        str(numof_cases),\
        str(numof_controls),\
        collection_name,\
        join(path_source_code_env, "sources", "qqplot_assoc_pval2CHISQ.r")))

    print "\n    another qq-plotpval from association results ...\n\n"
    os.system("/home/sukmb113/progs/R-2.12.1/bin/R --slave --args %s %s %s %s < %s" \
        %(file_PLINK_noxMHC,\
        str(numof_cases),\
        str(numof_controls),\
        collection_name,\
        join(path_source_code_env, "sources", "qqplotP.r")))

    print "\n    another qq-plotpval from association results ...\n\n"
    os.system("/home/sukmb113/progs/R-2.12.1/bin/R --slave --args %s %s %s %s < %s" \
        %(file_PLINK_noxMHC,\
        str(numof_cases),\
        str(numof_controls),\
        collection_name,\
        join(path_source_code_env, "sources", "qqplotP_v2.r")))

    # -------------------------------------- #
    # -- including xMHC: chr6, [25,34[ Mb -- #
    # -------------------------------------- #

    print "\n    qq-plot from association results ...\n\n"
    os.system("/home/sukmb113/progs/R-2.12.1/bin/R --slave --args %s %s %s %s < %s" \
        %(file_PLINK,\
        str(numof_cases),\
        str(numof_controls),\
        collection_name,\
        join(path_source_code_env, "sources", "qqplot_assoc_pval2CHISQ.r")))

    print "\n    another qq-plotpval from association results ...\n\n"
    os.system("/home/sukmb113/progs/R-2.12.1/bin/R --slave --args %s %s %s %s < %s" \
        %(file_PLINK,\
        str(numof_cases),\
        str(numof_controls),\
        collection_name,\
        join(path_source_code_env, "sources", "qqplotP.r")))

    print "\n    another qq-plotpval from association results ...\n\n"
    os.system("/home/sukmb113/progs/R-2.12.1/bin/R --slave --args %s %s %s %s < %s" \
        %(file_PLINK,\
        str(numof_cases),\
        str(numof_controls),\
        collection_name,\
        join(path_source_code_env, "sources", "qqplotP_v2.r")))

    if QQplot_SNPexcludeList_exists:

        excludeQQplot_variants(assoc_input=file_PLINK, \
                               assoc_output=file_PLINK + "_excludeRegions")

        excludeQQplot_variants(assoc_input=file_PLINK_noxMHC, \
                               assoc_output=file_PLINK_noxMHC + "_excludeRegions")

        # ---------------------------- #
        # -- filtered for non-xMHC  -- #
        # -- xMHC: chr6, [25,34[ Mb -- #
        # ---------------------------- #

        print "\n    qq-plot from association results ...\n\n"
        os.system("/home/sukmb113/progs/R-2.12.1/bin/R --slave --args %s %s %s %s < %s" \
            %(file_PLINK_noxMHC + "_excludeRegions",\
            str(numof_cases),\
            str(numof_controls),\
            collection_name,\
            join(path_source_code_env, "sources", "qqplot_assoc_pval2CHISQ.r")))

        print "\n    another qq-plotpval from association results ...\n\n"
        os.system("/home/sukmb113/progs/R-2.12.1/bin/R --slave --args %s %s %s %s < %s" \
            %(file_PLINK_noxMHC + "_excludeRegions",\
            str(numof_cases),\
            str(numof_controls),\
            collection_name,\
            join(path_source_code_env, "sources", "qqplotP.r")))

        print "\n    another qq-plotpval from association results ...\n\n"
        os.system("/home/sukmb113/progs/R-2.12.1/bin/R --slave --args %s %s %s %s < %s" \
            %(file_PLINK_noxMHC + "_excludeRegions",\
            str(numof_cases),\
            str(numof_controls),\
            collection_name,\
            join(path_source_code_env, "sources", "qqplotP_v2.r")))

        os.system("rm -f %s" %(file_PLINK_noxMHC + "_excludeRegions"))

        # -------------------------------------- #
        # -- including xMHC: chr6, [25,34[ Mb -- #
        # -------------------------------------- #

        print "\n    qq-plot from association results ...\n\n"
        os.system("/home/sukmb113/progs/R-2.12.1/bin/R --slave --args %s %s %s %s < %s" \
            %(file_PLINK + "_excludeRegions",\
            str(numof_cases),\
            str(numof_controls),\
            collection_name,\
            join(path_source_code_env, "sources", "qqplot_assoc_pval2CHISQ.r")))

        print "\n    another qq-plotpval from association results ...\n\n"
        os.system("/home/sukmb113/progs/R-2.12.1/bin/R --slave --args %s %s %s %s < %s" \
            %(file_PLINK + "_excludeRegions",\
            str(numof_cases),\
            str(numof_controls),\
            collection_name,\
            join(path_source_code_env, "sources", "qqplotP.r")))

        print "\n    another qq-plotpval from association results ...\n\n"
        os.system("/home/sukmb113/progs/R-2.12.1/bin/R --slave --args %s %s %s %s < %s" \
            %(file_PLINK + "_excludeRegions",\
            str(numof_cases),\
            str(numof_controls),\
            collection_name,\
            join(path_source_code_env, "sources", "qqplotP_v2.r")))




def generate_merge_file(file_PLINK, merge_file):
    """ generate merge file for all batches """

    try:
        fh_w  = file(merge_file, "w")
    except IOError, e:
        print e
        sys.exit(1)

    sep = " "
    # leave first data set out: chr1
    for i in xrange(2,23):

        bed = file_PLINK +".chr" +str(i)+ ".bed"
        bim = file_PLINK +".chr" +str(i)+ ".bim"
        fam = file_PLINK +".chr" +str(i)+ ".fam"

        fh_w.writelines(bed +sep+ bim +sep+ fam +"\n")

    fh_w.close()




def check_number_of_individuals_pre_and_post_imputation(file_PLINKdosage):
    """ """

    p = subprocess.Popen("wc -l " +disease_data_set_prefix_release_statistics + ".fam | cut -d \" \" -f 1", shell=True, stdout=subprocess.PIPE)
    count_1 = p.stdout.read().strip()
    p = subprocess.Popen("wc -l " +file_PLINKdosage +".rsq0.4.chr1-22.fam | cut -d \" \" -f 1", shell=True, stdout=subprocess.PIPE)
    count_2 = p.stdout.read().strip()
    if not (count_1 == count_2):
        print >> sys.stderr, "abort: different number of individuals in fam file pre versus post imputation."
        print >> sys.stderr, "    Expected " +str(count_1)+ " individuals in " +disease_data_set_prefix_release_statistics + ".fam"
        print >> sys.stderr, "    Observed " +str(count_2)+ " individuals in " +file_PLINKdosage +".rsq0.4.chr1-22.fam"
        sys.exit(1)




def DeFinetti_plots(file_PLINK_geno_imp):
    """ """

    print "\n    do Hardy-Weinberg test ...\n\n"

    os.system("${HOME}/progs/wdist/wdist_linux_x86_64_9Nov2016/plink --bfile %s --hardy --out %s_hardy --hwe 0.0 --threads 16" \
        %(file_PLINK_geno_imp +".rsq0.4.chr1-22",\
          file_PLINK_geno_imp +".rsq0.4.chr1-22"))

    print "\n    plot DeFinetti diagram from genotype counts ...\n\n"
    os.system("R --slave --args %s_hardy.hwe %s_controls_DeFinetti %s_cases_DeFinetti %s_cases_controls_DeFinetti < %s" \
        %(file_PLINK_geno_imp +".rsq0.4.chr1-22",\
        file_PLINK_geno_imp +".rsq0.4.chr1-22",\
        file_PLINK_geno_imp +".rsq0.4.chr1-22",\
        file_PLINK_geno_imp +".rsq0.4.chr1-22",\
        join(path_source_code_env, "sources", "DeFinetti_hardy.r")))

    os.system("${HOME}/progs/wdist/wdist_linux_x86_64_9Nov2016/plink --bfile %s --hardy --out %s_hardy --hwe 0.0 --threads 16" \
        %(file_PLINK_geno_imp +".rsq0.8.chr1-22",\
          file_PLINK_geno_imp +".rsq0.8.chr1-22"))

    print "\n    plot DeFinetti diagram from genotype counts ...\n\n"
    os.system("R --slave --args %s_hardy.hwe %s_controls_DeFinetti %s_cases_DeFinetti %s_cases_controls_DeFinetti < %s" \
        %(file_PLINK_geno_imp +".rsq0.8.chr1-22",\
        file_PLINK_geno_imp +".rsq0.8.chr1-22",\
        file_PLINK_geno_imp +".rsq0.8.chr1-22",\
        file_PLINK_geno_imp +".rsq0.8.chr1-22",\
        join(path_source_code_env, "sources", "DeFinetti_hardy.r")))




def addL95U95(assoc_output):
    """ add L95 and U95 calculated from PVAL and OR """

    # header = CHR  SNP BP  A1  A2  FRQ_A   FRQ_U   INFO    OR  SE  P
    os.system("awk '{ print $9, $11 }' %s > %s.P.OR.tmp1" %(assoc_output, assoc_output))
    os.system("/home/sukmb113/progs/R-2.12.1/bin/R --no-save --args %s.P.OR.tmp1 %s.P.OR.tmp2 < %s" \
        %(assoc_output,\
          assoc_output,\
          join(path_source_code_env, "sources", "pval2qchisq_2.r")))
    os.system("gawk '{ print $1\"\t\"$2 }' %s.P.OR.tmp2 > %s.P.OR.tmp3" \
              %(assoc_output, assoc_output))
    os.system("paste %s %s.P.OR.tmp3 > %s.P.OR.tmp4" %(assoc_output, assoc_output, assoc_output))
    os.system("mv %s.P.OR.tmp4 %s" %(assoc_output, assoc_output))
    os.system("rm -f %s.P.OR.tmp1 %s.P.OR.tmp2 %s.P.OR.tmp3" %(assoc_output, assoc_output, assoc_output))




def PLINK_run_plots(file_prefix_PLINK):
    """ run plots from PLINK association program """

    file_PLINKlogistic               = file_prefix_PLINK + ".genotyped"
    file_PLINKdosage                 = file_prefix_PLINK + ".imputed"
    file_PLINK_geno_imp              = file_prefix_PLINK + ".genotyped.imputed"

    file_PLINK_geno_imp_rsq0_4       = file_PLINK_geno_imp +"_chr1-22.assoc.rsq0.4"
    file_PLINK_geno_imp_rsq0_8       = file_PLINK_geno_imp +"_chr1-22.assoc.rsq0.8"

    ## I found out that PLINK behaves the following way:
    ## IF MAF<=1 in cases and controls  -> NAs in PLINK --dosage result file     -> no assoc results for SNVs with MAF<=1
    ## If MAF>1 in cases and controls   -> no NAs in PLINK --dosage result file  -> assoc results for SNVs with MAF>1
    ## THESE FILES ONLY CONTAIN MAF>1% SNPS from now on

    # ---------------------------------------------- ##
    # -- extract samples, genotyped variants only -- ##
    # ---------------------------------------------- ##
    os.system("gawk '{ print $1, $2, $6 }' %s > %s"\
              %(disease_data_set_prefix_release_statistics + ".fam",\
                join(PLINK_dir, basename(disease_data_set_prefix_release_statistics)) + ".fam.single-id.pheno"))

    os.system("/${HOME}/progs/wdist/wdist_linux_x86_64_9Nov2016/plink --bfile %s --threads 16 --remove %s --keep %s --pheno %s --make-bed --allow-no-sex --out %s"\
        %(disease_data_set_prefix_release,\
          disease_data_set_prefix_release + "_flag.relatives.txt",\
          disease_data_set_prefix_release_statistics + ".fam",\
          join(PLINK_dir, basename(disease_data_set_prefix_release_statistics)) + ".fam.single-id.pheno",\
          join(PLINK_dir, basename(disease_data_set_prefix_release_statistics))))

    # doubleID was used for conversion to vcf file before imputation
    os.system("gawk '{ print $1\"_\"$2, $1\"_\"$2, $6 }' %s > %s"\
              %(join(PLINK_dir, basename(disease_data_set_prefix_release_statistics)) + ".fam",\
                join(PLINK_dir, basename(disease_data_set_prefix_release_statistics)) + ".fam.double-id.pheno"))
    # update phenotype because gender and cases status was lost during imputation
    os.system("gawk '{ print $1\"_\"$2, $1\"_\"$2, $5 }' %s > %s"\
              %(join(PLINK_dir, basename(disease_data_set_prefix_release_statistics)) + ".fam",\
                join(PLINK_dir, basename(disease_data_set_prefix_release_statistics)) + ".fam.double-id.gender"))

    os.system("gawk '{ print $1\"_\"$2, $1\"_\"$2, $3, $4, $5, $6 }' %s > %s"\
              %(join(PLINK_dir, basename(disease_data_set_prefix_release_statistics)) + ".fam",\
                join(PLINK_dir, basename(disease_data_set_prefix_release_statistics)) + ".fam.tmp"))
    os.system("rm -f %s"\
              %(join(PLINK_dir, basename(disease_data_set_prefix_release_statistics)) + ".fam"))
    os.system("mv %s %s"\
              %(join(PLINK_dir, basename(disease_data_set_prefix_release_statistics)) + ".fam.tmp",\
                join(PLINK_dir, basename(disease_data_set_prefix_release_statistics)) + ".fam"))

    for i in xrange(1,23):

        ## -------------------------------------------------- ##
        ## -- Part 1.2 PLINK --dosage for imputed dat set  -- ##
        ## -------------------------------------------------- ##

        if (numofPCs == 1):
            os.system("/${HOME}/progs/wdist/wdist_linux_x86_64_9Nov2016/plink --threads 16 --fam %s --map %s --dosage %s skip0=2 skip1=0 skip2=1 format=3 case-control-freqs --covar %s --covar-name PC%s --allow-no-sex --ci 0.95 --out %s"\
                %(join(PLINK_dir, basename(disease_data_set_prefix_release_statistics)) + ".fam",\
                  join(Imputation_orig_dir, str(i)+"."+disease_data_set_suffix_release_imputed +".PLINKdosage.map"),\
                  join(Imputation_orig_dir, str(i)+"."+disease_data_set_suffix_release_imputed +".PLINKdosage.gz"),\
                  disease_data_set_prefix_release + ".dat.pca.evec",\
                  str(numofPCs),\
                  file_PLINKdosage +"_chr"+str(i)))
        else:
            os.system("/${HOME}/progs/wdist/wdist_linux_x86_64_9Nov2016/plink --threads 16 --fam %s --map %s --dosage %s skip0=2 skip1=0 skip2=1 format=3 case-control-freqs --covar %s --covar-name PC1-PC%s --allow-no-sex --ci 0.95 --out %s"\
                %(join(PLINK_dir, basename(disease_data_set_prefix_release_statistics)) + ".fam",\
                  join(Imputation_orig_dir, str(i)+"."+disease_data_set_suffix_release_imputed +".PLINKdosage.map"),\
                  join(Imputation_orig_dir, str(i)+"."+disease_data_set_suffix_release_imputed +".PLINKdosage.gz"),\
                  disease_data_set_prefix_release + ".dat.pca.evec",\
                  str(numofPCs),\
                  file_PLINKdosage +"_chr"+str(i)))

        if (i == 1):
            # remove NA values
            os.system("cat %s | gawk '{ OFS=\"\t\"; print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11 }' | grep -v \"NA\" > %s"\
                %(file_PLINKdosage +"_chr"+str(i)+".assoc.dosage",\
                  file_PLINKdosage +"_chr1-22.assoc.dosage"))
            os.system("cat %s > %s"\
                %(file_PLINKdosage +"_chr"+str(i)+".log",\
                  file_PLINKdosage +"_chr1-22.log"))
        else:
            # remove NA values
            os.system("tail -n +2 %s | gawk '{ OFS=\"\t\"; print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11 }' | grep -v \"NA\" >> %s"\
                %(file_PLINKdosage +"_chr"+str(i)+".assoc.dosage",\
                  file_PLINKdosage +"_chr1-22.assoc.dosage"))
            os.system("tail -n +2 %s >> %s"\
                %(file_PLINKdosage +"_chr"+str(i)+".log",\
                  file_PLINKdosage +"_chr1-22.log"))

        ## ------------------------------------------------------ ##
        ## -- Part 1.1 PLINK --logistic for genotyped dat set  -- ##
        ## ------------------------------------------------------ ##
        # -------------------------------------------------------------------------------------- #
        # -- SNPs NOT ONLY FROM HD REGIONS because all SNPs from Immunochip are analysed here -- #
        # -------------------------------------------------------------------------------------- #
        # -- Program Stats_Immunochip_RMW_HRC1_1_plots.py will restrict the SNP set (also only gentoyped SNPs not in imputation reference) to the HD regions only -- #
        # ---------------------------------------------------------------------------------------------------------------------------------------------------------- #

        # tmp1 --logistic
        if (numofPCs == 1):
            os.system("/${HOME}/progs/wdist/wdist_linux_x86_64_9Nov2016/plink --threads 16 --bfile %s --logistic hide-covar --covar %s --covar-name PC%s --allow-no-sex --ci 0.95 --chr %s --out %s"\
                %(join(PLINK_dir, basename(disease_data_set_prefix_release_statistics)),\
                  disease_data_set_prefix_release + ".dat.pca.evec",\
                  str(numofPCs),\
                  str(i),\
                  file_PLINKlogistic +"_chr"+str(i)+"_tmp1"))

        else:
            os.system("/${HOME}/progs/wdist/wdist_linux_x86_64_9Nov2016/plink --threads 16 --bfile %s --logistic hide-covar --covar %s --covar-name PC1-PC%s --allow-no-sex --ci 0.95 --chr %s --out %s"\
                %(join(PLINK_dir, basename(disease_data_set_prefix_release_statistics)),\
                  disease_data_set_prefix_release + ".dat.pca.evec",\
                  str(numofPCs),\
                  str(i),\
                  file_PLINKlogistic +"_chr"+str(i)+"_tmp1"))

        # tmp2 --assoc for additional columns A1, A2, F_A, F_U
        os.system("/${HOME}/progs/wdist/wdist_linux_x86_64_9Nov2016/plink --threads 16 --bfile %s --assoc --allow-no-sex --ci 0.95 --chr %s --out %s"\
            %(join(PLINK_dir, basename(disease_data_set_prefix_release_statistics)),\
              str(i),\
              file_PLINKlogistic +"_chr"+str(i)+"_tmp2"))

        # merge tmp1 and tmp2
        if (i == 1):
            # workaround for gawk/head paste command in 3 steps, because the following error else:
            # gawk: (FILENAME=- FNR=1123) Fatal: print to "Standardausgabe"
            # fehlgeschlagen (Datenuebergabe unterbrochen (broken pipe))
            # paste: Schreibfehler: Datenuebergabe unterbrochen (broken pipe)
            # paste: Schreibfehler
            os.system("paste %s %s | gawk '{ OFS=\"\t\"; print $1, $2, $3, $4, $19, $17, $18, \"INFO\", $7, $8, $12 }' > %s"\
                      %(file_PLINKlogistic +"_chr"+str(i)+"_tmp1.assoc.logistic",\
                       file_PLINKlogistic +"_chr"+str(i)+"_tmp2.assoc",\
                       file_PLINKlogistic +"_chr1-22.assoc.logistic.tmp"))
            os.system("head -n 1 %s > %s"\
                      %(file_PLINKlogistic +"_chr1-22.assoc.logistic.tmp",\
                       file_PLINKlogistic +"_chr1-22.assoc.logistic"))
            os.system("rm -f %s"\
                      %(file_PLINKlogistic +"_chr1-22.assoc.logistic.tmp"))
            ####os.system("gawk '{ OFS=\"\t\"; print $1, $2, $3, $4, $19, $17, $18, \"INFO\", $7, $8, $12 }' %s | head -n 1 > %s"\
            ####          %(file_PLINKlogistic +"_chr1-22.assoc.logistic.tmp",\
            ####           file_PLINKlogistic +"_chr1-22.assoc.logistic"))
            ####os.system("paste %s %s | gawk '{ OFS=\"\t\"; print $1, $2, $3, $4, $19, $17, $18, \"INFO\", $7, $8, $12 }' | head -n 1 > %s"\
            ####          %(file_PLINKlogistic +"_chr"+str(i)+"_tmp1.assoc.logistic",\
            ####           file_PLINKlogistic +"_chr"+str(i)+"_tmp2.assoc",\
            ####           file_PLINKlogistic +"_chr1-22.assoc.logistic"))
        os.system("paste %s %s | gawk '{ OFS=\"\t\"; print $1, $2, $3, $4, $19, $17, $18, \"genotyped\", $7, $8, $12 }' > %s"\
                  %(file_PLINKlogistic +"_chr"+str(i)+"_tmp1.assoc.logistic",\
                   file_PLINKlogistic +"_chr"+str(i)+"_tmp2.assoc",\
                   file_PLINKlogistic +"_chr1-22.assoc.logistic.tmp"))
        os.system("tail -n +2 %s >> %s"\
                  %(file_PLINKlogistic +"_chr1-22.assoc.logistic.tmp",\
                   file_PLINKlogistic +"_chr1-22.assoc.logistic"))
        os.system("rm -f %s"\
                  %(file_PLINKlogistic +"_chr1-22.assoc.logistic.tmp"))
        #####os.system("paste %s %s | gawk '{ OFS=\"\t\"; print $1, $2, $3, $4, $19, $17, $18, \"genotyped\", $7, $8, $12 }' | tail -n +2 >> %s"\
        #####          %(file_PLINKlogistic +"_chr"+str(i)+"_tmp1.assoc.logistic",\
        #####            file_PLINKlogistic +"_chr"+str(i)+"_tmp2.assoc",\
        #####            file_PLINKlogistic +"_chr1-22.assoc.logistic"))

    # merge files part 1.1 and 1.2 and remove low INFO score SNPs from PLINK --dosage results, stratified by RSQ (INFO column) filtering of imputed variants rsq0.8, rsq0.4
    extract_Rsq_variants(assoc_logistic_input=file_PLINKlogistic +"_chr1-22.assoc.logistic",\
                         assoc_dosage_input=file_PLINKdosage +"_chr1-22.assoc.dosage",\
                         assoc_merge_output_rsq0_4=file_PLINK_geno_imp_rsq0_4,\
                         assoc_merge_output_rsq0_8=file_PLINK_geno_imp_rsq0_8)

    # add confidence interval from ORs
    addL95U95(assoc_output=file_PLINK_geno_imp_rsq0_4)
    addL95U95(assoc_output=file_PLINK_geno_imp_rsq0_4+".locuszoom")
    addL95U95(assoc_output=file_PLINK_geno_imp_rsq0_8)
    addL95U95(assoc_output=file_PLINK_geno_imp_rsq0_8+".locuszoom")

    # ----------------------------- #
    # -- Q-Q and manhattan plots -- #
    # ----------------------------- #

    qqplots_manhattanplots(file_PLINK_geno_imp_rsq0_4=file_PLINK_geno_imp_rsq0_4,\
                           file_PLINK_geno_imp_rsq0_8=file_PLINK_geno_imp_rsq0_8)

    # --------------------------------------------------------------------- #
    # -- Part 1.3 convert dosages to PLINK genotypes for LD calculations -- #
    # --------------------------------------------------------------------- #
    # ------------------------------------------------------------------------------------------------ #
    # -- rsq0.4 - these files can be used for every filtered subsets (e.g. rsq0.8, null snps, etc.) -- #
    # ------------------------------------------------------------------------------------------------ #

    # for each chromosome, convert into PLINK bed/bim/fam
    for i in xrange(1,23):

        os.system("${HOME}/progs/wdist/wdist_linux_x86_64_9Nov2016/plink --vcf %s --double-id --remove %s --keep %s --exclude %s --out %s --allow-no-sex --make-bed --threads 16"\
            %(join(Imputation_orig_dir, str(i)+"."+disease_data_set_suffix_release_imputed +".gz"),\
              disease_data_set_prefix_release + "_flag.relatives.doubleID.txt",\
              disease_data_set_prefix_release_statistics + ".ped",\
              Multiallelic_SNPexcludeList,\
              file_PLINKdosage + ".rsq0.4.chr" +str(i)+ ".rs"))

        # create another version with CHR:POS SNPIDs instead of rsnumbers and "." SNPIDs
        os.system("cp %s.fam %s.fam" %(file_PLINKdosage +".rsq0.4.chr" +str(i)+ ".rs", file_PLINKdosage +".rsq0.4.chr" +str(i)+ ".tmp"))
        os.system("cp %s.bed %s.bed" %(file_PLINKdosage +".rsq0.4.chr" +str(i)+ ".rs", file_PLINKdosage +".rsq0.4.chr" +str(i)+ ".tmp"))
        os.system("cp %s.log %s.log" %(file_PLINKdosage +".rsq0.4.chr" +str(i)+ ".rs", file_PLINKdosage +".rsq0.4.chr" +str(i)+ ".tmp"))
        os.system("LC_NUMERIC=POSIX gawk -f %s -- %s > %s" \
            %(join(path_source_code_env, "sources", "awk_rs2CHRPOS_bimfiles.awk"),\
              file_PLINKdosage +".rsq0.4.chr" +str(i)+ ".rs.bim",\
              file_PLINKdosage +".rsq0.4.chr" +str(i)+ ".tmp.bim"))

        os.system("${HOME}/progs/wdist/wdist_linux_x86_64_9Nov2016/plink --bfile %s --exclude %s --out %s --allow-no-sex --make-bed --threads 16"\
            %(file_PLINKdosage +".rsq0.4.chr" +str(i)+ ".tmp",\
              Multiallelic_SNPexcludeList,\
              file_PLINKdosage + ".rsq0.4.chr" +str(i)))

        # remove tmp files
        os.system("rm -f %s.bim %s.bed %s.fam %s.log %s.nosex" \
            %(file_PLINKdosage + ".rsq0.4.chr" +str(i)+ ".rs", \
              file_PLINKdosage + ".rsq0.4.chr" +str(i)+ ".rs", \
              file_PLINKdosage + ".rsq0.4.chr" +str(i)+ ".rs", \
              file_PLINKdosage + ".rsq0.4.chr" +str(i)+ ".rs", \
              file_PLINKdosage + ".rsq0.4.chr" +str(i)+ ".rs"))
        os.system("rm -f %s.bim %s.bed %s.fam %s.log %s.nosex" \
            %(file_PLINKdosage + ".rsq0.4.chr" +str(i)+ ".tmp", \
              file_PLINKdosage + ".rsq0.4.chr" +str(i)+ ".tmp", \
              file_PLINKdosage + ".rsq0.4.chr" +str(i)+ ".tmp", \
              file_PLINKdosage + ".rsq0.4.chr" +str(i)+ ".tmp", \
              file_PLINKdosage + ".rsq0.4.chr" +str(i)+ ".tmp"))

    ## ---------------------------- ##
    ## - no filtering of variants - ##
    ## ---------------------------- ##
    merge_file = file_PLINKdosage + "_merge-list.txt"
    generate_merge_file(file_PLINK=file_PLINKdosage +".rsq0.4", merge_file=merge_file)
    os.system("${HOME}/progs/wdist/wdist_linux_x86_64_9Nov2016/plink --bfile %s --merge-list %s --make-bed --out %s --allow-no-sex --threads 16" \
        %(file_PLINKdosage +".rsq0.4.chr1",\
          merge_file,\
          file_PLINKdosage +".rsq0.4.chr1-22.tmp"))

    # --------------------------------------------------- #
    # --  QCed PLINK DATA SET from HRC imputation only -- #
    # --------------------------------------------------- #
    os.system("${HOME}/progs/wdist/wdist_linux_x86_64_9Nov2016/plink --bfile %s --pheno %s --update-sex %s --make-bed --out %s --allow-no-sex --threads 16" \
        %(file_PLINKdosage +".rsq0.4.chr1-22.tmp",\
          join(PLINK_dir, basename(disease_data_set_prefix_release_statistics)) + ".fam.double-id.pheno",\
          join(PLINK_dir, basename(disease_data_set_prefix_release_statistics)) + ".fam.double-id.gender",\
          file_PLINKdosage +".rsq0.4.chr1-22"))

    # remove tmp files
    for i in xrange(1,23):
        os.system("rm -f %s.bim %s.bed %s.fam %s.log %s.nosex" \
            %(file_PLINKdosage + ".rsq0.4.chr" +str(i), \
              file_PLINKdosage + ".rsq0.4.chr" +str(i), \
              file_PLINKdosage + ".rsq0.4.chr" +str(i), \
              file_PLINKdosage + ".rsq0.4.chr" +str(i), \
              file_PLINKdosage + ".rsq0.4.chr" +str(i)))
    os.system("rm -f %s.bim %s.bed %s.fam %s.log %s.nosex"\
        %(file_PLINKdosage +".rsq0.4.chr1-22.tmp",\
          file_PLINKdosage +".rsq0.4.chr1-22.tmp",\
          file_PLINKdosage +".rsq0.4.chr1-22.tmp",\
          file_PLINKdosage +".rsq0.4.chr1-22.tmp",\
          file_PLINKdosage +".rsq0.4.chr1-22.tmp"))

    os.system("gawk '{ print $2 }' %s | tail -n +2 > %s"\
        %(file_PLINK_geno_imp_rsq0_8,\
          file_PLINKdosage +".rsq0.8.rs.txt"))

    os.system("${HOME}/progs/wdist/wdist_linux_x86_64_9Nov2016/plink --bfile %s --extract %s --make-bed --out %s --allow-no-sex --threads 16" \
        %(file_PLINKdosage +".rsq0.4.chr1-22",\
          file_PLINKdosage +".rsq0.8.rs.txt",\
          file_PLINKdosage +".rsq0.8.chr1-22"))

    # --------------------------------------------------------- #
    # -- check number of individuals pre and post imputation -- #
    # --------------------------------------------------------- #
    check_number_of_individuals_pre_and_post_imputation(file_PLINKdosage)

    # ------------------------------------------------------------------------------- #
    # --  THIS IS THE FINAL QCed PLINK DATA SET from genotyped plus HRC imputation -- #
    # ------------------------------------------------------------------------------- #

    # the tmp version still has duplicates from genotyping (rs numbers) and imputation (chr.:...)
    os.system("${HOME}/progs/wdist/wdist_linux_x86_64_9Nov2016/plink --bfile %s --bmerge %s %s %s --make-bed --out %s --allow-no-sex --threads 16" \
        %(file_PLINKdosage +".rsq0.4.chr1-22",\
          join(PLINK_dir, basename(disease_data_set_prefix_release_statistics)) + ".bed",\
          join(PLINK_dir, basename(disease_data_set_prefix_release_statistics)) + ".bim",\
          join(PLINK_dir, basename(disease_data_set_prefix_release_statistics)) + ".fam",\
          file_PLINK_geno_imp +".rsq0.4.chr1-22.tmp"))
    os.system("${HOME}/progs/wdist/wdist_linux_x86_64_9Nov2016/plink --bfile %s --bmerge %s %s %s --make-bed --out %s --allow-no-sex --threads 16" \
        %(file_PLINKdosage +".rsq0.8.chr1-22",\
          join(PLINK_dir, basename(disease_data_set_prefix_release_statistics)) + ".bed",\
          join(PLINK_dir, basename(disease_data_set_prefix_release_statistics)) + ".bim",\
          join(PLINK_dir, basename(disease_data_set_prefix_release_statistics)) + ".fam",\
          file_PLINK_geno_imp +".rsq0.8.chr1-22.tmp"))

    # determine the duplicate SNPs from genotyping (rs numbers) and imputation (chr.:...)
    os.system("gawk '{ print $1\":\"$4 }' %s.bim | sort | uniq -d > %s"\
        %(file_PLINK_geno_imp +".rsq0.4.chr1-22.tmp",\
          file_PLINK_geno_imp +".rsq0.4.chr1-22.excluded.duplicates.imputed"))

    # tmp version without duplicate genotyped and imputed SNPs - exclude imputed SNPs which are already genotyped
    os.system("${HOME}/progs/wdist/wdist_linux_x86_64_9Nov2016/plink --bfile %s --exclude %s --make-bed --out %s --allow-no-sex --threads 16" \
        %(file_PLINK_geno_imp +".rsq0.4.chr1-22.tmp",\
          file_PLINK_geno_imp +".rsq0.4.chr1-22.excluded.duplicates.imputed",\
          file_PLINK_geno_imp +".rsq0.4.chr1-22.tmp2"))
    os.system("${HOME}/progs/wdist/wdist_linux_x86_64_9Nov2016/plink --bfile %s --exclude %s --make-bed --out %s --allow-no-sex --threads 16" \
        %(file_PLINK_geno_imp +".rsq0.8.chr1-22.tmp",\
          file_PLINK_geno_imp +".rsq0.4.chr1-22.excluded.duplicates.imputed",\
          file_PLINK_geno_imp +".rsq0.8.chr1-22.tmp2"))

    # determine monomorphic SNPs
    os.system("${HOME}/progs/wdist/wdist_linux_x86_64_9Nov2016/plink --bfile %s --freq --out %s --allow-no-sex --threads 16" \
        %(file_PLINK_geno_imp +".rsq0.4.chr1-22.tmp2",\
          file_PLINK_geno_imp +".rsq0.4.chr1-22.tmp2_freq"))
    os.system("gawk '{ if ($5 == 0) print $2 }' %s > %s"\
        %(file_PLINK_geno_imp +".rsq0.4.chr1-22.tmp2_freq.frq",\
          file_PLINK_geno_imp +".rsq0.4.chr1-22.excluded.monomorphic"))

    # final version without duplicate genotyped and imputed SNPs - exclude imputed SNPs which are already genotyped
    os.system("${HOME}/progs/wdist/wdist_linux_x86_64_9Nov2016/plink --bfile %s --exclude %s --make-bed --out %s --allow-no-sex --threads 16" \
        %(file_PLINK_geno_imp +".rsq0.4.chr1-22.tmp2",\
          file_PLINK_geno_imp +".rsq0.4.chr1-22.excluded.monomorphic",\
          file_PLINK_geno_imp +".rsq0.4.chr1-22.tmp3"))
    os.system("${HOME}/progs/wdist/wdist_linux_x86_64_9Nov2016/plink --bfile %s --not-chr 0 --make-bed --out %s --allow-no-sex --threads 16" \
        %(file_PLINK_geno_imp +".rsq0.4.chr1-22.tmp3",\
          file_PLINK_geno_imp +".rsq0.4.chr1-22"))
    os.system("${HOME}/progs/wdist/wdist_linux_x86_64_9Nov2016/plink --bfile %s --exclude %s --make-bed --out %s --allow-no-sex --threads 16" \
        %(file_PLINK_geno_imp +".rsq0.8.chr1-22.tmp2",\
          file_PLINK_geno_imp +".rsq0.4.chr1-22.excluded.monomorphic",\
          file_PLINK_geno_imp +".rsq0.8.chr1-22.tmp3"))
    os.system("${HOME}/progs/wdist/wdist_linux_x86_64_9Nov2016/plink --bfile %s --not-chr 0 --make-bed --out %s --allow-no-sex --threads 16" \
        %(file_PLINK_geno_imp +".rsq0.8.chr1-22.tmp3",\
          file_PLINK_geno_imp +".rsq0.8.chr1-22"))

    # another final version without duplicate genotyped and imputed SNPs - # exclude imputed SNPs which are already genotyped - but here chr:pos as SNPids for Locuszoom
    os.system("awk '{ print $1\":\"$4 }' %s.bim | sort | uniq -d > %s.bim.duplicates" \
        %(file_PLINK_geno_imp +".rsq0.4.chr1-22",\
          file_PLINK_geno_imp +".rsq0.4.chr1-22"))
    os.system("${HOME}/progs/wdist/wdist_linux_x86_64_9Nov2016/plink --bfile %s --make-bed --out %s.tmp --allow-no-sex --threads 16" \
        %(file_PLINK_geno_imp +".rsq0.4.chr1-22",\
          file_PLINK_geno_imp +".rsq0.4.chr1-22.locuszoom"))
    os.system("awk '{ print $1\"\t\"$1\":\"$4\"\t\"$3\"\t\"$4\"\t\"$5\"\t\"$6 }' %s.bim > %s.tmp.bim" \
        %(file_PLINK_geno_imp +".rsq0.4.chr1-22",\
          file_PLINK_geno_imp +".rsq0.4.chr1-22.locuszoom"))
    os.system("${HOME}/progs/wdist/wdist_linux_x86_64_9Nov2016/plink --bfile %s.tmp --exclude %s.bim.duplicates --make-bed --out %s --allow-no-sex --threads 16" \
        %(file_PLINK_geno_imp +".rsq0.4.chr1-22.locuszoom",\
          file_PLINK_geno_imp +".rsq0.4.chr1-22",\
          file_PLINK_geno_imp +".rsq0.4.chr1-22.locuszoom"))

    os.system("awk '{ print $1\":\"$4 }' %s.bim | sort | uniq -d > %s.bim.duplicates" \
        %(file_PLINK_geno_imp +".rsq0.8.chr1-22",\
          file_PLINK_geno_imp +".rsq0.8.chr1-22"))
    os.system("${HOME}/progs/wdist/wdist_linux_x86_64_9Nov2016/plink --bfile %s --make-bed --out %s.tmp --allow-no-sex --threads 16" \
        %(file_PLINK_geno_imp +".rsq0.8.chr1-22",\
          file_PLINK_geno_imp +".rsq0.8.chr1-22.locuszoom"))
    os.system("awk '{ print $1\"\t\"$1\":\"$4\"\t\"$3\"\t\"$4\"\t\"$5\"\t\"$6 }' %s.bim > %s.tmp.bim" \
        %(file_PLINK_geno_imp +".rsq0.8.chr1-22",\
          file_PLINK_geno_imp +".rsq0.8.chr1-22.locuszoom"))
    os.system("${HOME}/progs/wdist/wdist_linux_x86_64_9Nov2016/plink --bfile %s.tmp --exclude %s.bim.duplicates --make-bed --out %s --allow-no-sex --threads 16" \
        %(file_PLINK_geno_imp +".rsq0.8.chr1-22.locuszoom",\
          file_PLINK_geno_imp +".rsq0.8.chr1-22",\
          file_PLINK_geno_imp +".rsq0.8.chr1-22.locuszoom"))

    # remove tmp files
    os.system("rm -f %s.bim %s.bed %s.fam %s.log %s.nosex"\
        %(file_PLINK_geno_imp +".rsq0.4.chr1-22.tmp",\
          file_PLINK_geno_imp +".rsq0.4.chr1-22.tmp",\
          file_PLINK_geno_imp +".rsq0.4.chr1-22.tmp",\
          file_PLINK_geno_imp +".rsq0.4.chr1-22.tmp",\
          file_PLINK_geno_imp +".rsq0.4.chr1-22.tmp"))
    os.system("rm -f %s.bim %s.bed %s.fam %s.log %s.nosex"\
        %(file_PLINK_geno_imp +".rsq0.8.chr1-22.tmp",\
          file_PLINK_geno_imp +".rsq0.8.chr1-22.tmp",\
          file_PLINK_geno_imp +".rsq0.8.chr1-22.tmp",\
          file_PLINK_geno_imp +".rsq0.8.chr1-22.tmp",\
          file_PLINK_geno_imp +".rsq0.8.chr1-22.tmp"))
    os.system("rm -f %s.bim %s.bed %s.fam %s.log %s.nosex"\
        %(file_PLINK_geno_imp +".rsq0.4.chr1-22.tmp2",\
          file_PLINK_geno_imp +".rsq0.4.chr1-22.tmp2",\
          file_PLINK_geno_imp +".rsq0.4.chr1-22.tmp2",\
          file_PLINK_geno_imp +".rsq0.4.chr1-22.tmp2",\
          file_PLINK_geno_imp +".rsq0.4.chr1-22.tmp2"))
    os.system("rm -f %s.bim %s.bed %s.fam %s.log %s.nosex"\
        %(file_PLINK_geno_imp +".rsq0.8.chr1-22.tmp2",\
          file_PLINK_geno_imp +".rsq0.8.chr1-22.tmp2",\
          file_PLINK_geno_imp +".rsq0.8.chr1-22.tmp2",\
          file_PLINK_geno_imp +".rsq0.8.chr1-22.tmp2",\
          file_PLINK_geno_imp +".rsq0.8.chr1-22.tmp2"))
    os.system("rm -f %s.bim %s.bed %s.fam %s.log %s.nosex"\
        %(file_PLINK_geno_imp +".rsq0.4.chr1-22.locuszoom.tmp",\
          file_PLINK_geno_imp +".rsq0.4.chr1-22.locuszoom.tmp",\
          file_PLINK_geno_imp +".rsq0.4.chr1-22.locuszoom.tmp",\
          file_PLINK_geno_imp +".rsq0.4.chr1-22.locuszoom.tmp",\
          file_PLINK_geno_imp +".rsq0.4.chr1-22.locuszoom.tmp"))
    os.system("rm -f %s.bim %s.bed %s.fam %s.log %s.nosex"\
        %(file_PLINK_geno_imp +".rsq0.8.chr1-22.locuszoom.tmp",\
          file_PLINK_geno_imp +".rsq0.8.chr1-22.locuszoom.tmp",\
          file_PLINK_geno_imp +".rsq0.8.chr1-22.locuszoom.tmp",\
          file_PLINK_geno_imp +".rsq0.8.chr1-22.locuszoom.tmp",\
          file_PLINK_geno_imp +".rsq0.8.chr1-22.locuszoom.tmp"))
    os.system("rm -f %s.bim %s.bed %s.fam %s.log %s.nosex"\
        %(file_PLINK_geno_imp +".rsq0.4.chr1-22.tmp3",\
          file_PLINK_geno_imp +".rsq0.4.chr1-22.tmp3",\
          file_PLINK_geno_imp +".rsq0.4.chr1-22.tmp3",\
          file_PLINK_geno_imp +".rsq0.4.chr1-22.tmp3",\
          file_PLINK_geno_imp +".rsq0.4.chr1-22.tmp3"))
    os.system("rm -f %s.bim %s.bed %s.fam %s.log %s.nosex"\
        %(file_PLINK_geno_imp +".rsq0.8.chr1-22.tmp3",\
          file_PLINK_geno_imp +".rsq0.8.chr1-22.tmp3",\
          file_PLINK_geno_imp +".rsq0.8.chr1-22.tmp3",\
          file_PLINK_geno_imp +".rsq0.8.chr1-22.tmp3",\
          file_PLINK_geno_imp +".rsq0.8.chr1-22.tmp3"))

    # convert ot vcf format, bgzip and tabix to use this as reference LD data set for Locuszoom program
    os.system("${HOME}/progs/wdist/wdist_linux_x86_64_9Nov2016/plink --bfile %s --recode vcf --out %s --allow-no-sex --threads 16" \
        %(file_PLINK_geno_imp +".rsq0.4.chr1-22",\
          file_PLINK_geno_imp +".rsq0.4.chr1-22_vcf"))
    os.system("bgzip %s" \
        %(file_PLINK_geno_imp +".rsq0.4.chr1-22_vcf.vcf"))
    os.system("tabix -p vcf %s" \
        %(file_PLINK_geno_imp +".rsq0.4.chr1-22_vcf.vcf.gz"))
    os.system("${HOME}/progs/wdist/wdist_linux_x86_64_9Nov2016/plink --bfile %s --recode vcf --out %s --allow-no-sex --threads 16" \
        %(file_PLINK_geno_imp +".rsq0.8.chr1-22",\
          file_PLINK_geno_imp +".rsq0.8.chr1-22_vcf"))
    os.system("bgzip %s" \
        %(file_PLINK_geno_imp +".rsq0.8.chr1-22_vcf.vcf"))
    os.system("tabix -p vcf %s" \
        %(file_PLINK_geno_imp +".rsq0.8.chr1-22_vcf.vcf.gz"))

    #### if multi-allelic SNPs available, remove one SNV from duplicates, then merge again
    ####if os.path.isfile(file_PLINKdosage +".rsq0.4.chr1-22-merge.missnp"):
    ####    os.system("mv %s %s"\
    ####              %(file_PLINKdosage +".rsq0.4.chr1-22-merge.missnp",\
    ####              file_PLINKdosage +".rsq0.4.chr1-22-merge.missnp.multi-allelicSNVs.txt"))
    ####    os.system("${HOME}/progs/wdist/wdist_linux_x86_64_9Nov2016/plink --bfile %s --merge-list %s --make-bed --out %s --threads 16 --exclude %s" \
    ####        %(file_PLINKdosage +".rsq0.4.chr1",\
    ####          merge_file,\
    ####          file_PLINKdosage +".rsq0.4.chr1-22",\
    ####          file_PLINKdosage +".rsq0.4.chr1-22-merge.missnp.multi-allelicSNVs.txt"))

    # --------------------- #
    # -- DeFinetti plots -- #
    # --------------------- #
    DeFinetti_plots(file_PLINK_geno_imp=file_PLINK_geno_imp)




def PLINK_clumping(file_prefix_PLINK):
    """ run clumping from PLINK association program """

    file_PLINK_geno_imp              = file_prefix_PLINK + ".genotyped.imputed"
    file_PLINK_geno_imp_rsq0_4       = file_PLINK_geno_imp +"_chr1-22.assoc.rsq0.4.locuszoom"
    file_PLINK_geno_imp_rsq0_8       = file_PLINK_geno_imp +"_chr1-22.assoc.rsq0.8.locuszoom"

    # --------------------- #
    # -- clumping rsq0.4 -- #
    # --------------------- #
    print "\n    do clumping with PLINK ...\n\n"
    os.system("${HOME}/progs/wdist/wdist_linux_x86_64_9Nov2016/plink --bfile %s --out %s_clump --clump %s --clump-r2 %s --clump-p1 %s --clump-p2 %s --clump-kb %s --allow-no-sex --threads 16"\
        %(file_PLINK_geno_imp +".rsq0.4.chr1-22.locuszoom",\
        file_PLINK_geno_imp_rsq0_4,\
        file_PLINK_geno_imp_rsq0_4,\
        clumpr2,\
        clumpp1,\
        clumpp2,\
        clumpkb))

    print "\n    filter clumps with supporting snps ...\n\n"
    clump = Clump(clump_file=file_PLINK_geno_imp_rsq0_4 + "_clump.clumped",\
                  write_file=file_PLINK_geno_imp_rsq0_4 + "_clump.clumped_all")
    clump.write_snps_typed_noChr23() # these snps are all typed
    del clump
    clump = Clump(clump_file=file_PLINK_geno_imp_rsq0_4 + "_clump.clumped_all",\
                  write_file=file_PLINK_geno_imp_rsq0_4 + "_clump.clumped_groups")
    clump.write_supported_snps_from_altered_file()
    del clump

    # ---------------------------- #
    # -- filter for non-xMHC    -- #
    # -- xMHC: chr6, [25,34[ Mb -- #
    # ---------------------------- #

    os.system("gawk '{ if (!($1 == 6 && ($4 >= 25000000 && $4 < 34000000))) print }' %s_clump.clumped_all > %s_clump.clumped_all.noXMHC"\
              %(file_PLINK_geno_imp_rsq0_4, file_PLINK_geno_imp_rsq0_4))
    os.system("gawk '{ if (!($1 == 6 && ($4 >= 25000000 && $4 < 34000000))) print }' %s_clump.clumped_groups > %s_clump.clumped_groups.noXMHC"\
              %(file_PLINK_geno_imp_rsq0_4, file_PLINK_geno_imp_rsq0_4))

    # --------------------------------------------------------- #
    # -- clumping                                            -- #
    # -- rsq0.8                                              -- #
    # --------------------------------------------------------- #
    print "\n    do clumping with PLINK ...\n\n"

    os.system("${HOME}/progs/wdist/wdist_linux_x86_64_9Nov2016/plink --bfile %s --out %s_clump --clump %s --clump-r2 %s --clump-p1 %s --clump-p2 %s --clump-kb %s --allow-no-sex --threads 16"\
        %(file_PLINK_geno_imp +".rsq0.8.chr1-22.locuszoom",\
        file_PLINK_geno_imp_rsq0_8,\
        file_PLINK_geno_imp_rsq0_8,\
        clumpr2,\
        clumpp1,\
        clumpp2,\
        clumpkb))
    print "\n    filter clumps with supporting snps ...\n\n"
    clump = Clump(clump_file=file_PLINK_geno_imp_rsq0_8 + "_clump.clumped",\
                  write_file=file_PLINK_geno_imp_rsq0_8 + "_clump.clumped_all")
    clump.write_snps_typed_noChr23() # these snps are all typed
    del clump
    clump = Clump(clump_file=file_PLINK_geno_imp_rsq0_8 + "_clump.clumped_all",\
                  write_file=file_PLINK_geno_imp_rsq0_8 + "_clump.clumped_groups")
    clump.write_supported_snps_from_altered_file()
    del clump

    # ---------------------------- #
    # -- filter for non-xMHC    -- #
    # -- xMHC: chr6, [25,34[ Mb -- #
    # ---------------------------- #

    os.system("gawk '{ if (!($1 == 6 && ($4 >= 25000000 && $4 < 34000000))) print }' %s_clump.clumped_all > %s_clump.clumped_all.noXMHC"\
              %(file_PLINK_geno_imp_rsq0_8, file_PLINK_geno_imp_rsq0_8))
    os.system("gawk '{ if (!($1 == 6 && ($4 >= 25000000 && $4 < 34000000))) print }' %s_clump.clumped_groups > %s_clump.clumped_groups.noXMHC"\
              %(file_PLINK_geno_imp_rsq0_8, file_PLINK_geno_imp_rsq0_8))




def generateHitSpecfile(chromosome, pos, snp, GWAS_regions, file_HitSpec):
    """ """
    try:
        fh_w = open(file_HitSpec, "w")
    except IOError, e:
        print e
        sys.exit(1)

    sep = "\t"
    fh_w.writelines("snp" +sep+\
          "chr" +sep+\
          "start" +sep+\
          "end" +sep+\
          "flank" +sep+\
          "run" +sep+\
          "m2zargs\n")

    flank        = "500kb"
    run          = "yes"
    start        = "NA"
    end          = "NA"
    #phenos       = list[5]

    m2zargs      = "title=\""+ basename(disease_data_set_prefix_release_statistics) +" "+ snp +"\""

    # -- highlight known GWAS region if in plotting region -- #
    if QQplot_SNPexcludeList_exists:
        left  = int(pos) - 500000
        right = int(pos) + 500000
        for region in GWAS_regions:
            if chromosome == str(region[0]):
                # whole plotting region is known locu
                if int(region[1]) <= left and right <= int(region[2]):
                    hiStart      = str(left)
                    hiEnd        = str(right)
                    m2zargs      += " hiStart=\""+ hiStart + "\"" + " hiEnd=\""+ hiEnd + "\" hiColor=blue"
                # known locus intersects left plotting border
                elif int(region[1]) <= left and (left < int(region[2]) and int(region[2]) <= right):
                    hiStart      = str(left)
                    hiEnd        = str(region[2])
                    m2zargs      += " hiStart=\""+ hiStart + "\"" + " hiEnd=\""+ hiEnd + "\" hiColor=blue"
                # known locus intersects right plotting border
                elif (left <= int(region[1]) and int(region[1]) < right) and right <= int(region[2]):
                    hiStart      = str(region[1])
                    hiEnd        = str(right)
                    m2zargs      += " hiStart=\""+ hiStart + "\"" + " hiEnd=\""+ hiEnd + "\" hiColor=blue"
                # known locus within left and right plotting border
                elif (left <= int(region[1]) and int(region[1]) < right) and (left < int(region[2]) and int(region[2]) <= right):
                    hiStart      = str(region[1])
                    hiEnd        = str(region[2])
                    m2zargs      += " hiStart=\""+ hiStart + "\"" + " hiEnd=\""+ hiEnd + "\" hiColor=blue"

    m2zargs      += " showAnnot=T showRefsnpAnnot=T annotPch=\"21,24,24,25,22,22,8,7\""
    m2zargs      += " rfrows=4"
    m2zargs      += " refsnpTextSize=0.8"
    m2zargs      += " legendSize=0.4"
    #m2zargs      += " snpset=\"HapMap\""
    #m2zargs      += " metalRug=\"Immunochip\""
    #m2zargs      += " theme=\"publication\""

    if (snp == chromosome+":"+pos):
        fh_w.writelines("chr" +snp +sep+\
              chromosome +sep+\
              start +sep+\
              end +sep+\
              flank +sep+\
              run +sep+\
              m2zargs +"\n")
    elif (snp[0:3] == "imm" or snp[0:3] == "seq"):
        fh_w.writelines("chr" +chromosome+":"+pos +sep+\
              chromosome +sep+\
              start +sep+\
              end +sep+\
              flank +sep+\
              run +sep+\
              m2zargs +"\n")
    else:
        fh_w.writelines(snp +sep+\
              chromosome +sep+\
              start +sep+\
              end +sep+\
              flank +sep+\
              run +sep+\
              m2zargs +"\n")

    fh_w.close()




def generateCredibleSets(chromosome, pos, snp, pval, file_CredibleSets):
    """ """

    sep = "\t"
    #counter2color = ["black", "red", "green3", "blue", "cyan", "magenta", "yellow", "gray"]
    counter2color = ["black", "red", "green3", "blue", "cyan", "magenta", "yellow", "gray", "antiquewhite4", "aquamarine4", "azure4", "bisque4", "blueviolet", "brown4", "cadetblue4", "chartreuse", "chocolate4", "cornflowerblue", "cornsilk4", "cyan4", "darkgoldenrod4", "darkmagenta", "darkolivegreen4", "darkorange1", "darkorange4", "darkorchid", "darkorchid4", "darkred", "blue4", "aquamarine1", "brown2", "cornsilk2"]
    snp_counter = 0

    try:
        fh_w = open(file_CredibleSets, "w")
    except IOError, e:
        print e
        sys.exit(1)

    fh_w.writelines("snp" +sep+\
          "chr" +sep+\
          "pos" +sep+\
          "pp" +sep+\
          "p-value" +sep+\
          "group" +sep+\
          "color\n")

    snp_counter = 0
    pp  = "NA"

    if (snp == chromosome+":"+pos):
        fh_w.writelines("chr" +snp +sep+\
              chromosome +sep+\
              pos +sep+\
              pp +sep+\
              pval +sep+\
              snp +sep+\
              counter2color[snp_counter] + "\n")
    elif (snp[0:3] == "imm" or snp[0:3] == "seq"):
        fh_w.writelines("chr" +chromosome+":"+pos +sep+\
              chromosome +sep+\
              pos +sep+\
              pp +sep+\
              pval +sep+\
              snp +sep+\
              counter2color[snp_counter] + "\n")
    else:
        fh_w.writelines(snp +sep+\
              chromosome +sep+\
              pos +sep+\
              pp +sep+\
              pval +sep+\
              snp +sep+\
              counter2color[snp_counter] + "\n")

    ####LDvariants = sp2.replace("_typed(1)","").split(",")
    ####for variant in LDvariants:
    ####    fh_w.writelines(variant +sep+\
    ####          chromosome +sep+\
    ####          pos +sep+\
    ####          pp +sep+\
    ####          pval +sep+\
    ####          snp +" and LD variants"+ +sep+\
    ####          counter2color[snp_counter] + "\n")

    ####snp_counter += 1
    fh_w.close()




def generateDenoteMarker(chromosome, pos, snp, file_DenoteMarker):
    """ """

    sep = "\t"

    try:
        fh_w = open(file_DenoteMarker, "w")
    except IOError, e:
        print e
        sys.exit(1)

    fh_w.writelines("snp" +sep+\
          "chr" +sep+\
          "pos" +sep+\
          "string" +sep+\
          "color\n")

    if (snp == chromosome+":"+pos):
        fh_w.writelines("chr" +snp +sep+\
              chromosome +sep+\
              pos +sep+\
              "" +sep+\
              "black\n")
    elif (snp[0:3] == "imm" or snp[0:3] == "seq"):
        fh_w.writelines("chr" +chromosome+":"+pos +sep+\
              chromosome +sep+\
              pos +sep+\
              "" +sep+\
              "black\n")
    else:
        fh_w.writelines(snp +sep+\
              chromosome +sep+\
              pos +sep+\
              "" +sep+\
              "black\n")

    fh_w.close()




def run_locuszoom(chromosome, file_PLINK_assoc, file_HitSpec, file_CredibleSets, file_DenoteMarker, prefix_out, file_PLINK_ld):
    """ """

    # workaround: Because Locuszoom tries to display all loci in file_CredibleSets, I need a separate input for each locus
    try:
        fh_r1 = open(file_CredibleSets, "r")
        fh_r2 = open(file_HitSpec, "r")
        fh_r3 = open(file_DenoteMarker, "r")
        # write new tmp files
        fh_w1_tmp = open(file_CredibleSets +".tmp", "w")
        fh_w2_tmp = open(file_HitSpec +".tmp", "w")
        fh_w3_tmp = open(file_DenoteMarker +".tmp", "w")
    except IOError, e:
        print e
        sys.exit(1)

    # read header
    header_1 = fh_r1.readline().rstrip('\n')
    header_2 = fh_r2.readline().rstrip('\n')
    header_3 = fh_r3.readline().rstrip('\n')
    header = True
    line_1 = header_1
    while line_1:

        if header:
            # print header
            fh_w1_tmp.writelines(header_1 +"\n")
            fh_w2_tmp.writelines(header_2 +"\n")
            fh_w3_tmp.writelines(header_3 +"\n")

            # print first line
            line_1   = fh_r1.readline().rstrip('\n')
            line_2   = fh_r2.readline().rstrip('\n')
            line_3   = fh_r3.readline().rstrip('\n')
            fh_w1_tmp.writelines(line_1 +"\n")
            fh_w2_tmp.writelines(line_2 +"\n")
            fh_w3_tmp.writelines(line_3 +"\n")
            list = re.split("\s+",line_1)
            locus_nr_prev = list[5]
            header = False
        else:
            list = re.split("\s+",line_1)
            locus_nr = list[5]
            if locus_nr == locus_nr_prev:
                fh_w1_tmp.writelines(line_1 +"\n")
                fh_w2_tmp.writelines(line_2 +"\n")
                fh_w3_tmp.writelines(line_3 +"\n")
                locus_nr_prev = locus_nr
            else:
                # close current file and plot with Locuszoom
                fh_w1_tmp.close()
                fh_w2_tmp.close()
                fh_w3_tmp.close()
                os.system("locuszoom --metal=%s --markercol SNP --pvalcol P --ld %s --build hg19 --hitspec=%s fineMap=\"%s\" --denote-markers-file %s --gwas-cat whole-cat_significant-only --plotonly --no-date --prefix %s" %(file_PLINK_assoc, file_PLINK_ld, file_HitSpec +".tmp", file_CredibleSets +".tmp", file_DenoteMarker +".tmp", prefix_out))
                #os.system("locuszoom --metal=%s --markercol SNP --pvalcol P --ld-vcf %s --build hg19 --hitspec=%s fineMap=\"%s\" --denote-markers-file %s --gwas-cat whole-cat_significant-only --plotonly --no-date --prefix %s" %(file_PLINK_assoc, join(Imputation_orig_dir, str(chromosome)+"."+disease_data_set_suffix_release_imputed +".gz"), file_HitSpec +".tmp", file_CredibleSets +".tmp", file_DenoteMarker +".tmp", prefix_out))
                # vcf file tabixed from filtered PLINK data does not work, don't know why, but LD cannot be calculated by Locuszoom
                #os.system("locuszoom --metal=%s --markercol SNP --pvalcol P --ld-vcf %s --build hg19 --hitspec=%s fineMap=\"%s\" --denote-markers-file %s --gwas-cat whole-cat_significant-only --plotonly --no-date --prefix %s" %(file_PLINK_assoc, file_PLINK_rsq0_4, file_HitSpec +".tmp", file_CredibleSets +".tmp", file_DenoteMarker +".tmp", prefix_out))

                # write new tmp files
                try:
                    fh_w1_tmp = open(file_CredibleSets +".tmp", "w")
                    fh_w2_tmp = open(file_HitSpec +".tmp", "w")
                    fh_w3_tmp = open(file_DenoteMarker +".tmp", "w")
                except IOError, e:
                    print e
                    sys.exit(1)

                # print header
                fh_w1_tmp.writelines(header_1 +"\n")
                fh_w2_tmp.writelines(header_2 +"\n")
                fh_w3_tmp.writelines(header_3 +"\n")

                fh_w1_tmp.writelines(line_1 +"\n")
                fh_w2_tmp.writelines(line_2 +"\n")
                fh_w3_tmp.writelines(line_3 +"\n")

                list = re.split("\s+",line_1)
                locus_nr_prev = list[5]

        line_1   = fh_r1.readline().rstrip('\n')
        line_2   = fh_r2.readline().rstrip('\n')
        line_3   = fh_r3.readline().rstrip('\n')

    # close current file and plot with Locuszoom
    fh_w1_tmp.close()
    fh_w2_tmp.close()
    fh_w3_tmp.close()
    os.system("locuszoom --metal=%s --markercol SNP --pvalcol P --ld %s --build hg19 --hitspec=%s fineMap=\"%s\" --denote-markers-file %s --gwas-cat whole-cat_significant-only --plotonly --no-date --prefix %s" %(file_PLINK_assoc, file_PLINK_ld, file_HitSpec +".tmp", file_CredibleSets +".tmp", file_DenoteMarker +".tmp", prefix_out))
    #os.system("locuszoom --metal=%s --markercol SNP --pvalcol P --ld-vcf %s --build hg19 --hitspec=%s fineMap=\"%s\" --denote-markers-file %s --gwas-cat whole-cat_significant-only --plotonly --no-date --prefix %s" %(file_PLINK_assoc, join(Imputation_orig_dir, str(chromosome)+"."+disease_data_set_suffix_release_imputed +".gz"), file_HitSpec +".tmp", file_CredibleSets +".tmp", file_DenoteMarker +".tmp", prefix_out))
    # vcf file tabixed from filtered PLINK data does not work, don't know why, but LD cannot be calculated by Locuszoom
    #os.system("locuszoom --metal=%s --markercol SNP --pvalcol P --ld-vcf %s --build hg19 --hitspec=%s fineMap=\"%s\" --denote-markers-file %s --gwas-cat whole-cat_significant-only --plotonly --no-date --prefix %s" %(file_PLINK_assoc, file_PLINK_rsq0_4, file_HitSpec +".tmp", file_CredibleSets +".tmp", file_DenoteMarker +".tmp", prefix_out))
    os.system("rm -f %s" %(file_CredibleSets +".tmp"))
    os.system("rm -f %s" %(file_HitSpec +".tmp"))
    os.system("rm -f %s" %(file_DenoteMarker +".tmp"))
    fh_r1.close()
    fh_r2.close()
    fh_r3.close()




def extract_known_GWAS_regions():
    """ extract known GWAS region from file """

    try:
        fh_r = file(QQplot_SNPexcludeList, "r")
    except IOError, e:
      print e
      sys.exit(1)

    GWAS_regions = []
    line = fh_r.readline().rstrip('\n')
    while line:
        list = re.split("\s+",line)
        chr    = int(list[0])
        left   = int(list[1])
        right  = int(list[2])
        GWAS_regions.append((chr, left, right))
        line = fh_r.readline().rstrip('\n')

    fh_r.close()

    return GWAS_regions




def locuszoom_run(snp_colnr, file_suffix, file_PLINK_assoc, file_PLINK):
    """ do the plotting """

    # ----------------------------------------------------------- #
    # -- get numoftophits top snp_ids from clumped groups file -- #
    # ----------------------------------------------------------- #

    # this list has the correct order of SNPs how they appear in the clump file
    snp_tophits = []
    # clump rs : chr, rs_lead, pos, p_value, SP2
    snp_hash = {}
    #snp_tophits_sp2_hash = {}

    # if clumped groups or clumped all file
    if file_suffix == "_clump.clumped_all" or file_suffix == "_clump.clumped_groups" or \
       file_suffix == "_clump.clumped_all.noXMHC" or file_suffix == "_clump.clumped_groups.noXMHC":
        clump = Clump(clump_file=file_PLINK_assoc + file_suffix)
        clump.map()

        snp_tophits = clump.rs_get_typed_imputed_tophits(numoftophits=numoftophits)
        snp_hash = clump.rs_get_all_typed_imputed_tophits_hash(numoftophits=numoftophits)
        #snp_tophits_sp2_hash = clump.rs_sp2_get_typed_imputed_tophits_hash(numoftophits=numoftophits)

        clump.free_map() ; del clump

    # ------------------------------------------------------- #
    # -- create new output directory for association plots -- #
    # ------------------------------------------------------- #

    plot_dir = join(PLINK_dir,\
        "tmpfiles_locuszoom_" + basename(file_PLINK_assoc) + file_suffix)
    os.system("mkdir -p %s" %(plot_dir))

    # ---------------------------------------------------------------------- #
    # -- extract known GWAS region if file exists with known GWAS regions -- #
    # ---------------------------------------------------------------------- #
    GWAS_regions = []
    if QQplot_SNPexcludeList_exists:
        GWAS_regions = extract_known_GWAS_regions()

    # ----------------------- #
    # -- Locuszoom per SNP -- #
    # ----------------------- #
    string = ""
    for snp in snp_tophits:

        chromosome = snp_hash[snp][0]
        pos  = snp_hash[snp][2]
        pval = snp_hash[snp][3]
        #sp2  = snp_hash[snp][4]

        # generate hitspec file for Locuszoom plotting
        file_HitSpec = join(plot_dir, snp +".hitspec.csv")
        generateHitSpecfile(chromosome=chromosome, pos=pos, snp=snp, GWAS_regions=GWAS_regions, file_HitSpec=file_HitSpec)

        # generate credibleSets file for displaying different disease subsets
        file_CredibleSets = join(plot_dir, snp +".credibleSets.csv")
        #####generateCredibleSets(chromosome=chromosome, pos=pos, snp=snp, sp2=sp2, pval=pval, file_CredibleSets=file_CredibleSets)
        generateCredibleSets(chromosome=chromosome, pos=pos, snp=snp, pval=pval, file_CredibleSets=file_CredibleSets)

        # generate denoteMarker file for displaying SNV IDs in plot
        file_DenoteMarker = join(plot_dir, snp +".denoteMarker.csv")
        generateDenoteMarker(chromosome=chromosome, pos=pos, snp=snp, file_DenoteMarker=file_DenoteMarker)

        # calculate LD on my own for --ld option in Locuszoom
        os.system("${HOME}/progs/wdist/wdist_linux_x86_64_9Nov2016/plink --bfile %s --r2 --ld-window-r2 0.0 --ld-snp %s --ld-window 1000000 --ld-window-kb 1000 --out %s_ld --allow-no-sex" %(file_PLINK, snp, join(plot_dir, snp)))
        os.system("echo \"snp1 snp2 dprime rsquare\" > %s_ld.ld.csv" %(join(plot_dir, snp)))
        if (snp == chromosome+":"+pos):
            os.system("awk '{ print \"chr\"$3, $6, \"0\", $7 }' %s_ld.ld | tail -n +2 >> %s_ld.ld.csv" %(join(plot_dir, snp), join(plot_dir, snp)))
        elif (snp[0:3] == "imm" or snp[0:3] == "seq"):
            os.system("awk '{ print \"chr\"$1\":\"$2, $6, \"0\", $7 }' %s_ld.ld | tail -n +2 >> %s_ld.ld.csv" %(join(plot_dir, snp), join(plot_dir, snp)))
        else:
            os.system("awk '{ print $3, $6, \"0\", $7 }' %s_ld.ld | tail -n +2 >> %s_ld.ld.csv" %(join(plot_dir, snp), join(plot_dir, snp)))

        # run locuszoom
        run_locuszoom(chromosome, file_PLINK_assoc, file_HitSpec, file_CredibleSets, file_DenoteMarker, join(plot_dir, snp +".locuszoom"), join(plot_dir, snp+"_ld.ld.csv"))
        if (snp == chromosome+":"+pos):
            string += " %s" %(join(plot_dir, snp +".locuszoom_chr" +chromosome+ "_" +pos+ ".pdf"))
        elif (snp[0:3] == "imm" or snp[0:3] == "seq"):
            string += " %s" %(join(plot_dir, snp +".locuszoom_chr" +chromosome+ "_" +pos+ ".pdf"))
        else:
            string += " %s" %(join(plot_dir, snp +".locuszoom_" +snp+ ".pdf"))

    print "\n    merging single pdf regional association plots into one final pdf  ...\n\n"
    print "\n\n    OutputFile: "+ file_PLINK_assoc + file_suffix + ".pdf" +"  ...\n\n"
    print "\n\n    String    : "+ string +"  ...\n\n"
    os.system("gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite -sOutputFile=%s %s"\
              %(file_PLINK_assoc + file_suffix + ".pdf",\
                string))




######################################
#
# -- MAIN --
#
######################################


def run_PLINK(sample_with_params):
    """ run PLINK plots from association testing """

    # -------------- #
    # -- aliasing -- #
    # -------------- #

    aliasing(sample_with_params=sample_with_params)

    file_prefix_PLINK = join(PLINK_dir, \
                             basename(disease_data_set_prefix_release_statistics))


    ## ------------------------------------------- ##
    ## -- PLINK --dosage and --logistic MAF>=1% -- ##
    ## ------------------------------------------- ##

    # I found out that PLINK behaves the following way with --dosage
    # IF MAF<=1 in cases and controls  -> NAs in PLINK --dosage result file     -> no assoc results for SNVs with MAF<=1
    # If MAF>1 in cases and controls   -> no NAs in PLINK --dosage result file  -> assoc results for SNVs with MAF>1

    # PLINK --logistic for all genotyped variants (every MAF) which are not available in HRC1.1 reference
    # plus PLINK --dosage for imputed variants with MAF>=1% which are available in HRC1.1 reference and not genotyped on Ichip
    # imputed variants with MAF<1% (not genotyped on Ichip) must be analyzed and added with Stats_Immunochip_RMW_HRC1_1_plots.py

    print "\nrun plots of PLINK --dosage and --logistic program ...\n\n"
    PLINK_run_plots(file_prefix_PLINK=file_prefix_PLINK)

    print "\nrun clumping of PLINK --dosage and --logistic program ...\n\n"
    PLINK_clumping(file_prefix_PLINK=file_prefix_PLINK)

    ## -------------------------------- ##
    ## -- Regional association plots -- ##
    ## -------------------------------- ##

    print "\ndo the Locuszoom plotting of clumped files ...\n\n"

    ## --------------------------------------------------------- #
    ## -- PLINK                                               -- #
    ## -- only SNPs above MAF 0.01                            -- #
    ## -- without xMHC                                        -- #
    ## --                                                     -- #
    ## -- rsq0.4, rsq0.8                                      -- #
    ## --------------------------------------------------------- #

    file_PLINK_geno_imp              = file_prefix_PLINK + ".genotyped.imputed"
    file_PLINK_geno_imp_assoc_rsq0_4 = file_PLINK_geno_imp +"_chr1-22.assoc.rsq0.4.locuszoom"
    file_PLINK_geno_imp_assoc_rsq0_8 = file_PLINK_geno_imp +"_chr1-22.assoc.rsq0.8.locuszoom"
    ##### vcf file tabixed from filtered PLINK data does not work, don't know why, but LD cannot be calculated by Locuszoom
    file_PLINK_geno_imp_rsq0_4 = file_PLINK_geno_imp + ".rsq0.4.chr1-22.locuszoom"
    file_PLINK_geno_imp_rsq0_8 = file_PLINK_geno_imp + ".rsq0.8.chr1-22.locuszoom"
    ####locuszoom_run(snp_colnr=2, file_suffix="_clump.clumped_groups.noXMHC", file_PLINKdosage=file_PLINK_geno_imp_rsq0_4, file_PLINK_rsq0_4=file_PLINK_rsq0_4)

    # snp_colnr starting from 0
    ## clumps all (without LD support)
    locuszoom_run(snp_colnr=2, file_suffix="_clump.clumped_all.noXMHC", file_PLINK_assoc=file_PLINK_geno_imp_assoc_rsq0_4, file_PLINK=file_PLINK_geno_imp_rsq0_4)
    # these jobs take a lot of time, perhaps clumps all is sufficient
    ####locuszoom_run(snp_colnr=2, file_suffix="_clump.clumped_all.noXMHC", file_PLINK_assoc=file_PLINK_geno_imp_assoc_rsq0_8, file_PLINK=file_PLINK_geno_imp_rsq0_8)
    ## clumps groups
    ####locuszoom_run(snp_colnr=2, file_suffix="_clump.clumped_groups.noXMHC", file_PLINK_assoc=file_PLINK_geno_imp_assoc_rsq0_4, file_PLINK=file_PLINK_geno_imp_rsq0_4)
    ####locuszoom_run(snp_colnr=2, file_suffix="_clump.clumped_groups.noXMHC", file_PLINK_assoc=file_PLINK_geno_imp_assoc_rsq0_8, file_PLINK=file_PLINK_geno_imp_rsq0_8)
