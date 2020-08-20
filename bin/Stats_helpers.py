#!/usr/bin/env python

import sys
import re
import os
import decimal

from os.path import *


# import Ingos lib
#sys.path.append(os.path.join(os.path.dirname[0], "../../all_scripts"))
sys.path.append(os.environ['PYLIB_DIR'] + "/all_scripts")
sys.path.append(os.environ['PYLIB_DIR'] + "/lib")
#from all_common import Command

# import my lib
# sys.path.append(join(sys.path[0], "../lib"))
# sys.path.append(os.environ['PYLIB_DIR'] + "/lib")

from plink_classes import Clump
from eigenstrat_classes import PackedPed


def extractQQplot_null_variants(qqplot_null_variants, assoc_input, assoc_output):
    """ extract null variants from assoc file """

    try:
        fh_r = file(qqplot_null_variants, "r")
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


def extract_Rsq_variants(assoc_logistic_input, assoc_dosage_input, assoc_merge_output_rsq0_3, assoc_merge_output_rsq0_8):
    """ extract imputed variants based on Rsq value """

    try:
        fh_w1 = file(assoc_merge_output_rsq0_8, "w")
        fh_w2 = file(assoc_merge_output_rsq0_3, "w")
    except IOError, e:
        print e
        sys.exit(1)

    store_lines_rsq0_3 = []
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
        list = re.split("\s+", line)
        if list[0] == "":
            del list[0]
        if list[-1] == "":
            del list[-1]

        # rs numbers instead of chr:pos in column SNP
        snp2line[list[0] + ":" + list[2]] = line
        chr  = decimal.Decimal(list[0])
        pos  = decimal.Decimal(list[2])
        store_lines_rsq0_3.append((chr, pos, line))
        store_lines_rsq0_8.append((chr, pos, line))
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
    fh_w1.writelines(header + "\n")
    fh_w2.writelines(header + "\n")
    line = fh_r.readline().rstrip('\n')
    while line:
        list = re.split("\s+", line)
        if list[0] == "":
            del list[0]
        if list[-1] == "":
            del list[-1]
        # if not already genotyped --> add
        if not (list[0] + ":" + list[2]) in snp2line:
            chr  = decimal.Decimal(list[0])
            pos  = decimal.Decimal(list[2])
            store_lines_rsq0_3.append((chr, pos, line))
            store_lines_rsq0_8.append((chr, pos, line))
#            if float(list[7]) >= 0.3:
#                chr  = decimal.Decimal(list[0])
#                pos  = decimal.Decimal(list[2])
#                store_lines_rsq0_3.append((chr, pos, line))
#                if float(list[7]) >= 0.8:
#                    store_lines_rsq0_8.append((chr, pos, line))
        line = fh_r.readline().rstrip('\n')
    fh_r.close()

    # sort merged files by chr, pos
    s_rsq0_3 = sorted(store_lines_rsq0_3, key=lambda tupel: tupel[1])
    t_rsq0_3 = sorted(s_rsq0_3, key=lambda tupel: tupel[0], reverse=False)
    s_rsq0_8 = sorted(store_lines_rsq0_8, key=lambda tupel: tupel[1])
    t_rsq0_8 = sorted(s_rsq0_8, key=lambda tupel: tupel[0], reverse=False)

    for tupel in t_rsq0_3:
        line = tupel[2]
        fh_w2.writelines(line + "\n")

    for tupel in t_rsq0_8:
        line = tupel[2]
        fh_w1.writelines(line + "\n")

    fh_w1.close()
    fh_w2.close()

    # save another version with chr:pos as SNPids for Locuszoom
    try:
        fh_w1 = file(assoc_merge_output_rsq0_8 + ".locuszoom", "w")
        fh_w2 = file(assoc_merge_output_rsq0_3 + ".locuszoom", "w")
    except IOError, e:
        print e
        sys.exit(1)

    fh_w1.writelines(header + "\n")
    fh_w2.writelines(header + "\n")

    check_duplicates_rsq0_3 = {}
    check_duplicates_rsq0_8 = {}

    for tupel in t_rsq0_3:
        line = tupel[2]
        list = re.split("\s+", line)
        if not check_duplicates_rsq0_3.has_key(list[0] + ":" + list[2]):
            fh_w2.writelines(list[0] + "\t" +
                             list[0] + ":" + list[2] + "\t" +
                             list[2] + "\t" +
                             list[3] + "\t" +
                             list[4] + "\t" +
                             list[5] + "\t" +
                             list[6] + "\t" +
                             list[7] + "\t" +
                             list[8] + "\t" +
                             list[9] + "\t" +
                             list[10] + "\n")
            check_duplicates_rsq0_3[list[0] + ":" + list[2]] = True

    for tupel in t_rsq0_8:
        line = tupel[2]
        list = re.split("\s+", line)
        if not check_duplicates_rsq0_8.has_key(list[0] + ":" + list[2]):
            fh_w1.writelines(list[0] + "\t" +
                             list[0] + ":" + list[2] + "\t"
                             + list[2] + "\t"
                             + list[3] + "\t"
                             + list[4] + "\t"
                             + list[5] + "\t"
                             + list[6] + "\t"
                             + list[7] + "\t"
                             + list[8] + "\t"
                             + list[9] + "\t"
                             + list[10] + "\n")
            check_duplicates_rsq0_8[list[0] + ":" + list[2]] = True

    fh_w1.close()
    fh_w2.close()


def excludeQQplot_variants(assoc_input, assoc_output, snpexclude):
    """ exclude specific regions from assoc null file """

    try:
        fh_r = file(snpexclude, "r")
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


def qqplot_xMHC_noxMHC(file_PLINK, file_PLINK_noxMHC, numof_cases, numof_controls, qqplotp, qqplotp2, pval2chisq, snpexclude, collection_name, qqman, qqman2):
    """ extract variants based on Rsq value """

    # ---------------------------- #
    # -- filtered for non-xMHC  -- #
    # -- xMHC: chr6, [25,34[ Mb -- #
    # ---------------------------- #

    print "\n    qq-plot from association results ...\n\n"
    os.system("R --slave --args %s %s %s %s < %s" \
        %(file_PLINK_noxMHC,\
        str(numof_cases),\
        str(numof_controls),\
        collection_name,\
        pval2chisq))

    print "\n    another qq-plotpval from association results ...\n\n"
    os.system("R --slave --args %s %s %s %s < %s" \
        %(file_PLINK_noxMHC,\
        str(numof_cases),\
        str(numof_controls),\
        collection_name,\
        qqplotp))

    print "\n    another qq-plotpval from association results ...\n\n"
    os.system("R --slave --args %s %s %s %s < %s" \
        %(file_PLINK_noxMHC,\
        str(numof_cases),\
        str(numof_controls),\
        collection_name,\
        qqplotp2))

    # -------------------------------------- #
    # -- including xMHC: chr6, [25,34[ Mb -- #
    # -------------------------------------- #

    print "\n    qq-plot from association results ...\n\n"
    os.system("R --slave --args %s %s %s %s < %s" \
        %(file_PLINK,\
        str(numof_cases),\
        str(numof_controls),\
        collection_name,\
        pval2chisq))

    print "\n    another qq-plotpval from association results ...\n\n"
    os.system("R --slave --args %s %s %s %s < %s" \
        %(file_PLINK,\
        str(numof_cases),\
        str(numof_controls),\
        collection_name,\
        qqplotp))

    print "\n    another qq-plotpval from association results ...\n\n"
    os.system("R --slave --args %s %s %s %s < %s" \
        %(file_PLINK,\
        str(numof_cases),\
        str(numof_controls),\
        collection_name,\
        qqplotp2))

    if snpexclude:

        excludeQQplot_variants(assoc_input=file_PLINK, \
                               assoc_output=file_PLINK + "_excludeRegions",\
                               snpexclude=snpexclude)

        excludeQQplot_variants(assoc_input=file_PLINK_noxMHC, \
                               assoc_output=file_PLINK_noxMHC + "_excludeRegions",\
                               snpexclude=snpexclude)

        # ---------------------------- #
        # -- filtered for non-xMHC  -- #
        # -- xMHC: chr6, [25,34[ Mb -- #
        # ---------------------------- #

        print "\n    qq-plot from association results ...\n\n"
        os.system("R --slave --args %s %s %s %s < %s" \
            %(file_PLINK_noxMHC + "_excludeRegions",\
            str(numof_cases),\
            str(numof_controls),\
            collection_name,\
            pval2chisq))

        print "\n    another qq-plotpval from association results ...\n\n"
        os.system("R --slave --args %s %s %s %s < %s" \
            %(file_PLINK_noxMHC + "_excludeRegions",\
            str(numof_cases),\
            str(numof_controls),\
            collection_name,\
            qqplotp))

        print "\n    another qq-plotpval from association results ...\n\n"
        os.system("R --slave --args %s %s %s %s < %s" \
            %(file_PLINK_noxMHC + "_excludeRegions",\
            str(numof_cases),\
            str(numof_controls),\
            collection_name,\
            qqplotp2))

        os.system("rm -f %s" %(file_PLINK_noxMHC + "_excludeRegions"))

        # -------------------------------------- #
        # -- including xMHC: chr6, [25,34[ Mb -- #
        # -------------------------------------- #

        print "\n    qq-plot from association results ...\n\n"
        os.system("R --slave --args %s %s %s %s < %s" \
            %(file_PLINK + "_excludeRegions",\
            str(numof_cases),\
            str(numof_controls),\
            collection_name,\
            pval2chisq))

        print "\n    another qq-plotpval from association results ...\n\n"
        os.system("R --slave --args %s %s %s %s < %s" \
            %(file_PLINK + "_excludeRegions",\
            str(numof_cases),\
            str(numof_controls),\
            collection_name,\
            qqplotp))

        print "\n    another qq-plotpval from association results ...\n\n"
        os.system("R --slave --args %s %s %s %s < %s" \
            %(file_PLINK + "_excludeRegions",\
            str(numof_cases),\
            str(numof_controls),\
            collection_name,\
            qqplotp2))



def qqplot_null(file_PLINK_null, numof_cases, numof_controls, qqplotp, qqplotp2, pval2chisq, snpexclude, collection_name, qqman, qqman2):
    """ extract null variants """

    # use normal qqplot scripts instead of null qqplot scripts because usage of
    # ld pruned variants as null SNPs
    print "\n    qq-plotpval for null SNVs from association results ...\n\n"
    os.system("R --no-save --args %s %s %s %s %s < %s" \
        %(file_PLINK_null,\
        str(numof_cases),\
        str(numof_controls),\
        collection_name,\
        qqman,\
        qqplotp))

    print "\n    qq-plotpval for null SNVs from association results ...\n\n"
    os.system("R --no-save --args %s %s %s %s %s < %s" \
        %(file_PLINK_null,\
        str(numof_cases),\
        str(numof_controls),\
        collection_name,\
        qqman2,\
        qqplotp2))

    print "\n    qq-plot for null SNVs from association results ...\n\n"
    os.system("R --no-save --args %s %s %s %s < %s" \
        %(file_PLINK_null,\
        str(numof_cases),\
        str(numof_controls),\
        collection_name,\
        pval2chisq))

    if snpexclude:

        excludeQQplot_variants(assoc_input=file_PLINK_null, \
                               assoc_output=file_PLINK_null + "_excludeRegions", \
                               snpexclude=snpexclude)

        print "\n    qq-plotpval for null SNVs with specific regions excluded from association results ...\n\n"
        os.system("R --no-save --args %s %s %s %s %s < %s" \
            %(file_PLINK_null + "_excludeRegions",\
            str(numof_cases),\
            str(numof_controls),\
            collection_name,\
            qqman,\
            qqplotp))

        os.system("R --no-save --args %s %s %s %s %s < %s" \
            %(file_PLINK_null + "_excludeRegions",\
            str(numof_cases),\
            str(numof_controls),\
            collection_name,\
            qqman2,\
            qqplotp2))


# rsq4: PLINK basename, INFO-filtered
# rsq8: PLINK basename, INFO-filtered for different threshold (optional, will not be processed if rsq8=="")
# file_PLINK_geno_imp_rsq*: association test report that was prepared for use with LocusZoom
# clump*: floating-point parameters for clumping, see PLINK manual
def PLINK_clumping(rsq4, rsq8, file_PLINK_geno_imp_rsq0_4, file_PLINK_geno_imp_rsq0_8, clumpr2, clumpp1, clumpp2, clumpkb):
    """ run clumping from PLINK association program """

#    file_PLINK_geno_imp              = file_prefix_PLINK + ".genotyped.imputed"
#    file_PLINK_geno_imp_rsq0_4       = file_PLINK_geno_imp +"_chr1-22.assoc.rsq0.4.locuszoom"
#    file_PLINK_geno_imp_rsq0_8       = file_PLINK_geno_imp +"_chr1-22.assoc.rsq0.8.locuszoom"

    # --------------------- #
    # -- clumping rsq0.4 -- #
    # --------------------- #
    print "\n    do clumping with PLINK ...\n\n"
    os.system("plink --bfile %s --out %s_clump --clump %s --clump-r2 %s --clump-p1 %s --clump-p2 %s --clump-kb %s --allow-no-sex --threads 16"\
#        %(file_PLINK_geno_imp +".rsq0.4.chr1-22.locuszoom",\
        %(rsq4,\
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

    if file_PLINK_geno_imp_rsq0_8 == "":
        return None

    # --------------------------------------------------------- #
    # -- clumping                                            -- #
    # -- rsq0.8                                              -- #
    # --------------------------------------------------------- #
    print "\n    do clumping with PLINK ...\n\n"

    os.system("plink --bfile %s --out %s_clump --clump %s --clump-r2 %s --clump-p1 %s --clump-p2 %s --clump-kb %s --allow-no-sex --threads 16"\
#        %(file_PLINK_geno_imp +".rsq0.8.chr1-22.locuszoom",\
        %(rsq8,\
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

def generateHitSpecfile(chromosome, pos, snp, GWAS_regions, file_HitSpec, QQplot_SNPexcludeList, disease_data_set_prefix_release_statistics):
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

    QQplot_SNPexcludeList_exists = os.path.isfile(QQplot_SNPexcludeList)


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
                os.system("/opt/locuszoom-1.3/bin/locuszoom --metal=%s --markercol SNP --pvalcol P --ld %s --build hg19 --hitspec=%s fineMap=\"%s\" --denote-markers-file %s --gwas-cat whole-cat_significant-only --plotonly --no-date --prefix %s" %(file_PLINK_assoc, file_PLINK_ld, file_HitSpec +".tmp", file_CredibleSets +".tmp", file_DenoteMarker +".tmp", prefix_out))
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
    os.system("/opt/locuszoom-1.3/bin/locuszoom --metal=%s --markercol SNP --pvalcol P --ld %s --build hg19 --hitspec=%s fineMap=\"%s\" --denote-markers-file %s --gwas-cat whole-cat_significant-only --plotonly --no-date --prefix %s" %(file_PLINK_assoc, file_PLINK_ld, file_HitSpec +".tmp", file_CredibleSets +".tmp", file_DenoteMarker +".tmp", prefix_out))
    #os.system("locuszoom --metal=%s --markercol SNP --pvalcol P --ld-vcf %s --build hg19 --hitspec=%s fineMap=\"%s\" --denote-markers-file %s --gwas-cat whole-cat_significant-only --plotonly --no-date --prefix %s" %(file_PLINK_assoc, join(Imputation_orig_dir, str(chromosome)+"."+disease_data_set_suffix_release_imputed +".gz"), file_HitSpec +".tmp", file_CredibleSets +".tmp", file_DenoteMarker +".tmp", prefix_out))
    # vcf file tabixed from filtered PLINK data does not work, don't know why, but LD cannot be calculated by Locuszoom
    #os.system("locuszoom --metal=%s --markercol SNP --pvalcol P --ld-vcf %s --build hg19 --hitspec=%s fineMap=\"%s\" --denote-markers-file %s --gwas-cat whole-cat_significant-only --plotonly --no-date --prefix %s" %(file_PLINK_assoc, file_PLINK_rsq0_4, file_HitSpec +".tmp", file_CredibleSets +".tmp", file_DenoteMarker +".tmp", prefix_out))
    os.system("rm -f %s" %(file_CredibleSets +".tmp"))
    os.system("rm -f %s" %(file_HitSpec +".tmp"))
    os.system("rm -f %s" %(file_DenoteMarker +".tmp"))
    fh_r1.close()
    fh_r2.close()
    fh_r3.close()




def extract_known_GWAS_regions(QQplot_SNPexcludeList):
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




def locuszoom_run(snp_colnr, file_suffix, file_PLINK_assoc, file_PLINK, numoftophits, QQplot_SNPexcludeList, name):
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

    plot_dir = ""

    # ---------------------------------------------------------------------- #
    # -- extract known GWAS region if file exists with known GWAS regions -- #
    # ---------------------------------------------------------------------- #

    QQplot_SNPexcludeList_exists = os.path.isfile(QQplot_SNPexcludeList)

    GWAS_regions = []
    if QQplot_SNPexcludeList_exists:
        GWAS_regions = extract_known_GWAS_regions(QQplot_SNPexcludeList)

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
        generateHitSpecfile(chromosome=chromosome, pos=pos, snp=snp, GWAS_regions=GWAS_regions, file_HitSpec=file_HitSpec, QQplot_SNPexcludeList=QQplot_SNPexcludeList, disease_data_set_prefix_release_statistics=name)

        # generate credibleSets file for displaying different disease subsets
        file_CredibleSets = join(plot_dir, snp +".credibleSets.csv")
        #####generateCredibleSets(chromosome=chromosome, pos=pos, snp=snp, sp2=sp2, pval=pval, file_CredibleSets=file_CredibleSets)
        generateCredibleSets(chromosome=chromosome, pos=pos, snp=snp, pval=pval, file_CredibleSets=file_CredibleSets)

        # generate denoteMarker file for displaying SNV IDs in plot
        file_DenoteMarker = join(plot_dir, snp +".denoteMarker.csv")
        generateDenoteMarker(chromosome=chromosome, pos=pos, snp=snp, file_DenoteMarker=file_DenoteMarker)

        # calculate LD on my own for --ld option in Locuszoom
        os.system("plink --bfile %s --r2 --ld-window-r2 0.0 --ld-snp %s --ld-window 1000000 --ld-window-kb 1000 --out %s_ld --allow-no-sex" %(file_PLINK, snp, join(plot_dir, snp)))
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
