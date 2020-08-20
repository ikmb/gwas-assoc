#!/usr/bin/env python

import sys
import re
import os

from os.path import join, dirname
# import string

# import gzip
# import math
# import decimal
# import datetime
# from os import listdir
# import subprocess

# may also need some of these:

# import Ingos lib
# sys.path.append(os.path.join(os.path.dirname[0], "../../all_scripts"))
sys.path.append(os.environ['PYLIB_DIR'] + "/all_scripts")
sys.path.append(os.environ['PYLIB_DIR'] + "/lib")
# from all_common import *

# import my lib
# sys.path.append(join(sys.path[0], "../lib"))
# sys.path.append(os.environ['PYLIB_DIR'] + "/lib")

from plink_classes import Test_missing,Frq
from eigenstrat_classes import PcaLog


def determine_unknown_diagnosis(annotationfile, outfile, diagnoses):
    """ determine samples with unknown diagnosis """

    print "\n        check for samples with unknown diagnosis ..\n\n"

    try:
        fh   = file(annotationfile, "r")
        fh_w = file(outfile, "w")
    except IOError, e:
        print e
        sys.exit(1)

    # header line
    line = fh.readline().rstrip('\n')
    list = re.split("\s+", line)
    # delete empty elements
    if list[0] == "":
        del list[0]
    if list[-1] == "":
        del list[-1]
    assert(list[8] == "diagnosis")

    # body
    line = fh.readline().rstrip('\n')
    while line:

        list = re.split("\s+", line)
        # delete empty elements
        if list[0] == "":
            del list[0]
        if list[-1] == "":
            del list[-1]

        diag = list[8]
        if not (diag in diagnoses):
            fh_w.writelines(line + "\n")

        line = fh.readline().rstrip('\n')

    fh.close()
    fh_w.close()


def extract_QCsamples_annotationfile_relativesfile(fam, individuals_annotation_QCed, related_samples_file, related_samples_file_QCed, individuals_annotation, diagnoses):
    """ extract only QCed samples from original annotation file """

    if isinstance(diagnoses, basestring):
        diagnoses = diagnoses.split(',')

    individualIDs = {}

    try:
        fh_r  = file(fam, "r")
    except IOError, e:
        print e
        sys.exit(1)

    line = fh_r.readline().rstrip('\n')
    while line:

        list = re.split("\s+", line)
        if list[0] == "":
            del list[0]

        if not (list[1] in individualIDs):
            individualIDs[list[1]] = True

        line = fh_r.readline().rstrip('\n')

    fh_r.close()

    # -------------------------- #
    # -- scan annotation file -- #
    # -------------------------- #
    try:
        fh_r  = file(individuals_annotation, "r")
        fh_w  = file(individuals_annotation_QCed, "w")
    except IOError, e:
        print e
        sys.exit(1)

    # read header
    line = fh_r.readline().rstrip('\n')
    fh_w.writelines(line + "\n")

    # read body
    count_samples = 0
    line = fh_r.readline().rstrip('\n')
    while line:

        list = re.split("\s+", line)
        if list[0] == "":
            del list[0]
        if list[-1] == "":
            del list[-1]

        if list[1] in individualIDs:
            fh_w.writelines(line + "\n")
            count_samples += 1

        line = fh_r.readline().rstrip('\n')

    fh_r.close()
    fh_w.close()
    assert(len(individualIDs) == count_samples)

    # generate diagnosis individual annotation files
    for diag in diagnoses:
        try:
            fh_r  = file(individuals_annotation, "r")
            # flake8 complains that neigher join nor dirname are defined...
            fh_w  = file(join(dirname(individuals_annotation_QCed), diag + "_QCed_annotation.txt"), "w")
        except IOError, e:
            print e
            sys.exit(1)

        # read header
        line = fh_r.readline().rstrip('\n')
        fh_w.writelines(line + "\n")

        # read body
        count_samples = 0
        line = fh_r.readline().rstrip('\n')
        while line:

            list = re.split("\s+", line)
            if list[0] == "":
                del list[0]
            if list[-1] == "":
                del list[-1]

            if list[8] == diag and (list[1] in individualIDs):
                fh_w.writelines(line + "\n")
                count_samples += 1

            line = fh_r.readline().rstrip('\n')

        fh_r.close()
        fh_w.close()

    # --------------------------------- #
    # -- scan flagged relatives file -- #
    # --------------------------------- #
    try:
        fh_r  = file(related_samples_file, "r")
        fh_w  = file(related_samples_file_QCed, "w")
    except IOError, e:
        print e
        sys.exit(1)

    # read body
    line = fh_r.readline().rstrip('\n')
    while line:

        list = re.split("\s+", line)
        if list[0] == "":
            del list[0]
        if list[-1] == "":
            del list[-1]

        if list[1] in individualIDs:
            fh_w.writelines(line + "\n")

        line = fh_r.readline().rstrip('\n')

    fh_r.close()
    fh_w.close()


def generate_exclude_file_CON(HFresults_file, individuals_annotation_QCed, variant_exclude_file, batches_names, prefix_merged_SNPQCII, FDR_index_remove_variants, plotscript):
    """ determine variants failed even with the worst batch removed """

    # IMPORTANT: Here exclude rejected variants from "all batches" instead of "worst batch removed" ###

    # THIS INDEX MUST BE SET BY USER
    # '4' corresponds to FDR at 1e-5, see list thresholds_FDR below
    # '6' corresponds to FDR at 1e-7, see list thresholds_FDR below
    # FDR_index_remove_variants = 4

    # --------------------------------------------------------- #
    # -- (1) count number of batches in QCed annotation file -- #
    # --------------------------------------------------------- #

    batches_dict = {}
    batches_list = []

    try:
        fh_r  = file(individuals_annotation_QCed, "r")
    except IOError, e:
        print e
        sys.exit(1)

    # read header
    line = fh_r.readline().rstrip('\n')
    list = re.split("\s+", line)
    assert(list[6] == "batch")

    # read body
    line = fh_r.readline().rstrip('\n')
    while line:

        list = re.split("\s+", line)
        if list[0] == "":
            del list[0]
        if list[-1] == "":
            del list[-1]

        if not list[6] in batches_dict:
            batches_dict[list[6]] = True
            batches_list.append(list[6])

        line = fh_r.readline().rstrip('\n')

    fh_r.close()

    numof_batches = len(batches_list)
    print "\n    Detected the following batches for HF test (n=" + str(numof_batches) + "): ..."
    for i in xrange(numof_batches):
        print "\n        " + batches_list[i]

    # QCed list of batches
    batches_list_dict = {}
    for i in xrange(numof_batches):
        batches_list_dict[batches_list[i]] = True

    # original list of batches
    print "\n    Following batches were removed: ..."
#    for i in xrange(len(batches_names)):
#        if not batches_names[i] in batches_list_dict:
#            print "\n        " + batches_names[i]

    # ------------------------------------------------------------ #
    # -- (2) read HF P-values from results file to control FDR -- #
    # ------------------------------------------------------------ #

    # TODO change HF to HF
    # HF_Pval_vectors
    HF_Pval_vector_allbatches_ctrls            = []  # (HF_p-value, variant_id)
    HF_Pval_vector_worstbatchremoved_ctrls     = []  # (HF_p-value, variant_id)

    try:
        fh_r  = file(HFresults_file, "r")
    except IOError, e:
        print e
        sys.exit(1)

    line = fh_r.readline().rstrip('\n')
    while line:

        list = re.split("\s+", line)
        if list[0] == "":
            del list[0]
        if list[-1] == "":
            del list[-1]

        # four preceding columns (CHR, SNP, POS, A1)
        numof_pvalues = len(list[4:])
        # numof_pvales - #phenotypes(CON,CD,UC) /  #phenotypes(CON,CD,UC)
        # numof_pvales - #phenotypes(CON) /  #phenotypes(CON)
        if not (numof_batches == (numof_pvalues - 1)):
            print >> sys.stderr, "abort: problem with results from splitted HF files, probably problem with jobs submitted to HP cluster."
            print >> sys.stderr, "    Expected #numof_pvalues=" + str(numof_batches - 1) + " in " + prefix_merged_SNPQCII + "_HF.auto.R"
            print >> sys.stderr, "    Observed #numof_pvalues=" + str(numof_batches - 1) + " in " + prefix_merged_SNPQCII + "_HF.auto.R"
            print >> sys.stderr, "abort: #batches_annotation=" + str(numof_batches) + " != " + str((numof_pvalues - 1)) + "=#batches_HFtest"
            sys.exit(1)

        # ----------------------------------------------------------- #
        # -- here the order of controls and diseases is important! -- #
        # ----------------------------------------------------------- #

        # (1) Controls
        if list[4] != "NA":
            HF_entire_collection_ctrls = float(list[4])
        else:
            HF_entire_collection_ctrls = 1.0

        # --------------------------------------------- #
        # -- look at p-values from entire collection -- #
        # --------------------------------------------- #

        variant_id = list[1]
        HF_Pval_vector_allbatches_ctrls.append((HF_entire_collection_ctrls, variant_id))

        # --------------------- #
        # -- look at batches -- #
        # --------------------- #

        # (1) Controls
        HF_max_excludebatch = 0.0
        count_NA = 0

        # find the highest p-value (when "worst" batch removed) when
        # running entire ctrl collection with one batch removed at one time
        for i in xrange(5, 5 + numof_batches, 1):
            if list[i] != "NA":
                HF_entire_collection_exclude_particularbatch = float(list[i])
                if HF_entire_collection_exclude_particularbatch > HF_max_excludebatch:
                    HF_max_excludebatch = HF_entire_collection_exclude_particularbatch
            else:
                count_NA += 1

        if numof_batches == count_NA:
                HF_max_excludebatch = 1.0

        HF_Pval_vector_worstbatchremoved_ctrls.append((HF_max_excludebatch, variant_id))

        line = fh_r.readline().rstrip('\n')

    fh_r.close()

    # ------------------------------------------------------------------- #
    # -- sort p-value vectors by first element of tuples, i.e. p-value -- #
    # ------------------------------------------------------------------- #
    HF_Pval_vector_allbatches_ctrls.sort(reverse=False)
    HF_Pval_vector_worstbatchremoved_ctrls.sort(reverse=False)

    assert(len(HF_Pval_vector_allbatches_ctrls) == len(HF_Pval_vector_worstbatchremoved_ctrls))

    # ---------------------------------------------------------------- #
    # -- count #variant failed at FDR at q=1e-1,1e-2,1e-3,...,1e-10 -- #
    # ---------------------------------------------------------------- #
    thresholds_FDR = [1e-1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8, 1e-9, 1e-10]

    counts_rejected_FDR_allbatches_ctrls           = [0 for i in range(10)]  # set count to 0
    counts_rejected_FDR_worstbatchremoved_ctrls    = [0 for i in range(10)]  # set count to 0

    # #### count total number of removeVariants from ctrls/CD/UC worstbatchremoved for each FDR
    # ###counts_rejected_FDR_allbatches_ctrls_CD_UC_cases        = [ {} for i in range(10) ] # add dictionaries
    # ###counts_rejected_FDR_worstbatchremoved_ctrls_CD_UC_cases = [ {} for i in range(10) ] # add dictionaries

    # fill this vector with HF_Pvalues at FDR thresholds
    thresholds_Pvals_allbatches_ctrls           = [float(0) for i in range(10)]
    thresholds_Pvals_worstbatchremoved_ctrls    = [float(0) for i in range(10)]

    # ------------------------------------------------------------------------------- #
    # -- calculate FDR for different FDR thresholds (Benjamini and Hochberg, 1995) -- #
    # ------------------------------------------------------------------------------- #

    # (a) ctrls - for all batches and worstbatchremoved
    n = len(HF_Pval_vector_allbatches_ctrls)
    for j in xrange(len(thresholds_FDR)):

        break_i_loop_part1 = False
        break_i_loop_part2 = False

        for i in xrange(1, n + 1, 1):

            rank = i / float(n)
            threshold = rank * thresholds_FDR[j]
            if (not break_i_loop_part1) and (HF_Pval_vector_allbatches_ctrls[i - 1][0] > threshold):
                thresholds_Pvals_allbatches_ctrls[j] = HF_Pval_vector_allbatches_ctrls[i - 2][0]
                counts_rejected_FDR_allbatches_ctrls[j] = i - 1

                # add variantIDs to dict for counting total number of
                # removeVariants from ctrls/CD/UC allbatches for each FDR
                # for k in xrange(i - 1):
                #    #counts_rejected_FDR_allbatches_ctrls_CD_UC_cases[j][HF_Pval_vector_allbatches_ctrls[k][1]] = True
                #    pass

                break_i_loop_part1 = True

            if (not break_i_loop_part2) and (HF_Pval_vector_worstbatchremoved_ctrls[i - 1][0] > threshold):
                thresholds_Pvals_worstbatchremoved_ctrls[j] = HF_Pval_vector_worstbatchremoved_ctrls[i - 2][0]
                counts_rejected_FDR_worstbatchremoved_ctrls[j] = i - 1
                break_i_loop_part2 = True

                # add variantIDs to dict for counting total number of
                # removeVariants from ctrls/CD/UC worstbatchremoved for each FDR
                # for k in xrange(i - 1):
                #    #counts_rejected_FDR_worstbatchremoved_ctrls_CD_UC_cases[j][HF_Pval_vector_worstbatchremoved_ctrls[k][1]] = True
                #    pass

            if break_i_loop_part1 and break_i_loop_part2:
                break

    # ---------------------------------------------------------------------------------------------- #
    # -- (3) extract rejected variants for FDRs at threshold with index FDR_index_remove_variants -- #
    # ---------------------------------------------------------------------------------------------- #

    # IMPORTANT: Here exclude rejected variants from "all batches" instead of "worst batch removed" ###

    try:
        fh_r  = file(HFresults_file, "r")
        fh_w  = file(variant_exclude_file, "w")
    except IOError, e:
        print e
        sys.exit(1)

    # THIS INDEX MUST BE SET BY USER at the beginning of this function
    # '4' corresponds to FDR at 1e-5
    # use variable FDR_index_remove_variants = 4 (per default)

    # raus
    dict_test = {}

    count_removeVariants_worstbatchremoved = 0
    line = fh_r.readline().rstrip('\n')
    while line:

        list = re.split("\s+", line)
        if list[0] == "":
            del list[0]
        if list[-1] == "":
            del list[-1]

        # four preceding columns (CHR, SNP, POS, A1)
        numof_pvalues = len(list[4:])
        # numof_pvales - #phenotypes(CON,CD,UC) /  #phenotypes(CON,CD,UC)
        # numof_pvales - #phenotypes(CON) /  #phenotypes(CON)
        if not (numof_batches == (numof_pvalues - 1)):
            print >> sys.stderr, "abort: problem with results from splitted HF files, probably problem with jobs submitted to HP cluster."
            print >> sys.stderr, "    Expected #numof_pvalues=" + str(numof_batches - 1) + " in " + prefix_merged_SNPQCII + "_HF.auto.R"
            print >> sys.stderr, "    Observed #numof_pvalues=" + str(numof_batches - 1) + " in " + prefix_merged_SNPQCII + "_HF.auto.R"
            print >> sys.stderr, "abort: #batches_annotation=" + str(numof_batches) + " != " + str((numof_pvalues - 1)) + "=#batches_HFtest"
            sys.exit(1)

        removeVariant = False

        # ----------------------------------------------------------- #
        # -- here the order of controls and diseases is important! -- #
        # ----------------------------------------------------------- #

        # (1) Controls
        if list[4] != "NA":
            HF_entire_collection_ctrls = float(list[4])
        else:
            HF_entire_collection_ctrls = 1.0

        # (1) Controls
        if HF_entire_collection_ctrls <= thresholds_Pvals_allbatches_ctrls[FDR_index_remove_variants]:

            # IMPORTANT: Here exclude rejected variants from "all batches" instead of "worst batch removed" ###
            removeVariant = True

            HF_max_excludebatch = 0.0

            # find the highest p-value (when "worst" batch removed) when
            # running entire ctrl collection with one batch removed at one time
            for i in xrange(5, 5 + numof_batches, 1):
                if list[i] != "NA":
                    HF_entire_collection_exclude_particularbatch = float(list[i])
                    if HF_entire_collection_exclude_particularbatch > HF_max_excludebatch:
                        HF_max_excludebatch = HF_entire_collection_exclude_particularbatch

            # if batch with smallest HF pvalue is removed AND
            # HF_max_excludebatch still below equal thresholds_Pvals_worstbatchremoved_ctrls[FDR_index_remove_variants], then remove variant.
            if HF_max_excludebatch <= thresholds_Pvals_worstbatchremoved_ctrls[FDR_index_remove_variants]:
                removeVariant = True

        if removeVariant:
            dict_test[list[1]] = True
            count_removeVariants_worstbatchremoved += 1
            fh_w.writelines(list[1] + "\n")

        line = fh_r.readline().rstrip('\n')

    fh_r.close()
    fh_w.close()

    # assert(count_removeVariants_worstbatchremoved == len(counts_rejected_FDR_worstbatchremoved_ctrls_CD_UC_cases[FDR_index_remove_variants]))

    # ----------------------------------------------------------------------------------- #
    # -- write #rejected variants for FDRs at different thresholds for plotting with R -- #
    # -- HF across entire collection and worstbatchremoved                             -- #
    # ----------------------------------------------------------------------------------- #

    try:
        fh_FDR_w = file(HFresults_file + ".FDRthresholds.SNPQCII.1.txt", "w")
    except IOError, e:
        print e
        sys.exit(1)

    # write header
    fh_FDR_w.writelines("FDR\tFail_allbatches_ctrls\tHF_pval_allbatches_ctrls\tFail_worstbatchremoved_ctrls\tHF_pval_worstbatchremoved_ctrls\n")

    for i in xrange(len(thresholds_FDR)):
        fh_FDR_w.writelines("%s" % (str(thresholds_FDR[i])))

        fh_FDR_w.writelines("\t%s" % (str(counts_rejected_FDR_allbatches_ctrls[i])))
        fh_FDR_w.writelines("\t%s" % (str(thresholds_Pvals_allbatches_ctrls[i])))
        fh_FDR_w.writelines("\t%s" % (str(counts_rejected_FDR_worstbatchremoved_ctrls[i])))
        fh_FDR_w.writelines("\t%s\n" % (str(thresholds_Pvals_worstbatchremoved_ctrls[i])))

    fh_FDR_w.close()

    # IMPORTANT: Here exclude rejected variants from "all batches" instead of "worst batch removed" ###

    # plot results applying FDR thresholds
    os.system("R --slave --args %s %s < %s"
              % (HFresults_file + ".FDRthresholds.SNPQCII",
                 str(FDR_index_remove_variants + 1),
                 plotscript))


def determine_pca_outlier(log, fam_file, outlier_file):
    """ determine outlier run eigenstrat program """

    pcalog  = PcaLog(input_file=log)
    # outlier dictionary key=outlier, val=sigmage
    outlier = pcalog.get_outlier()
    del pcalog

    fam = {}
    try:
        fh_r = file(fam_file, "r")
    except IOError, e:
        print e
        sys.exit(1)
    line = fh_r.readline()
    while line:
        list = re.split("\s+",line)
        fam[list[1]] = list[0]
        line = fh_r.readline()
    fh_r.close()

    try:
        fh_w = file(outlier_file, "w")
    except IOError, e:
        print e
        sys.exit(1)

    for indiv in outlier:
        #fh_w.writelines(indiv + "\n")
        fh_w.writelines(fam[indiv] +"\t" +indiv+ "\n")
    fh_w.close()

def addbatchinfo_32PCAs(fam, individuals_annotation, evec_file, eval_file, new_evec_file, new_eval_file):
    """ add batch information to final evec file """

    try:
      fh1 = file(fam, "r")
    except IOError, e:
      print e
      sys.exit(1)

    id2fid = {}
    line = fh1.readline().rstrip('\n')
    while line:
        list     = re.split("\s+",line)
        fid     = list[0]
        indivID = list[1]
        id2fid[indivID] = fid
        line = fh1.readline().rstrip('\n')
    fh1.close()

    try:
      fh1 = file(individuals_annotation, "r")
    except IOError, e:
      print e
      sys.exit(1)

    id2batch = {}

    line = fh1.readline().rstrip('\n')
    line = fh1.readline().rstrip('\n')
    while line:
        list     = re.split("\s+",line)
        indivID = list[1]
        batch   = list[6]
        id2batch[indivID] = batch
        line = fh1.readline().rstrip('\n')
    fh1.close()

    try:
      fh2 = file(evec_file, "r")
      fh3 = file(new_evec_file, "w")
    except IOError, e:
      print e
      sys.exit(1)

    line = fh2.readline().rstrip('\n')
    fh3.writelines("FID\tIID" +line.replace("indivID","") + "\tbatch\n")
    line = fh2.readline().rstrip('\n')
    while line:
        list = re.split("\s+",line)
        id = list[0]
        fh3.writelines(id2fid[id] +"\t"+\
              list[0] +"\t"+\
              list[1] +"\t"+\
              list[2] +"\t"+\
              list[3] +"\t"+\
              list[4] +"\t"+\
              list[5] +"\t"+\
              list[6] +"\t"+\
              list[7] +"\t"+\
              list[8] +"\t"+\
              list[9] +"\t"+\
              list[10] +"\t"+\
              list[11] +"\t"+\
              list[12] +"\t"+\
              list[13] +"\t"+\
              list[14] +"\t"+\
              list[15] +"\t"+\
              list[16] +"\t"+\
              list[17] +"\t"+\
              list[18] +"\t"+\
              list[19] +"\t"+\
              list[20] +"\t"+\
              list[21] +"\t"+\
              list[22] +"\t"+\
              list[23] +"\t"+\
              list[24] +"\t"+\
              list[25] +"\t"+\
              list[26] +"\t"+\
              list[27] +"\t"+\
              list[28] +"\t"+\
              list[29] +"\t"+\
              list[30] +"\t"+\
              list[31] +"\t"+\
              list[32] +"\t"+\
              id2batch[id] +"\n")
        line = fh2.readline().rstrip('\n')

    fh2.close()
    fh3.close()

    os.system("cp %s %s" %(eval_file, new_eval_file))

# This function does now require the HFresults_file and the individuals annotation to contain batches for exactly one disease.
# Please filter accordingly.
def generate_exclude_file_for_diagnosis(HFresults_file, individuals_annotation_QCed, variant_exclude_file, prefix_merged_SNPQCII, FDR_index_remove_variants, plotscript):
    """ determine variants failed even with the worst batch removed """

    ## IMPORTANT: Here exclude rejected variants from "worst batch removed" ###

    # THIS INDEX MUST BE SET BY USER
    # '4' corresponds to FDR at 1e-5, see list thresholds_FDR below
    # '6' corresponds to FDR at 1e-7, see list thresholds_FDR below
    ####FDR_index_remove_variants = 4

    # --------------------------------------------------------- #
    # -- (1) count number of batches in QCed annotation file -- #
    # --------------------------------------------------------- #

    batches_dict = {}
    batches_list = []

    try:
        fh_r  = file(individuals_annotation_QCed, "r")
    except IOError, e:
        print e
        sys.exit(1)

    # read header
    line = fh_r.readline().rstrip('\n')
    list = re.split("\s+",line)
    assert(list[6] == "batch")
    assert(list[8] == "diagnosis")

    # read body
    line = fh_r.readline().rstrip('\n')
    while line:

        list = re.split("\s+",line)
        if list[0] == "":
            del list[0]
        if list[-1] == "":
            del list[-1]

        if not batches_dict.has_key(list[6]):
            batches_dict[list[6]] = True
            batches_list.append(list[6])

        line = fh_r.readline().rstrip('\n')

    fh_r.close()

    numof_batches = len(batches_list)
    print "\n    Detected the following batches for HF test (n=" +str(numof_batches)+ "): ..."
    for i in xrange(numof_batches):
        print "\n        " +batches_list[i]

    # QCed list of batches
    batches_list_dict = {}
    for i in xrange(numof_batches):
        batches_list_dict[batches_list[i]] = True

    # original list of batches
    print "\n    Following batches were removed: ..."
    for i in xrange(len(batches_list)):
        if not batches_list_dict.has_key(batches_list[i]):
            print "\n        " +batches_list[i]

    # ------------------------------------------------------------ #
    # -- (2) read HF P-values from results file to control FDR -- #
    # ------------------------------------------------------------ #

    # HF_Pval_vectors
    #HF_Pval_vector_allbatches_ctrls            = [] # (HF_p-value, variant_id)
    HF_Pval_vector_allbatches_cases         = [] # (HF_p-value, variant_id)

    #HF_Pval_vector_worstbatchremoved_ctrls     = [] # (HF_p-value, variant_id)
    HF_Pval_vector_worstbatchremoved_cases  = [] # (HF_p-value, variant_id)


    try:
        fh_r  = file(HFresults_file, "r")
    except IOError, e:
        print e
        sys.exit(1)

    print("Scanning samples...\n")
    line = fh_r.readline().rstrip('\n')
    while line:

        list = re.split("\s+",line)
        if list[0] == "":
            del list[0]
        if list[-1] == "":
            del list[-1]

        # four preceding columns (CHR, SNP, POS, A1)
        numof_pvalues = len(list[4:])
        # numof_pvales - #phenotypes(CON,CD,UC) /  #phenotypes(CON,CD,UC)
        if not (numof_batches == (numof_pvalues-2)/2):
            print >> sys.stderr, "abort: problem with results from splitted HF files, probably problem with jobs submitted to HP cluster."
            print >> sys.stderr, "    Expected #numof_pvalues=" +str(numof_batches*2-2)+ " in " +HFresults_file
            print >> sys.stderr, "    Observed #numof_pvalues=" +str(numof_batches*2-2)+ " in " +HFresults_file
            print >> sys.stderr, "abort: #batches_annotation=" +str(numof_batches)+ " != " +str((numof_pvalues-2)/2)+ "=#batches_HFtest"
            sys.exit(1)
        #assert(numof_batches == (numof_pvalues-3)/3)

        # ----------------------------------------------------------- #
        # -- here the order of controls and diseases is important! -- #
        # ----------------------------------------------------------- #

        # (1) Controls
        #if list[4] != "NA":
        #    HF_entire_collection_ctrls = float(list[4])
        #else:
        #    HF_entire_collection_ctrls = 1.0

        ##### (2) PS
        ####if list[5] != "NA":
        ####    HF_entire_collection_PS_cases = float(list[5])
        ####else:
        ####    HF_entire_collection_PS_cases = 1.0

        ##### (3) AS
        ####if list[6] != "NA":
        ####    HF_entire_collection_AS_cases = float(list[6])
        ####else:
        ####    HF_entire_collection_AS_cases = 1.0

        # (4) Cases
        if list[4] != "NA":
            HF_entire_collection_cases = float(list[4])
        else:
            HF_entire_collection_cases = 1.0

        # (5) UC
        #if list[6] != "NA":
        #    HF_entire_collection_UC_cases = float(list[6])
        #else:
        #    HF_entire_collection_UC_cases = 1.0

        ##### (6) PSC
        ####if list[9] != "NA":
        ####    HF_entire_collection_PSC_cases = float(list[9])
        ####else:
        ####    HF_entire_collection_PSC_cases = 1.0


        # --------------------------------------------- #
        # -- look at p-values from entire collection -- #
        # --------------------------------------------- #

        variant_id = list[1]
    #    HF_Pval_vector_allbatches_ctrls.append( (HF_entire_collection_ctrls, variant_id) )
        HF_Pval_vector_allbatches_cases.append( (HF_entire_collection_cases, variant_id) )


        # --------------------- #
        # -- look at batches -- #
        # --------------------- #

        # (1) Controls
    #    HF_max_excludebatch = 0.0
    #    count_NA = 0

    #    # find the highest p-value (when "worst" batch removed) when
    #    # running entire ctrl collection with one batch removed at one time
    #    for i in xrange(6, 6+numof_batches, 1):
    #        if list[i] != "NA":
    #            HF_entire_collection_exclude_particularbatch = float(list[i])
    #            if HF_entire_collection_exclude_particularbatch > HF_max_excludebatch:
    #                HF_max_excludebatch = HF_entire_collection_exclude_particularbatch
    #        else:
    #            count_NA += 1

    #    if numof_batches == count_NA:
    #            HF_max_excludebatch = 1.0

    #    HF_Pval_vector_worstbatchremoved_ctrls.append( (HF_max_excludebatch, variant_id) )

        # (2) Cases
        HF_max_excludebatch = 0.0
        count_NA = 0

        # find the highest p-value (when "worst" batch removed) when
        # running entire case collection with one batch removed at one time
        for i in xrange(6+numof_batches, 6+2*numof_batches, 1):
            if list[i] != "NA":
                HF_entire_collection_exclude_particularbatch = float(list[i])
                if HF_entire_collection_exclude_particularbatch > HF_max_excludebatch:
                    HF_max_excludebatch = HF_entire_collection_exclude_particularbatch
            else:
                count_NA += 1

        if numof_batches == count_NA:
                HF_max_excludebatch = 1.0

        HF_Pval_vector_worstbatchremoved_cases.append( (HF_max_excludebatch, variant_id) )


        line = fh_r.readline().rstrip('\n')

    fh_r.close()

    # ------------------------------------------------------------------- #
    # -- sort p-value vectors by first element of tuples, i.e. p-value -- #
    # ------------------------------------------------------------------- #
    #HF_Pval_vector_allbatches_ctrls.sort(reverse=False)
    HF_Pval_vector_allbatches_cases.sort(reverse=False)

    #HF_Pval_vector_worstbatchremoved_ctrls.sort(reverse=False)
    HF_Pval_vector_worstbatchremoved_cases.sort(reverse=False)


    #assert(len(HF_Pval_vector_allbatches_ctrls) == len(HF_Pval_vector_allbatches_cases))

    #assert(len(HF_Pval_vector_allbatches_ctrls) == len(HF_Pval_vector_worstbatchremoved_ctrls))
    #assert(len(HF_Pval_vector_allbatches_ctrls) == len(HF_Pval_vector_worstbatchremoved_cases))


    # ---------------------------------------------------------------- #
    # -- count #variant failed at FDR at q=1e-1,1e-2,1e-3,...,1e-10 -- #
    # ---------------------------------------------------------------- #
    thresholds_FDR = [ 1e-1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8, 1e-9, 1e-10 ]

    #counts_rejected_FDR_allbatches_ctrls           = [ 0 for i in range(10) ] # set count to 0
    counts_rejected_FDR_allbatches_cases        = [ 0 for i in range(10) ] # set count to 0

    #counts_rejected_FDR_worstbatchremoved_ctrls    = [ 0 for i in range(10) ] # set count to 0
    counts_rejected_FDR_worstbatchremoved_cases = [ 0 for i in range(10) ] # set count to 0


    # count total number of removeVariants from ctrls/cases worstbatchremoved for each FDR
    counts_rejected_FDR_allbatches_ctrls_cases        = [ {} for i in range(10) ] # add dictionaries
    counts_rejected_FDR_worstbatchremoved_ctrls_cases = [ {} for i in range(10) ] # add dictionaries

    # fill this vector with HF_Pvalues at FDR thresholds
    #thresholds_Pvals_allbatches_ctrls           = [ float(0) for i in range(10) ]
    thresholds_Pvals_allbatches_cases        = [ float(0) for i in range(10) ]
    #thresholds_Pvals_worstbatchremoved_ctrls    = [ float(0) for i in range(10) ]
    thresholds_Pvals_worstbatchremoved_cases = [ float(0) for i in range(10) ]

    print("Calculating FDR for FDR thresholds...\n")
    # ------------------------------------------------------------------------------- #
    # -- calculate FDR for different FDR thresholds (Benjamini and Hochberg, 1995) -- #
    # ------------------------------------------------------------------------------- #

    # (a) ctrls - for all batches and worstbatchremoved

    #print(HF_Pval_vector_worstbatchremoved_ctrls)
    # (b) cases - for all batches and worstbatchremoved
    n = len (HF_Pval_vector_allbatches_cases)
    for j in xrange(len(thresholds_FDR)):

        break_i_loop_part1 = False
        break_i_loop_part2 = False

        for i in xrange(1,n+1,1):

            rank = i/float(n)
            threshold = rank*thresholds_FDR[j]
            if (not break_i_loop_part1) and (HF_Pval_vector_allbatches_cases[i-1][0] > threshold):
                thresholds_Pvals_allbatches_cases[j] = HF_Pval_vector_allbatches_cases[i-2][0]
                counts_rejected_FDR_allbatches_cases[j] = i-1

                # add variantIDs to dict for counting total number of
                # removeVariants from ctrls/CD/UC allbatches for each FDR
                for k in xrange(i-1):
                    counts_rejected_FDR_allbatches_ctrls_cases[j][HF_Pval_vector_allbatches_cases[k][1]] = True

                break_i_loop_part1 = True

            if (not break_i_loop_part2) and (HF_Pval_vector_worstbatchremoved_cases[i-1][0] > threshold):
                thresholds_Pvals_worstbatchremoved_cases[j] = HF_Pval_vector_worstbatchremoved_cases[i-2][0]
                counts_rejected_FDR_worstbatchremoved_cases[j] = i-1
                break_i_loop_part2 = True

                # add variantIDs to dict for counting total number of
                # removeVariants from ctrls/CD/UC worstbatchremoved for each FDR
                for k in xrange(i-1):
                    counts_rejected_FDR_worstbatchremoved_ctrls_cases[j][HF_Pval_vector_worstbatchremoved_cases[k][1]] = True

            if break_i_loop_part1 and break_i_loop_part2:
                break

    print "Case sample rejections for all batches by FDR: "
    print(counts_rejected_FDR_allbatches_cases)
    print "Case sample rejections for all but the worst batch by FDR:"
    print(counts_rejected_FDR_worstbatchremoved_cases)

    #print(HF_Pval_vector_worstbatchremoved_cases)
    print "Extracting rejected variants...\n"
    # ---------------------------------------------------------------------------------------------- #
    # -- (3) extract rejected variants for FDRs at threshold with index FDR_index_remove_variants -- #
    # ---------------------------------------------------------------------------------------------- #

    ## IMPORTANT: Here exclude rejected variants from "worst batch removed" ###

    try:
        fh_r  = file(HFresults_file, "r")
        fh_w  = file(variant_exclude_file, "w")
    except IOError, e:
        print e
        sys.exit(1)

    # THIS INDEX MUST BE SET BY USER at the beginning of this function
    # '4' corresponds to FDR at 1e-5
    # use variable FDR_index_remove_variants = 4 (per default)

    # raus
    dict_test = {}

    count_removeVariants_worstbatchremoved = 0
    line = fh_r.readline().rstrip('\n')
    while line:

        list = re.split("\s+",line)
        if list[0] == "":
            del list[0]
        if list[-1] == "":
            del list[-1]

        # four preceding columns (CHR, SNP, POS, A1)
        numof_pvalues = len(list[4:])
        # numof_pvales - #phenotypes(CON,CD,UC) /  #phenotypes(CON,CD,UC)
        if not (numof_batches == (numof_pvalues-2)/2):
            print >> sys.stderr, "abort: problem with results from splitted HF files, probably problem with jobs submitted to HP cluster."
            print >> sys.stderr, "    Expected #numof_pvalues=" +str(numof_batches*2-2)+ " in " +prefix_merged_SNPQCII+ "_HF.auto.R"
            print >> sys.stderr, "    Observed #numof_pvalues=" +str(numof_batches*2-2)+ " in " +prefix_merged_SNPQCII+ "_HF.auto.R"
            print >> sys.stderr, "abort: #batches_annotation=" +str(numof_batches)+ " != " +str((numof_pvalues-2)/2)+ "=#batches_HFtest"
            sys.exit(1)
        #assert(numof_batches == (numof_pvalues-3)/3)

        removeVariant = False

        # ----------------------------------------------------------- #
        # -- here the order of controls and diseases is important! -- #
        # ----------------------------------------------------------- #

        # (1) Controls
    #    if list[4] != "NA":
    #        HF_entire_collection_ctrls = float(list[4])
    #    else:
    #        HF_entire_collection_ctrls = 1.0
        # (4) cases
        if list[5] != "NA":
            HF_entire_collection_cases = float(list[4])
        else:
            HF_entire_collection_cases = 1.0

        # (1) Controls
    #   if HF_entire_collection_ctrls <= thresholds_Pvals_allbatches_ctrls[FDR_index_remove_variants]:

    #        HF_max_excludebatch = 0.0

    #        # find the highest p-value (when "worst" batch removed) when
    #        # running entire ctrl collection with one batch removed at one time
    #        for i in xrange(6, 6+numof_batches, 1):
    #            if list[i] != "NA":
    #                HF_entire_collection_exclude_particularbatch = float(list[i])
    #                if HF_entire_collection_exclude_particularbatch > HF_max_excludebatch:
    #                    HF_max_excludebatch = HF_entire_collection_exclude_particularbatch

    #        # if batch with smallest HF pvalue is removed AND
    #        # HF_max_excludebatch still below equal thresholds_Pvals_worstbatchremoved_ctrls[FDR_index_remove_variants], then remove variant.
    #        if HF_max_excludebatch <= thresholds_Pvals_worstbatchremoved_ctrls[FDR_index_remove_variants]:
    #            removeVariant = True

        # (2) cases
        if HF_entire_collection_cases <= thresholds_Pvals_allbatches_cases[FDR_index_remove_variants]:

            HF_max_excludebatch = 0.0

            # find the highest p-value (when "worst" batch removed) when
            # running entire ctrl collection with one batch removed at one time
            for i in xrange(6+numof_batches, 6+2*numof_batches, 1):
                if list[i] != "NA":
                    HF_entire_collection_exclude_particularbatch = float(list[i])
                    if HF_entire_collection_exclude_particularbatch > HF_max_excludebatch:
                        HF_max_excludebatch = HF_entire_collection_exclude_particularbatch

            # if batch with smallest HF pvalue is removed AND
            # HF_max_excludebatch still below equal thresholds_Pvals_worstbatchremoved_CD_cases[FDR_index_remove_variants], then remove variant.
            if HF_max_excludebatch <= thresholds_Pvals_worstbatchremoved_cases[FDR_index_remove_variants]:
                removeVariant = True

        if removeVariant:
            dict_test[list[1]] = True
            count_removeVariants_worstbatchremoved += 1
            fh_w.writelines(list[1] +"\n")

        line = fh_r.readline().rstrip('\n')

    fh_r.close()
    fh_w.close()

    #print("count_removeVariants_worstbatchremoved=%s counts_rejected_FDR_worstbatchremoved_ctrls_cases[FDR_index_remove_variants]=%s" %(count_removeVariants_worstbatchremoved, len(counts_rejected_FDR_worstbatchremoved_ctrls_cases[FDR_index_remove_variants])))
    #print(counts_rejected_FDR_worstbatchremoved_ctrls_cases)
    # schlaegt fehl, sollte es aber nicht! ...oder?
    #assert(count_removeVariants_worstbatchremoved == len(counts_rejected_FDR_worstbatchremoved_ctrls_cases[FDR_index_remove_variants]))
    print("Write rejected variants\n")
    # ----------------------------------------------------------------------------------- #
    # -- write #rejected variants for FDRs at different thresholds for plotting with R -- #
    # -- HF across entire collection and worstbatchremoved                             -- #
    # ----------------------------------------------------------------------------------- #

    try:
        fh_FDR_w = file(HFresults_file + ".FDRthresholds.SNPQCII.1.txt", "w")
    except IOError, e:
        print e
        sys.exit(1)

    # write header
    fh_FDR_w.writelines("FDR\tFail_allbatches_cases\tHF_pval_allbatches_cases\tFail_worstbatchremoved_cases\tHF_pval_worstbatchremoved_cases\n")

    for i in xrange(len(thresholds_FDR)):
        fh_FDR_w.writelines("%s" %(str(thresholds_FDR[i])))

        #fh_FDR_w.writelines("\t%s" %(str(len(counts_rejected_FDR_allbatches_ctrls_cases[i]))))
        #fh_FDR_w.writelines("\t%s" %(str(len(counts_rejected_FDR_worstbatchremoved_ctrls_cases[i]))))

        #fh_FDR_w.writelines("\t%s" %(str(counts_rejected_FDR_allbatches_ctrls[i])))
        #fh_FDR_w.writelines("\t%s" %(str(thresholds_Pvals_allbatches_ctrls[i])))
        #fh_FDR_w.writelines("\t%s" %(str(counts_rejected_FDR_worstbatchremoved_ctrls[i])))
        #fh_FDR_w.writelines("\t%s" %(str(thresholds_Pvals_worstbatchremoved_ctrls[i])))

        fh_FDR_w.writelines("\t%s" %(str(counts_rejected_FDR_allbatches_cases[i])))
        fh_FDR_w.writelines("\t%s" %(str(thresholds_Pvals_allbatches_cases[i])))
        fh_FDR_w.writelines("\t%s" %(str(counts_rejected_FDR_worstbatchremoved_cases[i])))
        fh_FDR_w.writelines("\t%s\n" %(str(thresholds_Pvals_worstbatchremoved_cases[i])))

    fh_FDR_w.close()
    print("plot\n")

    # plot results applying FDR thresholds
    os.system("R --slave --args %s %s < %s" \
        %(HFresults_file + ".FDRthresholds.SNPQCII",\
        str(FDR_index_remove_variants+1), \
        plotscript))
