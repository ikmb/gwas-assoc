#!/usr/bin/env python

import sys
import re
import os


def generate_exclude_file1_file2(HWEresults_file, batches_list, draw_script, all_output, perbatch_output, allbatches_output, FDR_index_remove_variants):
    """generate exclude file 1: From HWE calculations across the entire
    collection, remove variants for which HWE fails even if the worst batch
    removed (i.e. even if we remove the batch with the smallest p_value (below
    the HWE threshold), then the variant fails HWE for the entire collection
    without the particular batch

    generate exclude file 2: Remove variants (from HWE calculations by batch)
    failed HWE in >1 batches.
    """

    # THIS INDEX MUST BE SET BY USER
    # '4' corresponds to FDR at 1e-5, see list thresholds_FDR below

    # FDR_index_remove_variants = 4

    # -------------------------------------------------------- #
    # -- read HWE P-values from results file to control FDR -- #
    # -------------------------------------------------------- #

    numof_batches = len(batches_list)

    # HWE_Pval_vectors
    HWE_Pval_vector_allbatches        = []  # (HWE_p-value, variant_id)
    HWE_Pval_vector_worstbatchremoved = []  # (HWE_p-value, variant_id)

    # HWE_Pval_vector for each batch
    # 2D array
    HWE_Pval_vector_perbatch = []
    for b in xrange(numof_batches):
        HWE_Pval_vector_perbatch.append([])

    try:
        fh_r  = file(HWEresults_file, "r")
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

        numof_pvalues = (len(list[4:]) - 2) / 4
        assert(numof_batches == numof_pvalues)

        # ------------------------------- #
        # -- look at entire collection -- #
        # ------------------------------- #

        HWE_entire_collection = float(list[4])
        variant_id = list[1]
        HWE_Pval_vector_allbatches.append((HWE_entire_collection, variant_id))

        # --------------------- #
        # -- look at batches -- #
        # --------------------- #

        # find the ("worst") batch with the lowest p-value in HWE_particularbatch

        HWE_min_excludebatch = 1.0
        HWE_min_batch        = 1.0
        batch_index = 0
        for i in xrange(5, 5 + 2 * numof_batches, 2):
            HWE_particularbatch                           = float(list[i])
            HWE_entire_collection_exclude_particularbatch = float(list[i + 1])
            if HWE_particularbatch < HWE_min_batch:
                HWE_min_batch        = HWE_particularbatch
                HWE_min_excludebatch = HWE_entire_collection_exclude_particularbatch

            HWE_Pval_vector_perbatch[batch_index].append((HWE_particularbatch, variant_id))
            batch_index += 1

        HWE_Pval_vector_worstbatchremoved.append((HWE_min_excludebatch, variant_id))

        line = fh_r.readline().rstrip('\n')

    fh_r.close()

    # ------------------------------------------------------------------- #
    # -- sort p-value vectors by first element of tuples, i.e. p-value -- #
    # ------------------------------------------------------------------- #
    HWE_Pval_vector_allbatches.sort(reverse=False)
    HWE_Pval_vector_worstbatchremoved.sort(reverse=False)
    assert(len(HWE_Pval_vector_allbatches) == len(HWE_Pval_vector_worstbatchremoved))

    for b in xrange(numof_batches):
        HWE_Pval_vector_perbatch[b].sort(reverse=False)
        assert(len(HWE_Pval_vector_allbatches) == len(HWE_Pval_vector_perbatch[b]))

    # ---------------------------------------------------------------- #
    # -- count #variant failed at FDR at q=1e-1,1e-2,1e-3,...,1e-10 -- #
    # ---------------------------------------------------------------- #
    thresholds_FDR = [1e-1, 1e-2, 1e-3, 1e-4, 1e-5, 1e-6, 1e-7, 1e-8, 1e-9, 1e-10]

    counts_rejected_FDR_allbatches        = [0 for i in range(10)]  # set count to 0
    counts_rejected_FDR_worstbatchremoved = [0 for i in range(10)]  # set count to 0
    assert(len(counts_rejected_FDR_allbatches) == len(thresholds_FDR))

    # fill this vector with HWE_Pvalues at FDR thresholds
    thresholds_Pvals_allbatches        = [float(0) for i in range(10)]
    thresholds_Pvals_worstbatchremoved = [float(0) for i in range(10)]

    # 2D array, per batch
    counts_rejected_FDR_perbatch        = []
    for b in xrange(numof_batches):
        counts_rejected_FDR_perbatch.append([0 for i in range(10)])  # set count to 0

    # ------------------------------------------------------------------------------- #
    # -- calculate FDR for different FDR thresholds (Benjamini and Hochberg, 1995) -- #
    # ------------------------------------------------------------------------------- #

    # (a) for all batches and worstbatchremoved
    n = len(HWE_Pval_vector_allbatches)
    for j in xrange(len(thresholds_FDR)):

        break_i_loop_part1 = False
        break_i_loop_part2 = False

        for i in xrange(1, n + 1, 1):

            rank = i / float(n)
            threshold = rank * thresholds_FDR[j]
            if (not break_i_loop_part1) and (HWE_Pval_vector_allbatches[i - 1][0] > threshold):
                thresholds_Pvals_allbatches[j] = HWE_Pval_vector_allbatches[i - 2][0]
                counts_rejected_FDR_allbatches[j] = i - 1
                break_i_loop_part1 = True
            if (not break_i_loop_part2) and (HWE_Pval_vector_worstbatchremoved[i - 1][0] > threshold):
                thresholds_Pvals_worstbatchremoved[j] = HWE_Pval_vector_worstbatchremoved[i - 2][0]
                counts_rejected_FDR_worstbatchremoved[j] = i - 1
                break_i_loop_part2 = True
            if break_i_loop_part1 and break_i_loop_part2:
                break

    # (b) for each batch itself
    for j in xrange(len(thresholds_FDR)):
        for i in xrange(1, n + 1, 1):
            for b in xrange(numof_batches):
                rank = i / float(n)
                threshold = rank * thresholds_FDR[j]
                if HWE_Pval_vector_perbatch[b][i - 1][0] > threshold:
                    counts_rejected_FDR_perbatch[b][j] = i - 1
                    break

    # ------------ #
    # -- file 1 -- #
    # ------------ #

    # extract rejected variants for FDRs at threshold with index
    # FDR_index_remove_variants

    # THIS INDEX MUST BE SET BY USER at the beginning of this function
    # '4' corresponds to FDR at 1e-5
    # use variable FDR_index_remove_variants = 4 (per default)

    try:
        fh_worst_batch_removed_w  = file(all_output, "w")
        fh_allbatches_w = file(allbatches_output, "w")
    except IOError, e:
        print e
        sys.exit(1)

    sep = "\t"
    header = "Variant" + sep + "P_HWE\n"
    fh_worst_batch_removed_w.writelines(header)
    fh_allbatches_w.writelines(header)

    for i in xrange(counts_rejected_FDR_worstbatchremoved[FDR_index_remove_variants]):
        fh_worst_batch_removed_w.writelines(HWE_Pval_vector_worstbatchremoved[i][1] + sep +
                                            str(HWE_Pval_vector_worstbatchremoved[i][0]) + "\n")

    for i in xrange(counts_rejected_FDR_allbatches[FDR_index_remove_variants]):
        fh_allbatches_w.writelines(HWE_Pval_vector_allbatches[i][1] + sep +
                                            str(HWE_Pval_vector_allbatches[i][0]) + "\n")

    fh_worst_batch_removed_w.close()
    fh_allbatches_w.close()

    # ----------------------------------------------------------------------------------- #
    # -- write #rejected variants for FDRs at different thresholds for plotting with R -- #
    # -- HWE across entire collection and worstbatchremoved                            -- #
    # ----------------------------------------------------------------------------------- #

    try:
        fh_FDR_w = file(all_output + ".FDRthresholds.SNPQCI.1.txt", "w")
    except IOError, e:
        print e
        sys.exit(1)

    # write header
    fh_FDR_w.writelines("FDR\tFail_allbatches\tHWE_pval_allbatches\tFail_worstbatchremoved\tHWE_pval_worstbatchremoved\n")

    for i in xrange(len(thresholds_FDR)):
        fh_FDR_w.writelines("%s" % (str(thresholds_FDR[i])))
        fh_FDR_w.writelines("\t%s" % (str(counts_rejected_FDR_allbatches[i])))
        fh_FDR_w.writelines("\t%s" % (str(thresholds_Pvals_allbatches[i])))
        fh_FDR_w.writelines("\t%s" % (str(counts_rejected_FDR_worstbatchremoved[i])))
        fh_FDR_w.writelines("\t%s\n" % (str(thresholds_Pvals_worstbatchremoved[i])))

    fh_FDR_w.close()

    # ------------ #
    # -- file 2 -- #
    # ------------ #

    # ---------------------------------------------------------------------------------------------------------- #
    # -- for each batch: extract rejected variants for FDRs at threshold with index FDR_index_remove_variants -- #
    # ---------------------------------------------------------------------------------------------------------- #

    # THIS INDEX MUST BE SET BY USER at the beginning of this function
    # '4' corresponds to FDR at 1e-5
    # use again variable FDR_index_remove_variants = 4

    try:
        fh_failed2plusbatches_w  = file(perbatch_output, "w")
    except IOError, e:
        print e
        sys.exit(1)

    variant_exclude_dict_failed1plusbatches = {}
    variant_exclude_dict_failed2plusbatches = {}
    for b in xrange(numof_batches):
        for j in xrange(counts_rejected_FDR_perbatch[b][FDR_index_remove_variants]):
            if not HWE_Pval_vector_perbatch[b][j][1] in variant_exclude_dict_failed1plusbatches:
                variant_exclude_dict_failed1plusbatches[HWE_Pval_vector_perbatch[b][j][1]] = True
            elif HWE_Pval_vector_perbatch[b][j][1] in variant_exclude_dict_failed1plusbatches:
                variant_exclude_dict_failed2plusbatches[HWE_Pval_vector_perbatch[b][j][1]] = True

    for variant in variant_exclude_dict_failed2plusbatches.keys():
        fh_failed2plusbatches_w.writelines(variant + "\n")
    fh_failed2plusbatches_w.close()

    # ----------------------------------------------------------------------------------- #
    # -- write #rejected variants for FDRs at different thresholds for plotting with R -- #
    # -- HWE for each batch                                                            -- #
    # ----------------------------------------------------------------------------------- #

    try:
        fh_FDR_w  = file(perbatch_output + ".FDRthresholds.SNPQCI.2.txt", "w")
    except IOError, e:
        print e
        sys.exit(1)

    # write header
    fh_FDR_w.writelines("FDR\tFail_1plusbatches\tFail_2plusbatches\n")

    for i in xrange(len(thresholds_FDR)):

        variant_exclude_dict_failed1plusbatches = {}
        variant_exclude_dict_failed2plusbatches = {}
        for b in xrange(numof_batches):
            for j in xrange(counts_rejected_FDR_perbatch[b][i]):
                if not HWE_Pval_vector_perbatch[b][j][1] in variant_exclude_dict_failed1plusbatches:
                    variant_exclude_dict_failed1plusbatches[HWE_Pval_vector_perbatch[b][j][1]] = True
                elif HWE_Pval_vector_perbatch[b][j][1] in variant_exclude_dict_failed1plusbatches:
                    variant_exclude_dict_failed2plusbatches[HWE_Pval_vector_perbatch[b][j][1]] = True

        fh_FDR_w.writelines("%s" % (str(thresholds_FDR[i])))
        fh_FDR_w.writelines("\t%s" % (str(len(variant_exclude_dict_failed1plusbatches))))
        fh_FDR_w.writelines("\t%s\n" % (str(len(variant_exclude_dict_failed2plusbatches))))
        # do not write table of corresponding HWE_Pvalues thresholds for FDR
        # thresholds, because every batch has a different corresponding HWE_p-value

    fh_FDR_w.close()

    # plot results applying FDR thresholds
    os.system("R --slave --args %s %s %s < %s"
              % (all_output + ".FDRthresholds.SNPQCI.1.txt", perbatch_output + ".FDRthresholds.SNPQCI.2.txt",
                 str(FDR_index_remove_variants + 1),
                 draw_script))

    # TODO: include SNP_QCI_draw_FDR_CON_PS_AS_CD_UC_PSC.r as a special case


# Main
if __name__ == "__main__":

    # check args
    if len(sys.argv) != 8:
        print "Usage: " + sys.argv[0] + " <hwe-file> <individuals-annotation> <fdr-draw-script.R> <excluded-all> <excluded-perbatch> <excluded-allbatches> <fdr-index-remove-variants>\n"
        print "\twhere:\n"
        print "\t<hwe-file> HWE result file written by the Plink/Rserve HWE calculation template\n"
        print "\t<individuals-annotation> individuals annotation file that contains batches information\n"
        print "\t<fdr-draw-script.R> path to R script that draws the FDR values\n"
        print "\t<excluded-all> output file for variants that fail whole-collection HWE even if the worst batch is removed\n"
        print "\t<excluded-perbatch> output file for variants that fail HWE in at least two batches\n"
        print "\t..."
        sys.exit(1)

    # collect batches
    batches_dict = {}
    batches_list = []

    try:
        individuals_fh = file(sys.argv[2], "r")
    except IOError, e:
        print e
        sys.exit(1)

    line = individuals_fh.readline().rstrip('\n')
    list = re.split("\s+", line)
    # assert annotation file validity
    assert(list[6] == "batch")
    assert(list[8] == "diagnosis")

    line = individuals_fh.readline().rstrip('\n')
    while line:
        list = re.split("\s+", line)
        if list[0] == "":
            del list[0]
        if list[-1] == "":
            del list[-1]
        if list[8] != "Control":
            line = individuals_fh.readline().rstrip('\n')
            continue

        if not list[6] in batches_dict:
            batches_dict[list[6]] = True
            batches_list.append(list[6])

        line = individuals_fh.readline().rstrip('\n')
    individuals_fh.close()

    # generate exclusion files
    generate_exclude_file1_file2(sys.argv[1], batches_list, sys.argv[3], sys.argv[4], sys.argv[5], sys.argv[6], int(sys.argv[7]))
