#!/usr/bin/env python

import sys
import re
import os

from os.path import *
#import string

# import gzip
# import math
# import decimal
# import datetime
# from os import listdir
# import subprocess

# may also need some of these:

# import Ingos lib
#sys.path.append(os.path.join(os.path.dirname[0], "../../all_scripts"))
sys.path.append(os.environ['PYLIB_DIR'] + "/all_scripts")
sys.path.append(os.environ['PYLIB_DIR'] + "/lib")
from all_common import Command

# import my lib
# sys.path.append(join(sys.path[0], "../lib"))
# sys.path.append(os.environ['PYLIB_DIR'] + "/lib")

#from plink_classes import PackedPed
from eigenstrat_classes import PackedPed


def extract_QCsamples_from_pc_file(pc_file_old, pc_file_new, fam_file):
    """ extract only QCed samples from original pca file """

    individualIDs = {}

    # ------------------- #
    # -- scan fam_file -- #
    # ------------------- #
    try:
        fh_r  = file(fam_file, "r")
    except IOError, e:
        print e
        sys.exit(1)

    line = fh_r.readline().rstrip('\n')
    while line:

        list = re.split("\s+",line)
        if list[0] == "":
            del list[0]

        if not list[1] in individualIDs:
            individualIDs[list[1]] = True

        line = fh_r.readline().rstrip('\n')

    fh_r.close()

    # --------------------------- #
    # -- scan original pc file -- #
    # --------------------------- #
    try:
        fh_r  = file(pc_file_old, "r")
        fh_w  = file(pc_file_new, "w")
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


def extract_QCsamples_from_annotationfile(pc_file, individuals_annotation, individuals_annotation_QCed):
    """ extract only QCed samples from original annotation file """

    individualIDs = {}

    # ------------------ #
    # -- scan pc_file -- #
    # ------------------ #
    try:
        fh_r  = file(pc_file, "r")
    except IOError, e:
        print e
        sys.exit(1)

    # skip header
    line = fh_r.readline().rstrip('\n')

    line = fh_r.readline().rstrip('\n')
    while line:

        list = re.split("\s+", line)
        if list[0] == "":
            del list[0]

        if not list[1] in individualIDs:
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
    print (len(individualIDs))
    print (count_samples)
    assert(len(individualIDs) == count_samples)


def determine_pca_outlier(log, outlier_EIGENSTRAT_file):
    """ determine outlier run eigenstrat program """

    pcalog  = PcaLog(input_file=log)
    # outlier dictionary key=outlier, val=sigmage
    outlier = pcalog.get_outlier()
    del pcalog

    try:
        fh_w = file(outlier_EIGENSTRAT_file, "w")
    except IOError, e:
        print e
        sys.exit(1)

    for indiv in outlier:
        #fh_w.writelines(indiv + "\n")
        fh_w.writelines("0\t" +indiv+ "\n")
    fh_w.close()

def detect_duplicates_related_individuals_by_pihat(new_plink_pruned,  # currently unised
                                                   new_plink_IBS,     # expects ".genome" file, no basenames
                                                   new_plink_miss,    # expects ".imiss" file, no basenames
                                                   duplicate_threshold, # string or float
                                                   relative_threshold,
                                                   plot_script,       #join(path_source_code_env, "sources", "IBD-plot-genomefile.r")
                                                   dest_relatives,    # target filename for relative flags
                                                   dest_duplicates,    # target filename for duplicate list
                                                   batch_mode_par,    # batch mode
                                                   ):
    """ flag individuals exceeding maxmimal ibd """
    batch_mode = batch_mode_par == "True"

    print "\n    detect related and duplicated individuals ...\n\n"

    max_ibd_threshold_duplicates = float(duplicate_threshold)
    max_ibd_threshold_relatives = float(relative_threshold)

    comment_pattern = re.compile("^#.*$")
    blankline_pattern = re.compile("^\s*$")

    ## on local node
    #new_plink_pruned2 = ""
    #if batch_mode:
    #    new_plink_pruned2 = join(tmpdir, basename(new_plink_pruned))
    #    new_plink_IBS2    = join(tmpdir, basename(new_plink_IBS))
    #    new_plink_miss2   = join(tmpdir, basename(new_plink_miss))
    #else:
    #    new_plink_pruned2 = new_plink_pruned
    #    new_plink_IBS2    = new_plink_IBS
    #    new_plink_miss2   = new_plink_miss
    new_plink_pruned2 = new_plink_pruned
    new_plink_IBS2    = new_plink_IBS
    new_plink_miss2   = new_plink_miss

    # don't run this R script with more than 80,000 samples, it will crash
    # IBD plot - preQC
    os.system("gawk '{ print $7, $8 }' %s > %s.Z0.Z1" \
        %(new_plink_IBS2, \
          new_plink_IBS2))

    os.system("R --slave --args %s.Z0.Z1 < %s" \
        %(new_plink_IBS2, plot_script))

    if batch_mode:
        os.system("cp %s.IBD-plot.png %s.IBD-plot.png" %(new_plink_IBS2,new_plink_IBS))

    # genome file
    try:
        fh_gen = file(new_plink_IBS2, "r")
    except IOError, e:
        print e
        sys.exit(1)

    dict_aboveth = {} # store IDs above threshold
    dict_aboveth_ids_counts = {} # store counts of IDs above threshold

    line = fh_gen.readline().replace("\n", "") # skip first line
    line = fh_gen.readline().replace("\n", "")
    while line:

      # skip comment lines that start with "#"
      if(comment_pattern.search(line)):
        line = fh_gen.readline().replace("\n", "")
        continue
      # skip blank lines
      if(blankline_pattern.search(line)):
        line = fh_gen.readline().replace("\n", "")
        continue

      list = re.split("\s+",line)

      if list[0] == "":
          del list[0]

      FID1 = list[0]
      IID1 = list[1]
      FID2 = list[2]
      IID2 = list[3]

      key = FID1 + IID1 + FID2 + IID2
      pihat = float(list[9])

      # store values above IBD (PI_HAT=proportion IBD) threshold, often used
      # thresholds are:
      # PI_HAT = 0.25 : second-degree r elatives
      # PI_HAT = 0.5  : first-degree relatives
      # PI_HAT > 0.98 : duplicates or identical twins
      if( pihat >= max_ibd_threshold_relatives ):
        dict_aboveth[key] = (FID1 , IID1 , FID2 , IID2, pihat)

        if dict_aboveth_ids_counts.has_key(FID1 + IID1):
            dict_aboveth_ids_counts[FID1 + IID1] += 1
        else:
            dict_aboveth_ids_counts[FID1 + IID1] = 1

        if dict_aboveth_ids_counts.has_key(FID2 + IID2):
            dict_aboveth_ids_counts[FID2 + IID2] += 1
        else:
            dict_aboveth_ids_counts[FID2 + IID2] = 1

      line = fh_gen.readline().replace("\n", "")

    fh_gen.close()

    # imiss file
    try:
        fh_imiss = file(new_plink_miss2, "r")
    except IOError, e:
        print e
        sys.exit(1)

    dict_imiss = {}
    line = fh_imiss.readline().replace("\n", "")
    line = fh_imiss.readline().replace("\n", "") # skip first line
    while line:

        if(comment_pattern.search(line)):
          line = fh_imiss.readline().replace("\n", "")
          continue
        if(blankline_pattern.search(line)):
          line = fh_imiss.readline().replace("\n", "")
          continue

        list = re.split("\s+",line)

        if list[0] == "":
            del list[0]

        FID1 = list[0]
        IID1 = list[1]
        F_MISS = list[5]

        key = FID1 + IID1
        dict_imiss[key] = F_MISS

        line = fh_imiss.readline().replace("\n", "")

    fh_imiss.close()

    remove_dict_relatives  = {}
    remove_dict_duplicates = {}
    sep = " "

    thekeys = dict_aboveth.iterkeys()
    for k in thekeys:
        indiv1 = dict_aboveth[k][0] + dict_aboveth[k][1]
        indiv2 = dict_aboveth[k][2] + dict_aboveth[k][3]
        pihat = dict_aboveth[k][4]
        FMISS1 = dict_imiss[indiv1]
        FMISS2 = dict_imiss[indiv2]

        # first criterium: number of counts
        indiv1_counts = dict_aboveth_ids_counts[indiv1]
        indiv2_counts = dict_aboveth_ids_counts[indiv2]

        if indiv1_counts > indiv2_counts:

            # duplicates
            if( pihat >= max_ibd_threshold_duplicates ):

                remove_dict_duplicates[dict_aboveth[k][0] +sep+\
                    dict_aboveth[k][1]] = pihat

            # relatives
            else:

                remove_dict_relatives[dict_aboveth[k][0] +sep+\
                    dict_aboveth[k][1]] = pihat


        elif indiv1_counts < indiv2_counts:

            # duplicates
            if( pihat >= max_ibd_threshold_duplicates ):

                remove_dict_duplicates[dict_aboveth[k][2] +sep+\
                    dict_aboveth[k][3]] = pihat

            # relatives
            else:

                remove_dict_relatives[dict_aboveth[k][2] +sep+\
                    dict_aboveth[k][3]] = pihat

        # same number of counts
        else:

            # second criterium: missingness
            if(FMISS1 > FMISS2):

                # duplicates
                if( pihat >= max_ibd_threshold_duplicates ):

                    remove_dict_duplicates[dict_aboveth[k][0] +sep+\
                        dict_aboveth[k][1]] = pihat

                # relatives
                else:

                    remove_dict_relatives[dict_aboveth[k][0] +sep+\
                        dict_aboveth[k][1]] = pihat

            else:

                # duplicates
                if( pihat >= max_ibd_threshold_duplicates ):

                    remove_dict_duplicates[dict_aboveth[k][2] +sep+\
                        dict_aboveth[k][3]] = pihat

                # relatives
                else:

                    remove_dict_relatives[dict_aboveth[k][2] +sep+\
                        dict_aboveth[k][3]] = pihat

    try:
        fh_w1 = file(dest_duplicates, "w")
        fh_w2 = file(dest_relatives, "w")
    except IOError, e:
        print e
        sys.exit(1)

    # duplicates
    thekeys = remove_dict_duplicates.iterkeys()
    for k in thekeys:
        # print fam_id, indiv_id, ibs-pihat
        fh_w1.writelines(k +sep+ str(remove_dict_duplicates[k]) +"\n")

    # relatives
    thekeys = remove_dict_relatives.iterkeys()
    for k in thekeys:
        # print fam_id, indiv_id, ibs-pihat
        fh_w2.writelines(k +sep+ str(remove_dict_relatives[k]) +"\n")

    fh_w1.close()
    fh_w2.close()

    #if batch_mode:
    #    os.system("cp %s_remove.duplicates.txt %s_remove.duplicates.txt"\
    #    %(new_plink_pruned2,new_plink_pruned))
    #    os.system("cp %s_flag.relatives.txt %s_flag.relatives.txt"\
    #    %(new_plink_pruned2,new_plink_pruned))


def merge__new_plink_collection_pruned__1kG(new_plink_pruned, new_plink_pruned_1kG, PCA_SNPexcludeList, preQCIMDS_1kG):
    """ merge pruned data set with samples from 1kG for PCA """

    # convert SNP identifier to chr:pos
    os.system("cp %s.bed %s.chrpos.bed" %(new_plink_pruned, new_plink_pruned))
    os.system("cp %s.fam %s.chrpos.fam" %(new_plink_pruned, new_plink_pruned))
    os.system("gawk '{ print $1, $1\":\"$4, $3, $4, $5, $6 }' %s.bim > %s.chrpos.bim" \
        %(new_plink_pruned, new_plink_pruned))
    os.system("gawk '{ print $2 }' %s.chrpos.bim > %s.chrpos.bim.txt" \
        %(new_plink_pruned, new_plink_pruned))

    # calculate overlap in SNPs
    # remove some variants if specified in file PCA_SNPexcludeList, SNPIDs in chr:pos format required
    if os.path.isfile(PCA_SNPexcludeList):

        cmd = Command("plink --memory 10000 --bfile %s.chrpos --extract %s.bim.chrpos --exclude %s --make-bed --out %s --allow-no-sex" \
                   %(new_plink_pruned,\
                     preQCIMDS_1kG,\
                     PCA_SNPexcludeList,\
                     new_plink_pruned + "_tmp") )
        cmd.run() ; del cmd
        cmd = Command("plink --nmemory 10000 --bfile %s --extract %s.chrpos.bim.txt --exclude %s --make-bed --out %s --allow-no-sex" \
                   %(preQCIMDS_1kG,\
                     new_plink_pruned,\
                     PCA_SNPexcludeList,\
                     os.path.basename(preQCIMDS_1kG + "_tmp")) )
        cmd.run() ; del cmd
    else:
        cmd = Command("plink --memory 10000 --bfile %s.chrpos --extract %s.bim.chrpos --make-bed --out %s --allow-no-sex" \
                   %(new_plink_pruned,\
                     preQCIMDS_1kG,\
                     new_plink_pruned + "_tmp") )
        cmd.run() ; del cmd
        cmd = Command("plink --memory 10000 --bfile %s --extract %s.chrpos.bim.txt --make-bed --out %s --allow-no-sex" \
                   %(preQCIMDS_1kG,\
                     new_plink_pruned,\
                     os.path.basename(preQCIMDS_1kG + "_tmp")) )
        cmd.run() ; del cmd

    # merge
    cmd = Command("plink --memory 10000 --bfile %s --bmerge %s.bed %s.bim %s.fam --make-bed --out %s --allow-no-sex" \
               %(new_plink_pruned + "_tmp",\
                 os.path.basename(preQCIMDS_1kG + "_tmp"),\
                 os.path.basename(preQCIMDS_1kG + "_tmp"),\
                 os.path.basename(preQCIMDS_1kG + "_tmp"),\
                 new_plink_pruned_1kG) )
    cmd.run() ; del cmd

    # remove tmp file
    os.system("rm -f %s.* %s.* %s.*" %(new_plink_pruned + ".chrpos", new_plink_pruned + "_tmp", basename(preQCIMDS_1kG + "_tmp")))



def addbatchinfo_10PCs(evec_file, eval_file, new_evec_file, new_eval_file, individuals_annotation, preQCIMDS_1kG_sample):
    """ add batch information to final evec file """

    id2batch = {}

    try:
        fh1 = file(individuals_annotation, "r")
    except IOError, e:
        print e
        sys.exit(1)

    line = fh1.readline().rstrip('\n')
    line = fh1.readline().rstrip('\n')
    while line:
        list     = re.split("\s+", line)
        indivID = list[1]
        batch   = list[6]
        id2batch[indivID] = batch
        line = fh1.readline().rstrip('\n')
    fh1.close()

    try:
        fh1 = file(preQCIMDS_1kG_sample, "r")
    except IOError, e:
        print e
        sys.exit(1)

    line = fh1.readline().rstrip('\n')
    line = fh1.readline().rstrip('\n')
    while line:
        list     = re.split("\s+", line)
        indivID = list[0]
        batch   = list[1]
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
    fh3.writelines(line + "\tbatch\n")
    line = fh2.readline().rstrip('\n')
    while line:
        list = re.split("\s+", line)
        id = list[1]
        if id in id2batch:
            fh3.writelines(list[0] + "\t" +
                       list[1] + "\t" +
                       list[2] + "\t" +
                       list[3] + "\t" +
                       list[4] + "\t" +
                       list[5] + "\t" +
                       list[6] + "\t" +
                       list[7] + "\t" +
                       list[8] + "\t" +
                       list[9] + "\t" +
                       list[10] + "\t" +
                       list[11] + "\t" +
                       id2batch[id] + "\n")
        line = fh2.readline().rstrip('\n')

    fh2.close()
    fh3.close()

    os.system("cp %s %s" % (eval_file, new_eval_file))


def addphenoinfo_10PCs(evec_file, eval_file, new_evec_file, new_eval_file, individuals_annotation, preQCIMDS_1kG_sample):
    """ add batch information to final evec file """

    id2pheno = {}

    try:
        fh1 = file(individuals_annotation, "r")
    except IOError, e:
        print e
        sys.exit(1)

    line = fh1.readline().rstrip('\n')
    line = fh1.readline().rstrip('\n')
    while line:
        list     = re.split("\s+", line)
        indivID = list[1]
        pheno   = list[8]
        id2pheno[indivID] = pheno
        line = fh1.readline().rstrip('\n')
    fh1.close()

    try:
        fh1 = file(preQCIMDS_1kG_sample, "r")
    except IOError, e:
        print e
        sys.exit(1)

    line = fh1.readline().rstrip('\n')
    line = fh1.readline().rstrip('\n')
    while line:
        list     = re.split("\s+", line)
        indivID = list[0]
        pheno   = list[1]
        id2pheno[indivID] = pheno
        line = fh1.readline().rstrip('\n')
    fh1.close()

    try:
        fh2 = file(evec_file, "r")
        fh3 = file(new_evec_file, "w")
    except IOError, e:
        print e
        sys.exit(1)

    line = fh2.readline().rstrip('\n')
    fh3.writelines(line + "\tbatch\n")
    line = fh2.readline().rstrip('\n')
    while line:
        list = re.split("\s+", line)
        id = list[1]
        if id in id2pheno:
            fh3.writelines(list[0] + "\t" +
                       list[1] + "\t" +
                       list[2] + "\t" +
                       list[3] + "\t" +
                       list[4] + "\t" +
                       list[5] + "\t" +
                       list[6] + "\t" +
                       list[7] + "\t" +
                       list[8] + "\t" +
                       list[9] + "\t" +
                       list[10] + "\t" +
                       list[11] + "\t" +
                       id2pheno[id] + "\n")

        line = fh2.readline().rstrip('\n')

    fh2.close()
    fh3.close()

    os.system("cp %s %s" % (eval_file, new_eval_file))


def addcountryinfo_10PCs(evec_file, eval_file, new_evec_file, new_eval_file, ind_annotation_file, preQCIMDS_1kG_sample):
    """ add country information to final evec file """

    id2country = {}

    try:
        fh1 = file(ind_annotation_file, "r")
    except IOError, e:
        print e
        sys.exit(1)

    line = fh1.readline().rstrip('\n')
    line = fh1.readline().rstrip('\n')
    while line:
        list     = re.split("\s+", line)
        indivID = list[1]
        country   = list[9]
        id2country[indivID] = country
        line = fh1.readline().rstrip('\n')
    fh1.close()

    try:
        fh1 = file(preQCIMDS_1kG_sample, "r")
    except IOError, e:
        print e
        sys.exit(1)

    line = fh1.readline().rstrip('\n')
    line = fh1.readline().rstrip('\n')
    while line:
        list     = re.split("\s+", line)
        indivID = list[0]
        batch   = list[1]
        id2country[indivID] = batch
        line = fh1.readline().rstrip('\n')
    fh1.close()

    try:
        fh2 = file(evec_file, "r")
        fh3 = file(new_evec_file, "w")
    except IOError, e:
        print e
        sys.exit(1)

    line = fh2.readline().rstrip('\n')
    fh3.writelines(line + "\tcountry\n")
    line = fh2.readline().rstrip('\n')
    while line:
        list = re.split("\s+", line)
        id = list[1]
        if id in id2country:
            fh3.writelines(list[0] + "\t" +
                       list[1] + "\t" +
                       list[2] + "\t" +
                       list[3] + "\t" +
                       list[4] + "\t" +
                       list[5] + "\t" +
                       list[6] + "\t" +
                       list[7] + "\t" +
                       list[8] + "\t" +
                       list[9] + "\t" +
                       list[10] + "\t" +
                       list[11] + "\t" +
                       id2country[id] + "\n")
        line = fh2.readline().rstrip('\n')

    fh2.close()
    fh3.close()

    os.system("cp %s %s" % (eval_file, new_eval_file))



def pca_convert(plink, eigenstrat_parameter_file, annotation_file):
    """ convert PLINK file data set to eigenstrat format """

    # ----------------------------- #
    # - generate parameter file m - #
    # ----------------------------- #

    packedped = PackedPed(write_file=eigenstrat_parameter_file)

    packedped.set_input_PLINK_binary(
        bed=plink + ".bed",
        bim=plink + ".bim",
        fam=plink + ".fam")

    packedped.write_par_file()
    del packedped

    with open(eigenstrat_parameter_file, "a") as parfile:
        parfile.write("allowdups: YES\n")

    # ------------------------ #
    # - run convertf program - #
    # ------------------------ #

    cmd = Command("convertf -p %s" % (eigenstrat_parameter_file))
    cmd.run()
    del cmd

    os.system("mv %s.ind %s.ind.orig" % (plink, plink))

    # read individualIDs and HAPMAP info from from hapmap2 fam file
    try:
        fh_anno = file(annotation_file, "r")
    except IOError, e:
        print e
        sys.exit(1)

    individuals2batch_id = {}
    # skip header
    line = fh_anno.readline().replace("\n", "")
    line = fh_anno.readline().replace("\n", "")
    while line:

        list = re.split("\s+", line)
        IID = list[1]
        batch_id = list[6]
        individuals2batch_id[IID] = batch_id

        line = fh_anno.readline().replace("\n", "")

    fh_anno.close()

    # re-write ind file with info on HapMap samples and batch_info
    try:
        fh_ind     = file(plink + ".ind.orig", "r")
        fh_ind_new = file(plink + ".ind", "w")
    except IOError, e:
        print e
        sys.exit(1)

    batches = []
    batches_dict = {}

    # no header line
    line = fh_ind.readline().replace("\n", "")
    while line:

        list = re.split("\s+", line)
        if list[0] == "":
            del list[0]

        # change info last column from "Case/Control" to batch_id
        if list[0] in individuals2batch_id:

            batch_id = individuals2batch_id[list[0]]
            if batch_id in batches_dict:
                batches.append(batch_id)
                batches_dict[batch_id] = True
            if list[-1] == "Case":
                line = line.replace("Case", batch_id)
            elif list[-1] == "Control":
                line = line.replace("Control", batch_id)
            # nothing to replace
            else:
                print >> sys.stderr, "\n    warning: could not replace case/control status for sample " + list[0] + " by batch_id in file pca.evec file " + plink + ".pca.evec ...\n\n"
            fh_ind_new.writelines(line + "\n")

        # nothing to replace
        else:
            print >> sys.stderr, "\n    warning: could not find sample " + list[0] + " in annotation file " + annotation_file + " ...\n\n"
            fh_ind_new.writelines(line + "\n")

        line = fh_ind.readline().replace("\n", "")

    fh_ind.close()
    fh_ind_new.close()
    del batches_dict


pca_main_program = "smartpca.perl.DE"
#pca_main_program = "smartpca.perl"

def pca_run(plink, sigmathreshold, projection_on_populations, numof_pc, numof_threads, draw_evec, draw_without_projection):
    """ run eigenstrat program """

    # ------------------------ #
    # - run eigenstrat program - #
    # ------------------------ #

    plink_pca = plink + "_" + str(numof_pc) + "PC"

    teststring = "%s -i %s.eigenstratgeno -a %s.snp -b %s.ind -k %s -o %s.pca -p %s.plot -e %s.eval -l %s.log -m 5 -t %s -s %s -w %s -f %s -g %s.snpweights" \
                 % (pca_main_program,
                    plink,
                    plink,
                    plink,
                    numof_pc,
                    plink_pca,
                    plink_pca,
                    plink_pca,
                    plink_pca,
                    numof_pc,
                    sigmathreshold,
                    projection_on_populations,
                    numof_threads,
                    plink_pca)
    print >> sys.stderr, teststring
    cmd = Command("%s -i %s.eigenstratgeno -a %s.snp -b %s.ind -k %s -o %s.pca -p %s.plot -e %s.eval -l %s.log -m 5 -t %s -s %s -w %s -f 1 -g %s.snpweights"
                  % (pca_main_program,
                     plink,
                     plink,
                     plink,
                     numof_pc,
                     plink_pca,
                     plink_pca,
                     plink_pca,
                     plink_pca,
                     numof_pc,
                     sigmathreshold,
                     projection_on_populations,
                     plink_pca))
    cmd.run()
    del cmd

    # draw first two PCs
    os.system("R --slave --args %s < %s" % (plink_pca, draw_evec))

    # read which batches (HapMap) were used for projection
    projection_batches = {}
    try:
        fh_proj     = file(projection_on_populations, "r")
    except IOError, e:
        print e
        sys.exit(1)

    line = fh_proj.readline().rstrip('\n')
    while line:

        list = re.split("\s+", line)
        if list[0] == "":
            del list[0]

        projection_batches[list[0]] = list[0]
        line = fh_proj.readline().rstrip('\n')

    fh_proj.close()

    # re-write pca.evec file without projection samples (HapMap samples)
    try:
        fh_pcaevec     = file(plink_pca + ".pca.evec", "r")
        fh_pcaevec_new = file(plink_pca + ".withoutProjection.pca.evec", "w")
    except IOError, e:
        print e
        sys.exit(1)

    # skip header line
    line = fh_pcaevec.readline().rstrip('\n')

    line = fh_pcaevec.readline().rstrip('\n')
    while line:

        list = re.split("\s+", line)
        if list[0] == "":
            del list[0]
        if list[-1] == "":
            del list[-1]

        if not list[-1] in projection_batches:
            for i in xrange(len(list)):
                if i == 0:
                    fh_pcaevec_new.writelines(list[i])
                else:
                    fh_pcaevec_new.writelines("\t" + list[i])
            fh_pcaevec_new.writelines("\n")

        line = fh_pcaevec.readline().rstrip('\n')

    fh_pcaevec.close()
    fh_pcaevec_new.close()

    # draw first two PCs without HapMap samples
    os.system("R --slave --args %s %s < %s" % (plink_pca, plink_pca + ".withoutProjection", draw_without_projection))


def write_snps_autosomes_noLDRegions_noATandGC_noIndels(bim, outfile):
    """  write only autosomal snps, remove SNPs from high LD regions (also MHC),
    remove A/T and C/G SNPs, remove Indels """

    print "\n        remove SNPs from high LD regions ..\n\n"
    print "\n        remove A/T and C/G SNPs ...\n\n"
    print "\n        remove insertions/deletions ...\n\n"

    try:
        bim = file(bim, "r")
        out = file(outfile, "w")
    except IOError, e:
        print e
        sys.exit(1)

    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'D': 'D', 'I': 'I', '-': '-', 'N': 'N'}
    indels = {'D': 'D', 'I': 'I'}

    line = bim.readline().replace("\n", "")
    while line:

        list = re.split("\s+", line)
        chr  = int(list[0])
        pos  = int(list[3])
        a1   = list[4].upper()
        a2   = list[5].upper()
        # exclude non-autosomes
        if 0 < chr and chr < 23:
            # exclude xMHC SNPs AND exclude A/T and C/G SNPs AND exclude D/I SNPs
            if (not ((1 == chr and (48000000 <= pos and pos < 52000000)) or
                     (2 == chr and (86000000 <= pos and pos < 100500000)) or
                     (2 == chr and (134500000 <= pos and pos < 138000000)) or
                     (2 == chr and (183000000 <= pos and pos < 183000000)) or
                     (3 == chr and (47500000 <= pos and pos < 50000000)) or
                     (3 == chr and (83500000 <= pos and pos < 87000000)) or
                     (3 == chr and (89000000 <= pos and pos < 97500000)) or
                     (5 == chr and (44500000 <= pos and pos < 50500000)) or
                     (5 == chr and (98000000 <= pos and pos < 100500000)) or
                     (5 == chr and (129000000 <= pos and pos < 132000000)) or
                     (5 == chr and (135500000 <= pos and pos < 138500000)) or
                     (6 == chr and (25500000 <= pos and pos < 33500000)) or
                     (6 == chr and (57000000 <= pos and pos < 64000000)) or
                     (6 == chr and (140000000 <= pos and pos < 142500000)) or
                     (7 == chr and (55000000 <= pos and pos < 66000000)) or
                     (8 == chr and (8000000 <= pos and pos < 12000000)) or
                     (8 == chr and (43000000 <= pos and pos < 50000000)) or
                     (8 == chr and (112000000 <= pos and pos < 115000000)) or
                     (10 == chr and (37000000 <= pos and pos < 43000000)) or
                     (11 == chr and (46000000 <= pos and pos < 57000000)) or
                     (11 == chr and (87500000 <= pos and pos < 90500000)) or
                     (12 == chr and (33000000 <= pos and pos < 40000000)) or
                     (12 == chr and (109500000 <= pos and pos < 112000000)) or
                     (20 == chr and (32000000 <= pos and pos < 34500000))))\
                    and (len(a1) == 1) and (len(a2) == 1) \
                    and (a1 != complement[a2]) \
                    and (not (a1 in indels or a2 in indels)):

                # write variants for inclusion
                out.writelines("%s\n" % (list[1]))

        line = bim.readline().replace("\n", "")

    bim.close()
    out.close()


def write_snps_autosomes_noLDRegions_noIndels(bim, outfile):
    """  write only autosomal snps, remove SNPs from high LD regions (also MHC),
    remove Indels """

    print "\n        remove SNPs from high LD regions ..\n\n"
    print "\n        remove insertions/deletions ...\n\n"

    try:
        bim = file(bim, "r")
        out = file(outfile, "w")
    except IOError, e:
        print e
        sys.exit(1)

    complement = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'D': 'D', 'I': 'I'}
    indels = {'D': 'D', 'I': 'I'}

    line = bim.readline().replace("\n", "")
    while line:

        list = re.split("\s+", line)
        chr  = int(list[0])
        pos  = int(list[3])
        a1   = list[4].upper()
        a2   = list[5].upper()
        # exclude non-autosomes
        if 0 < chr and chr < 23:
            # exclude xMHC SNPs AND exclude A/T and C/G SNPs AND exclude D/I SNPs
            if (not ((1 == chr and (48000000 <= pos and pos < 52000000)) or
                     (2 == chr and (86000000 <= pos and pos < 100500000)) or
                     (2 == chr and (134500000 <= pos and pos < 138000000)) or
                     (2 == chr and (183000000 <= pos and pos < 183000000)) or
                     (3 == chr and (47500000 <= pos and pos < 50000000)) or
                     (3 == chr and (83500000 <= pos and pos < 87000000)) or
                     (3 == chr and (89000000 <= pos and pos < 97500000)) or
                     (5 == chr and (44500000 <= pos and pos < 50500000)) or
                     (5 == chr and (98000000 <= pos and pos < 100500000)) or
                     (5 == chr and (129000000 <= pos and pos < 132000000)) or
                     (5 == chr and (135500000 <= pos and pos < 138500000)) or
                     (6 == chr and (25500000 <= pos and pos < 33500000)) or
                     (6 == chr and (57000000 <= pos and pos < 64000000)) or
                     (6 == chr and (140000000 <= pos and pos < 142500000)) or
                     (7 == chr and (55000000 <= pos and pos < 66000000)) or
                     (8 == chr and (8000000 <= pos and pos < 12000000)) or
                     (8 == chr and (43000000 <= pos and pos < 50000000)) or
                     (8 == chr and (112000000 <= pos and pos < 115000000)) or
                     (10 == chr and (37000000 <= pos and pos < 43000000)) or
                     (11 == chr and (46000000 <= pos and pos < 57000000)) or
                     (11 == chr and (87500000 <= pos and pos < 90500000)) or
                     (12 == chr and (33000000 <= pos and pos < 40000000)) or
                     (12 == chr and (109500000 <= pos and pos < 112000000)) or
                     (20 == chr and (32000000 <= pos and pos < 34500000))))\
                    and (len(a1) == 1) and (len(a2) == 1) \
                    and (not (a1 in indels or a2 in indels)):

                # write variants for inclusion
                out.writelines("%s\n" % (list[1]))

        line = bim.readline().replace("\n", "")

    bim.close()
    out.close()
