#!/usr/bin/env python

import sys
import re
import os

from os.path import *
import string
import re
import gzip
import math
import decimal
import datetime
from os import listdir
import subprocess

# may also need some of these:

# import Ingos lib
#sys.path.append(join(sys.path[0], "../../all_scripts"))
sys.path.append(os.environ['PYLIB_DIR'] + "/all_scripts")
from all_common import *

# import my lib
#sys.path.append(join(sys.path[0], "../lib"))
sys.path.append(os.environ['PYLIB_DIR'] + "/lib")

from plink_classes import *
from eigenstrat_classes import *

def pca_convert(plink, eigenstrat_parameter_file, annotation_file):
    """ convert PLINK file data set to eigenstrat format """
    
    # ----------------------------- #
    # - generate parameter file m - #
    # ----------------------------- #
    
    packedped = PackedPed(write_file=eigenstrat_parameter_file)
            
    packedped.set_input_PLINK_binary(
        bed=plink + ".bed",\
        bim=plink + ".bim",\
        fam=plink + ".fam")
    
    packedped.write_par_file() ; del packedped
            
    # ------------------------ #
    # - run convertf program - #
    # ------------------------ #

    cmd = Command( "convertf -p %s" \
        %(eigenstrat_parameter_file) )
    cmd.run() ; del cmd

    os.system("mv %s.ind %s.ind.orig" \
        %(plink, plink) )

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
    
        list = re.split("\s+",line)
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
    
        list = re.split("\s+",line)
        if list[0] == "":
            del list[0]

        # change info last column from "Case/Control" to batch_id
        if individuals2batch_id.has_key(list[0]):
            
            batch_id = individuals2batch_id[list[0]]
            if not batches_dict.has_key(batch_id):
                batches.append(batch_id)
                batches_dict[batch_id] = True
            if list[-1] == "Case":
                line = line.replace("Case", batch_id)
            elif list[-1] == "Control":
                line = line.replace("Control", batch_id)
            # nothing to replace
            else:
                print >> sys.stderr, "\n    warning: could not replace case/control status for sample " +list[0]+ " by batch_id in file pca.evec file " +plink_pca + ".pca.evec ...\n\n"
            fh_ind_new.writelines(line +"\n")

        # nothing to replace
        else:
            print >> sys.stderr, "\n    warning: could not found sample " +list[0]+ " in annotation file " +individuals_annotation_cases_controls_hapmap2+ " ...\n\n"
            fh_ind_new.writelines(line +"\n")

        line = fh_ind.readline().replace("\n", "")
   
    fh_ind.close()
    fh_ind_new.close()
    del batches_dict


# Main
if __name__ == "__main__":

    # check args
    if len(sys.argv) != 4:
        print "Usage: " + sys.argv[0] + " <input plink basename> <eigenstrat parameter file> <annotations>\n"
        sys.exit(1)
        
    pca_convert(sys.argv[1], sys.argv[2], sys.argv[3])
	
