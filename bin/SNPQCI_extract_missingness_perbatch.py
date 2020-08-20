#!/usr/bin/env python

import sys
import re

def generate_exclude_file4(Missingness_file, threshold_par, outfile, batches_list):
    numof_batches = len(batches_list)

    threshold = float(threshold_par)

    try:
        fh_r  = file(Missingness_file, "r")
        fh_w  = file(outfile, "w")
    except IOError, e:
        print e
        sys.exit(1)

    print "Missingness: " + Missingness_file
    print "Threshold:   ", threshold
    print "Outfile:     " + outfile
    print "Batches:    ", batches_list

    # skip header
    line = fh_r.readline().rstrip('\n')
    line = fh_r.readline().rstrip('\n')
    while line:

        list = re.split("\s+", line)
        if list[0] == "":
            del list[0]
        if list[-1] == "":
            del list[-1]

        variant_excluded = False
        variant_id = list[1]
        for i in xrange(numof_batches):

            variant_id_tmp = list[1]
            assert(variant_id == variant_id_tmp)

            missingness = float(list[6])
            #print "missingness: ", missingness, "threshold: ", threshold
            if missingness > threshold and (not variant_excluded):
                fh_w.writelines(list[1] + "\n")
                variant_excluded = True
            line = fh_r.readline().rstrip('\n')

    fh_r.close()
    fh_w.close()


if __name__ == "__main__":
    if len(sys.argv) != 5:
        print "Usage: " + sys.argv[0] + " <missingness-file> <threshold> <annotations> <outfile>\n"
        sys.exit(1)

    # collect batches
    batches_dict = {}
    batches_list = []

    try:
        individuals_fh = file(sys.argv[3], "r")
    except IOError, e:
        print e
        sys.exit(1)

    line = individuals_fh.readline().rstrip('\n')
    list = re.split("\s+", line)
    # assert annotation file validity
    assert(list[6] == "batch")

    line = individuals_fh.readline().rstrip('\n')
    while line:
        list = re.split("\s+", line)
        if list[0] == "":
            del list[0]
        if list[-1] == "":
            del list[-1]

        if not list[6] in batches_dict:
            batches_dict[list[6]] = True
            batches_list.append(list[6])

        line = individuals_fh.readline().rstrip('\n')
    individuals_fh.close()

    generate_exclude_file4(sys.argv[1], sys.argv[2], sys.argv[4], batches_list)
