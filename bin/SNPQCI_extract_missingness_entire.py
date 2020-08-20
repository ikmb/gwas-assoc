#!/usr/bin/env python

import sys
import re

def generate_exclude_file3(Missingness_file, threshold_par, outfile):
    """ generate exclude file 3: Remove variants with certain missingness across the entire collection """

    print "Missingness: ", Missingness_file
    print "Threshold: ", threshold_par
    print "Outfile: ", outfile

    threshold = float(threshold_par)

    try:
        fh_r  = file(Missingness_file, "r")
        fh_w  = file(outfile, "w")
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
        if list[-1] == "":
            del list[-1]

        missingness = float(list[4])
        if missingness > threshold:
            fh_w.writelines(list[1] + "\n")

        line = fh_r.readline().rstrip('\n')

    fh_r.close()
    fh_w.close()


if __name__ == "__main__":
    if len(sys.argv) != 4:
        print "Usage: " + sys.argv[0] + " <missingness-file> <threshold> <results>\n"
        sys.exit(1)

    generate_exclude_file3(sys.argv[1], sys.argv[2], sys.argv[3])
