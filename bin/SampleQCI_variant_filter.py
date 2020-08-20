#!/usr/bin/env python

import sys
import re
import os

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
    
    complement = {'A':'T', 'T':'A', 'C':'G', 'G':'C', 'D':'D', 'I':'I'}
    indels = {'D':'D', 'I':'I'}
    
    line = bim.readline().replace("\n","")
    while line:

        list = re.split("\s+",line)
        chr  = int(list[0])
        pos  = int(list[3])
        a1   = list[4]
        a2   = list[5]
        # exclude non-autosomes
        if 0 < chr and chr < 23:
            # exclude xMHC SNPs AND exclude A/T and C/G SNPs AND exclude D/I SNPs
            if (not ( (1 == chr and (48000000 <= pos and pos < 52000000)) or \
                     (2 == chr and (86000000 <= pos and pos < 100500000)) or \
                     (2 == chr and (134500000 <= pos and pos < 138000000)) or \
                     (2 == chr and (183000000 <= pos and pos < 183000000)) or \
                     (3 == chr and (47500000 <= pos and pos < 50000000)) or \
                     (3 == chr and (83500000 <= pos and pos < 87000000)) or \
                     (3 == chr and (89000000 <= pos and pos < 97500000)) or \
                     (5 == chr and (44500000 <= pos and pos < 50500000)) or \
                     (5 == chr and (98000000 <= pos and pos < 100500000)) or \
                     (5 == chr and (129000000 <= pos and pos < 132000000)) or \
                     (5 == chr and (135500000 <= pos and pos < 138500000)) or \
                     (6 == chr and (25500000 <= pos and pos < 33500000)) or \
                     (6 == chr and (57000000 <= pos and pos < 64000000)) or \
                     (6 == chr and (140000000 <= pos and pos < 142500000)) or \
                     (7 == chr and (55000000 <= pos and pos < 66000000)) or \
                     (8 == chr and (8000000 <= pos and pos < 12000000)) or \
                     (8 == chr and (43000000 <= pos and pos < 50000000)) or \
                     (8 == chr and (112000000 <= pos and pos < 115000000)) or \
                     (10 == chr and (37000000 <= pos and pos < 43000000)) or \
                     (11 == chr and (46000000 <= pos and pos < 57000000)) or \
                     (11 == chr and (87500000 <= pos and pos < 90500000)) or \
                     (12 == chr and (33000000 <= pos and pos < 40000000)) or \
                     (12 == chr and (109500000 <= pos and pos < 112000000)) or \
                     (20 == chr and (32000000 <= pos and pos < 34500000)) )) \
                and (a1 != complement[a2]) \
                and (not (indels.has_key(a1) or indels.has_key(a2))): 
                
                # write variants for inclusion
                out.writelines("%s\n" %(list[1]))

        line = bim.readline().replace("\n","")

    bim.close()
    out.close()


# Main
if __name__ == "__main__":

    # check args
    if len(sys.argv) != 3:
        print "Usage: " + sys.argv[0] + " <bim-file> <target-file>\n"
        print "\twhere:\n"
        print "\t<bim-file> BIM input\n"
        print "\t<target-file> list of variants that should be filtered\n"
        sys.exit(1)
        
    write_snps_autosomes_noLDRegions_noATandGC_noIndels(sys.argv[1], sys.argv[2])
