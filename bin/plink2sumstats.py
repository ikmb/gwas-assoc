#!/usr/bin/python

#standard modules
import sys
#from string import *
import re
import math
import decimal
import scipy
from scipy.stats.distributions import chi2
from scipy.stats import norm

usage = "usage: " + sys.argv[0] +\
        " PLINK-log assoc_file"

numargs = len(sys.argv)
if numargs != 3:
    print(usage)
    sys.exit(1)



# read number of cases and controls
try:
    fh1 = open(sys.argv[1], "r")
except IOError as e:
    print(e)
    sys.exit(1)

special_line_pattern = re.compile("^Among remaining phenotypes, .* are cases and .* are controls.*$")
line = fh1.readline()
while line:
    if(special_line_pattern.search(line)):
        list = re.split("\s+",line)
        num_cases    = int(list[3])
        num_ctrls = int(list[7])
        break
    line = fh1.readline()
fh1.close()

num_samples = str(int(num_cases) + int(num_ctrls))




# read variants
try:
    fh1 = open(sys.argv[2], "r")
except IOError as e:
    print(e)
    sys.exit(1)

line = fh1.readline().replace("\n", "")
# print new header
print("CHR" +"\t"+ "BP" +"\t"+ "SNP" +"\t"+ "A1" +"\t"+ "A2"+"\t"+ "P" +"\t"+ "OR" +"\t"+ "BETA" +"\t"+ "SE" +"\t"+ "N" +"\t"+ "CHISQ" +"\t"+ "Z" +"\t"+ "SOURCE" +"\t"+ "FRQ_A_A1" +"\t"+ "FRQ_U_A1")# +"\t"+ "RefPanel_A1" +"\t"+ "RefPanelAF_A1" +"\t"+ "INFO" +"\t"+ "Genotyped" +"\t"+ "RSID")
line = fh1.readline().replace("\n", "")
while line:

    list = re.split("\s+",line)
    chromosome = int(list[0])
    SNPID = list[1]
    pos = int(list[2])
    INFO = list[7]

    if not (list[10] == "NA"):
        A1   = list[3]
        A2   = list[4]
        FRQ_A_A1 = list[5] # IMPORTANT: here this is FREQ from affected samples
        FRQ_U_A1 = list[6] # IMPORTANT: here this is FREQ from unaffected samples
        pval = list[10]
        if pval == "0":
            pval = "1e-324"
        chisq = chi2.isf(float(pval), 1) # https://www.biostars.org/p/261698/
        OR = float(list[8])
        source = "PLINK"
        BETA = round(math.log(OR),4)
        SE = list[9]
        zscore = float(BETA)/float(SE)

        print(str(chromosome) +"\t"+ str(pos) +"\t"+ SNPID +"\t"+ A1 +"\t"+ A2 +"\t"+ pval +"\t"+ str(OR) +"\t"+ str(BETA) +"\t"+ str(SE) +"\t"+ num_samples + "\t" + str(chisq) + "\t" + str(zscore) + "\t" + source + "\t" + str(FRQ_A_A1) + "\t" + str(FRQ_U_A1))
    line = fh1.readline().replace("\n", "")

fh1.close()

sys.exit(0)
