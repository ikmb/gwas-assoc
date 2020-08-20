import sys
import os
import re

pca_main_program = "smartpca.perl.DE"


# may also need some of these:

# import Ingos lib
#sys.path.append(join(sys.path[0], "../../all_scripts"))
#ys.path.append(os.environ['PYLIB_DIR'] + "/all_scripts")
from all_common import Command

# import my lib
#sys.path.append(join(sys.path[0], "../lib"))
#ys.path.append(os.environ['PYLIB_DIR'] + "/lib")

#rom plink_classes import *
#rom eigenstrat_classes import *


def pca_run(plink, sigmathreshold, projection_on_populations, numof_pc, numof_threads=1):
    """ run eigenstrat program """

    # ------------------------ #
    # - run eigenstrat program - #
    # ------------------------ #

    plink_pca = plink + "_" + numof_pc + "PC"

    teststring = "%s -i %s.eigenstratgeno -a %s.snp -b %s.ind -k %s -o %s.pca -p %s.plot -e %s.eval -l %s.log -m 5 -t %s -s %s -w %s -f %s -g %s.snpweights" \
                    %(pca_main_program,
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
    cmd = Command( "%s -i %s.eigenstratgeno -a %s.snp -b %s.ind -k %s -o %s.pca -p %s.plot -e %s.eval -l %s.log -m 5 -t %s -s %s -w %s -f %s -g %s.snpweights" \
                    %(pca_main_program,
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
                      plink_pca))
    cmd.run(); del cmd

    # draw first two PCs
    os.system("R --slave --args %s < %s" %(plink_pca, "draw_evec_EIGENSTRAT.r"))

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
    os.system("R --slave --args %s %s < %s" %(plink_pca, plink_pca + ".withoutProjection", "draw_evec_withoutProjection.r"))


# Main
if __name__ == "__main__":

    # check args
    if len(sys.argv) < 5:
        print "Usage: " + sys.argv[0] + " <input plink basename> <sigma threshold> <projection-on-populations> <number of PCs> [<number of threads=1>]\n"
        sys.exit(1)

    if len(sys.argv) == 5:
        numof_threads = 1
    else:
        numof_threads = sys.argv[5]

    pca_run(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], numof_threads)

