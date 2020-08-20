import sys
import re

# may also need some of these:

# import Ingos lib
# sys.path.append(join(sys.path[0], "../../all_scripts"))
# sys.path.append(os.environ['PYLIB_DIR'] + "/all_scripts")
# from all_common import Command

# import my lib
# sys.path.append(join(sys.path[0], "../lib"))
# sys.path.append(os.environ['PYLIB_DIR'] + "/lib")

# from plink_classes import *
# from eigenstrat_classes import *

def remove_projection_samples(plink_pca, projection_on_populations):
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


# Main
if __name__ == "__main__":

    # check args
    if len(sys.argv) < 3:
        print "Usage: " + sys.argv[0] + " <input plink basename> <projection-on-populations>\n"
        sys.exit(1)

    remove_projection_samples(sys.argv[1], sys.argv[2])
