##################################################
###  Genotype-phenotype overlap (GPO):         ###
###  2009, David Ellinghaus                    ###
###  d.ellinghaus(at)ikmb(dot)uni-kiel(dot)de  ###
##################################################

import os
import os.path
from os.path import *
import sys
import re
import gzip
import string

#RPy
#import rpy2.robjects as robjects
#import statistics

# import my lib
import plink_classes
import error


class PackedPed:
    """ class PackedPed implements tasks on PackedPed files from EIGENSTRAT software """

    def __init__(self, write_file="",\
                 outputformat="EIGENSTRAT",\
                 familynames="NO"):
        """ init """

        self._write_file      = write_file    # default name: par.PACKEDPED.EIGENSTRAT
        self._genotypename    = ""            # input PLINK bed file 
        self._snpname         = ""            # input PLINK bim file
        self._indivname       = ""            # input PLINK fam file
        self._outputformat    = outputformat  # default out format: EIGENSTRAT  
        self._genotypeoutname = ""            # output name eigenstratgeno file 
        self._snpoutname      = ""            # output name eigenstratsnp file
        self._indivoutname    = ""            # outout name eigenstratind file
        self._familynames     = familynames   # default family names: NO


    def set_input_PLINK_binary(self, bed, bim, fam):
        """ input files are in PLINK binary format """
       
        # determine prefix from PLINK files
        a = bed.split(".") ; a.pop() ; 
        prefix = string.join(a,".") ; del a

        self._genotypename = bed
        self._indivname    = fam
        self._snpname      = bim
        self._genotypeoutname   = prefix + ".eigenstratgeno" 
        self._snpoutname        = prefix + ".snp"
        self._indivoutname      = prefix + ".ind"


    def write_par_file(self):
        """  write par file """

        try:
            out = file(self._write_file, "w")
        except IOError, e:
            print e
            sys.exit(1)

        out.writelines("")        

        out.writelines("genotypename:\t%s\n" %(self._genotypename))
        out.writelines("snpname:\t%s\n" %(self._snpname))
        out.writelines("indivname:\t%s\n" %(self._indivname))
        out.writelines("outputformat:\t%s\n" %(self._outputformat))
        out.writelines("genotypeoutname:\t%s\n" %(self._genotypeoutname))
        out.writelines("snpoutname:\t%s\n" %(self._snpoutname))
        out.writelines("indivoutname:\t%s\n" %(self._indivoutname))
        out.writelines("familynames:\t%s\n" %(self._familynames))

        out.close()


class PcaEvec:
    """ class PcaEvec implements tasks on *.pca.evec files from EIGENSTRAT software """

    def __init__(self, input_file="", write_file="",\
                 tmpdir="", input2_file="", remove_outlier={},\
                 R_dir=""):
        """ init """

        self._input_file      = input_file      # input file *.pca.evac 
        self._input2_file     = input2_file     # second input file *.pca.evac
        self._write_file      = write_file      # write file for output
        self._tmpdir          = tmpdir          # tmpdir for intermediate files 
        self._R_dir           = R_dir           # R_dir for R scripts 
                                               
        self._sample_dim      = {}              # key=sample, val=listofdimensions
        self._samples_list    = []              # val=sample in order of file
        self._samples_dict    = {}              # key=sample, val=True
        self._numofdim        = 0               # number of dimensions
        self._outlier         = remove_outlier  # outlier to remove


    def map(self):
        """ map input_file into main memory """

        try:
            input  = file(self._input_file, "r")
        except IOError, e:
            print e
            sys.exit(1)

        header_pattern = re.compile("^.*#eigvals.*$")
        line = input.readline().replace("\n","")
        if not header_pattern.search(line):
            print >> sys.stderr, "error: wrong header in file \"" + self._input_file + "\"" 
            sys.exit(1)

        line = input.readline().replace("\n","")
        counter = 0
        while line:

            list = re.split("\s+",line)

            # delete first element if empty
            if list[0] == "":
                del list[0]

            sample_id = list[0]
            self._sample_dim[sample_id] = list[1:-1] 
            self._samples_list.append(sample_id)
            self._samples_dict[sample_id] = True
            if counter == 0:
                self._numofdim = len(list[1:-1])
            elif self._numofdim != len(list[1:-1]):
                print >> sys.stderr, "error: different number of dimensions in in file \"" +\
                         self._input_file + "\"" 
                print >> sys.stderr, str(self._numofdim) +" != "+\
                                     str(len(list[1:-1]))
                sys.exit(1)
            
            counter += 1
            line = input.readline().replace("\n","")

        input.close() 


    def compare_evecs(self):
        """ compare eigenvectors of input_file and input2_file """

        try:
            input2 = file(self._input2_file, "r")
            out    = file(self._write_file, "w")
        except IOError, e:
            print e
            sys.exit(1)

        samples_overlap_list  = []  # val=sample
        sample_dim    = {}          # key=sample, val=listofdimensions

        # read matrix 2
        header_pattern = re.compile("^.*#eigvals.*$")
        line = input2.readline().replace("\n","")
        if not header_pattern.search(line):
            print >> sys.stderr, "error: wrong header in file \"" + self._input_file + "\"" 
            sys.exit(1)

        line = input2.readline().replace("\n","")
        while line:

            list = re.split("\s+",line)
            
            # delete first element if empty
            if list[0] == "":
                del list[0]
            
            sample_id = list[0]
            if self._samples_dict.has_key(sample_id):
                sample_dim[sample_id] = list[1:-1] 
                samples_overlap_list.append(sample_id)

            line = input2.readline().replace("\n","")

        input2.close() 
        
        # ------------------------------ #
        # - dimensions of input file 1 - #
        # ------------------------------ #

        # list with self._numofdim lists in it
        dimensions_input1 = []
        dimensions_input1_samples = []
        for i in xrange(self._numofdim):
            dimensions_input1.append([])
        
        for sample_id in samples_overlap_list:

            # if outlier then ignore outlier
            if self._outlier.has_key(sample_id):
                pass
            else:
                # check for same number of dimensions in first input file
                if self._numofdim != len(self._sample_dim[sample_id]):
                    print >> sys.stderr, "error: different number of dimensions in file \"" +\
                             self._input1_file + "\"" 
                    print >> sys.stderr, str(self._numofdim) +" != "+\
                                         str(len(self._sample_dim[sample_id]))
                    sys.exit(1)
                
                # fill list with self._numofdim lists
                for i in xrange(len(self._sample_dim[sample_id])):
                    dimensions_input1[i].append(float(self._sample_dim[sample_id][i]))
                dimensions_input1_samples.append(sample_id)
        
        # ------------------------------ #
        # - dimensions of input file 2 - #
        # ------------------------------ #

        # list with self._numofdim lists in it
        dimensions_input2 = []
        dimensions_input2_samples = []
        for i in xrange(self._numofdim):
            dimensions_input2.append([])
        
        for sample_id in samples_overlap_list:

            # if outlier then ignore outlier
            if self._outlier.has_key(sample_id):
                pass
            else:
                # check for same number of dimensions in first input file
                if self._numofdim != len(sample_dim[sample_id]):
                    print >> sys.stderr, "error: different number of dimensions in file \"" +\
                             self._input2_file + "\"" 
                    print >> sys.stderr, str(self._numofdim) +" != "+\
                                         str(len(self._sample_dim[sample_id]))
                    sys.exit(1)
                
                # fill list with self._numofdim lists
                for i in xrange(len(sample_dim[sample_id])):
                    dimensions_input2[i].append(float(sample_dim[sample_id][i]))
                dimensions_input2_samples.append(sample_id)
             
        # ------------------------------------------------------------------ #
        # - calc correlation pearson for each dimension in file1 and file2 - #
        # ------------------------------------------------------------------ #
       
        assert(dimensions_input1_samples == dimensions_input2_samples)
        dimensions_correlation = []
        assert(len(dimensions_input1) == len(dimensions_input2))

        # write header
        for i in xrange(len(dimensions_input1)):
            if i == 0:
                out.writelines("dim" + str(i+1))
            else:
                out.writelines("\tdim" + str(i+1))
        out.writelines("\n")

        # write body
        for i in xrange(len(dimensions_input1)):
            assert(len(dimensions_input1[i]) == len(dimensions_input2[i]))
            #print dimensions_input1[i]
            #print dimensions_input2[i]
            dimensions_correlation.append(\
                                   statistics.correlation(\
                                       dimensions_input1[i],\
                                       dimensions_input2[i],\
                                       method="Pearson"))
            if i == 0:
                out.writelines(str(dimensions_correlation[i]))
            else:
                out.writelines("\t"+ str(dimensions_correlation[i]))
        out.writelines("\n")
        out.close()


class PcaLog:
    """ class PcaLog implements tasks on *.log files from EIGENSTRAT software """

    def __init__(self, input_file=""):
        """ init """

        self._input_file      = input_file    # input file *.log
        self._outlier         = {}            # key=outlier, val=sigmage


    def get_outlier(self):
        """ map outlier into outlier list """

        comment_pattern = re.compile("^#.*$")
        blankline_pattern = re.compile("^\s*$")
    
        try:
            input  = file(self._input_file, "r")
        except IOError, e:
            print e
            sys.exit(1)

        outlier_pattern = re.compile("^REMOVED.*$")
        
        line = input.readline()
        while line:

            # skip comment lines that start with "#"
            if(comment_pattern.search(line)):
                line = input.readline()
                continue
            # skip blank lines
            if(blankline_pattern.search(line)):
                line = input.readline()
                continue

            if outlier_pattern.search(line):
                
                list = re.split("\s+",line)
                outlier = list[2]
                sigmage = list[8]
                self._outlier[outlier] = sigmage
        
            line = input.readline()

        input.close() 
        return self._outlier.copy()
