import os
import os.path
import sys
import re
import gzip

class Gff_snps:
    """ class Gff_snps implements tasks on GFF files with snps """
  
    def __init__(self, gff_file, write_file=None, write_file_status_append=False):
        """ init """
        self._gff_file = gff_file       # name of gff file
        self._write_file = write_file    # name of output file
        self._write_file_status_append = write_file_status_append  # status of write file: new or append
        self._snps_order = []      # order of SNPs from gff file
        # only novel snps
        self._novel = {}          # gff snp : chr, genetic_distance, bp_position, allele1, allele2
        # novel and/or known snps
        self._novel_known = {}  # gff rs : chr, genetic_distance, bp_position, allele1, allele2 
        self._snps_order_perchr = []      # order of SNPs from gff file per chr, list of chr-lists
        # for each chr a list, index 0 -> chr1, index 1 -> chr2, etc.
        for i in xrange(26):
            self._snps_order_perchr.append([])

    def write_file_set(self, write_file):
        """ set name of output file """
        self._write_file = write_file 

    def write_file_set_append(self):
        """ switch to append mode """
        self._write_file_status_append = True
 
    def map_refAllele(self, novel=False):
        """ map GFF_snps file into memory """
        try:
            fh = file(self._gff_file,"r")
        except IOError, e:
            print e
            sys.exit(1)

        line = fh.readline()
        while line:

            list  = re.split("\s+",line)
          
            chr   = list[0].replace("chr","")
            if chr == "X":
                chr = "23"
            elif chr == "Y":
                chr = "24"
            elif chr == "XY":
                chr = "25"
            elif chr == "M":
                chr = "26"
            assert(list[2] == list[3])
            pos   = list[3]  
            feature_list = list[6:]
            if feature_list[-1] == "":
                feature_list.pop()
        
            # possible features
            feature_ref              = re.compile("^ref=")
            refAllele                = None

            for feature in feature_list:
                if feature_ref.search(feature):
                    refAllele        = feature.split("ref=")[1].split(";")[0]

            key = chr +"->"+ pos
            self._snps_order.append(key)
            # novel only 
            if novel: 
              self._novel[key]       = refAllele 
            # novel and/or known
            else:
              self._novel_known[key] = refAllele
            
            line = fh.readline()

        fh.close()

    def map_source_list(self, novel=False):
        """ map GFF_snps file into memory """
        try:
            fh = file(self._gff_file,"r")
        except IOError, e:
            print e
            sys.exit(1)

        line = fh.readline()
        while line:

            list  = re.split("\s+",line)
          
            chr   = list[0].replace("chr","")
            if chr == "X":
                chr = "23"
            elif chr == "Y":
                chr = "24"
            elif chr == "XY":
                chr = "25"
            elif chr == "M":
                chr = "26"
            assert(list[2] == list[3])
            pos   = list[3]  
            feature_list = list[6:]
            if feature_list[-1] == "":
                feature_list.pop()
        
            # possible features
            feature_sources          = re.compile("^sources=")
            sources_list             = []

            for feature in feature_list:
                if feature_sources.search(feature):
                    sources_list     = feature.split("sources=")[1].split(";")[0].split(",")

            key = chr +"->"+ pos
            self._snps_order.append(key)

            # novel only 
            if novel: 
                # if additional sources
                if self._novel.has_key(key):
                    for source in sources_list:
                        self._novel[key].append(source)
                # if new sources
                else:
                    self._novel[key]       = sources_list

            # novel and/or known
            else:
                # if additional sources
                if self._novel_known.has_key(key):
                    for source in sources_list:
                        self._novel_known[key].append(source)
                # if new sources
                else:
                    self._novel_known[key]       = sources_list
            
            line = fh.readline()

        fh.close()
    
    def free_map(self):
        """ free mapped Gff_snp file from memory """
        self._snps_order       = []
        self._novel            = {}
        self._novel_known      = {}

    def novel_known_get_feature_hash(self):
        """ get the feature for novel and/or known snps from gff file in hash """
        return self._novel_known
