import os
import os.path
import sys
import re
import gzip
from copy import deepcopy

class Mlma:
  """ class Mlma implements tasks on mlma output files """
  
  def __init__(self, mlma_file, write_file=None, write_file_status_append=False):
    """ init """
    self._mlma_file = mlma_file    # name of mlma file
    self._write_file = write_file    # name of output file
    self._write_file_status_append = write_file_status_append  # status of write file: new or append
    self._snps_order = []      # order of SNPs from mlma file
    # only typed snps not imputed
    self._typed = {}
    # typed and/or imputed snps
    self._typed_imputed = {}
    
    self._snps_order_perchr = []      # order of SNPs from mlma file per chr, list of chr-lists
    # for each chr a list, index 0 -> chr0, index 1 -> chr1, ..., index 26 -> chr26
    for i in xrange(27):
      self._snps_order_perchr.append([])

  def map_mlma(self, typed=False):
    """ map p-values from mlma file into memory """
    try:
      fh = file(self._mlma_file,"r")
    except IOError, e:
      print e
      sys.exit(1)
    
    typed_pattern = re.compile("^.*_typed.*$")
    imputed_pattern = re.compile("^.*_imputed.*$")

    line = fh.readline()
    while line:
      
      list = re.split("\s+",line)

      # delete empty elements
      if list[0] == "":
          del list[0]
      if list[-1] == "":
          del list[-1]

      rs   = str(list[1])  
      self._snps_order.append(rs.replace("_typed","").replace("_imputed",""))
      
      # typed only 
      if typed: 
        self._typed[rs.replace("_typed","").replace("_imputed","")] = (list[8], "typed")
      # typed and/or imputed 
      else:
        if typed_pattern.search(rs):
            self._typed_imputed[rs.replace("_typed","")] = (list[8], "typed")
        else:
            self._typed_imputed[rs.replace("_imputed","")] = (list[8], "imputed")

      line = fh.readline()
    
    fh.close()

  def free_map(self):
    """ free mapped mlma file from memory """
    self._snps_order       = []
    self._typed            = {}
    self._typed_imputed    = {}

  def rs_get_typed_hash(self):
    """ get the typed rs numbers from file in hash """
    return deepcopy(self._typed)
