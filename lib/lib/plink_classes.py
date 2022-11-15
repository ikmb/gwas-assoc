import os
import os.path
import sys
import re
import gzip
from copy import deepcopy
import random

class Fam:
  """ class Fam implements tasks on FAM files from PLINK """
  
  def __init__(self, fam_file, write_file=None, seed=None, numof_shuffle=None, write_file_status_append=False):
    """ init """
    self._fam_file = fam_file      # name of fam file
    self._write_file = write_file  # name of output file
    self._write_file_status_append = write_file_status_append  # status of write file: new or append
    self._samples_order = []  # order of samples from fam file
    self._cases    = {}       # fam_id->ind_id : pat_id, mat_id, sex, pheno
    self._controls = {}       # fam_id->ind_id : pat_id, mat_id, sex, pheno
    self._unknown  = {}       # missing affection status, fam_id->ind_id : pat_id, mat_id, sex, pheno
    self._seed = seed 
    self._numof_shuffle= int(numof_shuffle)

  def write_file_set(self, write_file):
    """ set name of output file """
    self._write_file = write_file 

  def write_file_set_append(self):
    """ switch to append mode """
    self._write_file_status_append = True
  
  def map(self):
    """ map FAM file into memory """
    try:
      fh = file(self._fam_file,"r")
    except IOError, e:
      print e
      sys.exit(1)

    line = fh.readline().replace("\n","")
    while line:
      list = re.split("\s+",line)
      id = str(list[0]) + "->" + str(list[1]) 
      self._samples_order.append(id)
      # cases
      if str(list[5]) == "1":
        self._controls[id]    = (list[2], list[3], list[4], list[5])
      # controls
      elif str(list[5]) == "2":
        self._cases[id] = (list[2], list[3], list[4], list[5])
      # unkown affection status
      else:
        self._unknown[id]  = (list[2], list[3], list[4], list[5])
      line = fh.readline().replace("\n","")
    fh.close()

  def free_map(self):
    """ free mapped FAM file from memory """
    self._samples_order = []
    self._cases    = []
    self._controls = []
    self._unknown  = []

  def permutate_case_control(self):
    """ map FAM file into memory, permutate (shuffle randomly using specified random_number)
        case/control status, a write permutated fam file """
    try:
      fh   = file(self._fam_file,"r")
      fh_w = file(self._write_file,"w")
    except IOError, e:
      print e
      sys.exit(1)

    phenotypes = []
    line = fh.readline().replace("\n","")
    while line:
      list = re.split("\s+",line)
      phenotypes.append(list[5])
      line = fh.readline().replace("\n","")

    # shuffle x times, x times with new random_number
    random.seed(self._seed)
    for i in xrange(self._numof_shuffle):
        # generate random number, float between 0 and 1, excluding 1
        random_number = random.random()  
        print >> sys.stderr, "perm="+str(i+1), "random_number="+str(random_number), phenotypes[1:20]
        random.shuffle(phenotypes, lambda: random_number)

    fh.seek(0)
    line = fh.readline().replace("\n","")
    counter = 0
    while line:
      list = re.split("\s+",line)
      fh_w.writelines(list[0] +" "+ list[1] +" "+ list[2] +" "+ list[3] +" "+ list[4] +" "+ phenotypes[counter] + "\n")
      counter += 1
      line = fh.readline().replace("\n","")

    assert(len(phenotypes) == counter)

    fh.close()
    fh_w.close()

  def write_file_cases_controls_even_random(self):
    """ write even selection cases and controls from fam file to file """
    if self._write_file_status_append:
      try:
        fh = file(self._write_file,"a")
      except IOError, e:
        print e
        sys.exit(1)
    else:
      try:
	    fh = file(self._write_file,"w")
      except IOError, e:
        print e
        sys.exit(1)

    for id in self._samples_order:
      if self._controls.has_key(id):
        fam_id, ind_id =  self.split_id(id)
        fh.writelines(fam_id + " " + ind_id + " ")
        for i in xrange(0,len(self._controls[id])-1):
          fh.writelines(self._controls[id][i] + " ")
        fh.writelines(self._controls[id][len(self._controls[id])-1] + "\n")

    cases_random_indices = random.sample(xrange(len(self._cases)), len(self._controls))
    case_counter = 0 
    for id in self._samples_order:
      if self._cases.has_key(id):
        if case_counter in cases_random_indices:
          fam_id, ind_id =  self.split_id(id)
          fh.writelines(fam_id + " " + ind_id + " ")
          for i in xrange(0,len(self._cases[id])-1):
            fh.writelines(self._cases[id][i] + " ")
          fh.writelines(self._cases[id][len(self._cases[id])-1] + "\n")
        case_counter += 1
    
    fh.close()

  def split_id(self, id):
    """ split fam_id->ind_id into fam_id, ind_id """
    return id.split("->")

  def write_file_orig(self):
    """ write original fam file to file """
    if self._write_file_status_append:
      try:
        fh = file(self._write_file,"a")
      except IOError, e:
        print e
        sys.exit(1)
    else:
      try:
        fh = file(self._write_file,"w")
      except IOError, e:
        print e
        sys.exit(1)

    for id in self._samples_order:
      alias_dict = None
      # find out which affection status dictionary
      if self._cases.has_key(id):
        alias_dict = self._cases
      elif self._controls.has_key(id):
        alias_dict = self._controls
      else: 
        alias_dict = self._unknown

      fam_id, ind_id =  self.split_id(id)
      fh.writelines(fam_id + " " + ind_id + " ")
      for i in xrange(0,len(alias_dict[id])-1):
        fh.writelines(alias_dict[id][i] + " ")
      fh.writelines(alias_dict[id][len(alias_dict[id])-1] + "\n")
    
    fh.close()

  def write_file_cases(self):
    """ write cases from fam file to file """
    if self._write_file_status_append:
      try:
        fh = file(self._write_file,"a")
      except IOError, e:
        print e
        sys.exit(1)
    else:
      try:
	fh = file(self._write_file,"w")
      except IOError, e:
        print e
        sys.exit(1)

    for id in self._samples_order:
      if self._cases.has_key(id):
        fam_id, ind_id =  self.split_id(id)
        fh.writelines(fam_id + " " + ind_id + " ")
        for i in xrange(0,len(self._cases[id])-1):
          fh.writelines(self._cases[id][i] + " ")
        fh.writelines(self._cases[id][len(self._cases[id])-1] + "\n")
    
    fh.close()

  def write_file_cases_sorted(self, sort_file):
    """ write sorted cases from fam file to file """
    if self._write_file_status_append:
      try:
        fh = file(self._write_file,"a")
        fh2 = gzip.open(sort_file, "rb")
      except IOError, e:
        print e
        sys.exit(1)
    else:
      try:
        fh = file(self._write_file,"w")
        fh2 = gzip.open(sort_file, "rb")
      except IOError, e:
        print e
        sys.exit(1)

    ids_from_sort_file = []
    line = fh2.readline()
    while line:
      ids_from_sort_file.append(re.split("\s+",line)[0])
      line = fh2.readline()
     
    for id in ids_from_sort_file:
      if self._cases.has_key(id):
        fam_id, ind_id =  self.split_id(id)
        fh.writelines(fam_id + " " + ind_id + " ")
        for i in xrange(0,len(self._cases[id])-1):
          fh.writelines(self._cases[id][i] + " ")
        fh.writelines(self._cases[id][len(self._cases[id])-1] + "\n")
    
    fh.close()
    fh2.close()

  def write_file_controls(self):
    """ write controls from fam file to file """
    if self._write_file_status_append:
      try:
        fh = file(self._write_file,"a")
      except IOError, e:
        print e
        sys.exit(1)
    else:
      try:
	    fh = file(self._write_file,"w")
      except IOError, e:
        print e
        sys.exit(1)

    for id in self._samples_order:
      if self._controls.has_key(id):
        fam_id, ind_id =  self.split_id(id)
        fh.writelines(fam_id + " " + ind_id + " ")
        for i in xrange(0,len(self._controls[id])-1):
          fh.writelines(self._controls[id][i] + " ")
        fh.writelines(self._controls[id][len(self._controls[id])-1] + "\n")
    
    fh.close()

  def write_file_controls_sorted(self, sort_file):
    """ write sorted controls from fam file to file """
    if self._write_file_status_append:
      try:
        fh = file(self._write_file,"a")
        fh2 = gzip.open(sort_file, "rb")
      except IOError, e:
        print e
        sys.exit(1)
    else:
      try:
        fh = file(self._write_file,"w")
        fh2 = gzip.open(sort_file, "rb")
      except IOError, e:
        print e
        sys.exit(1)

    ids_from_sort_file = []
    line = fh2.readline()
    while line:
      ids_from_sort_file.append(re.split("\s+",line)[0])
      line = fh2.readline()
 
    for id in ids_from_sort_file:
      if self._controls.has_key(id):
        fam_id, ind_id =  self.split_id(id)
        fh.writelines(fam_id + " " + ind_id + " ")
        for i in xrange(0,len(self._controls[id])-1):
          fh.writelines(self._controls[id][i] + " ")
        fh.writelines(self._controls[id][len(self._controls[id])-1] + "\n")
    
    fh.close()
    fh2.close()


  def write_file_cases_controls_sorted(self, sort_file):
    """ write sorted cases and controls from fam file to file """
    if self._write_file_status_append:
      try:
        fh = file(self._write_file,"a")
        fh2 = gzip.open(sort_file, "rb")
      except IOError, e:
        print e
        sys.exit(1)
    else:
      try:
        fh = file(self._write_file,"w")
        fh2 = gzip.open(sort_file, "rb")
      except IOError, e:
        print e
        sys.exit(1)

    ids_from_sort_file = []
    line = fh2.readline()
    while line:
      ids_from_sort_file.append(re.split("\s+",line)[0])
      line = fh2.readline()
 
    for id in ids_from_sort_file:
      if self._controls.has_key(id):
        fam_id, ind_id =  self.split_id(id)
        fh.writelines(fam_id + " " + ind_id + " ")
        for i in xrange(0,len(self._controls[id])-1):
          fh.writelines(self._controls[id][i] + " ")
        fh.writelines(self._controls[id][len(self._controls[id])-1] + "\n")
      elif self._cases.has_key(id):
        fam_id, ind_id =  self.split_id(id)
        fh.writelines(fam_id + " " + ind_id + " ")
        for i in xrange(0,len(self._cases[id])-1):
          fh.writelines(self._cases[id][i] + " ")
        fh.writelines(self._cases[id][len(self._cases[id])-1] + "\n")
  
    fh.close()
    fh2.close()

  def write_file_controls_sorted_list(self, sort_list):
    """ write sorted controls from fam file to file """
    if self._write_file_status_append:
      try:
        fh = file(self._write_file,"a")
      except IOError, e:
        print e
        sys.exit(1)
    else:
      try:
        fh = file(self._write_file,"w")
      except IOError, e:
        print e
        sys.exit(1)

    for id in sort_list:
      if self._controls.has_key(id):
        fam_id, ind_id =  self.split_id(id)
        fh.writelines(fam_id + " " + ind_id + " ")
        for i in xrange(0,len(self._controls[id])-1):
          fh.writelines(self._controls[id][i] + " ")
        fh.writelines(self._controls[id][len(self._controls[id])-1] + "\n")
    
    fh.close()

  def print_std(self):
    """ print fam file to stdout """
    for id in self._samples_order:
      fam_id, ind_id =  self.split_id(id)
      print fam_id, ind_id,
      for i in xrange(0,len(self._samples[id])-1):
	print self._samples[id][i],
      print self._samples[id][len(self._samples[id])-1] 
  
  def ids_cases_get(self):
    """ get the ids from cases in original order """
    ids = []
    for id in self._samples_order:
      if self._cases.has_key(id):
        ids.append(id) 
    return ids

  def ids_cases_hash_get(self):
    """ get the ids from cases in hash """
    ids = {}
    for id in self._samples_order:
      if self._cases.has_key(id):
        ids[id] = True 
    return ids.copy()

  def ids_cases_iid_hash_get(self):
    """ get the individual ids from cases in hash """
    ids = {}
    for id in self._samples_order:
      if self._cases.has_key(id):
        ids[id.split("->")[1]] = True 
    return ids.copy()

  def ids_controls_get(self):
    """ get the ids from controls in original order """
    ids = []
    for id in self._samples_order:
      if self._controls.has_key(id):
        ids.append(id) 
    return ids
  
  def ids_controls_hash_get(self):
    """ get the ids from controls in hash """
    ids = {}
    for id in self._samples_order:
      if self._controls.has_key(id):
        ids[id] = True
    return ids.copy()

  def ids_controls_iid_hash_get(self):
    """ get the individuals ids from controls in hash """
    ids = {}
    for id in self._samples_order:
      if self._controls.has_key(id):
        ids[id.split("->")[1]] = True
    return ids.copy()

  def ids_samples_get(self):
   """ get the ids from alle samples (cases, controls, unknown)
       in original order """
   return self._samples_order

  def indiv_ids_samples_get(self):
    """ get the individual ids from alle samples (cases, controls, unknown)
        as dictionary """
    iids = {}
    for id in self._samples_order:
      iids[id.split("->")[1]] = True
    return iids.copy()

  def fam_file_name_get(self):
    """ get the name of the fam_file """
    return self._fam_file

  def get_ranges_cases(self):
    """ get ranges cases """
    try:
      fh = file(self._fam_file,"r")
    except IOError, e:
      print e
      sys.exit(1)

    line = fh.readline().replace("\n","")
    ranges_cases = []
    current_case = False
    # start counting by 1
    indiv_counter = 1
    while line:
      list = re.split("\s+",line)
      # case
      if str(list[5]) == "2" and not current_case:
        ranges_cases.append(indiv_counter)
        current_case = True
      elif str(list[5]) == "2" and current_case:
        pass
      elif str(list[5]) != "2" and not current_case:
        pass
      else:
        ranges_cases.append(indiv_counter)
        current_case = False

      indiv_counter += 1
      line = fh.readline().replace("\n","")

    if len(ranges_cases)%2 == 1:
        ranges_cases.append(indiv_counter)

    fh.close()
    return ranges_cases[:]

  def get_ranges_controls(self):
    """ get ranges controls """
    try:
      fh = file(self._fam_file,"r")
    except IOError, e:
      print e
      sys.exit(1)

    line = fh.readline().replace("\n","")
    ranges_controls = []
    current_control = False
    # start counting by 1
    indiv_counter = 1
    while line:
      list = re.split("\s+",line)
      # control
      if str(list[5]) == "1" and not current_control:
        ranges_controls.append(indiv_counter)
        current_control = True
      elif str(list[5]) == "1" and current_control:
        pass
      elif str(list[5]) != "1" and not current_control:
        pass
      else:
        ranges_controls.append(indiv_counter)
        current_control = False

      indiv_counter += 1
      line = fh.readline().replace("\n","")

    if len(ranges_controls)%2 == 1:
        ranges_controls.append(indiv_counter)

    fh.close()
    return ranges_controls[:]




class Bim:
  """ class Bim implements tasks on BIM files from PLINK """
  
  def __init__(self, bim_file, write_file=None, write_file_status_append=False):
    """ init """
    self._bim_file = bim_file       # name of bim file
    self._write_file = write_file    # name of output file
    self._write_file_status_append = write_file_status_append  # status of write file: new or append
    self._snps_order = []      # order of SNPs from bim file
    # only typed snps not imputed
    self._typed = {}          # bim rs : chr, genetic_distance, bp_position, allele1, allele2
    # typed and/or imputed snps
    self._typed_imputed = {}  # bim rs : chr, genetic_distance, bp_position, allele1, allele2 
    
    self._snps_order_perchr = []      # order of SNPs from bim file per chr, list of chr-lists
    # for each chr a list, index 0 -> chr0, index 1 -> chr1, ..., index 26 -> chr26
    for i in xrange(27):
      self._snps_order_perchr.append([])

  def write_file_set(self, write_file):
    """ set name of output file """
    self._write_file = write_file 

  def write_file_set_append(self):
    """ switch to append mode """
    self._write_file_status_append = True
 
  def map(self, typed=False):
    """ map BIM file into memory """
    try:
      fh = file(self._bim_file,"r")
    except IOError, e:
      print e
      sys.exit(1)

    line = fh.readline()
    while line:
      list = re.split("\s+",line)
      rs   = str(list[1])  
      self._snps_order.append(rs)
      # typed only 
      if typed: 
        self._typed[rs]         = (list[0], list[2], list[3], list[4], list[5])
      # typed and/or imputed 
      else:
        self._typed_imputed[rs] = (list[0], list[2], list[3], list[4], list[5])
      line = fh.readline()
    fh.close()

  def free_map(self):
    """ free mapped BIM file from memory """
    self._snps_order       = []
    self._typed            = {}
    self._typed_imputed    = {}

  def map_perchr(self, typed=False):
    """ map BIM file into memory per chr """
    try:
      fh = file(self._bim_file,"r")
    except IOError, e:
      print e
      sys.exit(1)

    line = fh.readline()
    list = re.split("\s+",line)
    chr_counter = int(list[0])
    while line:

      list = re.split("\s+",line)
      if chr_counter != int(list[0]):
        chr_counter = int(list[0])

      rs   = str(list[1])  
      self._snps_order_perchr[chr_counter].append(rs)
      # typed only 
      if typed: 
        self._typed[rs]         = (list[0], list[2], list[3], list[4], list[5])
      # typed and/or imputed 
      else:
        self._typed_imputed[rs] = (list[0], list[2], list[3], list[4], list[5])
      line = fh.readline()

    fh.close()

  def free_map_perchr(self):
    """ free mapped BIM file from memory """
    self._snps_order       = []
    self._typed            = []
    self._typed_imputed    = []
    for i in xrange(27):
      self._snps_order_perchr[i] = []

  def rs_get_typed(self):
    """ get the typed rs numbers from bim file in original order """
    rs_list = []
    for rs in self._snps_order: 
      if self._typed.has_key(rs):
        rs_list.append(rs)    
    return rs_list

  def rs_perchr_get(self, chr_list=range(1,23)):
    """ get the typed rs numbers from bim file in original order """
    snp_list = []
    for chr in chr_list:
      for snp in self._snps_order_perchr[chr]:
        snp_list.append(snp)
    return snp_list

  def rs_get_typed_hash(self):
    """ get the typed rs numbers from bim file in hash """
    rs_hash = {}
    keys = self._typed.iterkeys()  
    for k in keys:
      rs_hash[k] = 1 
    assert(len(rs_hash) == len(self._typed))
    return rs_hash

  def rs_get_typed_imputed(self):
    """ get the typed and/or imputed rs numbers from bim file in original order """
    rs_list = []
    for rs in self._snps_order: 
      if self._typed_imputed.has_key(rs):
        rs_list.append(rs)    
    return rs_list

  def rs_get_typed_imputed_hash(self):
    """ get the typed and/or imputed rs numbers from bim file in hash """
    rs_hash = {}
    keys = self._typed_imputed.iterkeys()  
    for k in keys:
      rs_hash[k] = 1 
    assert(len(rs_hash) == len(self._typed_imputed))
    return rs_hash

  def snps_get_typed(self):
    """ get the typed snps with information from bim file in original order """
    snp_list = [] 
    for rs in self._snps_order: 
      if self._typed.has_key(rs):
        snp_list.append( (self._typed[rs][0],rs,self._typed[rs][1], \
	  self._typed[rs][2],self._typed[rs][3],self._typed[rs][4]) )
    return snp_list

  def snps_get_typed_hash(self):
    """ get the typed snps with all information from bim file in hash """
    return self._typed

  def snps_get_typed_imputed_hash(self):
    """ get the typed and imputed snps with all information from bim file in hash """
    return self._typed_imputed

  def get_rs_typed_imputed_nearest_neighbour(self, chr, pos, rs_pvalue_type_hash):
    """ get nearest neighbour for pos, must be also in rs_pvalue_type_hash"""
    chr = int(chr)
    pos = int(pos)
    pos_lower   = 0
    rs_lower    = ""
    pos_upper   = int(self._typed_imputed[self._snps_order_perchr[chr][-1]][2])
    rs_upper    = "" 

    for i in xrange(len(self._snps_order_perchr[chr])):
      rs = self._snps_order_perchr[chr][i]
      if rs_pvalue_type_hash.has_key(rs):
        if pos >= int(self._typed_imputed[rs][2]):
          pos_lower   = int(self._typed_imputed[rs][2])
          rs_lower    = rs
        else:
          pos_upper   = int(self._typed_imputed[rs][2])
          rs_upper    = rs
          break

    # boundary condition left
    if rs_lower == "":
      return rs_upper
    # boundary condition right
    elif rs_upper == "":
      return rs_lower
    # nearest neighbour
    elif (pos - pos_lower) <= (pos_upper - pos):
      return rs_lower
    else:
      return rs_upper


class Map:
  """ class Map implements tasks on MAP files from PLINK """
  
  def __init__(self, map_file, write_file=None, write_file_status_append=False):
    """ init """
    self._map_file = map_file       # name of map file
    self._write_file = write_file    # name of output file
    self._write_file_status_append = write_file_status_append  # status of write file: new or append
    self._snps_order = []      # order of SNPs from map file
    # only typed snps not imputed
    self._typed = {}          # map rs : chr, genetic_distance, bp_position, allele1, allele2
    # typed and/or imputed snps
    self._typed_imputed = {}  # map rs : chr, genetic_distance, bp_position, allele1, allele2 

  def write_file_set(self, write_file):
    """ set name of output file """
    self._write_file = write_file 

  def write_file_set_append(self):
    """ switch to append mode """
    self._write_file_status_append = True
 
  def map(self, typed=False):
    """ map MAP file into memory """
    try:
      fh = file(self._map_file,"r")
    except IOError, e:
      print e
      sys.exit(1)

    line = fh.readline()
    while line:
      list = re.split("\s+",line)
      rs   = str(list[1])  
      self._snps_order.append(rs)
      # typed only 
      if typed: 
        self._typed[rs]         = (list[0], list[2], list[3])
      # typed and/or imputed 
      else:
        self._typed_imputed[rs] = (list[0], list[2], list[3])
      line = fh.readline()
    fh.close()

  def free_map(self):
    """ free mapped MAP file from memory """
    self._snps_order       = []
    self._typed            = []
    self._typed_imputed    = []

  def write_map_increasing_ids_novelsnps(self, typed=False):
    """ write map file with increasing id starting with 1 
        for novel "snp" snps """
    try:
      fh = file(self._write_file,"w")
    except IOError, e:
      print e
      sys.exit(1)

    list_pointer = None
    if typed:
      list_pointer = self._typed
    else:
      list_pointer = self._typed_imputed
    
    snp_pattern = re.compile("^snp.*$")
    snp_counter = 1

    snp_list = [] 
    for snp in self._snps_order: 
      if list_pointer.has_key(snp):
        if snp_pattern.search(snp) and int(snp.split("snp")[1].replace("x", "")) >= 1:
          snp_list.append( (list_pointer[snp][0],\
                            "snp" + str(snp_counter),\
                            list_pointer[snp][1],\
                            list_pointer[snp][2]) )
          snp_counter += 1
        else:
          snp_list.append( (list_pointer[snp][0],\
                            snp,\
                            list_pointer[snp][1],\
                            list_pointer[snp][2]) )


    assert(len(self._snps_order) == len(snp_list))
    for tuple in snp_list:
      fh.writelines(tuple[0] +"\t"+ tuple[1] +"\t"+ tuple[2] +"\t"+ tuple[3] +"\n")
    fh.close()

  def rs_get_typed(self):
    """ get the typed rs numbers from map file in original order """
    rs_list = []
    for rs in self._snps_order: 
      if self._typed.has_key(rs):
        rs_list.append(rs)    
    return rs_list

  def rs_get_typed_hash(self):
    """ get the typed rs numbers from map file in hash """
    rs_hash = {}
    keys = self._typed.iterkeys()  
    for k in keys:
      rs_hash[k] = 1 
    assert(len(rs_hash) == len(self._typed))
    return rs_hash

  def rs_get_typed_hash_rsfunction(self):
    """ get the typed rs numbers from map file in hash with their function """
    rs_hash = {}
    keys = self._typed.iterkeys()  
    for k in keys:
      k_tmp = k.replace("_typed","")
      k_tmp = k_tmp.replace("_imputed","")
      rs_hash[k_tmp] = k
    assert(len(rs_hash) == len(self._typed))
    return rs_hash

  def rs_get_typed_position_hash(self):
    """ get the typed rs numbers from map file in hash with position values """
    rs_hash = {}
    keys = self._typed.iterkeys()  
    for k in keys:
      # key=rs, value=chr->position
      rs_hash[k] = self._typed[k][0] +"->"+ self._typed[k][2]
    assert(len(rs_hash) == len(self._typed))
    return rs_hash

  def position_get_typed_rs_hash(self):
    """ get the typed position numbers from map file in hash with rs values """
    pos_hash = {}
    keys = self._typed.iterkeys()  
    for k in keys:
      # key=chr->position, value=rs
      pos_hash[self._typed[k][0] +"->"+ self._typed[k][2]] = k
    assert(len(pos_hash) == len(self._typed))
    return pos_hash

  def rs_get_typed_imputed(self):
    """ get the typed and/or imputed rs numbers from map file in original order """
    rs_list = []
    for rs in self._snps_order: 
      if self._typed_imputed.has_key(rs):
        rs_list.append(rs)    
    return rs_list

  def rs_get_typed_imputed_hash(self):
    """ get the typed and/or imputed rs numbers from map file in hash """
    rs_hash = {}
    keys = self._typed_imputed.iterkeys()  
    for k in keys:
      rs_hash[k] = 1 
    assert(len(rs_hash) == len(self._typed_imputed))
    return rs_hash

  def snps_get_typed(self):
    """ get the typed snps with information from map file in original order """
    snp_list = [] 
    for rs in self._snps_order: 
      if self._typed.has_key(rs):
        snp_list.append( (self._typed[rs][0],rs,self._typed[rs][1], \
                          self._typed[rs][2]) )
    return snp_list

  def snps_original_get_typed(self):
    """ get the typed snps with all information (original) from map file in original order """
    snp_list = [] 
    for rs in self._snps_order: 
      if self._typed.has_key(rs):
        snp_list.append( (self._typed[rs][0], rs, self._typed[rs][1], \
                          self._typed[rs][2]) )
    return snp_list

  def snps_get_typed_imputed_hash(self):
    """ get the typed and imputed snps with all information from map file in hash """
    return self._typed_imputed



class Ped:
  """ class Ped implements tasks on PED files from PLINK """
  
  def __init__(self, ped_file, write_file=None, write_file_status_append=False):
    """ init """
    self._ped_file = ped_file      # name of ped file
    self._write_file = write_file  # name of output file
    self._write_file_status_append = write_file_status_append  # status of write file: new or append
    self._samples_order = []  # order of samples from ped file
    self._cases    = {}       # ped_id->ind_id : pat_id, mat_id, sex, pheno, alleles ...
    self._controls = {}       # ped_id->ind_id : pat_id, mat_id, sex, pheno, alleles ...
    self._unknown  = {}       # missing affection status, ped_id->ind_id : pat_id, mat_id, sex, pheno, alleles ...

  def write_file_set(self, write_file):
    """ set name of output file """
    self._write_file = write_file 

  def write_file_set_append(self):
    """ switch to append mode """
    self._write_file_status_append = True
 
  def map(self):
    """ map PED file into memory """
    try:
      fh = file(self._ped_file,"r")
    except IOError, e:
      print e
      sys.exit(1)

    line = fh.readline().replace("\n","")
    while line:
      list = re.split("\s+",line)
      id = str(list[0]) + "->" + str(list[1]) 
      self._samples_order.append(id)
      # cases
      if str(list[5]) == "1":
	self._controls[id]    = (list[2:])
      # controls
      elif str(list[5]) == "2":
	self._cases[id] = (list[2:])
      # unkown affection status
      else:
	self._unknown[id]  = (list[2:])
      line = fh.readline().replace("\n","")
    fh.close()

  def free_map(self):
    """ free mapped PED file from memory """
    self._samples_order = []
    self._cases    = []
    self._controls = []
    self._unknown  = []

  def write_new_mldose_cases(self, ids_list, new_snp_list_first_allele=[]):
    """ write new mldose file with allele dosages for samples from ids_list """
    if self._write_file_status_append:
      try:
        out = gzip.open(self._write_file, "a")
      except IOError, e:
        print e
        sys.exit(1)
    else:
      try:
        out = gzip.open(self._write_file, "w")
      except IOError, e:
        print e
        sys.exit(1)
  
    for id in ids_list:
      if self._cases.has_key(id):
	out.writelines("%s ML_DOSE" %(id))
        
	snp_counter = 0
	for i in xrange(4,len(self._cases[id]), 2):
          dosage = 0.0
	  dosage_allele = new_snp_list_first_allele[snp_counter]
	  if dosage_allele == self._cases[id][i]:
	    dosage = dosage + 1.0
	  if dosage_allele == self._cases[id][i+1]:
	    dosage = dosage + 1.0
	  snp_counter = snp_counter + 1
	  out.writelines(" %s" %(str(dosage)))
	out.writelines("\n")

    out.close()

  def write_new_mldose_controls(self, ids_list, new_snp_list_first_allele):
    """ write new mldose file with allele dosages for samples from ids_list """
    if self._write_file_status_append:
      try:
        out = gzip.open(self._write_file, "a")
      except IOError, e:
        print e
        sys.exit(1)
    else:
      try:
        out = gzip.open(self._write_file, "w")
      except IOError, e:
        print e
        sys.exit(1)
  
    for id in ids_list:
      if self._controls.has_key(id):
	out.writelines("%s ML_DOSE" %(id))
        
	snp_counter = 0
	for i in xrange(4,len(self._controls[id]), 2):
          dosage = 0.0
	  dosage_allele = new_snp_list_first_allele[snp_counter]
	  if dosage_allele == self._controls[id][i]:
	    dosage = dosage + 1.0
	  if dosage_allele == self._controls[id][i+1]:
	    dosage = dosage + 1.0
	  snp_counter = snp_counter + 1
	  out.writelines(" %s" %(str(dosage)))
	out.writelines("\n")

    out.close()

  def write_new_mldose_cases_controls(self, ids_list, new_snp_list_first_allele=[]):
    """ write new mldose file with allele dosages for samples from ids_list """
    if self._write_file_status_append:
      try:
        out = gzip.open(self._write_file, "a")
      except IOError, e:
        print e
        sys.exit(1)
    else:
      try:
        out = gzip.open(self._write_file, "w")
      except IOError, e:
        print e
        sys.exit(1)
  
    for id in ids_list:

      cases_controls_unknown = None 
      
      if self._cases.has_key(id):
        cases_controls_unknown = self._cases
      elif self._controls.has_key(id):
        cases_controls_unknown = self._controls
      elif self._unknown.has_key(id):
        cases_controls_unknown = self._unknown

      assert(cases_controls_unknown != None)
        
      out.writelines("%s ML_DOSE" %(id))
      snp_counter = 0
      for i in xrange(4,len(cases_controls_unknown[id]), 2):
        dosage = 0.0
        dosage_allele = new_snp_list_first_allele[snp_counter]
        if dosage_allele == cases_controls_unknown[id][i]:
          dosage = dosage + 1.0
        if dosage_allele == cases_controls_unknown[id][i+1]:
          dosage = dosage + 1.0
        snp_counter = snp_counter + 1
        out.writelines(" %s" %(str(dosage)))
      out.writelines("\n")

    out.close()


class Clump:
  """ class Clump implements tasks on CLUMP files from PLINK """
  
  def __init__(self, clump_file, write_file=None):
    """ init """
    self._clump_file = clump_file    # name of clump file
    self._write_file = write_file    # name of output file
    self._snps_order = []            # order of SNPs from clump file
    self._lines      = []            # order of lines from clump file
    # only typed snps not imputed
    self._typed = {}          # clump rs : chr, rs_lead, pos, p_value, SP2
    # typed and/or imputed snps
    self._typed_imputed = {}  # clump rs : chr, rs_lead, pos, p_value, SP2

  def map(self, typed=False):
    """ map CLUMP file into memory """
    try:
      fh = file(self._clump_file,"r")
    except IOError, e:
      print e
      sys.exit(1)

    # check header line
    line = fh.readline()
    blankline_pattern = re.compile("^\s*$")
    header_pattern = re.compile("^CHR\s+F.*$")
    if not header_pattern.search(line):
      print >> sys.stderr, "error: wrong header in file \"" + self._clump_file + "\"" 
      sys.exit(1)
   
    # map whole file
    underscore_pattern = re.compile("^.*_.*$")
    line = fh.readline()
    while line:
      if blankline_pattern.search(line):
          line = fh.readline()
          continue
   
      list = re.split("\s+",line)

      # delete empty elements
      if list[0] == "":
          del list[0]
      if list[-1] == "":
          del list[-1]

      if str(list[0]) != "NA":

          chr     = str(list[0])
          rs_lead = ""
          if underscore_pattern.search(list[2]):
            rs_lead = list[2].replace("_typed","")
            rs_lead = rs_lead.replace("_imputed","")
          else:
            rs_lead = list[2].replace("_typed","")
            rs_lead = rs_lead.replace("_imputed","")
          self._snps_order.append(rs_lead)
          # typed only 
          if typed: 
            self._typed[rs_lead]         = (chr, rs_lead, list[3], list[4], list[11])
          # typed and/or imputed 
          else:
            self._typed_imputed[rs_lead] = (chr, rs_lead, list[3], list[4], list[11])

      line = fh.readline()

    fh.close()

  def free_map(self):
    """ free mapped CLUMP file from memory """
    self._snps_order       = []
    self._typed            = []
    self._typed_imputed    = []

  def rs_get_typed_imputed_tophits(self, numoftophits=100):
    """ get the typed and/or imputed rs numbers from clump file in original order """
    rs_list = []
    hit_counter = 0
    for rs in self._snps_order: 
      if self._typed_imputed.has_key(rs):
        rs_list.append(rs)    
        hit_counter = hit_counter + 1
      if hit_counter >= int(numoftophits):
        break
    return rs_list

  def rs_sp2_get_typed_imputed_tophits_hash(self, numoftophits=100):
    """ get the typed and/or imputed supported snps (SP2) from rs numbers 
        from clump file in hash """
    rs_sp2_hash = {}
    hit_counter = 0
    underscore_pattern = re.compile("^.*_.*$")
    for rs in self._snps_order: 

      if self._typed_imputed.has_key(rs):

        if self._typed_imputed[rs][4] == "NONE":
            
          rs_sp2_hash[rs] = ["NA"]

        else:
          
          # get supported snps for rs in list_sp2
          list_sp2_tmp = self._typed_imputed[rs][4].split(",")
          list_sp2 = []
          for snp in list_sp2_tmp:
            if underscore_pattern.search(snp):
              snp = snp.replace("_typed(1)","").replace("_imputed(1)","")
              list_sp2.append(snp)
            else:
              list_sp2.append(snp)
        
          # clone list
          rs_sp2_hash[rs] = list_sp2[:]

        hit_counter = hit_counter + 1
      if hit_counter >= int(numoftophits):
        break
    return rs_sp2_hash.copy()

  def rs_pos_get_typed_imputed_tophits_hash(self, numoftophits=100):
    """ get the typed and/or imputed positions from rs numbers 
        from clump file in hash """
    rs_pos_hash = {}
    hit_counter = 0
    for rs in self._snps_order: 

      if self._typed_imputed.has_key(rs):

        # get chr->position for rs
        rs_pos_hash[rs] = self._typed_imputed[rs][0] + "->" + self._typed_imputed[rs][2]

        hit_counter = hit_counter + 1
      if hit_counter >= int(numoftophits):
        break
    return rs_pos_hash.copy()

  def sp2_rs_get_typed_imputed_tophits_hash(self, numoftophits=100):
    """ get the typed and/or imputed supported snps (SP2) from rs numbers 
        from clump file in hash """
    sp2_rs_hash = {}
    hit_counter = 0
    for rs in self._snps_order: 

      if self._typed_imputed.has_key(rs):

        # get supported snps for rs in list_sp2
        list_sp2_tmp = self._typed_imputed[rs][4].split(",")
        list_sp2 = []
        for snp in list_sp2_tmp:
          snp = snp.replace("_typed(1)","").replace("_imputed(1)","")
          list_sp2.append(snp)

        # clone list
        for sp in list_sp2:
            sp2_rs_hash[sp] = rs 

        hit_counter = hit_counter + 1
      if hit_counter >= int(numoftophits):
        break
    return sp2_rs_hash.copy()

  def rs_get_all_typed_imputed_tophits_hash(self, numoftophits=100):
    """ get all info from rs numbers 
        from clump file in hash """
    rs_hash = {}
    hit_counter = 0
    for rs in self._snps_order: 

      if self._typed_imputed.has_key(rs):

        # clump rs : chr, rs_lead, pos, p_value, SP2
        rs_hash[rs] = self._typed_imputed[rs]

        hit_counter = hit_counter + 1
      if hit_counter >= int(numoftophits):
        break
    return rs_hash.copy()

  def get_clump_linenumber(self, leadsnp):
      """ get line number/rank of lead snp in clump file """
      try:
        return self._snps_order.index(leadsnp) + 1
      except ValueError, e:
        return "NA" 

  def write_only_check_snps_from_clump_group(self, check_snps, gene_file):
    """ write only check snps from clumped file ot write_file """
    try:
      fh  = file(self._clump_file,"r")
      fh2 = file(self._write_file,"w")
      fh3 = file(gene_file,"r")
    except IOError, e:
      print e
      sys.exit(1)
   
    # read genes from gene_file
    genes = {}
    line = fh3.readline()
    while line:
      list = re.split("\s+",line)
      chr = list[0] 
      if not genes.has_key(chr):
        genes[chr] = []
      genes[chr].append((list[1], list[2], list[3], list[4], list[5]))
      line = fh3.readline()
    fh3.close()  

    # check header line
    line = fh.readline()
    header_pattern = re.compile("^CHR F.*$")
    if not header_pattern.search(line):
      print >> sys.stderr, "error: wrong header in file \"" + self._clump_file + "\"" 
      sys.exit(1)
    fh2.writelines("rank_in_clumped_groups genes num_of_hits hits " + line)
 
    blankline_pattern = re.compile("^\s*$")
    line = fh.readline()
    rank_counter = 0
    while line:
      
      if blankline_pattern.search(line):
        print >> sys.stderr, "error: blank line in file \"" + self._clump_file + "\"" 
        sys.exit(1)
      
      rank_counter += 1
      list = re.split("\s+",line)
 
      chr = list[0]
      pos = list[3]

      # check for genes at chr,pos
      gene_names = []
      if genes.has_key(chr):
        for i in xrange(len(genes[chr])):
          gene_start = genes[chr][i][0]
          gene_stop  = genes[chr][i][1]
          if gene_start <= pos and pos <= gene_stop:
            gene_names.append(genes[chr][i][4])
      if len(gene_names) == 0:
        gene_names.append("NA")

      # generate list of rs numbers 
      rs_list     = []
      rs_lead             = list[2].replace("_imputed","").replace("_typed","")
      rs_support_list_tmp = list[11].split(",")
      rs_list.append(rs_lead)
      for elem in rs_support_list_tmp:
          # it should be no None values in there
          if elem != "None":
              rs_list.append(elem.replace("_imputed(1)","").replace("_typed(1)",""))

      hit_counter = 0
      list = []
      for rs in rs_list:
          if check_snps.has_key(rs):
              if hit_counter == 0:
                list.append(rs)
                hit_counter += 1
              else: 
                list.append(",")
                list.append(rs)
                hit_counter += 1
      
      # print only if hit from check snps
      if hit_counter >= 1:
        fh2.writelines(str(rank_counter))
        for i in xrange(len(gene_names)):
          if i == 0:
            fh2.writelines(" " + gene_names[i]) 
          else:
            fh2.writelines("," + gene_names[i]) 
        fh2.writelines(" " + str(hit_counter) + " ") 
        for elem in list:
          fh2.writelines(str(elem))
        fh2.writelines(" " + line)

      line = fh.readline()

    fh.close()
    fh2.close()  

  def write_file_set(self, write_file):
    """ set name of output file """
    self._write_file = write_file 

  def write_snps_typed(self):
    """ map Clump file into memory """
    try:
      fh  = file(self._clump_file,"r")
      fh2 = file(self._write_file,"w")
    except IOError, e:
      print e
      sys.exit(1)
   
    # check header line
    line = fh.readline()
    header_pattern = re.compile("^\s*CHR.*$")
    if not header_pattern.search(line):
      print >> sys.stderr, "error: wrong header in file \"" + self._clump_file + "\"" 
      sys.exit(1)
    list = re.split("\s+",line)
    for i in xrange(1, len(list)):
      if i == 1:
        fh2.writelines(list[i])
      else:
        fh2.writelines(" " + list[i])
    fh2.writelines("\n")
    
    blankline_pattern = re.compile("^\s*$")
    # parse only clumped SNPs with supported SNPs
    line = fh.readline()
    while line:
      if not blankline_pattern.search(line):
        
        list = re.split("\s+",line)

        # categorize snps: typed or imputed
        rs_lead             = list[3]
        rs_support_list_tmp = list[12].split(",")
        rs_support_list     = []
 
        # lead snp
        rs_lead = rs_lead + "_typed"

        # if supported snps not exist
        if rs_support_list_tmp[0] == "NONE" or rs_support_list_tmp[0] == "NA":
          rs_support_list.append(rs_support_list_tmp[0])
        # if supported snps exist
        else:
          for snp in rs_support_list_tmp:
            prefix, suffix = snp.split("(")
            snp_new = ""
            snp_new = prefix + "_typed(" + suffix
            rs_support_list.append(snp_new)

        # write altered line
        fh2.writelines(list[1] + " " + list[2] + " " + rs_lead + " " +\
                       list[4] + " " + list[5] + " " + list[6] + " " +\
                       list[7] + " " + list[8] + " " + list[9] + " " +\
                       list[10]+ " " + list[11]+ " " )
        for i in xrange(len(rs_support_list)):
          if i == 0:
            fh2.writelines(rs_support_list[i])
          else:
            fh2.writelines("," + rs_support_list[i])
        fh2.writelines("\n")

      line = fh.readline()

    fh.close()
    fh2.close()  

  def write_snps_typed_noChr23(self):
    """ map Clump file into memory """
    try:
      fh  = file(self._clump_file,"r")
      fh2 = file(self._write_file,"w")
    except IOError, e:
      print e
      sys.exit(1)
   
    # check header line
    line = fh.readline()
    header_pattern = re.compile("^\s*CHR.*$")
    if not header_pattern.search(line):
      print >> sys.stderr, "error: wrong header in file \"" + self._clump_file + "\"" 
      sys.exit(1)
    list = re.split("\s+",line)
    for i in xrange(1, len(list)):
      if i == 1:
        fh2.writelines(list[i])
      else:
        fh2.writelines(" " + list[i])
    fh2.writelines("\n")
    
    blankline_pattern = re.compile("^\s*$")
    # parse only clumped SNPs with supported SNPs
    line = fh.readline()
    while line:
      if not blankline_pattern.search(line):
        
        list = re.split("\s+",line)

        chr                 = list[1]
        if chr == "23" or chr == "24" or chr == "25" or chr == "26":
          line = fh.readline()
          continue
        
        # categorize snps: typed or imputed
        rs_lead             = list[3]
        rs_support_list_tmp = list[12].split(",")
        rs_support_list     = []
 
        # lead snp
        rs_lead = rs_lead + "_typed"

        # if supported snps not exist
        if rs_support_list_tmp[0] == "NONE" or rs_support_list_tmp[0] == "NA":
          rs_support_list.append(rs_support_list_tmp[0])
        # if supported snps exist
        else:
          for snp in rs_support_list_tmp:
            prefix, suffix = snp.split("(")
            snp_new = ""
            snp_new = prefix + "_typed(" + suffix
            rs_support_list.append(snp_new)

        # write altered line
        fh2.writelines(list[1] + " " + list[2] + " " + rs_lead + " " +\
                       list[4] + " " + list[5] + " " + list[6] + " " +\
                       list[7] + " " + list[8] + " " + list[9] + " " +\
                       list[10]+ " " + list[11]+ " " )
        for i in xrange(len(rs_support_list)):
          if i == 0:
            fh2.writelines(rs_support_list[i])
          else:
            fh2.writelines("," + rs_support_list[i])
        fh2.writelines("\n")

      line = fh.readline()

    fh.close()
    fh2.close()  
    
  def write_genes_with_clump_group(self, vicinity_kb, gene_file):
    """ write only check snps from clumped file ot write_file """
    try:
      fh  = file(self._clump_file,"r")
      fh2 = file(self._write_file,"w")
      fh3 = file(gene_file,"r")
    except IOError, e:
      print e
      sys.exit(1)
  
    vicinity = 1000 * int(vicinity_kb) # in bp

    # read genes from gene_file
    genes = {}
    line = fh3.readline()
    while line:
      list = re.split("\s+",line)
      chr = list[0] 
      if not genes.has_key(chr):
        genes[chr] = []
      genes[chr].append((list[1], list[2], list[3], list[4], list[5]))
      line = fh3.readline()
    fh3.close()  

    # check header line
    line = fh.readline()
    header_pattern = re.compile("^CHR F.*$")
    if not header_pattern.search(line):
      print >> sys.stderr, "error: wrong header in file \"" + self._clump_file + "\"" 
      sys.exit(1)
    fh2.writelines("genes within_genes pos_from_gene_start pos_from_gene_end " + line)
 
    blankline_pattern = re.compile("^\s*$")
    line = fh.readline()
    while line:
      
      if blankline_pattern.search(line):
        print >> sys.stderr, "error: blank line in file \"" + self._clump_file + "\"" 
        sys.exit(1)
      
      list = re.split("\s+",line)
 
      chr = list[0]
      pos = int(list[3])

      # check for genes at chr,pos
      gene_names          = []
      within_genes        = []
      pos_from_gene_start = []
      pos_from_gene_end   = []
      if genes.has_key(chr):
        for i in xrange(len(genes[chr])):
          gene_range_start = int(genes[chr][i][0])
          gene_range_end   = int(genes[chr][i][1])
          gene_start = gene_range_start + vicinity
          gene_end   = gene_range_end   - vicinity
          if gene_range_start <= pos and pos <= gene_range_end:
            gene_names.append(genes[chr][i][4])
            if gene_start <= pos and pos <= gene_end:
              within_genes.append("y")
            else:
              within_genes.append("n")
            pos_from_gene_start.append(pos - gene_start)
            pos_from_gene_end.append(pos - gene_end)
      if len(gene_names) == 0:
        gene_names.append("NA")
        within_genes.append("NA")
        pos_from_gene_start.append("NA")
        pos_from_gene_end.append("NA")

      # add column genes, pos_from_gene_start  
      for i in xrange(len(gene_names)):
        if i == 0:
          fh2.writelines(gene_names[i]) 
        else:
          fh2.writelines("," + gene_names[i]) 
      fh2.writelines(" ") 
      for i in xrange(len(within_genes)):
        if i == 0:
          fh2.writelines(" " + within_genes[i]) 
        else:
          fh2.writelines("," + within_genes[i]) 
      fh2.writelines(" ") 
      for i in xrange(len(pos_from_gene_start)):
        if i == 0:
          fh2.writelines(" " + str(pos_from_gene_start[i])) 
        else:
          fh2.writelines("," + str(pos_from_gene_start[i])) 
      fh2.writelines(" ") 
      for i in xrange(len(pos_from_gene_end)):
        if i == 0:
          fh2.writelines(" " + str(pos_from_gene_end[i]))
        else:
          fh2.writelines("," + str(pos_from_gene_end[i]))
      fh2.writelines(" " + line)

      line = fh.readline()

    fh.close()
    fh2.close()  

  def write_snps_typed_imputed(self, typed_snps):
    """ map Clump file into memory """
    try:
      fh  = file(self._clump_file,"r")
      fh2 = file(self._write_file,"w")
    except IOError, e:
      print e
      sys.exit(1)
   
    # check header line
    line = fh.readline()
    header_pattern = re.compile("^\s*CHR.*$")
    if not header_pattern.search(line):
      print >> sys.stderr, "error: wrong header in file \"" + self._clump_file + "\"" 
      sys.exit(1)
    
    list = re.split("\s+",line)
    # delete empty elements
    if list[0] == "":
      del list[0]
    if list[-1] == "":
      del list[-1]

    for i in xrange(len(list)):
      if i == 0:
        fh2.writelines(list[i])
      else:
        fh2.writelines(" " + list[i])
    fh2.writelines("\n")
    
    blankline_pattern = re.compile("^\s*$")
    # parse only clumped SNPs with supported SNPs
    line = fh.readline()
    while line:
      if not blankline_pattern.search(line):
        
        list = re.split("\s+",line)

        # delete empty elements
        if list[0] == "":
          del list[0]
        if list[-1] == "":
          del list[-1]

        # categorize snps: typed or imputed
        rs_lead             = list[2]
        rs_support_list_tmp = list[11].split(",")
        rs_support_list     = []
 
        # lead snp
        if typed_snps.has_key(rs_lead):
          rs_lead = rs_lead + "_typed"
        else:
          rs_lead = rs_lead + "_imputed"

        # if supported snps not exist
        if rs_support_list_tmp[0] == "NONE" or rs_support_list_tmp[0] == "NA":
          rs_support_list.append(rs_support_list_tmp[0])
        # if supported snps exist
        else:
          for snp in rs_support_list_tmp:
            prefix, suffix = snp.split("(")
            snp_new = ""
            if typed_snps.has_key(prefix):
              snp_new = prefix + "_typed(" + suffix
            else:
              snp_new = prefix + "_imputed(" + suffix
            rs_support_list.append(snp_new)

        # write altered line
        fh2.writelines(list[0] + " " + list[1] + " " + rs_lead + " " +\
                       list[3] + " " + list[4] + " " + list[5] + " " +\
                       list[6] + " " + list[7] + " " + list[8] + " " +\
                       list[9]+ " " + list[10]+ " " )
        for i in xrange(len(rs_support_list)):
          if i == 0:
            fh2.writelines(rs_support_list[i])
          else:
            fh2.writelines("," + rs_support_list[i])
        fh2.writelines("\n")

      line = fh.readline()

    fh.close()
    fh2.close()  
    
  def write_supported_snps_from_altered_file(self):
    """ write supported snps from altered snp filemap Clump file into memory """
    try:
      fh  = file(self._clump_file,"r")
      fh2 = file(self._write_file,"w")
    except IOError, e:
      print e
      sys.exit(1)
   
    # check header line
    line = fh.readline().replace("\n","")
    header_pattern = re.compile("^CHR.*$")
    if not header_pattern.search(line):
      print >> sys.stderr, "error: wrong header in file \"" + self._clump_file + "\"" 
      sys.exit(1)
    fh2.writelines(line + "num_of_SP2\n")
    
    blankline_pattern = re.compile("^\s*$")

    # parse only clumped SNPs with supported SNPs
    line = fh.readline().replace("\n","")
    while line:
      if not blankline_pattern.search(line):
        list = re.split("\s+",line)
        # if snp has support
        if list[11] != "NA" and list[11] != "NONE":
          num_of_support = len(list[11].split(","))
          fh2.writelines(line + " " + str(num_of_support) + "\n")
      line = fh.readline().replace("\n","")

    fh.close()
    fh2.close()


class Ld:
  """ class Ld implements tasks on ld (linkage disequilibirum) files from PLINK """
  
  def __init__(self, ld_file, write_file=None):
    """ init """
    self._ld_file = ld_file          # name of ld file
    self._snp_a = None               # snp 1 for ld pairs
    self._snps_b_order = []            # order of SNPs snp_b from ld file
    self._snps_b = {}                  # snp_b -> r2 

  def map(self):
    """ map ld file into memory """
    try:
      fh = file(self._ld_file,"r")
    except IOError, e:
      print e
      sys.exit(1)

    # check header line
    line = fh.readline()
    header_pattern = re.compile("^ CHR_A.*$")
    if not header_pattern.search(line):
      print >> sys.stderr, "error: wrong header in file \"" + self._ld_file + "\"" 
      sys.exit(1)
   
    # map second line
    line = fh.readline()
    # if file not empty
    if line:
      list = re.split("\s+",line)
      snp_a = str(list[3])  
      self._snp_a = snp_a
      snp_b = str(list[6])  
      r2    = str(list[7])  
      self._snps_b_order.append(snp_b)
      self._snps_b[snp_b] = r2 
      
      # map whole file
      line = fh.readline()
      while line:
        list = re.split("\s+",line)
        snp_a = str(list[3])  
        snp_b = str(list[6])  
        r2    = str(list[7])  
    
        if self._snp_a != snp_a:
            print >> sys.stderr, "error: SNP_A is not always the same in file \"" + self._ld_file + "\"" 
            sys.exit(1)
    
        self._snps_b_order.append(snp_b)
        self._snps_b[snp_b] = r2 
        
        line = fh.readline()
    
    fh.close()

  def free_map(self):
    """ free ld file from memory """
    self._snp_a = None 
    self._snps_b_order = []
    self._snps_b = {}

  def get_snp_r2_hash(self):
    """ get snp -> r2 in hash """ 
    return self._snps_b

  def get_snp_r2_list(self):
    """ get snp list in original order """ 
    return self._snps_b_order



class Assoc:
  """ class Assoc implements tasks on Assoc files from PLINK """
  
  def __init__(self, assoc_file, write_file=None, write_file_status_append=False):
    """ init """
    self._assoc_file = assoc_file    # name of assoc file
    self._write_file = write_file    # name of output file
    self._write_file_status_append = write_file_status_append  # status of write file: new or append
    self._snps_order = []      # order of SNPs from assoc file
    # only typed snps not imputed
    self._typed = {}
    # typed and/or imputed snps
    self._typed_imputed = {}
    
    self._snps_order_perchr = []      # order of SNPs from assoc file per chr, list of chr-lists
    # for each chr a list, index 0 -> chr0, index 1 -> chr1, ..., index 26 -> chr26
    for i in xrange(27):
      self._snps_order_perchr.append([])

  def write_file_set(self, write_file):
    """ set name of output file """
    self._write_file = write_file 

  def write_file_set_append(self):
    """ switch to append mode """
    self._write_file_status_append = True
 
  def map(self, typed=False):
    """ map p-values from ASSOC file into memory """
    try:
      fh = file(self._assoc_file,"r")
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

  def map_assoc_logistic(self, typed=False):
    """ map p-values from ASSOC file into memory """
    try:
      fh = file(self._assoc_file,"r")
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
        self._typed[rs.replace("_typed","").replace("_imputed","")] = (list[11], "typed")
      # typed and/or imputed 
      else:
        if typed_pattern.search(rs):
            self._typed_imputed[rs.replace("_typed","")] = (list[11], "typed")
        else:
            self._typed_imputed[rs.replace("_imputed","")] = (list[11], "imputed")

      line = fh.readline()
    
    fh.close()

  def map_assoc_dosage(self, typed=False):
    """ map p-values from ASSOC dosage file into memory """
    try:
      fh = file(self._assoc_file,"r")
    except IOError, e:
      print e
      sys.exit(1)
    
    typed_pattern = re.compile("^.*_typed.*$")
    imputed_pattern = re.compile("^.*_imputed.*$")

    # skip header
    line = fh.readline()
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
        self._typed[rs.replace("_typed","").replace("_imputed","")] = (list[9], "typed")
      # typed and/or imputed 
      else:
        if typed_pattern.search(rs):
            self._typed_imputed[rs.replace("_typed","")] = (list[9], "typed")
        else:
            self._typed_imputed[rs.replace("_imputed","")] = (list[9], "imputed")

      line = fh.readline()
    
    fh.close()

  def map_assoc_OR_dosage(self, typed=False):
    """ map p-values and ORs from ASSOC dosage file into memory """
    try:
      fh = file(self._assoc_file,"r")
    except IOError, e:
      print e
      sys.exit(1)
    
    typed_pattern = re.compile("^.*_typed.*$")
    imputed_pattern = re.compile("^.*_imputed.*$")

    # skip header
    line = fh.readline()
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
        self._typed[rs.replace("_typed","").replace("_imputed","")] = (list[9], "typed", list[7])
      # typed and/or imputed 
      else:
        if typed_pattern.search(rs):
            self._typed_imputed[rs.replace("_typed","")] = (list[9], "typed", list[7])
        else:
            self._typed_imputed[rs.replace("_imputed","")] = (list[9], "imputed", list[7])

      line = fh.readline()
    
    fh.close()

  def map_assoc_meta(self, typed=False):
    """ map p-values from ASSOC meta-analysis file into memory """
    try:
      fh = file(self._assoc_file,"r")
    except IOError, e:
      print e
      sys.exit(1)
    
    typed_pattern = re.compile("^.*_typed.*$")
    imputed_pattern = re.compile("^.*_imputed.*$")

    # skip header
    line = fh.readline()
    line = fh.readline()
    while line:
      
      list = re.split("\s+",line)

      # delete empty elements
      if list[0] == "":
          del list[0]
      if list[-1] == "":
          del list[-1]

      rs   = str(list[2])  
      self._snps_order.append(rs.replace("_typed","").replace("_imputed",""))
      
      # typed only 
      if typed: 
        self._typed[rs.replace("_typed","").replace("_imputed","")] = (list[6], "typed")
      # typed and/or imputed 
      else:
        if typed_pattern.search(rs):
            self._typed_imputed[rs.replace("_typed","")] = (list[6], "typed")
        else:
            self._typed_imputed[rs.replace("_imputed","")] = (list[6], "imputed")

      line = fh.readline()
    
    fh.close()

  def free_map(self):
    """ free mapped ASSOC file from memory """
    self._snps_order       = []
    self._typed            = {}
    self._typed_imputed    = {}

  def rs_get_typed_hash(self):
    """ get the typed rs numbers from file in hash """
    return deepcopy(self._typed)

  def rs_get_typed_imputed_hash(self):
    """ get the typed imputed rs numbers from file in hash """
    return deepcopy(self._typed_imputed)

  def rs_get_typed_imputed_OR_hash(self):
    """ get the typed imputed rs numbers and ORs from file in hash """
    return deepcopy(self._typed_imputed)

  def rs_get(self):
    """ get the rs numbers in original order """
    return self._snps_order[:]


class Frq:
  """ class Frq implements tasks on frq files from PLINK """
 
  def __init__(self, frq_file, write_monomorphic_file=None):
    """ init """
    self._frq_file = frq_file
    self._write_monomorphic_file = write_monomorphic_file
    self._monomorphic_variants = []

  def map_monomorphic(self):
    """ map monomorphic variants into memory """
    try:
      fh = file(self._frq_file,"r")
    except IOError, e:
      print e
      sys.exit(1)
    
    # header line
    line = fh.readline().rstrip('\n')
    list = re.split("\s+",line)
    # delete empty elements
    if list[0] == "":
        del list[0]
    if list[-1] == "":
        del list[-1]
    assert(list[1] == "SNP")
    assert(list[4] == "MAF")
    
    # body
    line = fh.readline().rstrip('\n')
    while line:
      
      list = re.split("\s+",line)
      # delete empty elements
      if list[0] == "":
          del list[0]
      if list[-1] == "":
          del list[-1]

      rs   = str(list[1])  
      maf  = str(list[4])  
      
      if maf == "0":
        self._monomorphic_variants.append(rs)
      
      line = fh.readline().rstrip('\n')
    
    fh.close()

  def get_monomorphic_variants(self):
    """ get the monomorphic variants """
    return self._monomorphic_variants[:]

  def free_map(self):
    """ free mapped frq file from memory """
    self._monomorphic_variants = []

  def write_monomorphic_variants_file(self):
    """ write the monomorphic variants to file """

    try:
      fh   = file(self._frq_file,"r")
      fh_w = file(self._write_monomorphic_file,"w")
    except IOError, e:
      print e
      sys.exit(1)
    
    # header line
    line = fh.readline().rstrip('\n')
    list = re.split("\s+",line)
    # delete empty elements
    if list[0] == "":
        del list[0]
    if list[-1] == "":
        del list[-1]
    assert(list[1] == "SNP")
    assert(list[4] == "MAF")
    
    # body
    line = fh.readline().rstrip('\n')
    while line:
      
      list = re.split("\s+",line)
      # delete empty elements
      if list[0] == "":
          del list[0]
      if list[-1] == "":
          del list[-1]

      rs   = str(list[1])  
      maf  = str(list[4])  
      
      if maf == "0":
        fh_w.writelines(rs +"\n")
      
      line = fh.readline().rstrip('\n')
    
    fh.close()
    fh_w.close()

class FrqStrat:
  """ class FrqStrat implements tasks on frq.strat files from PLINK """
 
  def __init__(self, frq_file, write_nearly_monomorphic_file=None, maf_thresh=None):
    """ init """
    self._frq_file = frq_file
    self._write_monomorphic_file = write_nearly_monomorphic_file
    self._maf_thresh = float(maf_thresh)
    self._monomorphic_variants = []

  def map_monomorphic(self):
    """ map monomorphic variants into memory """
    try:
      fh = file(self._frq_file,"r")
    except IOError, e:
      print e
      sys.exit(1)
    
    # header line
    line = fh.readline().rstrip('\n')
    list = re.split("\s+",line)
    # delete empty elements
    if list[0] == "":
        del list[0]
    if list[-1] == "":
        del list[-1]
    assert(list[1] == "SNP")
    assert(list[4] == "MAF")
    
    # body
    line = fh.readline().rstrip('\n')
    while line:
      
      list = re.split("\s+",line)
      # delete empty elements
      if list[0] == "":
          del list[0]
      if list[-1] == "":
          del list[-1]

      rs   = str(list[1])  
      maf  = str(list[4])  
      
      if maf == "0":
        self._monomorphic_variants.append(rs)
      
      line = fh.readline().rstrip('\n')
    
    fh.close()

  def get_monomorphic_variants(self):
    """ get the monomorphic variants """
    return self._monomorphic_variants[:]

  def free_map(self):
    """ free mapped frq file from memory """
    self._monomorphic_variants = []

  def write_nearly_monomorphic_variants_PS_AS_IBD_PSC_file(self):
    """ write nearly monomorphic variants to file """

    try:
      fh   = file(self._frq_file,"r")
      fh_w = file(self._write_monomorphic_file,"w")
    except IOError, e:
      print e
      sys.exit(1)
    
    # header line
    line = fh.readline().rstrip('\n')
    list = re.split("\s+",line)
    # delete empty elements
    if list[0] == "":
        del list[0]
    if list[-1] == "":
        del list[-1]
    assert(list[1] == "SNP")
    assert(list[5] == "MAF")
    
    # body
    line = fh.readline().rstrip('\n')
    while line:
      
      # AS
      list = re.split("\s+",line)
      # delete empty elements
      if list[0] == "":
          del list[0]
      if list[-1] == "":
          del list[-1]
      maf_AS = float(list[5])

      # CD
      line = fh.readline().rstrip('\n')
      list = re.split("\s+",line)
      # delete empty elements
      if list[0] == "":
          del list[0]
      if list[-1] == "":
          del list[-1]
      maf_CD = float(list[5])

      # Control
      line = fh.readline().rstrip('\n')
      list = re.split("\s+",line)
      # delete empty elements
      if list[0] == "":
          del list[0]
      if list[-1] == "":
          del list[-1]
      maf_CON = float(list[5])

      # PS
      line = fh.readline().rstrip('\n')
      list = re.split("\s+",line)
      # delete empty elements
      if list[0] == "":
          del list[0]
      if list[-1] == "":
          del list[-1]
      maf_PS = float(list[5])

      # PSC
      line = fh.readline().rstrip('\n')
      list = re.split("\s+",line)
      # delete empty elements
      if list[0] == "":
          del list[0]
      if list[-1] == "":
          del list[-1]
      maf_PSC = float(list[5])

      # UC
      line = fh.readline().rstrip('\n')
      list = re.split("\s+",line)
      # delete empty elements
      if list[0] == "":
          del list[0]
      if list[-1] == "":
          del list[-1]
      maf_UC = float(list[5])

      if maf_AS < self._maf_thresh or \
         maf_CD < self._maf_thresh or \
         maf_CON < self._maf_thresh or \
         maf_PS < self._maf_thresh or \
         maf_PSC < self._maf_thresh or \
         maf_UC < self._maf_thresh :
          fh_w.writelines(list[1] + "\n")
      
      line = fh.readline().rstrip('\n')
    
    fh.close()
    fh_w.close()

class Test_missing:
  """ class Test_missing implements tasks on missing files from PLINK """

  def __init__(self, missing_file, write_file=None, threshold=1e-50):
    """ init """
    self._missing_file = missing_file
    self._write_file   = write_file
    self._threshold    = threshold

  def write_variants_file(self):
    """ write variants to file """

    try:
      fh   = file(self._missing_file,"r")
      fh_w = file(self._write_file,"w")
    except IOError, e:
      print e
      sys.exit(1)
    
    # header line
    line = fh.readline().rstrip('\n')
    list = re.split("\s+",line)
    # delete empty elements
    if list[0] == "":
        del list[0]
    if list[-1] == "":
        del list[-1]
    assert(list[1] == "SNP")
    assert(list[4] == "P")
    
    # body
    line = fh.readline().rstrip('\n')
    while line:
      
      list = re.split("\s+",line)
      # delete empty elements
      if list[0] == "":
          del list[0]
      if list[-1] == "":
          del list[-1]

      rs = str(list[1])  
      p  = float(list[4])  
      
      if p < float(self._threshold):
        fh_w.writelines(rs +"\n")
      
      line = fh.readline().rstrip('\n')
    
    fh.close()
    fh_w.close()
