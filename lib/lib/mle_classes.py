import os
import os.path
import sys
import re
import gzip
import decimal

## -- several mlinfo and mldose files can be put into an object of Mlinfo -- #
## -- or Mldose                                                           -- #

class Mlinfo:
  """ class Mlinfo implements tasks on Mlinfo files from MACH software """
  
  def __init__(self, mlinfo_file_cases, write_file=None, rsq=0.0, post_prob=0.0):
    """ init with cases """
    self._mlinfo_files_cases    = []  # names of mlinfo files from cases
    self._mlinfo_files_controls = []  # names of mlinfo files from controls
    self._mlinfo_files          = []  # names of mlinfo files from cases and controls
    self._mlinfo_files_cases.append(mlinfo_file_cases) 
    self._mlinfo_files.append(mlinfo_file_cases) 
    self._write_file = write_file     # name of output file
    self._write_file_filtered_snps = None 
    self._rsq = float(rsq)             # chosen rsquare threshold after imputation
    self._post_prob = float(post_prob) # chosen post_prob threshold after imputation
    self._map_snps = [] 

  def map_snps(self, file_number=0):
    """ map mlinfo file into memory """
   
    mlinfo_files_numof = len(self._mlinfo_files)
    if mlinfo_files_numof <= file_number:
      return

    try:
      fh  = gzip.open(self._mlinfo_files[file_number], "rb")
    except IOError, e:
      print e
      sys.exit(1)
    
    comment_pattern = re.compile("^#.*$")
    header_pattern = re.compile("^SNP.*$")
    
    # skip comment lines that start with "#" and header line
    line = fh.readline()
    while comment_pattern.search(line):
      line = fh.readline()
    
    if not header_pattern.search(line):
      print >> sys.stderr, "error: no header in mlinfo file." 
      sys.exit(1)
    line = fh.readline() #skip header
    while line:
      self._map_snps.append(re.split("\s+",line)[0])
      line = fh.readline()

    fh.close()

  def snps_get(self):
    """ get mapped snps from mlinfo file in order they appeared in file """
    return self._map_snps

  def snps_get_hash(self):
    """ get mapped snps from mlinfo file im hash """
    snp_hash = {}
    for snp in self._map_snps:
      snp_hash[snp] = 1 
    return snp_hash

  def free_snps(self):
    """ free mapped snps from mlinfo file """
    self._map_snps = []

  def rsq_set(self, rsq):
    """ set rsq for filtering """
    self._rsq = float(rsq)

  def post_prob_set(self, post_prob):
    """ set rsq for filtering """
    self._post_prob = float(post_prob)

  def write_file_set(self, write_file):
    """ set name of output file """
    self._write_file = write_file 

  def write_file_set_filtered_SNPs(self, write_file_filtered_snps):
    """ set name of output file for filtered SNPs """
    self._write_file_filtered_snps = write_file_filtered_snps 

  def mlinfo_file_add_controls(self, mlinfo_file_controls):
    """ add another mlinfo file from controls to mlinfo file list """
    self._mlinfo_files_controls.append(mlinfo_file_controls)
    self._mlinfo_files.append(mlinfo_file_controls)

  def mlinfo_file_add_cases(self, mlinfo_file_cases):
    """ add another mlinfo file from cases to mlinfo file list """
    self._mlinfo_files_cases.append(mlinfo_file_cases)
    self._mlinfo_files.append(mlinfo_file_cases)

  def check_same_snps_same_order_same_alleles(self):
    """ check for same snps and same snp order and for same minor alleles in
        all mlinfo file """
    files_numof = len(self._mlinfo_files)
    if files_numof < 1:
      return

    # first mlinfo file
    rs_previous = [] # rs numbers
    a1_previous = [] # allele 1
    a2_previous = [] # allele 2
    try:
      fh = gzip.open(self._mlinfo_files[0], "rb")
    except IOError, e:
      print e
      sys.exit(1)
    line = fh.readline()
    while line:
      rs_previous.append(re.split("\s+",line)[0])
      a1_previous.append(re.split("\s+",line)[1])
      a2_previous.append(re.split("\s+",line)[2])
      line = fh.readline()
    fh.close()

    # next mlinfo files
    for i in xrange(1,files_numof):
      rs_next = []
      a1_next = []
      a2_next = []
      try:
        fh = gzip.open(self._mlinfo_files[i], "rb")
      except IOError, e:
        print e
        sys.exit(1)
      line = fh.readline()
      while line:
        rs_next.append(re.split("\s+",line)[0])
        a1_next.append(re.split("\s+",line)[1])
        a2_next.append(re.split("\s+",line)[2])
        line = fh.readline()
      fh.close()
      if rs_previous != rs_next or a1_previous != a1_previous \
	  or a2_previous != a2_previous:
	print >> sys.stderr, "error: mlinfo files differ." 
        sys.exit(1)
      rs_previous = rs_next
      a1_previous = a1_next
      a2_previous = a2_next

  def paste_mlinfo_files(self):
    """ paste all mlinfo files into one file """

    files_numof = len(self._mlinfo_files)
    if files_numof < 1:
      return
    
    fh = [] # list of mlinfo filehandles
    for i in xrange(files_numof):
      fh.append(gzip.open(self._mlinfo_files[i], "rb"))
    out = gzip.open(self._write_file, "w")

    line = fh[0].readline().replace("\n","")
    while line:
      out.writelines(line + " ") 
      for i in xrange(1,files_numof-1):
        line = fh[i].readline().replace("\n","")
        out.writelines(line + " ") 
      line = fh[files_numof-1].readline().replace("\n","")
      out.writelines(line) 
      out.writelines("\n")
      line = fh[0].readline().replace("\n","")
    
    for i in xrange(files_numof):
      fh[i].close()
    out.close()

  def select_filtered_snps_rsq_postprob(self, file_number=0):
    """ filter snps by rsq and average posterior probability for cases mlinfo file """
    
    mlinfo_files_numof = len(self._mlinfo_files)
    if mlinfo_files_numof <= file_number:
      return

    try:
      fh  = gzip.open(self._mlinfo_files[file_number], "rb")
      out = gzip.open(self._write_file, "w")
    except IOError, e:
      print e
      sys.exit(1)
    
    #### skip comment lines that start with "#" and header line
    #### write comment lines
    ###comment_pattern = re.compile("^#.*$")
    ###while comment_pattern.search(line):
    ###  out.writelines(line) 
    ###  line = fh.readline()
    #### filter snps for r2 and postprob
    ###out.writelines("# filtered Rsq     >=%s\n" %(str(self._rsq))) 
    ###out.writelines("# filtered Quality >=%s\n" %(str(self._post_prob))) 
    
    # read from header how many mlinfo have been merge
    # and which columns to check for r2 and postprob
    # write header
    header_pattern = re.compile("^SNP.*$")
    line = fh.readline()
    if not header_pattern.search(line):
      print >> sys.stderr, "error: no header in merged mlinfo file." 
      sys.exit(1)
    mlinfo_files_numof = len(re.findall('SNP', line))
    post_prob_columns = []
    r2_columns = []
    for i in xrange(1,len(re.findall("SNP", line))+1):
       post_prob_columns.append(i*7-2)
       r2_columns.append(i*7-1)
    assert(mlinfo_files_numof     == len(r2_columns))
    assert(len(post_prob_columns) == len(r2_columns))
    out.writelines(line) 

    # write to out file
    line = fh.readline()
    while line:
      list = re.split("\s+",line)
      keep_snp = True
      for i in xrange(mlinfo_files_numof):
	if   float(list[post_prob_columns[i]]) < self._post_prob:
          keep_snp = False
	  break
	elif float(list[r2_columns[i]]) < self._rsq:  
          keep_snp = False
	  break
      if keep_snp:
        out.writelines(line) 
      line = fh.readline()

    fh.close()
    out.close()

  def snp_list_first_allele_get(self, new_snp_list=[]):
    """ get list of first alleles from new_snp_list """
    # copy list
    new_snp_list_first_allele = new_snp_list[:]
    
    for i in xrange(len(new_snp_list)):
      # check allele order
      a1 = ""
      #if cmp(new_snp_list[i][4],new_snp_list[i][5]) == 0:
      #  print >> sys.stderr, "error: monoallelic snps not allowed, rs=" +\
      #  str(new_snp_list[i][1])
      #  sys.exit(1)
      #elif cmp(new_snp_list[i][4],new_snp_list[i][5]) > 0:
      if cmp(new_snp_list[i][4],new_snp_list[i][5]) > 0:
	# switch alleles
        a1 = new_snp_list[i][5] 
        new_snp_list_first_allele[i] = a1
      else:
	# do not switch alleles
        a1 = new_snp_list[i][4] 
        new_snp_list_first_allele[i] = a1
    
    return new_snp_list_first_allele
 

  def write_new_mlinfo(self, freq1_only_typed_snps=[], new_snp_list=[]):
    """ write new mlinfo file from new_snp_list """
    try:
      out = gzip.open(self._write_file, "w")
    except IOError, e:
      print e
      sys.exit(1)
    assert(len(new_snp_list) == len(freq1_only_typed_snps))
    header_line = "SNP\tAl1\tAl2\tFreq1\tMAF\tQuality\tRsq\n"

    out.writelines(header_line) 
    for i in xrange(len(new_snp_list)):
     
      # check allele order
      a1 = ""
      a2 = "" 
      if cmp(new_snp_list[i][4],new_snp_list[i][5]) == 0:
        print >> sys.stderr, "error: monoallelic snps not allowed, rs=" +\
	    str(new_snp_list[i][1])
        sys.exit(1)
      elif cmp(new_snp_list[i][4],new_snp_list[i][5]) > 0:
	# switch alleles
        a1 = new_snp_list[i][5] 
	a2 = new_snp_list[i][4] 
      else:
	# do not switch alleles
        a1 = new_snp_list[i][4] 
	a2 = new_snp_list[i][5] 

      freq1 = float('%.4f' %freq1_only_typed_snps[i])
      if freq1_only_typed_snps[i] > 0.5: 
        maf = 1.0 - freq1
      else:
        maf = freq1_only_typed_snps[i]
      out.writelines(new_snp_list[i][1] + "\t" + a1 + "\t" + a2 +\
          "\t%.4f\t%.4f\t1.0\t1.0\n" %(freq1, maf))

    out.close()

  def snps_al1_get(self, file_number=0):
    """ return list: (snp, allele1) """
  
    snps_allele1_list = []

    mlinfo_files_numof = len(self._mlinfo_files)
    if mlinfo_files_numof <= file_number:
      return

    try:
      fh  = gzip.open(self._mlinfo_files[file_number], "rb")
    except IOError, e:
      print e
      sys.exit(1)
    
    comment_pattern = re.compile("^#.*$")
    header_pattern = re.compile("^SNP.*$")
    
    # skip comment lines that start with "#" and header line
    line = fh.readline()
    while comment_pattern.search(line):
      line = fh.readline()
    
    # scan header line
    if not header_pattern.search(line):
      print >> sys.stderr, "error: no header in mlinfo file." 
      sys.exit(1)
    line = fh.readline() #skip header
    
    # scan whole mlinfo file line by line
    while line:
      list = re.split("\s+",line)
      snp = list[0]
      al1 = list[1]
      snps_allele1_list.append( (snp, al1) )
      line = fh.readline()

    fh.close()
    return snps_allele1_list


class Mldose:
  """ class Mldose implements tasks on Mldose files from MACH software """
  
  def __init__(self, mldose_file_cases, write_file=None, write_file_ids=None):
    """ init with mldose file from cases"""
    self._mldose_files_cases    = [] # names of mldose files from cases
    self._mldose_files_controls = [] # names of mldose files from controls
    self._mldose_files          = [] # names of mldose files from cases and controls
    self._mldose_files_cases.append(mldose_file_cases) 
    self._mldose_files.append(mldose_file_cases) 
    self._write_file     = write_file     # name of output file
    self._write_file_ids = write_file_ids # name of output file for ids
    self._samples_controls = []      # control sample list
    self._samples_cases    = []      # cases sample list

  def write_file_set(self, write_file):
    """ set name of output file """
    self._write_file = write_file 

  def write_file_set_ids(self, write_file_ids):
    """ set name of output file """
    self._write_file_ids = write_file_ids 

  def mldose_file_add_cases(self, mldose_file_cases):
    """ add another mldose file from cases to mldose file list """
    self._mldose_files_cases.append(mldose_file_cases)
    self._mldose_files.append(mldose_file_cases)

  def mldose_file_add_controls(self, mldose_file_controls):
    """ add another mldose file from controls to mldose file list """
    self._mldose_files_controls.append(mldose_file_controls)
    self._mldose_files.append(mldose_file_controls)

  def mldose_files_cases_get(self):
    """ get name of cases mldose file """
    return self._mldose_files_cases

  def mldose_files_get(self):
    """ get name of mldose file """
    return self._mldose_files 

  def ids_all_get(self, file_number=0):
    """ select all ids from mldose file and return ids as list """
    mldose_files_numof = len(self._mldose_files)
    if mldose_files_numof <= file_number:
      return

    try:
      fh = gzip.open(self._mldose_files[file_number], "rb")
    except IOError, e:
      print e
      sys.exit(1)
  
    ids_list = []
    line = fh.readline()
    while line:
      id = re.split("\s+",line)[0]
      ids_list.append(id)
      line = fh.readline()
     
    fh.close()
    return ids_list

  def select_cases_from_mldose(self, cases=[], file_number=0):
    """ select cases from mldose file(s) and write to output file """
    mldose_files_cases_numof = len(self._mldose_files_cases)
    if mldose_files_cases_numof <= file_number:
      return
 
    try:
      fh = gzip.open(self._mldose_files_cases[file_number], "rb")
      out = gzip.open(self._write_file, "w")
      if self._write_file_ids:
        out2 = gzip.open(self._write_file_ids, "w")
    except IOError, e:
      print e
      sys.exit(1)
    
    line = fh.readline()
    while line:
      id = re.split("\s+",line)[0]
      try:
        index = cases.index(id)
      except ValueError:
	# skip sample in mldose file
        line = fh.readline()
	continue
      out.writelines(line) 
      if self._write_file_ids:
        out2.writelines(id + "\n") 
      del cases[index]
      line = fh.readline()
    
    fh.close()
    out.close()
    if self._write_file_ids:
      out2.close()

  def select_cases_from_mldose_and_append(self, cases=[], file_number=0):
    """ select cases from mldose file(s) and write to output file
        file_number starts with 0 """
    
    mldose_files_cases_numof = len(self._mldose_files_cases)
    if mldose_files_cases_numof <= file_number:
      return
   
    try:
      fh  = gzip.open(self._mldose_files_cases[file_number], "rb")
      out = gzip.open(self._write_file, "a")
      if self._write_file_ids:
        out2 = gzip.open(self._write_file_ids, "a")
    except IOError, e:
      print e
      sys.exit(1)
    
    line = fh.readline()
    while line:
      id = re.split("\s+",line)[0]
      try:
        index = cases.index(id)
      except ValueError:
	# skip sample in mldose file
        line = fh.readline()
	continue
      out.writelines(line) 
      if self._write_file_ids:
        out2.writelines(id + "\n") 
      del cases[index]
      line = fh.readline()

    fh.close()
    out.close()
    if self._write_file_ids:
      out2.close()
 
  def select_controls_from_mldose(self, controls=[], file_number=0):
    """ select controls from mldose file(s) and write to output file """
    mldose_files_controls_numof = len(self._mldose_files_controls)
    if mldose_files_controls_numof <= file_number:
      return
  
    try:
      fh = gzip.open(self._mldose_files_controls[file_number], "rb")
      out = gzip.open(self._write_file, "w")
      if self._write_file_ids:
        out2 = gzip.open(self._write_file_ids, "w")
    except IOError, e:
      print e
      sys.exit(1)
    
    line = fh.readline()
    while line:
      id = re.split("\s+",line)[0]
      try:
        index = controls.index(id)
      except ValueError:
	# skip sample in mldose file
        line = fh.readline()
	continue
      out.writelines(line) 
      if self._write_file_ids:
        out2.writelines(id + "\n") 
      del controls[index]
      line = fh.readline()
    
    fh.close()
    out.close()
    if self._write_file_ids:
      out2.close()
 
  def select_controls_from_mldose_and_append(self, controls=[], file_number=0):
    """ select controls from mldose file(s) and write to output file
        file_number starts with 0 """
    
    mldose_files_controls_numof = len(self._mldose_files_controls)
    if mldose_files_controls_numof <= file_number:
      return
   
    try:
      fh  = gzip.open(self._mldose_files_controls[file_number], "rb")
      out = gzip.open(self._write_file, "a")
      if self._write_file_ids:
        out2 = gzip.open(self._write_file_ids, "a")
    except IOError, e:
      print e
      sys.exit(1)
    
    line = fh.readline()
    while line:
      id = re.split("\s+",line)[0]
      try:
        index = controls.index(id)
      except ValueError:
	# skip sample in mldose file
        line = fh.readline()
	continue
      out.writelines(line) 
      if self._write_file_ids:
        out2.writelines(id + "\n") 
      del controls[index]
      line = fh.readline()

    fh.close()
    out.close()
    if self._write_file_ids:
      out2.close()

  def select_ids_cases_controls_from_mldose(self, cases_controls=[], file_number=0):
    """ select cases from mldose file(s) and write to output file """
    mldose_files_cases_numof = len(self._mldose_files_cases)
    if mldose_files_cases_numof <= file_number:
      return
 
    try:
      fh = gzip.open(self._mldose_files_cases[file_number], "rb")
      out2 = gzip.open(self._write_file_ids, "w")
    except IOError, e:
      print e
      sys.exit(1)
    
    line = fh.readline()
    while line:
      id = re.split("\s+",line)[0]
      try:
        index = cases_controls.index(id)
      except ValueError:
        # skip sample in mldose file
        line = fh.readline()
        continue
      if self._write_file_ids:
        out2.writelines(id + "\n") 
      del cases_controls[index]
      line = fh.readline()
    
    fh.close()
    out2.close()


  def select_freq1_from_mldose(self, file_number=0):
    """ select freq1 in order of mldose file and return
        note that control ids are utilized for extracting cases
        select freq1 from controls only in order of mldose file and return """
    mldose_files_numof = len(self._mldose_files)
    if mldose_files_numof <= file_number:
      return
  
    try:
      fh = gzip.open(self._mldose_files[file_number], "rb")
    except IOError, e:
      print e
      sys.exit(1)
   
    freq1_snps = []
    sample_counter = 0
   
    # first line, generate list
    line = fh.readline().replace("\n","")
    dosages = re.split("\s+",line)[2:]

    # if last element is empty then remove it
    if dosages[-1] == "":
      dosages.pop()

    dosages_numof = len(dosages)
    for i in xrange(dosages_numof):
      # generate list of dosages
      freq1_snps.append(float(dosages[i]))
    sample_counter = sample_counter + 1

    line = fh.readline().replace("\n","")
    # sum all dosages for allele 1
    while line:
      dosages = re.split("\s+",line)[2:]

      # if last element is empty then remove it
      if dosages[-1] == "":
        dosages.pop()
 
      if len(dosages) != dosages_numof:
        print >> sys.stderr, "error: different num of columns in mldose file."
        sys.exit(1)

      for i in xrange(dosages_numof):
        freq1_snps[i] = freq1_snps[i] + float(dosages[i])
      
      sample_counter = sample_counter + 1
      line = fh.readline().replace("\n","")
    
    fh.close()
  
    # calc freq1 for allele 1
    for i in xrange(dosages_numof):
      freq1_snps[i] = freq1_snps[i] / (2*float(sample_counter))
    
    return freq1_snps

  def select_freq1_from_mldose_all_cases_controls(self, ids_controls_list, file_number=0):
    """ select freq1 in order of mldose file and return
        note that control ids are utilized for extracting cases
        select freq1 from controls only in order of mldose file and return """
    mldose_files_numof = len(self._mldose_files)
    if mldose_files_numof <= file_number:
      return
  
    try:
      fh = gzip.open(self._mldose_files[file_number], "rb")
    except IOError, e:
      print e
      sys.exit(1)
  
    ids_controls_hash = {}
    for id in ids_controls_list:
        ids_controls_hash[id] = True

    freq1_snps_all      = []
    freq1_snps_cases    = []
    freq1_snps_controls = []
    sample_counter_all = 0
    sample_counter_ca = 0
    sample_counter_co = 0
    dosages_numof = 0 

    line = fh.readline().replace("\n","")
    # sum all dosages for allele 1
    while line:
      list = re.split("\s+",line)
      id = list[0]
      dosages = list[2:]

      # if last element is empty then remove it
      if dosages[-1] == "":
        dosages.pop()
      
      # find first case
      if sample_counter_all == 0:

        # count only if case
        if not ids_controls_hash.has_key(id):
          dosages_numof = len(dosages)
          for i in xrange(dosages_numof):
            # generate list of dosages
            freq1_snps_cases.append(float(dosages[i]))
            freq1_snps_controls.append(float(0.0))
            freq1_snps_all.append(float(dosages[i]))
          sample_counter_ca = sample_counter_ca + 1
        
        # count only if control
        elif ids_controls_hash.has_key(id):
          dosages_numof = len(dosages)
          for i in xrange(dosages_numof):
            # generate list of dosages
            freq1_snps_cases.append(float(0.0))
            freq1_snps_controls.append(float(dosages[i]))
            freq1_snps_all.append(float(dosages[i]))
          sample_counter_co = sample_counter_co + 1
       
        # error: neither case nor control
        else:
          print >> sys.stderr, "error: sample id neither case nor control."
          sys.exit(1)
          
        sample_counter_all = sample_counter_all + 1
        line = fh.readline().replace("\n","")
 
      # parse next cases/controls
      else:

        if len(dosages) != dosages_numof:
          print >> sys.stderr, "error: different num of columns in mldose file."
          sys.exit(1)

        # count only if case
        if not ids_controls_hash.has_key(id):
          for i in xrange(dosages_numof):
            freq1_snps_cases[i] = freq1_snps_cases[i] + float(dosages[i])
            # count all samples
            freq1_snps_all[i] = freq1_snps_all[i] + float(dosages[i])
          sample_counter_ca = sample_counter_ca + 1
          sample_counter_all = sample_counter_all + 1
        
        # count only if controls
        elif ids_controls_hash.has_key(id):
          for i in xrange(dosages_numof):
            freq1_snps_controls[i] = freq1_snps_controls[i] + float(dosages[i])
            # count all samples
            freq1_snps_all[i] = freq1_snps_all[i] + float(dosages[i])
          sample_counter_co = sample_counter_co + 1
          sample_counter_all = sample_counter_all + 1
      
        # error: neither case nor control
        else:
          print >> sys.stderr, "error: sample id neither case nor control."
          sys.exit(1)

        line = fh.readline().replace("\n","")
  
    fh.close()
  
    # calc freq1 for allele 1
    for i in xrange(dosages_numof):
      freq1_snps_cases[i]    = freq1_snps_cases[i]    / (2*float(sample_counter_ca))
      freq1_snps_controls[i] = freq1_snps_controls[i] / (2*float(sample_counter_co))
      freq1_snps_all[i]      = freq1_snps_all[i]      / (2*float(sample_counter_all))
    
    return freq1_snps_all, freq1_snps_cases, freq1_snps_controls


def compare(a,b):
    """ compare decimals """
    if float(a) < float(b):
        return -1
    elif a == b:
        return 0
    else:
        return 1


def col_ten(t):
    """ sort 10th column """
    return t[9]


class Mlids:
  """ class Mlids implements tasks on Mlids files from MACH software """
  
  def __init__(self, mlids_file_cases, write_file=None):
    """ init with cases """
    self._mlids_files_cases    = []  # names of mlids files from cases
    self._mlids_files_controls = []  # names of mlids files from controls
    self._mlids_files          = []  # names of mlids files from cases and controls
    self._mlids_files_cases.append(mlids_file_cases) 
    self._mlids_files.append(mlids_file_cases) 
    self._write_file = write_file     # name of output file
    self._write_file_filtered_snps = None 
    self._map_indivs = []             # sample names

  def map_indivs(self, file_number=0):
    """ map mlids file into memory """
   
    mlids_files_numof = len(self._mlids_files)
    if mlids_files_numof <= file_number:
      return

    try:
      fh  = gzip.open(self._mlids_files[file_number], "rb")
    except IOError, e:
      print e
      sys.exit(1)
    
    # skip comment lines that start with "#" and header line
    comment_pattern = re.compile("^#.*$")
    line = fh.readline()
    while comment_pattern.search(line):
      line = fh.readline()
    
    while line:
      self._map_indivs.append(re.split("\s+",line)[0])
      line = fh.readline()

    fh.close()

  def indivs_numof_get(self):
    """ get number of mapped indivs from mlids file """
    return len(self._map_indivs)

  def free_indivs(self):
    """ free mapped indivs from mlids file """
    self._map_indivs = []



class Mach2dat:
    """ class Mach2dat implements tasks on results files from MACH2DAT software """
    
    def __init__(self, mach2dat_file, write_file=None, mlinfo_filter=None):
        """ init """
        self._mach2dat_files          = []  # names of mach2dat files
        self._mach2dat_files.append(mach2dat_file) 
        self._write_file = write_file       # name of output file
        self._mlinfo_filter = mlinfo_filter # name of mlinfo file of filtered snps
        self._map_snps_header = []
        self._map_snps = {}

    def write_file_set(self, write_file):
        """ set name of output file """
        self._write_file = write_file 

    def mlinfo_filter_set(self, mlinfo_filter):
        """ set name of mlinfo file of filtered snps """
        self._mlinfo_filter = mlinfo_filter

    def map_snps(self, file_number=0):
        """ map mach2dat file into memory """
        
        mach2dat_files_numof = len(self._mach2dat_files)
        if mach2dat_files_numof <= file_number:
          return
        
        try:
          fh  = file(self._mach2dat_files[file_number], "r")
        except IOError, e:
          print e
          sys.exit(1)
        
        header_pattern = re.compile("^TRAIT.*$")
        
        # skip lines until header line
        line = fh.readline()
        while not header_pattern.search(line):
          line = fh.readline()
        
        # write new header
        header_columns = re.split("\s+",line)[0:3] + re.split("\s+",line)[5:12]
        self._map_snps_header = header_columns[:]

        # read all snps
        line = fh.readline().replace("\n","")
        while line:
        
            list = re.split("\s+",line)
            # delete empty elements
            if list[-1] == "":
                del list[-1]
            # workaround because of MACH2DAT bug:
            # If affy ids in output, then col2 and 3 are accidently merged (one
            # less column)

            # no workaround
            if len(list) == 12:
                if list[5] != "NA":
                    self._map_snps[list[1]] = (list[0], list[1], list[2],\
                                           decimal.Decimal(list[5]),\
                                           decimal.Decimal(list[6]),\
                                           decimal.Decimal(list[7]),\
                                           decimal.Decimal(list[8]),\
                                           decimal.Decimal(list[9]),\
                                           decimal.Decimal(list[10]),\
                                           decimal.Decimal(list[11]))
            
            # workaround, col 2 and col 3 merged and has to be splitted
            elif len(list) == 11:
                if list[4] != "NA" and list[5] != "NA":
                    split_list = list[1].split(",")
                    marker  = split_list[0][:-1] 
                    alleles = split_list[0][-1:] + "," + split_list[1]

                    self._map_snps[marker] = (list[0], marker, alleles,\
                                           decimal.Decimal(list[4]),\
                                           decimal.Decimal(list[5]),\
                                           decimal.Decimal(list[6]),\
                                           decimal.Decimal(list[7]),\
                                           decimal.Decimal(list[8]),\
                                           decimal.Decimal(list[9]),\
                                           decimal.Decimal(list[10]))
            line = fh.readline().replace("\n","")

        fh.close()

    def free_snps(self):
        """ free mapping of snps """
        self._map_snps_header = []
        self._map_snps = {}

    def write_snps_filtered_sorted(self, snps_rs_al1_al1freq, maf_threshold, file_number=0):
        """ write snps filtered by mlinfo_filter file and snps
            sorted by p-value with additionally al1 and al1_freq
            and al1_freq for cases and controls separately """
        
        if self._mlinfo_filter == None:
            print >> sys.stderr, "error: no mlinfo_filter file specified." 
            sys.exit(1)
        try:
            fh  = file(self._write_file, "w")
        except IOError, e:
            print e
            sys.exit(1)
 
        mlinfo = Mlinfo(mlinfo_file_cases=self._mlinfo_filter)
        mlinfo.map_snps()
        snps_filtered = mlinfo.snps_get()
        mlinfo.free_snps() ; del mlinfo

        mach2dat_filtered = []
        for snp in snps_filtered:
            if self._map_snps.has_key(snp):
                mach2dat_filtered.append(self._map_snps[snp])
        # sort by pvalues
        mach2dat_filtered.sort(cmp=compare, key=col_ten) # sort 10th column
   
        # print new header
        for i in xrange(len(self._map_snps_header)):
            if i == 0:
                fh.writelines("%s" %(self._map_snps_header[i]))
            else:
                fh.writelines(" %s" %(self._map_snps_header[i]))
        fh.writelines(" AL1 FREQ1 MAF FREQ1_CA FREQ1_CO\n")

        # print the rest
        for line in mach2dat_filtered:

            # round 4 decimals
            al1_freq    = decimal.Decimal(str(round(snps_rs_al1_al1freq[line[1]][1], 4)))
            al1_freq_ca = decimal.Decimal(str(round(snps_rs_al1_al1freq[line[1]][2], 4)))
            al1_freq_co = decimal.Decimal(str(round(snps_rs_al1_al1freq[line[1]][3], 4)))
            
            # determine MAF for SNP
            maf = decimal.Decimal(str(0.0))
            if al1_freq >= decimal.Decimal(str(0.5)):
                maf = 1 - al1_freq 
            else:
                maf = al1_freq 
            
            # filter for maf frequency
            if ( al1_freq >= decimal.Decimal(str(maf_threshold)) and \
                 al1_freq <= decimal.Decimal(str(1-maf_threshold)) ):
                for i in xrange(len(line)):
                    if i == 0:
                        fh.writelines("%s"  %(line[i]))
                    else:
                        fh.writelines(" %s" %(line[i]))
                al1      = snps_rs_al1_al1freq[line[1]][0]
                fh.writelines( " %s %s %s %s %s\n" %(al1, al1_freq, maf, \
                                                     al1_freq_ca, al1_freq_co) )
        
        fh.close()

    def write_snps_sorted(self, snps_rs_al1_al1freq, maf_threshold, file_number=0):
        """ write snps filtered by mlinfo_filter file and snps
            sorted by p-value with additionally al1 and al1_freq"""
        
        try:
            fh  = file(self._write_file, "w")
        except IOError, e:
            print e
            sys.exit(1)
 
        mach2dat = []
        keys = self._map_snps.iterkeys()
        for snp in keys:
            mach2dat.append(self._map_snps[snp])
        
        # sort by pvalues
        mach2dat.sort(cmp=compare, key=col_ten) # sort 10th column
   
        # print new header
        for i in xrange(len(self._map_snps_header)):
            if i == 0:
                fh.writelines("%s" %(self._map_snps_header[i]))
            else:
                fh.writelines(" %s" %(self._map_snps_header[i]))
        fh.writelines(" AL1 FREQ1 MAF FREQ1_CA FREQ1_CO\n")

        # print the rest
        for line in mach2dat:
            
            # round 4 decimals
            al1_freq    = decimal.Decimal(str(round(snps_rs_al1_al1freq[line[1]][1], 4)))
            al1_freq_ca = decimal.Decimal(str(round(snps_rs_al1_al1freq[line[1]][2], 4)))
            al1_freq_co = decimal.Decimal(str(round(snps_rs_al1_al1freq[line[1]][3], 4)))
            
            # determine MAF for SNP
            maf = decimal.Decimal(str(0.0))
            if al1_freq >= decimal.Decimal(str(0.5)):
                maf = 1 - al1_freq 
            else:
                maf = al1_freq 
     
            # filter for maf frequency
            if (al1_freq >= decimal.Decimal(str(maf_threshold)) and \
                al1_freq <= decimal.Decimal(str(1-maf_threshold)) ):
                for i in xrange(len(line)):
                    if i == 0:
                        fh.writelines("%s"  %(line[i]))
                    else:
                        fh.writelines(" %s" %(line[i]))
                al1      = snps_rs_al1_al1freq[line[1]][0]
                fh.writelines( " %s %s %s %s %s\n" %(al1, al1_freq, maf, \
                                                    al1_freq_ca, al1_freq_co) )
        
        fh.close()

class Mach2dat_merge:
    """ class Mach2dat_merge implements tasks on filtered/sorted 
        results files from MACH2DAT software """
    
    def __init__(self, write_file=None, write_file_assoc=None):
        """ init """
        self._mach2dat_files          = []  # names of mach2dat files
        self._write_file = write_file       # name of output file
        self._write_file_assoc = write_file_assoc  # name of assoc output file
        self._write_file_combined_short = None # name of combined short output file
        self._map_snps_header = []
        self._map_snps = []   # list of dictionaries
        self._map_snps_lists = []   # list of lists

    def write_file_set(self, write_file):
        """ set name of output file """
        self._write_file = write_file 

    def write_file_assoc_set(self, write_file_assoc):
        """ set name of output assoc file """
        self._write_file_assoc = write_file_assoc 

    def write_file_combined_short_set(self, write_file_combined_short):
        """ set name of combined short output file """
        self._write_file_combined_short = write_file_combined_short 

    def mach2dat_file_add(self, mach2dat_file):
        """ add another mach2dat file from cases to mach2dat file list """
        self._mach2dat_files.append(mach2dat_file)

    def map_snps_filter_snps(self, snp_filter_dict, file_number=0):
        """ map mach2dat file into memory and filter by using snp_filter_dict """
        
        mach2dat_files_numof = len(self._mach2dat_files)
        if mach2dat_files_numof <= file_number:
            return
        
        try:
          fh  = file(self._mach2dat_files[file_number], "r")
        except IOError, e:
          print e
          sys.exit(1)
        
        header_pattern = re.compile("^TRAIT.*$")
        
        # skip lines until header line
        line = fh.readline().replace("\n","")
        while not header_pattern.search(line):
          line = fh.readline().replace("\n","")
        
        # write new header
        self._map_snps_header.append(re.split("\s+",line)[:])

        # read all snps
        line = fh.readline()
        map_snps = {}
        while line:
            list = re.split("\s+",line)
            if list[3] != "NA" and (not snp_filter_dict.has_key(list[1])):
                map_snps[list[1]] = (list[0], list[1], list[2],\
                                           decimal.Decimal(list[3]),\
                                           decimal.Decimal(list[4]),\
                                           decimal.Decimal(list[5]),\
                                           decimal.Decimal(list[6]),\
                                           decimal.Decimal(list[7]),\
                                           decimal.Decimal(list[8]),\
                                           decimal.Decimal(list[9]),\
                                           list[10],\
                                           decimal.Decimal(list[11]),\
                                           decimal.Decimal(list[12]),\
                                           decimal.Decimal(list[13]),\
                                           decimal.Decimal(list[14]))
            line = fh.readline()
                
        self._map_snps.append(map_snps.copy())

        fh.close()

    def map_snps(self, file_number=0):
        """ map mach2dat file into memory """
        
        mach2dat_files_numof = len(self._mach2dat_files)
        if mach2dat_files_numof <= file_number:
            return
        
        try:
          fh  = file(self._mach2dat_files[file_number], "r")
        except IOError, e:
          print e
          sys.exit(1)
        
        header_pattern = re.compile("^TRAIT.*$")
        
        # skip lines until header line
        line = fh.readline().replace("\n","")
        while not header_pattern.search(line):
          line = fh.readline().replace("\n","")
        
        # write new header
        self._map_snps_header.append(re.split("\s+",line)[:])

        # read all snps
        line = fh.readline()
        map_snps = {}
        while line:
            list = re.split("\s+",line)
            if list[3] != "NA":
                map_snps[list[1]] = (list[0], list[1], list[2],\
                                           decimal.Decimal(list[3]),\
                                           decimal.Decimal(list[4]),\
                                           decimal.Decimal(list[5]),\
                                           decimal.Decimal(list[6]),\
                                           decimal.Decimal(list[7]),\
                                           decimal.Decimal(list[8]),\
                                           decimal.Decimal(list[9]),\
                                           list[10],\
                                           decimal.Decimal(list[11]),\
                                           decimal.Decimal(list[12]),\
                                           decimal.Decimal(list[13]),\
                                           decimal.Decimal(list[14]))
            line = fh.readline()
                
        self._map_snps.append(map_snps.copy())

        fh.close()

    def map_snps_list(self, file_number=0):
        """ map mach2dat file into memory but into list instead
            of hash """
        
        mach2dat_files_numof = len(self._mach2dat_files)
        if mach2dat_files_numof <= file_number:
            return
        
        try:
          fh  = file(self._mach2dat_files[file_number], "r")
        except IOError, e:
          print e
          sys.exit(1)
        
        header_pattern = re.compile("^TRAIT.*$")
        
        # skip lines until header line
        line = fh.readline().replace("\n","")
        while not header_pattern.search(line):
          line = fh.readline().replace("\n","")
        
        # write new header
        self._map_snps_header.append(re.split("\s+",line)[:])

        # read all snps
        line = fh.readline()
        map_snps = []
        while line:
            list = re.split("\s+",line)
            if list[3] != "NA":
                map_snps.append( (list[0], list[1], list[2],\
                                           decimal.Decimal(list[3]),\
                                           decimal.Decimal(list[4]),\
                                           decimal.Decimal(list[5]),\
                                           decimal.Decimal(list[6]),\
                                           decimal.Decimal(list[7]),\
                                           decimal.Decimal(list[8]),\
                                           decimal.Decimal(list[9]),\
                                           list[10],\
                                           decimal.Decimal(list[11]),\
                                           decimal.Decimal(list[12]),\
                                           decimal.Decimal(list[13]),\
                                           decimal.Decimal(list[14])) )
            line = fh.readline()
               
        # clone list
        self._map_snps_lists.append(map_snps[:])

        fh.close()

    def free_snps(self):
        """ free mapping of snps """
        self._map_snps_header = []
        self._map_snps = []
        self._map_snps_lists = []

    def write_mach2dat_merged_sorted(self, file_number=0):
        """ write all snps to merged file, sorted """ 
        
        try:
            fh   = file(self._write_file, "w")
            fh2  = file(self._write_file_assoc, "w")
        except IOError, e:
            print e
            sys.exit(1)

        mach2dat = []
        for j in xrange(len(self._map_snps)):
             keys = self._map_snps[j].iterkeys()
             for snp in keys:
                 mach2dat.append(self._map_snps[j][snp])
             
        # sort by pvalues
        mach2dat.sort(cmp=compare, key=col_ten) # sort 10th column
  
        # print new header from first mach2dat result file
        for i in xrange(len(self._map_snps_header[0])):
            if i == 0:
                fh.writelines("%s" %(self._map_snps_header[0][i]))
            else:
                fh.writelines(" %s" %(self._map_snps_header[0][i]))
        fh.writelines("\n")
        # print new header for new assoc file with pvalues from MACH2DAT 
        fh2.writelines(" SNP P\n")
        
        # print the rest
        for line in mach2dat:
            for i in xrange(len(line)):
                if i == 0:
                    fh.writelines("%s"  %(line[i]))
                else:
                    fh.writelines(" %s" %(line[i]))
            fh.writelines("\n")
            fh2.writelines(" %s %s\n" %(line[1], line[9]))

        fh.close()
        fh2.close()
    
    def write_mach2dat_merged_sorted_no_assoc_file(self, file_number=0):
        """ write all snps to merged file, sorted, do not write assoc file """ 
        
        try:
            fh   = file(self._write_file, "w")
        except IOError, e:
            print e
            sys.exit(1)

        mach2dat = []
        for j in xrange(len(self._map_snps)):
             keys = self._map_snps[j].iterkeys()
             for snp in keys:
                 mach2dat.append(self._map_snps[j][snp])
             
        # sort by pvalues
        mach2dat.sort(cmp=compare, key=col_ten) # sort 10th column
  
        # print new header from first mach2dat result file
        for i in xrange(len(self._map_snps_header[0])):
            if i == 0:
                fh.writelines("%s" %(self._map_snps_header[0][i]))
            else:
                fh.writelines(" %s" %(self._map_snps_header[0][i]))
        fh.writelines("\n")
        
        # print the rest
        for line in mach2dat:
            for i in xrange(len(line)):
                if i == 0:
                    fh.writelines("%s"  %(line[i]))
                else:
                    fh.writelines(" %s" %(line[i]))
            fh.writelines("\n")

        fh.close()

    def write_mach2dat_typed_imputed_info(self, typed_snps, file_number=0):
        """ write all snps with imputed or typed info """ 
        
        try:
            fh   = file(self._write_file, "w")
        except IOError, e:
            print e
            sys.exit(1)

        mach2dat = []
        for j in xrange(len(self._map_snps)):
             keys = self._map_snps[j].iterkeys()
             for snp in keys:
                 mach2dat.append(self._map_snps[j][snp])
             
        # sort by pvalues
        mach2dat.sort(cmp=compare, key=col_ten) # sort 10th column
  
        # print new header from first mach2dat result file
        for i in xrange(len(self._map_snps_header[0])):
            if i == 0:
                fh.writelines("%s" %(self._map_snps_header[0][i]))
            else:
                fh.writelines(" %s" %(self._map_snps_header[0][i]))
        fh.writelines("\n")
        
        # print the rest in order how it was mapped
        for line in mach2dat:
            for i in xrange(len(line)):
                if i == 0:
                    fh.writelines("%s"  %(line[i]))
                elif i == 1:
                    if typed_snps.has_key(line[i]):
                      fh.writelines(" %s_typed" %(line[i]))
                    else:
                      fh.writelines(" %s_imputed" %(line[i]))
                else:
                    fh.writelines(" %s" %(line[i]))
            fh.writelines("\n")

        fh.close()
   
    def typed_imputed_info2pos(self, ref_file, pos_write_file):
        """ add chr and pos to typed imputed info """
        
        try:
            fh       = file(self._write_file, "r")
            fh_ref   = file(ref_file, "r")
            fh_pos   = file(pos_write_file, "w")
        except IOError, e:
            print e
            sys.exit(1)
    
        dict_SNPs = {} # store SNPs with line

        comment_pattern = re.compile("^#.*$")
        blankline_pattern = re.compile("^\s*$")

        # typed imputed info file
        # print header
        line = fh.readline().replace("\n", "")
        fh_pos.writelines("CHR POS " + line + "\n")
        line = fh.readline().replace("\n", "")
        while line:

            # skip comment lines that start with "#"
            if(comment_pattern.search(line)):
                line = fh.readline().replace("\n", "")
                continue
            # skip blank lines
            if(blankline_pattern.search(line)):
                line = fh.readline().replace("\n", "")
                continue
          
            list = re.split("\s+",line)
            snp  = list[1].split("_")[0]
            dict_SNPs[snp] = line
            
            line = fh.readline().replace("\n", "")
        
        fh.close()

        # hapmap reference bim file 
        sep = " "
        line = fh_ref.readline().replace("\n", "")
        while line:

            # skip comment lines that start with "#"
            if(comment_pattern.search(line)):
                line = fh_ref.readline().replace("\n", "")
                continue
            # skip blank lines
            if(blankline_pattern.search(line)):
                line = fh_ref.readline().replace("\n", "")
                continue
           
            list = re.split("\s+",line)
            snp  = list[1]
            if dict_SNPs.has_key(snp):
                fh_pos.writelines(list[0] +sep+\
                                  list[3] +sep+\
                                  dict_SNPs[snp] +"\n")
           
            line = fh_ref.readline().replace("\n", "")

        fh_ref.close()
        fh_pos.close()


    def get_rs_pvalue_type_hash(self, file_number=0):
        """ get rs, pvalue, type in hash """ 

        rs_pvalue_type_hash = {}
        mach2dat = []
        for j in xrange(len(self._map_snps)):
             keys = self._map_snps[j].iterkeys()
             for snp in keys:
                 mach2dat.append(self._map_snps[j][snp])
 
        #for snp in mach2dat:
            #rs, type =  mach2dat[snp][1].split("_")
            #rs_pvalue_type_hash[rs] = (mach2dat[snp][9], type)
        for i in xrange(len(mach2dat)):
            rs, type =  mach2dat[i][1].split("_")
            rs_pvalue_type_hash[rs] = (mach2dat[i][9], type)

        return rs_pvalue_type_hash

    def get_best_pvalue(self, file_number=0):
        """ get best_pvalue as Decimal.decimal object """ 

        best_pvalue = 1.0
        mach2dat = []
        for j in xrange(len(self._map_snps)):
             keys = self._map_snps[j].iterkeys()
             for snp in keys:
                 mach2dat.append(self._map_snps[j][snp])
 
        for i in xrange(len(mach2dat)):
            if i == 0:
                best_pvalue = mach2dat[i][9]
            elif mach2dat[i][9] < best_pvalue:
                best_pvalue = mach2dat[i][9]

        return best_pvalue 

    def write_metal_format_from_map_snps_list(self,\
            numof_cases_controls, hapmap_rs_hash, file_number=0):
        """ write metal format """ 
        
        try:
            fh   = file(self._write_file, "w")
        except IOError, e:
            print e
            sys.exit(1)

        mach2dat = []
        for j in xrange(len(self._map_snps_lists)):
             for snp_line in self._map_snps_lists[j]:
                 mach2dat.append(snp_line)
             
        # print new header from first mach2dat result file
        fh.writelines("TRAIT MARKER CHR POS AL1 AL2 FREQ1 MAF SAMPLE_SIZE EFFECT1 OR STDERR LRCHISQ LRPVAL\n")
        
        # print the rest in order how it was mapped
        for tuple in mach2dat:
            
            trait   = tuple[0]
            snp     = tuple[1]
            chr     = str(hapmap_rs_hash[snp][0])
            pos     = str(hapmap_rs_hash[snp][2])
            a1, a2  = tuple[2].split(",")
            freq1   = tuple[11]
            maf     = tuple[12]
            sample_size = str(numof_cases_controls)
            effect1 = tuple[3]
            odds_ratio  = tuple[4]
            stderr  = tuple[5]
            lrchisq = tuple[8]
            lrpval  = tuple[9]
           
            assert(a1 == tuple[10])

            fh.writelines( "%s %s %s %s %s %s %s %s %s %s %s %s %s %s\n" \
                           %(trait, snp, chr, pos, a1, a2, freq1, maf, \
                             sample_size, effect1, odds_ratio, stderr, \
                             lrchisq, lrpval) )

        fh.close()

    def write_mach2dat_typed_imputed_info_combined_without_map(self, disease_name,\
                                                               disease_name_comparison_list,\
                                                               rs_tophits_hash, file_number_disease=0):
        """ write combined outfile without mapping files into memory,
            that means files were not mapped previous to this function call 
            only print combined output for rs_tophits_hash 
            DO NOT MAP snps previously into memory """ 

        # check if valid file_number_disease
        mach2dat_files_numof = len(self._mach2dat_files)
        if mach2dat_files_numof <= file_number_disease:
            return
       
        # ----------------------------------------------------------- # 
        # -- (1) read rs_tophits snps from MACH2dat merged outfile -- #
        # ----------------------------------------------------------- # 
        
        try:
            fh  = file(self._mach2dat_files[file_number_disease], "r")
        except IOError, e:
            print e
            sys.exit(1)
        
        header_pattern = re.compile("^TRAIT.*$")
        # skip lines until header line
        line = fh.readline().replace("\n","")
        while not header_pattern.search(line):
            line = fh.readline().replace("\n","")
        
        # store header
        self._map_snps_header.append(re.split("\s+",line)[:])

        # read all snps
        line = fh.readline()
        map_snps_rs_tophits_disease = []
        while line:

            list = re.split("\s+",line)
            rs   =  list[1].split("_")[0]
            
            if rs_tophits_hash.has_key(rs) and list[3] != "NA":
                map_snps_rs_tophits_disease.append( (list[0], list[1], list[2],\
                                                     list[3],\
                                                     list[4],\
                                                     list[5],\
                                                     list[6],\
                                                     list[7],\
                                                     list[8],\
                                                     list[9],\
                                                     list[10],\
                                                     list[11],\
                                                     list[12],\
                                                     list[13],\
                                                     list[14]) )
            line = fh.readline()
               
        fh.close()

        # ------------------------------------------------------------------ #
        # -- (2) read rs_tophits snps (if available) from MACH2DAT merged -- #
        # -- outfile(s) for comparison                                    -- #
        # ------------------------------------------------------------------ #
 
        map_snps_rs_tophits_compare_list = [] # list of hashes for each
                                              # comparison file

        # for all comparison file(s) (without disease file)
        for i in xrange(mach2dat_files_numof):

            # skip disease file
            if i == file_number_disease:
               continue 

            # try to open disease file for comparison
            try:
                fh = file(self._mach2dat_files[i], "r")
            except IOError, e:
                print e
                sys.exit(1)

            header_pattern = re.compile("^TRAIT.*$")
            # skip lines until header line
            line = fh.readline().replace("\n","")
            while not header_pattern.search(line):
                line = fh.readline().replace("\n","")
            
            # store additional header
            self._map_snps_header.append(re.split("\s+",line)[:])
        
            # read all snps
            line = fh.readline()
            map_snps_rs_tophits_hash = {}
            while line:

                list = re.split("\s+",line)
                rs   =  list[1].split("_")[0]
                
                if rs_tophits_hash.has_key(rs) and list[3] != "NA":
                    map_snps_rs_tophits_hash[rs] = ( list[0], list[1], list[2],\
                                                     list[3],\
                                                     list[4],\
                                                     list[5],\
                                                     list[6],\
                                                     list[7],\
                                                     list[8],\
                                                     list[9],\
                                                     list[10],\
                                                     list[11],\
                                                     list[12],\
                                                     list[13],\
                                                     list[14] )
                line = fh.readline()
                   
            # copy hash 
            map_snps_rs_tophits_compare_list.append(map_snps_rs_tophits_hash.copy())
            fh.close()

        # ---------------------------- #
        # -- print comparison table -- #
        # ---------------------------- #

        try:
            fh   = file(self._write_file, "w")
            fh2   = file(self._write_file_combined_short, "w")
        except IOError, e:
            print e
            sys.exit(1)

        short_columns_disease = [1,2,4,9,10,11,13,14] # count from 0
        short_columns_compare_disease = [4,9,10,11,13,14]
 
        # -- print new header  -- #
        # print columns from disease mach2dat result file
        for i in xrange(len(self._map_snps_header[0])):
            
            # short version
            if i in short_columns_disease:
                if i == short_columns_disease[0]:
                    fh2.writelines("%s_%s" %(self._map_snps_header[0][i],\
                                             disease_name))
                else:
                    fh2.writelines("\t%s_%s" %(self._map_snps_header[0][i],\
                                               disease_name))

            # long version
            if i == 0:
                fh.writelines("%s_%s" %(self._map_snps_header[0][i], disease_name))
            else:
                fh.writelines("\t%s_%s" %(self._map_snps_header[0][i], disease_name))
        # print additional columns from additional comparison file(s)
        # for each comparison file
        for j in xrange(1, len(self._map_snps_header)):
            
            for i in xrange(len(self._map_snps_header[j])):
          
                # short version
                if i in short_columns_compare_disease:
                    fh2.writelines("\t%s_%s" %(self._map_snps_header[j][i],\
                                               disease_name_comparison_list[j-1]))
                
                # long version
                fh.writelines("\t%s_%s" %(self._map_snps_header[j][i],\
                                               disease_name_comparison_list[j-1]))
        fh.writelines("\n")
        fh2.writelines("\n")

        # -- print rs_tophits lines -- #
        for line in map_snps_rs_tophits_disease:

            # rs of line
            rs   =  line[1].split("_")[0]

            # print columns from disease file
            for j in xrange(len(line)):
            
                # short version
                if j in short_columns_disease:
                    if j == short_columns_disease[0]:
                        fh2.writelines(line[j])
                    else:
                        fh2.writelines("\t" + line[j])
                
                # long version
                if j == 0:
                    fh.writelines(line[j])
                else:
                    fh.writelines("\t" + line[j])

            # print columns for each comparison file
            for map_snps_rs_tophits_hash in map_snps_rs_tophits_compare_list:

                # if snp is also in comparison file
                if map_snps_rs_tophits_hash.has_key(rs):
                    assert(len(line) == len(map_snps_rs_tophits_hash[rs]))
                    for j in xrange(len(map_snps_rs_tophits_hash[rs])):
                        
                        # short version
                        if j in short_columns_compare_disease:
                            fh2.writelines("\t%s" %((map_snps_rs_tophits_hash[rs][j])))
                        
                        # long version
                        fh.writelines("\t%s" %((map_snps_rs_tophits_hash[rs][j])))

                # if snp is not in comparison file
                else:
                    for j in xrange(len(line)):

                        # short version
                        if j in short_columns_compare_disease:
                            fh2.writelines("\t---")

                        # long version
                        fh.writelines("\t---")
        
            fh.writelines("\n")
            fh2.writelines("\n")

        fh.close()
        fh2.close()


    def write_mach2dat_typed_imputed_info_combined_without_map_columns_grouped(self, disease_name,\
                                                               disease_name_comparison_list,\
                                                               rs_tophits_hash, file_number_disease=0):
        """ write combined outfile without mapping files into memory,
            that means files were not mapped previous to this function call 
            only print combined output for rs_tophits_hash 
            DO NOT MAP snps previously into memory """ 

        # check if valid file_number_disease
        mach2dat_files_numof = len(self._mach2dat_files)
        if mach2dat_files_numof <= file_number_disease:
            return
       
        # ----------------------------------------------------------- # 
        # -- (1) read rs_tophits snps from MACH2dat merged outfile -- #
        # ----------------------------------------------------------- # 
        
        try:
            fh  = file(self._mach2dat_files[file_number_disease], "r")
        except IOError, e:
            print e
            sys.exit(1)
        
        header_pattern = re.compile("^TRAIT.*$")
        # skip lines until header line
        line = fh.readline().replace("\n","")
        while not header_pattern.search(line):
            line = fh.readline().replace("\n","")
        
        # store header
        self._map_snps_header.append(re.split("\s+",line)[:])

        # read all snps
        line = fh.readline()
        map_snps_rs_tophits_disease = []
        while line:

            list = re.split("\s+",line)
            rs   =  list[1].split("_")[0]
            
            if rs_tophits_hash.has_key(rs) and list[3] != "NA":
                map_snps_rs_tophits_disease.append( (list[0], list[1], list[2],\
                                                     list[3],\
                                                     list[4],\
                                                     list[5],\
                                                     list[6],\
                                                     list[7],\
                                                     list[8],\
                                                     list[9],\
                                                     list[10],\
                                                     list[11],\
                                                     list[12],\
                                                     list[13],\
                                                     list[14]) )
            line = fh.readline()
               
        fh.close()

        # ------------------------------------------------------------------ #
        # -- (2) read rs_tophits snps (if available) from MACH2DAT merged -- #
        # -- outfile(s) for comparison                                    -- #
        # ------------------------------------------------------------------ #
 
        map_snps_rs_tophits_compare_list = [] # list of hashes for each
                                              # comparison file

        # for all comparison file(s) (without disease file)
        for i in xrange(mach2dat_files_numof):

            # skip disease file
            if i == file_number_disease:
               continue 

            # try to open disease file for comparison
            try:
                fh = file(self._mach2dat_files[i], "r")
            except IOError, e:
                print e
                sys.exit(1)

            header_pattern = re.compile("^TRAIT.*$")
            # skip lines until header line
            line = fh.readline().replace("\n","")
            while not header_pattern.search(line):
                line = fh.readline().replace("\n","")
            
            # store additional header
            self._map_snps_header.append(re.split("\s+",line)[:])
        
            # read all snps
            line = fh.readline()
            map_snps_rs_tophits_hash = {}
            while line:

                list = re.split("\s+",line)
                rs   =  list[1].split("_")[0]
                
                if rs_tophits_hash.has_key(rs) and list[3] != "NA":
                    map_snps_rs_tophits_hash[rs] = ( list[0], list[1], list[2],\
                                                     list[3],\
                                                     list[4],\
                                                     list[5],\
                                                     list[6],\
                                                     list[7],\
                                                     list[8],\
                                                     list[9],\
                                                     list[10],\
                                                     list[11],\
                                                     list[12],\
                                                     list[13],\
                                                     list[14] )
                line = fh.readline()
                   
            # copy hash 
            map_snps_rs_tophits_compare_list.append(map_snps_rs_tophits_hash.copy())
            fh.close()

        # ---------------------------- #
        # -- print comparison table -- #
        # ---------------------------- #

        try:
            fh   = file(self._write_file, "w")
            fh2   = file(self._write_file_combined_short, "w")
        except IOError, e:
            print e
            sys.exit(1)

        long_columns = range(len(self._map_snps_header[0]))  # count from 0
        short_columns = [1,2,4,9,10,11,13,14] # count from 0
        short_columns_compare = [4,9,10,11,13,14] # count from 0
 
        for i in long_columns:
    
            # -- print new header  -- #
            # print columns from disease mach2dat result file
            # short version
            if i in short_columns:
                if i == short_columns[0]:
                    fh2.writelines("%s_%s" %(self._map_snps_header[0][i],\
                                             disease_name))
                else:
                    fh2.writelines("\t%s_%s" %(self._map_snps_header[0][i],\
                                               disease_name))
        
            # long version
            if i == 0:
                fh.writelines("%s_%s" %(self._map_snps_header[0][i], disease_name))
            else:
                fh.writelines("\t%s_%s" %(self._map_snps_header[0][i], disease_name))
            
            # print additional columns from additional comparison file(s)
            # for each comparison file
            for j in xrange(1, len(self._map_snps_header)):
                
                # short version
                if i in short_columns_compare:
                    fh2.writelines("\t%s_%s" %(self._map_snps_header[j][i],\
                                               disease_name_comparison_list[j-1]))
                # long version
                fh.writelines("\t%s_%s" %(self._map_snps_header[j][i],\
                                                   disease_name_comparison_list[j-1]))
        fh.writelines("\n")
        fh2.writelines("\n")

        # -- print rs_tophits lines -- #
        for line in map_snps_rs_tophits_disease:

            # rs of line
            rs   =  line[1].split("_")[0]

            for i in long_columns:
        
                # short version
                if i in short_columns:
                    if i == short_columns[0]:
                        fh2.writelines(line[i])
                    else:
                        fh2.writelines("\t" + line[i])
                
                # long version
                if i == 0:
                    fh.writelines(line[i])
                else:
                    fh.writelines("\t" + line[i])
            
                # print columns for each comparison file
                for map_snps_rs_tophits_hash in map_snps_rs_tophits_compare_list:
            
                    # if snp is also in comparison file
                    if map_snps_rs_tophits_hash.has_key(rs):
                        
                        assert(len(line) == len(map_snps_rs_tophits_hash[rs]))
                        # short version
                        if i in short_columns_compare:
                            fh2.writelines("\t%s" %((map_snps_rs_tophits_hash[rs][i])))
                        
                        # long version
                        fh.writelines("\t%s" %((map_snps_rs_tophits_hash[rs][i])))
            
                    # if snp is not in comparison file
                    else:
                        
                        # short version
                        if i in short_columns_compare:
                            fh2.writelines("\t---")
            
                        # long version
                        fh.writelines("\t---")
        
            fh.writelines("\n")
            fh2.writelines("\n")

        fh.close()
        fh2.close()
