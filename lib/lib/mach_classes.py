import os
import os.path
import sys
import re
import gzip

class Dat:
  """ class Dat implements tasks on Dat files from MACH software
      pedigree file organization taken from QTDT format:
      http://www.sph.umich.edu/csg/abecasis/qtdt/docs/data.html  """
  
  def __init__(self, write_file=None, affection_status_switch=True):
    """ init """
    self._affection_status           = 'A' 
    self._affection_status_switch = affection_status_switch
    self._write_file = write_file    # name of output file

  def write_file_set(self, write_file):
    """ set name of output file """
    self._write_file = write_file 

  def affection_status_switch_on(self):
    """ switch on affection status """
    self._affection_status_switch = True
   
  def affection_status_switch_on(self):
    """ switch off affection status """
    self._affection_status_switch = False 
 
  def write_file(self):
    """ write dat information to write_file """
    if self._write_file == None:
      return

    try:
      out = file(self._write_file, "w")
    except IOError, e:
      print e
      sys.exit(1)
    out.writelines("A cases") 
    out.close()
