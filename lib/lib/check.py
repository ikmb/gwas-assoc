import os
import os.path

class Check:
  """ class Check implements several check for path and file """

  def path_exists(self, path):
    """ check if path exists """
    return os.path.exists(path) 

  def file_exists(self, file):
    """ check if file exists """
    return os.access(file) 

