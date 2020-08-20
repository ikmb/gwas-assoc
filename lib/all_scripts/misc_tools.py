import string
from subprocess import call

import csv
import sys
from os.path import *


####################################################################
#
#  COMMON FUNCTIONS
#
####################################################################

def read_table_from_csv(file_name, header_starts=["##hdr", "## hdr"], delimiter="\t"):
		
	file = open(file_name ) 
	reader = csv.reader(file, delimiter=delimiter)
						
	result_list = []
	column_names = None
												
	for row in reader:
		if row != []:
			
			col1 = row[0].strip()
			if col1 in header_starts:		
				# a header line
				column_names = []
				for cn in row[1:]:
					temp = cn.strip()
					if temp == "":
						abort("Found an 'empty' column name in the header")
					column_names.append(temp)
												
			elif col1[0] != "#":
				
				if column_names == None:
					abort("Found a data line without any preceeding header line.")
				
				if len(row) > len(column_names):
					abort("Data line contains more entries than header line column names.")
				
				# a data line
				temp_row = []
				for r in row:
					temp_row.append(r.strip())
					
				result_list.append( dict (zip ( column_names, temp_row) ) )				
			else:
				pass # ignoring a comment 
	
	file.close()
	return result_list



def all_members(aClass):
    members = {}
    bases = list(aClass.__bases__)
    bases.reverse()
    for base in bases:
        members.update(all_members(base))
    members.update(vars(aClass))
    return members


__item_line_string__=""
def item_line(s="__item_line_none_value__"):
	global __item_line_string__
	
	if not s == "__item_line_none_value__":
		__item_line_string__ +=str(s).strip() + "	"
	else:		
		if not __item_line_string__ == "":
			print __item_line_string__
			__item_line_string__ = ""

def print_underlined(s, c="-", d=1):
	print " " * d + s
	print c * (len(s) + d*2)


def print_header(s):

	ss = str(s).split("\n")
	textlen = 0
	for s in ss:
		textlen = max(textlen, len(s.strip()))
	
	print
	print
	print "#" * (textlen + 8)
	print "#" + " " * (textlen+6) + "#"
	for s in ss:
		print "#   " + s + " " * (textlen - len(s) + 1) + "  #"
	print "#" + " " * (textlen+6) + "#"
	print "#" * (textlen + 8)
	print


def print_sub_header(s):
	print
	text = "   " + str(s).strip() + "   "
	textlen = len(text)
	print "#" * (textlen + 2)
	print "#" + text 
	print "#" * (textlen + 2)
	print


def abort(s = "Undefinied reason !", where=None):

	if where:
		print
		print "Aborting script in Class:", str(where.__class__)
		print s
		print
	else:
		print
		print "Aborting script:"
		print s
		print
	sys.exit (1)


####################################################################
#
#   Classe for easy shell execution
#
####################################################################


class Command:

	def __init__(self, cmd):
		self.__cmd = str(cmd)
		self.__args = []


	def add(self, arg):
		self.adds("",arg)


	def adds(self, switch, arg):

		a = string.strip(str(arg))
		if " " in a:
			a = '"' +a + '"'
		if self.__args == None:
			self.__args = [switch, a]
		else:
			self.__args.append( str (switch) + " " + a)

	def addj(self, *parts):
		self.add( join(*parts))


	def set_cmd(self,cmd):
		self.__cmd = str(cmd)


	def get(self):
		return str(self.__cmd) + "  " + string.join(join (self.__args), "  ")


	def reset(self):
		self.__cmd = None
		self.__args = None


	def run(self, ignore_return = False):
		if self.__cmd == None:
			abort("Shell Command must contain at least a command to run.")

		ret = call(self.get(), shell=True)

		if ret and not ignore_return:
			abort ("Error while running following shell command:\n\n" + self.get())

		return ret
