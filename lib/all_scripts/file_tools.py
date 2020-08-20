import os
import os.path

from misc_tools import *

def ensure_dir_exists_for_file(path_to_file):
	"""
		checks if dir of a file (full path)
		exists and if not creates it                
	"""
	ap = os.path.dirname(os.path.abspath(os.path.expanduser(path_to_file)))

	if os.path.exists(ap):                    
		if not os.path.isdir(ap):
			abort(ap + " exists but is not a directory.")
	else:
		os.makedirs(ap)  


def ensure_dir_exists(path_to_dir):
	"""
		checks if a dir exists and if not creates it		
	"""
	
	ap = os.path.abspath(os.path.expanduser(path_to_dir))	
	if os.path.exists(ap):
		
		if not os.path.isdir(ap):
			abort(ap + " exists but is not a directory.")
	else:
		os.makedirs(ap)		
		
	
	

def read_sets_of_lines_from_file(file_name, line_ranges, converter = None, verbose=True):
	"""
		Opens and reads a file and extracts sets of lines specified by a list of start/stop pairs.

		file_name:	Name of the file to be read from
		line_ranges:	A list of pairs of start and stop lines. These are integers starting at 0 for the first line.
		converter:	A function (e. g. 'int()') which will be applied to each line before being put into the resulting list

		Returns a list containing lists of values. The order is accordant to the order of the line_ranges

		ASSUMPTION: The regions do must not overlap
	"""


	#
	# Create a local, sorted version of the line ranges, that
	# preserves the information about the order
	#
	temp_line_ranges =[]
	line_range_no = 0
	for lr in line_ranges:
		temp_line_ranges.append( [lr[0], lr[1], line_range_no])
		
		line_range_no += 1
	temp_line_ranges.sort()


	#
	# Now, go through the file
	#
	result_list = [None] * len(line_ranges)
	line_range_iter = iter(temp_line_ranges)

	act_line_range = line_range_iter.next()
	act_start = act_line_range[0]
	act_stop = act_line_range[1]
	act_result_range = []

	in_range = False
	line_no = 0
	file = open(file_name)

	if verbose:
		print
		print "Reading", len(line_ranges), "region(s) from this file:"
		print file_name
		print

	for line in file:

		if in_range:
			#
			# Do the reading
			# (if the last pos has been reached
			#
			if converter:
				act_result_range.append(converter(line))
				#print converter(line)
			else:
				act_result_range.append(line)


			if line_no ==act_stop:
				#
				# Stoping the reading
				#
				if verbose:
					print "	Stop reading at ", line_no
					print
				in_range = False
				result_list[act_line_range[2]] = act_result_range

				# Prepare for the next range
				try:
					act_line_range = line_range_iter.next()
				except StopIteration:
					break
				act_start = act_line_range[0]
				act_stop = act_line_range[1]
				act_result_range = []


		else:
			if line_no ==act_start:
				#	
				# Starting the reading
				#
				if verbose:
					print "	Start reading at", line_no
				in_range = True

				if converter:
					act_result_range.append(converter(line))
				else:
					act_result_range.append(line)

		line_no+=1


	#
	# It's done
	#
	file.close()
	return result_list


