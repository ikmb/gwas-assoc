import os
import os.path

from misc_tools import *
from file_tools import *


class job_script:
	"""
		Creates a (PBS) job list name.
		If 'file_name' does not start with 'job.' or '.sh' the name
		will be changed accordingly	
	"""
	
	def __init__(self, file_name, queue="exc_s", script_l_options="select=1:ncpus=1:exc_s=true"):
		
		new_basename = basename(file_name) 
	
		if new_basename[0:3] != "job.":
			new_basename = "job." + new_basename
		
		if new_basename[-3:] != ".sh":
			new_basename += ".sh"
			
		self.script_name = join(dirname(file_name), new_basename)		
		ensure_dir_exists_for_file(self.script_name)
		
		fd = open(self.script_name, "w")
		fd.close()		

		self.command_lines=[]
		self.is_finished = False

		self.command_line_temp = None
		self.command_line_temp_ident = ""

		#
		# The automatically created PBS parameters
		#
		self.script_output_name = self.script_name + ".output"		
		self.script_job_name = basename(file_name)[0:15]
		self.script_queue = queue
		self.script_l_options = script_l_options
								
	def get_script_name(self):
		return self.script_name

	def get_script_output_name(self):
		return self.script_output_name

	def job_name(self, name):
		if len(name) > 15:
			self.script_job_name = name[0:15]
		else:
			self.script_job_name = name

	def ac(self, *command_components):
		"""
			Add components of a long command.
			The command line will be concluded when no further command component
			is handed over to this method. All components added at once are put
			into one line with a tab in front of them.
		"""
		if len (command_components) == 0:
			if self.command_line_temp:
				self.command_lines.append(self.command_line_temp)				
				#self.add(self.command_line_temp)
				self.command_line_temp = None
				self.command_line_temp_ident = ""
				

		else:
			if self.command_line_temp:
				self.command_line_temp += "   \\\n	" + self.command_line_temp_ident
			else:
				# determine the ident				
				first_component = command_components[0]
				len1 = len(first_component)
				len2 = len(first_component.lstrip())
				whitespace_length = len1 - len2
				if whitespace_length > 0:
					self.command_line_temp_ident = first_component[0:whitespace_length]				
				self.command_line_temp = self.command_line_temp_ident
				
			for cc in command_components:
					add = str(cc).strip() + " "
					if add[0] == "-":
						add = " " + add
					self.command_line_temp += add


	def add(self, *command_line_components):
		self.ac()
		temp = ""
		if len(command_line_components) > 0:		
			for clc in command_line_components:
				if temp != "":
					temp +=	" " + clc.strip()
				else:					
					temp = clc
		self.command_lines.append(temp)


	def _check(self):
		self.add("_check")

	def submit(self, depend_afterok=None):
		if not self.is_finished:
			self.finish()
			
		if depend_afterok:
			options = " -W depend=afterok:" + depend_afterok
		else:
			options = ""
				
		
		if call( "qsub --unchanged " + options + " " + self.script_name, shell=True):
			abort("Error while submitting PBS script:\n" + self.script_name)


		
	def finish(self, verbose=True):
		
		self.ac()
		
		if not self.is_finished:
			fd = open(self.script_name, "w")
			
			def p(s=""):			
				print >> fd, s
			
			#
			# print header
			#
			p("#!/bin/bash")
			
			p()
			p("#######################################################################################")
			p("# An automatically created PBS job")
			p("#######################################################################################")
			p()
				
			p("#PBS -N " + self.script_job_name)
			p("#PBS -o " + self.script_output_name)
			p("#PBS -j oe")
			p("#PBS -q " + self.script_queue)
			p("#PBS -l " + self.script_l_options)
			p()
			
			p("#######################################################################################")
			if verbose:
				p("# Checking and logging functions")
			else:
				p("# Checking functions")
			p("#######################################################################################")
			p()
			if verbose:
				p("function _start_log {	echo '##### [JOB-ID]' $PBS_JOBID")			
				p("			echo '### [JOB-HOST]' `cat $PBS_NODEFILE`")			
				p("			echo '## [JOB-START]' `date '+%Y-%m-%d %k:%M:%S'`")				
				p("}")
				p()
				p("function _stop_log { 	if [ -z $1 ]  ")
				p("				then  EXIT_CODE=$?")
				p("				else  EXIT_CODE=$1")
				p("			fi")
				p("			echo '### [JOB-STOP]' `date '+%Y-%m-%d %k:%M:%S'`")
				p("			exit $EXIT_CODE")
				p("}")
				p()
				p("function _check { 	EXIT_CODE=$?")	
				p("			if [ $EXIT_CODE -ne 0 ]  ")
				p("				then  echo ERROR!") 
				p("				      _stop_log $EXIT_CODE")
				p('				else  if [ -n "$*" ]; then echo $*; fi')								
				p("				      echo '## [JOB-CHECK]' `date '+%Y-%m-%d %k:%M:%S'`")
				p("			fi")				
				p("}")
			else:
				p("function _check { EXIT_CODE=$?; if [ $EXIT_CODE -ne 0 ] ; then exit $EXIT_CODE; fi }")
	
			p()
			if verbose:
				p("_start_log")
	
			#
			# print command lines
			#		
			
			p("#######################################################################################")
			p("# THE ACTUAL COMMANDS")
			p("#######################################################################################")
			p()
			
			for cl in self.command_lines:
				p(cl)		
			p()
			
			#
			# print footer
			#
			p("#######################################################################################")
			p("# end of ACTUAL COMMANDS")
			p("#######################################################################################")
			if verbose:
				p("_stop_log")
				
			fd.close()		
			self.is_finished = True



class job_submission_list:
	
	def __init__(self, file_name):
	
		new_basename = basename(file_name) 
	
		if new_basename[0:3] != "job.list.":
			new_basename = "job.list." + new_basename
		
		if new_basename[-4:] != ".txt":
			new_basename += ".txt"
				
		self.file_name = join(dirname(file_name), new_basename)
		ensure_dir_exists_for_file(self.file_name)
			
		# create the job list file		
		#self.file_name = file_name		
		ensure_dir_exists_for_file(self.file_name)
		fd = open(self.file_name, "w")
		fd.close()
				
		# no jobs so far
		self.job_no = 0
		self.job_no_last_depend_all = 1
		
			
	def get_list_name(self):
		return self.file_name
		
		
	def add(self, pbs_script):
		
		if isinstance(pbs_script, job_script):
			pbs_script_name = pbs_script.get_script_name()
			pbs_script.finish()
		else:
			pbs_script_name = pbs_script
		
		fd = open(self.file_name, "a")
		self.job_no += 1
		print >> fd, str(self.job_no) + "	" + pbs_script_name
		fd.close()
		
		
	def add_depend_all(self,pbs_script):
	
		if isinstance(pbs_script, job_script):
			pbs_script_name = pbs_script.get_script_name()
			pbs_script.finish()
		else:
			pbs_script_name = pbs_script
		
		depend_list = ""
		self.job_no += 1

		if self.job_no == 1:
			depend_list = ""
		else:		
		
			for n in range(self.job_no_last_depend_all, self.job_no):
				if depend_list == "":
					depend_list = str(n)
				else:
					depend_list += "," + str(n)
					
		self.job_no_last_depend_all = self.job_no
		
		fd = open(self.file_name, "a")
		
		print >> fd, str(self.job_no) + "	" + pbs_script_name + "	" + depend_list
			
			
	def submit(self, depend_afterok=None):
	
		if not self.job_no_last_depend_all == self.job_no:
			
			dirn = os.path.dirname(self.file_name)
			basn = os.path.basename(self.file_name)
			
			job_script_basename = "job.TRAK_" + basn[9:-3] + "sh" 
			job_script_basename = "TRAK_" + basn[9:-4]  
			js = job_script(join(dirn, job_script_basename))
			js.add("#")
			js.add("# Doing nothing - except being the last (tracking) job of the job list")
			js.add("#")
			js.finish(verbose=False)
			self.add_depend_all(js)
			
	
		options = ""
		if depend_afterok:			
			for s in depend_afterok.split(":"):
		 	 	options += " " + s
		 
		if call( "job_list_submitter " + self.file_name + options, shell=True):
		 	abort("Error while submitting PBS job list:\n" + self.file_name)			
