import os
import sys
from os.path import *

import ConfigParser

from optparse import *
from copy import deepcopy

from operator import itemgetter

from misc_tools import *
from file_tools import *
from job_submission import *



####################################################################
#
#  Main classes and functions
#
####################################################################

ALL_POSSIBLE_SWITCHES = []

class Parameter_Set:
	def __init__(self, dict_of_params):		
		self.__p__ = {}		
		for d in dict_of_params:
			self.__p__[d] = deepcopy(dict_of_params[d])

	#
	# simple getters
	# 
	def __getitem__(self, key):
		self.__getattr__(key)

	
	def __getattr__(self, key):

		#
		# Any special prefix?
		#
		if key[0:4] == "SET_":
			k = key[4:]
			if self.__p__[k] == None:
				return False
			else:
				return True
		else:
			try:
				return self.__p__[key]
			except KeyError:
				return None

	
def get_parameter_sets(parameters,options, description, switch_file):
	parse_possible_switches(switch_file)
	
	aop = all_option_parsing(parameters, options, description, "")
	return_list = []
	for a in aop:
		return_list.append(Parameter_Set(aop[a]))

	return return_list


def parse_possible_switches(file):

	global ALL_POSSIBLE_SWITCHES

	try:
		config = ConfigParser.SafeConfigParser()
		config.read(realpath(file))

		#
		# Build list of all "Possible_Switch"es
		#
		for sec in config.sections():
			dic = {}
			dic.update(config.items(sec))

			p = Possible_Switch(sec, eval(dic['value_as_list']), dic['optparse_type'], dic['optparse_metavar'], dic['optparse_help'], eval(dic['default_value']), eval(dic['optparse_choices']), str(dic['short_version']), eval(dic['shell_variable']), eval(dic['ignore_params_file']))

			ALL_POSSIBLE_SWITCHES.append(p)

	except ConfigParser.InterpolationError:
		abort("Error in configuration file while trying to substitute params.\n" + str(sys.exc_info()[1]))

	except ConfigParser.ParsingError:
		abort("Error in configuration file.\n" + str(sys.exc_info()[1]))

		
	
class Possible_Switch:

	def __init__(self, switch_name, value_as_list, optparse_type, optparse_metavar, optparse_help, default_value, optparse_choices, short_version, shell_variable, ignore_params_file):

		# TODO: checking for params
		self.switch_name = switch_name
		self.value_as_list = value_as_list
		self.optparse_type = optparse_type
		self.optparse_metavar = optparse_metavar
		self.optparse_help = optparse_help
		self.optparse_choices = optparse_choices
		self.short_version = short_version
		self.shell_variable = shell_variable
		self.default_value = default_value
		#if default_value:
		#	print "switch:", switch_name, "has a default value of", default_value
		self.ignore_params_file = ignore_params_file


####################################################################
#
#  Other functions
#
####################################################################


def merge_opts_and_parameter_files(p_file, defaults_from_command_line, defaults_from_shell_variables):
	"""
	PARAMS:
		- parameter file to read		
	RETURNS
		- A dict of dicts merged parameters for all samples
	"""

	# If there is a config file in the home directory, get the "first set of defaults" from there
	try:
		first_defaults = {}

		

		config = ConfigParser.SafeConfigParser(defaults_from_shell_variables)
		config.optionxform = str
		config.read(p_file)	

		

		dfcl = eval(str(defaults_from_command_line))
				
		
		#
		# A "false" from the command line means a --option has not been set!
		# So change the value to a None
		#
		for d in dfcl:
			if type(dfcl[d]) == type(False):
				if dfcl[d] == False:
					dfcl[d] = None
		
		result_dict = {}
		
		def special_sort(x,y):
			try:
				xx = float (x)
				try:
					yy = float (y)					
					if xx < yy:
						return -1
					elif xx > yy:
						return 1
					else:
						return 0
				except ValueError:
					return -1
			except ValueError:
				try:
					yy = float(y)
					return 1
				except ValueError:					
					if x < y:
						return -1
					elif x > y:
						return 1
					else:
						return 0

		sections = config.sections()
		sections.sort(special_sort)

		for sec in sections:

			dic = {}
			dic.update(config.items(sec))

			for key in dfcl:							

				value = dfcl[key]

				if not value == None:
					dic[key.lower()]=value				
					
			result_dict[sec]= dic
	
		return result_dict, config.defaults()


	
	except ConfigParser.InterpolationError:
		abort("Error in configuration file while trying to substitute params.\n" + str(sys.exc_info()[1]))

	except ConfigParser.ParsingError:
		abort("Error in configuration file.\n" + str(sys.exc_info()[1]))

	
def parse_all_parameter_sets(all_pars, req_pars, req_samples, instance_wide_defaults):
	"""
		RETURN:
			A dict of parameter sets. Maybe this sets are modifed by other sets
			(e.g. with "COMMAND  init_set1,mod_set1,mod_set2,... ")
	"""
	#
	# What are requested and what are modifying sets?
	# Create set_modifcations!
	#
	global_modifying_sets = []
	set_modifications = []
	
	for r in req_samples:
		parts = r.split(",")
		head = parts[0]

		if len(parts) == 1:
			set_modifications.append([head])
		else:
			tail = parts[1:]
	
			if head == "":
				global_modifying_sets.extend(tail)
			else:
				set_modifications.append(parts)

	#
	# If there are globally modifying sets, 
	# append them to existens sets...
	#
	if len(global_modifying_sets) > 0:
		for sm in set_modifications:
			sm.extend(global_modifying_sets)
	
	#
	# Check if all requested and modifying sets exists
	#
	unfound_requested_sets = []
	unfound_modifying_sets = []
	for sm in set_modifications:
		if not sm[0] in all_pars:
			if not sm[0] in unfound_requested_sets:
				unfound_requested_sets.append(sm[0])
			
		if len(sm) > 1:
			for mod_set in sm[1:]:
				if not mod_set in all_pars:
					if not mod_set in unfound_modifying_sets:					
						unfound_modifying_sets.append(mod_set)
			
	abort_string = ""
	
	if len(unfound_requested_sets):
		abort_string = "\nUnable to find these requested sets: " + str(unfound_requested_sets)
		
	if len(unfound_modifying_sets):
		abort_string += "\nUnable to find these modifying sets: " + str(unfound_modifying_sets)
		
	if not abort_string == "":
		abort(abort_string)
	

	#
	# Now, go through the requested sets
	#

	# The resultig dict
	result={}

	for sm in set_modifications:
		
		# head and tail of each set_modification
		head = sm[0]
		tail = []
		if len(sm) > 1:
			tail = sm[1:]
		
				
		#
		# the name of the resulting parameter set
		#
		result_set_name = head
		for par_set_name in tail:		
				result_set_name += "," + par_set_name
		
		#
		# Start with the head
		#			
		result[result_set_name] = deepcopy(all_pars[head])
		
		#
		# Now go through the modifying sets of the tail
		#
		for par_set_name in tail:
			
			par_set = all_pars[par_set_name]
			
			for par_set_key in par_set:
				
				if par_set_key == "force" or par_set_key in instance_wide_defaults:
					pass
				else:
					par_set_value = par_set[par_set_key]
					result[result_set_name][par_set_key] = par_set_value				
	return result	


	

def all_option_parsing(par_parameters,par_options, p_description, p_default_param_file ):

	"""
		PARAMS:
			 par_parameters:	list of required list of params
			 par_options:	list of optional parameters
			 
			p_description:	description of the calling script for optparse
			p_default_param_file:	name of the default parameter file to look for			
		RETURN:
			TODO
	"""

	p_parameters = []
	for pp in par_parameters:
		if not pp in p_parameters:
			p_parameters.append(pp)
	
	p_options = []
	for op in par_options:
		if not op in p_options:
			p_options.append(op)
	
	parameters_and_options = p_parameters + p_options	

	
	#
	# prepare and start the parser
	#
	script_name = basename(sys.argv[0])

	#
	# A dict with all possible switches
	#
	aps = {}
	for s in ALL_POSSIBLE_SWITCHES:
		aps[s.switch_name] = s

	default_description = "The parameter sets will be read from PARAM_SRC and can also be set/overridden using command parameters and options. PARAM_SRC can be a file or directory, where a file named '" + p_default_param_file + "' is expected."
	parser = OptionParser(script_name + "  [options]  PARAM_SRC  [ set1, set2, set3 ... ]", description=p_description + "\t " + default_description, formatter=TitledHelpFormatter())

	def add_switch_to_option_group(s, g):

		help_temp = s.optparse_help		

		#if not help_temp[-1] in ".:,;()":
		#	help_temp += "."

		if s.optparse_choices:
			help_temp += " CHOICES: " + str(s.optparse_choices[0])
			for c in s.optparse_choices[1:]:
				help_temp += "," + str(c)
			
		if s.default_value:
			help_temp = help_temp + " \t\tDEFAULT: " + str(s.default_value)

		

		if s.optparse_type == "boolean":
			g.add_option(s.short_version, "--" + s.switch_name, action="store_true", default=False, help=help_temp)
		elif s.optparse_type == "explict_boolean":			
			g.add_option(s.short_version, "--" + s.switch_name, metavar="[Yes|No]", help=help_temp, choices=["Yes","No"])
		else:
			#print "switch_name", s.switch_name
			#g.add_option(s.short_version, "--" + s.switch_name, type=s.optparse_type, metavar=s.optparse_metavar, help=help_temp, choices=s.optparse_choices, default=s.default_value)
			g.add_option(s.short_version, "--" + s.switch_name, type=s.optparse_type, metavar=s.optparse_metavar, help=help_temp, choices=s.optparse_choices)


	#
	# Are there switches requested as parameters?
	#
	if len(p_parameters) > 0:

		# sort & unqify req. params

		p_parameters = sorted(list(set(p_parameters)))

		group_params = OptionGroup(parser, "PARAMETERS (required)", "The required parameters are usually they are defined by the parameter file specified by PARAM_SRC, but they can be overridden using options on the command line:")
		group_params_any_switches = False

		group_shell_params = OptionGroup(parser, "Shell PARAMETERS (required)", "These params are read from shell variables - if existing. The are overridden by params read from parameter file, which can themselves can be overruled using them as command line options.")
		group_shell_params_any_switches = False

		for xxx in p_parameters:
			po = aps[xxx]

			if po.shell_variable:
				add_switch_to_option_group(po, group_shell_params)
				group_shell_params_any_switches = True
			else:
				add_switch_to_option_group(po, group_params)
				group_params_any_switches = True

		if group_params_any_switches:
			parser.add_option_group(group_params)

		if group_shell_params_any_switches:
			parser.add_option_group(group_shell_params)


	#
	# Are there switches requested as options?
	#
	parser.add_option("--force",  action="store_true", help="In case of no given args forces to start with the standard options without asking.")
	if len(p_options) > 0:

		group_shell_options = OptionGroup(parser, "Shell OPTIONS" ,"Options to be read from shell variables.")
		group_shell_options_any_switches = False

		for xxx in p_options:
			po = aps[xxx]

			if po.shell_variable:
				add_switch_to_option_group(po, group_shell_options)
				group_shell_options_any_switches = True
			else:
				#add_switch_to_option_group(po, group_options)
				add_switch_to_option_group(po, parser)
		
		if group_shell_options_any_switches:
			parser.add_option_group(group_shell_options)

	#
	# Now the magic moment: Parse the args!
	#	
	opts, args = parser.parse_args()
	
	#
	# Is any option set that prevents reading the config file with param sets?
	# (Ist ein bisschen durch die Brust ins Auge...)
	#
	cmd_line_dict = eval(str(opts))
	
	for pao in parameters_and_options:
			
		if aps[pao].ignore_params_file:
			# Ok, this a switch that could prevent reading the file.
			# Is it set?

			if cmd_line_dict[pao]:				

				#
				# Build a result list just with the args and without reading
				# from a parameters sets files
				#
				d = {}
				for pao in parameters_and_options:
					d[pao] = cmd_line_dict[pao]
					

				cmd_line_params = {"DUMMY" : d}

				return cmd_line_params

	#
	# Really start with all standard args (if nothing is specified)?
	#

	# Check for args beginning with ,
	modifying_sets_as_args = False
	requested_sets_as_args = False
	for a in args:
		if a.startswith (","):
			modifying_sets_as_args = True
		else:
			requested_sets_as_args =True

	if opts.force and requested_sets_as_args:
		abort("Either specify sets of params or use the '--force' option to choose all 'standard' sets of parameters.")

	if not requested_sets_as_args and not opts.force:
		print
		print "No sets of parameters are specified."
		if modifying_sets_as_args:
			print "(Modifying parameters sets are given.)"
		print
		s = raw_input("Do you really want to start with all 'standard' sets of parameters? [y/n]:")

		if s.lower() not in ("y", "yes"):
			abort("For help start with option '-h")

	#
	# Define the TARGET parameter file name
	#

	# Without any args, look in the local directory for a config file
	# If the first arg is nothing readable (like a sample number ;-)
	# try also the local directory
	if len(args) == 0 or not exists(args[0]):
		args.insert(0, ".")

	parameter_file=""	

	if isdir(args[0]):
		temp = join(args[0], p_default_param_file)
		
		if isfile(temp):
			parameter_file=temp					
		else:
			parser.error("No default parameter file '" + p_default_param_file + "' found in " + args[0])
	else:
		if isfile(args[0]):
			parameter_file=args[0]			
		else:
			parser.error("No default parameter named " + args[0] + "found.")


	#
	# If necessary try to read defaults from shell variables
	#
	defaults_from_shell_variables = {}

	for x in parameters_and_options :
		pao = aps[x]		
		if pao.shell_variable:
			try:
				defaults_from_shell_variables[x.lower()] = os.environ[x]
			except KeyError:
				0 

	#
	# Read all params from config file and merge/overide them with the comman line options
	#
	all_params, instance_wide_defaults = merge_opts_and_parameter_files(parameter_file, opts, defaults_from_shell_variables)

	req_samples = []
	
	try:
		force_sets = instance_wide_defaults["opp.force-sets"].split(",")
	except KeyError:		
		force_sets = []

		

	requested_sets_present = False

	
	if requested_sets_as_args:
		for i in args[1:]:
			req_samples.append(i)
			requested_sets_present = True	
	else:
		for ap in all_params:
			
			try:
				to_be_ignored = eval(all_params[ap]["ignore_this_sample"])			
			except KeyError:
				to_be_ignored = True
			
			
			if ap in force_sets:
				to_be_ignored = False
			
			if not to_be_ignored:
				req_samples.append(ap)
				requested_sets_present = True	
			
		for i in args[1:]:
			if i.startswith(","):
				req_samples.append(i)

	if not requested_sets_present:
		if len(req_samples) == 0:
			abort("No sets specified (either in file or command line)")
		else:
			abort("Only modifying sets specified (either in file or command line)")
			
	parameter_sets = parse_all_parameter_sets(all_params, p_parameters, req_samples, instance_wide_defaults)
	

	#
	# Final checks for all sets of parameters
	#
	abort_missing_string =""
	abort_choices_string = ""
	
	for ps_name in parameter_sets:


		ps = parameter_sets[ps_name]
		
		#
		# 1. Are all requested parameters present?
		#
		for rp in p_parameters:

			# Is there a not (yet) set parameter?
			if rp not in ps:
				
				# If the missing parameter is a possible shell variable ...				
				if aps[rp].shell_variable:
					# ... get it from the environment
					ps[rp] = os.environ[rp]
	
	
				# Is the parameter still not set?
				if rp not in ps:
					# The use the switch default ...
					if aps[rp].default_value:
						ps[rp] = aps[rp].default_value
					else:
						# ... or write a "missing note"
						abort_missing_string += "	" + ps_name + ":	'" + rp + "'\n"
	
		#
		# 2. Do the params comply to any given choices?
		#
		
		for pao in parameters_and_options :
			if aps[pao].optparse_choices:
				if not str(ps[pao]) in aps[pao].optparse_choices:
					abort_choices_string += "	" + ps_name + ":	'" + pao + "' is '" + str(ps[pao]) + " 'instead one of: " + str(aps[pao].optparse_choices) + "\n"
	
	
	
		#
		# 3. Translate explicit boolean in True/False
		#
	
		for pao in parameters_and_options :
			if aps[pao].optparse_type == "explict_boolean":
				if ps[pao] == "Yes":
					ps[pao] = True
				else:
					ps[pao] = False
	
	
		#
		# 4. In case of boolean transfer a non-boolean value into False
		#
		for pao in parameters_and_options :
			if aps[pao].optparse_type == "boolean":
				
				try:
					if type(ps[pao]) != type(True):				
						temp = str(ps[pao]).strip()
						if temp.lower() == "true" or temp == "":
							ps[pao] = True
						elif temp.lower() == "false":
							ps[pao] = False
						else:
							abort_choices_string += "Option --" + pao + " in parameter file is set to: '" + temp + "'\n	(Must either be 'False' or 'True' (or empty, which equals 'True'))"  
				except KeyError: # an optional boolean switch, that is not set
					ps[pao] = False				
				
	
		#
		# 5. If requested split values into lists
		#
		for pao in parameters_and_options:
			if aps[pao].value_as_list:
	
				try:
					actual_type = type(ps[pao])
					if not actual_type == type([]):
						ps[pao] = ps[pao].split(",")
				except:
					pass
	
	
	abort_string = ""
	if abort_missing_string != "":
		abort_string += "Missing parameters ...\n" + abort_missing_string
	
	if abort_choices_string != "":
		abort_string += "Parameter not one of the possible choices ...\n" + abort_choices_string
	
	if abort_string != "":
		abort(abort_string)
	
	return parameter_sets
