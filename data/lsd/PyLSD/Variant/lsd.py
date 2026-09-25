"""driver for LSD with Variable Molecular Formula (VMF) and Variable Status Atoms (VSA)"""

import sys
import os
import defaults

class LSD:
	"""All about LSD"""
	def __init__(self, default, fn):
		"""Creation of an LSD instance from an LSD input file of from standard input"""
# remember defaults/os-dependent file-related constants
		self.defaults = default
# file handle to LSD input data
		fh = sys.stdin if fn == None else open(fn)
# remember where input data come from, for the solution listing
		self.source = "standard input" if fh == sys.stdin else "file: " + fn
# save all lines for futur use at the beginning of the solution listing
		self.originals = fh.readlines()
# done with the input file
		fh.close()
# import linesets module
		import linesets
# generate LSD data sets without VMF ambiguities
		self.linesets = linesets.LineSets(self.originals)
# generate LSD data sets without ambiguities
		self.linesets.run()
# prepare test for success of the generation test
		self.notEmpty = len(self.linesets.alllsddata) != 0
# print experimental chemical shift data files, if needed
		self.linesets.cs.write_shift_files(self.defaults)

	def generate_files(self):
		"""Create input files for LSD, using rootname as root name."""
# common part of all lsd input and output file names
		rootname = self.defaults.rootname
# folder in which all lsd input and output files will be created
		datafolder = self.defaults.datafolder
# list of LSD data sets
		alllsddata = self.linesets.alllsddata
# number of LSD input file to create
		numfiles = len(alllsddata)
# how many digital places for the number of LSD input files?
		numdigits = len(str(numfiles))
# file index format with possible left-padding with zeros
		numformat = "%0" + str(numdigits) + "d"
# initial list of LSD input file paths
		self.inputpaths = []
# initial list of LSD output file paths
		self.outputpaths = []
# scan through input file indexes and correspond sets of MULT lines
		for fileindex, lsddata in zip(range(1, numfiles+1), alllsddata):
# first part of current LSD input file name
			filenamehead = rootname + numformat % fileindex
# assembles LSD input file name
			inputname = filenamehead + self.defaults.lsd_extension
# assembles LSD input file path
			inputpath = os.path.join(datafolder, inputname)
# open new LSD input file name
			fh = open(inputpath, "w+")
# write the single string of all the initials lines, once replacements were carried out
			fh.write(lsddata)
# close the newly created LSD input file
			fh.close()
# update the list of the input file paths
			self.inputpaths.append(inputpath)
# create corresponding LSD output file name
			outputname = filenamehead + ".sol"
# create corresponding LSD output file path
			outputpath = os.path.join(datafolder, outputname)
# update the list of the output file paths
			self.outputpaths.append(outputpath)
# remember the name of the single final output file with all the solutions, no extension
		self.resultpath = os.path.join(datafolder, rootname + "0")
# final solution output file, .sol extension
		self.solpath = self.resultpath + ".sol"

	def run_files(self):
		"""run all the LSD files that were created by generate_files()"""
# remember directory with LSD binaries (lsd, outlsd, genpos)
		self.lsdbin = self.defaults.lsdbin
# remember path to the LSD executable
		self.lsdpath = os.path.normpath(self.lsdbin + "/lsd ")
# get null device for stderr/stdout LSD catching
		nulldev = self.defaults.nulldev
# redirect part of the LSD command
		redirect = " >" + nulldev + " 2>&1"
# initial list of solution numbers
		self.localnsols = [0] * len(self.inputpaths)
# solution numbers are given by solncounter once an LSD run is finished
		self.solncounter = os.path.normpath(self.linesets.getSolutionCounter(self.defaults.solncounter))
# if stopfile exists, stop everything
		self.stopfile = os.path.normpath(self.linesets.getStopFile(self.defaults.stopfile))
# not stopped yet
		self.stopped = False
# scan through LSD input file names
		for isol, inputpath in enumerate(self.inputpaths):
# to be stopped?
			if os.path.exists(self.stopfile):
# remember answer
				self.stopped = True
# exit from loop
				break
# assembles lsd command
			command = self.lsdpath + inputpath + redirect
			print (command)
# run LSD (all that business for that!!)
			code = os.system(command)
# get LSD exit code 
			code = self.syscode(code)
# -1 means error in input file and 0 if no solution. get number of solutions
			if code > 0:
# open solution counter
				fh = open(self.solncounter)
# get the number of solutions for the current LSD run
				localnsol = int(fh.readlines()[0].rstrip("\n"))
# done with the solution counter
				fh.close()
# in case of error or no solution
			else:
# no solution
				localnsol = 0
			print (localnsol, "solution" + ("s" if localnsol > 1 else ""))
# update list of number of solutions
			self.localnsols[isol] = localnsol
# one all LSD input files has been run, remove the solution counter file
		if os.path.isfile(self.solncounter):
			os.unlink(self.solncounter)

	def syscode(self, code):
		"""decode value from os.system() for Unix platforms"""
		if not self.defaults.isWin:
# get most significant byte
			code = code / 256
# get code if negative (255 is -1 for true)
			if code > 127:
				code -= 256
		return code

	def collect_results(self):
		"""collect the solutions into a single file, same rootname but 0 as index"""
# solutions so far
		self.nsol = 0
# open the global solution file
		self.resulth = open(self.solpath, "w+")
# header line, where do data come from?
		self.resulth.write("# From " + self.source + "\n")
# copy of initial LSD input file (the one with VSAs)
		for sourceline in self.originals:
			self.resulth.write(sourceline)
# end of preamble
		self.resulth.write("#\n")
# scan through LSD output file paths and corresponding number of solution
		for outputpath, localnsol in zip(self.outputpaths, self.localnsols):
# no time to loose with with empty LSD solution files
			if localnsol:
# report from current non empty LSD solution file
				self.collect_solfile(outputpath, localnsol)
# clean end of global solution file if at least one solution was found
		if self.nsol != 0:
# end of file marker is a 0
			self.resulth.write("0\n")
# close the global solution file
		self.resulth.close()

	def collect_solfile(self, outputpath, localnsol):
		"""transfer a solution file named outputpath, having localnsol solutions, to the global solution file"""
		print ("Collect solutions for %s (%d)" % (outputpath, localnsol))
# for sure, there will be one more solution
		self.nsol += 1
# write the magic word OUTLSD before the very first global solution
		if self.nsol == 1:
			self.resulth.write("OUTLSD\n")
# open for reading in the current LSD output file
		fh = open(outputpath)
# all the lines up to OUTLSD are part of the header. useless.
		line = ""
		while line != "OUTLSD\n":
			line = fh.readline()
# new index of the first solution
		iinit = self.nsol
# new index of the last solution, +1
		ifinal = iinit + localnsol
# scan through LSD output file for the current solution
		for isol in range(iinit, ifinal):
# solution header line: number of atoms, old solution index (useless)
			line = fh.readline()
# get new solution header line
			l = line.rstrip("\n")
			nat = l.split()[0]
			newline = nat + " " + str(isol) + "\n"
# write new solution header line in global solution file
			self.resulth.write(newline)
# scan through atom lines in LSD output file, copy them to global output file
			for iat in range(int(nat)):
				line = fh.readline()
				self.resulth.write(line)
# done with the current LSD outpul file
		fh.close()
# update the number of solutions (remember ifinal is excluded)
		self.nsol = ifinal - 1

	def postproc(self):
		"""ranking and graphical display of solutions"""
		print()
		if self.nsol == 0:
			print ("No solution. Handling of Variable Status Atoms was OK.")
			return
# flag is true if ranking was achieved
		self.rankedsolutions = False
# is a solution ranking possible?

		# if self.linesets.cs.filenames:
		# 	self.rank()
# path to solution file, either ranked or not
		solpath = self.rankedsolpath if self.rankedsolutions else self.solpath
# path to the OUTLSD executable
		outlsdpath = os.path.normpath(self.lsdbin + "/outlsd")
# final coordinate file, .coo extension
		coopath = self.resultpath + ".coo"
# outlsd command, for coo file
		outlsdcommandcoo = outlsdpath + " 6 < " + solpath + " > " + coopath
		print ("Running:",  outlsdcommandcoo)
		code = os.system(outlsdcommandcoo)
		code = self.syscode(code)
		if code != 0:
			print( "Failed. (%s)" % code)
			return
# final coordinate file, .sdf extension
		# sdfpath = self.resultpath + ".sdf"
		sdfpath = self.resultpath + ".smiles"
# outlsd command, for SD file
		# outlsdcommandsdf = outlsdpath + " 7 < " + solpath + " > " + sdfpath
		outlsdcommandsdf = outlsdpath + " 5 < " + solpath + " > " + sdfpath
		print ("Running:",  outlsdcommandsdf)
		code = os.system(outlsdcommandsdf)
		code = self.syscode(code)
		if code != 0:
			print ("Failed. (%s)" % code)
			return
# path to the GENPOS executable
		genpospath = os.path.normpath(self.lsdbin + "/genpos")
# # final postscript file, .ps extension
# 		pspath = self.resultpath + ".ps"
# # genpos command
# 		genposcommand = genpospath + " < " + coopath + " > " + pspath
# 		print ("Running:", genposcommand)
# 		os.system(genposcommand)
# 		code = self.syscode(code)
# 		if code != 0:
# 			print ("Failed. (%s)" % code)
# 			return
# # path to the postscript viewer
# 		viewerpath = self.defaults.viewerpath
# # viewer command
# 		viewercommand = viewerpath + ' ' + pspath
# 		print ("Running:",  viewercommand)
# 		os.system(viewercommand)
# message about solution ranking done or not
		print ("With" if self.rankedsolutions else "No", "solution ranking")

	def rank(self):
		print ("Ranking is possible", )
# for which element? Only 1 element for the moment, sorry
		pilotfilenames = self.linesets.cs.filenames
		print(pilotfilenames)
		element = list(pilotfilenames.keys())[0]
		pilotfilename = pilotfilenames[element]
		print ("for", element, "as piloted by", pilotfilename)
# path to the OUTLSD executable
		outlsdpath = os.path.join(self.lsdbin, "outlsd")
# directory in which predictors live
		predspath = self.defaults.predictorspath
# path to the SD file that will be used for chemical shift prediction
		sdfpath = os.path.join(predspath, self.defaults.rootname + '0' + ".sdf")
# get the SD file by means of genpos. The coordinates are useless -> new option for outlsd: sdf without coordinates
		outlsdcommandsdf = outlsdpath + " 10 < " + self.solpath + " > " + sdfpath
		print ("Running:",  outlsdcommandsdf)
# run outlsd, catch return code.
		code = os.system(outlsdcommandsdf)
		code = self.syscode(code)
# check return code
		if code != 0:
			print ("Failed. (%s)" %code)
			return
# build prediction command
# command name
		# command = self.defaults.predictorsbat[element]
		command = self.defaults.predictors[element]
# rootname for predictor
		predroot = os.path.join(predspath, self.defaults.rootname)
# prediction command
		predcommand = command + " " + predroot
		print ("predict with:", predcommand)
# run outlsd, catch return code. Creates a rootname_D.sdf and a rootname_R.txt file in predspath
		code = os.system(predcommand)
		code = self.syscode(code)
# check return code
		if code != 0:
			print ("Failed. (%s)" % code)
			return
# get shift order as a list of solution indexes
		shiftorder = []
# ranking result filename
		predresultname = predroot + "R.txt"
		fhresult = open(predresultname)
# read first three lines (about SD, LSD numbers and experimental chemical shift
		line = fhresult.readline()
		line = fhresult.readline()
		line = fhresult.readline()
# read useful lines
		while 1:
			line = fhresult.readline()
			if not line:
				break
			pos = line.find(' ')
			rank = int(line[0:pos])
			shiftorder.append(rank)
		fhresult.close()
# open solution file for ranking
		solfh = open(self.solpath)
		self.rankedsolpath = os.path.join(self.defaults.datafolder, self.defaults.rootname + "Ranked.sol")
		rankedsolfh = open(self.rankedsolpath, "w+")
# copy first line of header
		line = solfh.readline()
		rankedsolfh.write(line)
# copy following lines up to the # line
		while 1:
			line = solfh.readline()
			rankedsolfh.write(line)
			if len(line) and line[0] == "#":
				break
# copy the outlsd line
		line = solfh.readline()
		rankedsolfh.write(line)
# solutions in an array of strings
		solutions = []
# loop over solutions in solfh
		while 1:
# get solution header line
			line = solfh.readline()
# solution header line contains 0 if the end of the .sol file is reached
			if not line or line[0] == '0':
				break
# get number of atoms, equal to the number of atom lines to read
			natoms = int(line[0:line.find(' ')])
# one solution as a string, starts with the solution header line
			onesol = line
# collect solution lines
			for i in range(natoms):
				line = solfh.readline()
				onesol += line
# store one solution
			solutions.append(onesol)
# finished with the solution file
		solfh.close()
# initialize solution counter for ranked solutions
		isol = 1
# scan solution ranks
		for rank in shiftorder:
# get solution from rank
			onesol = solutions[rank]
# solution as array of lines
			lines = onesol.splitlines(True)
# get first line
			line = lines[0]
# get pieces of first line
			pieces = line.split()
# replace the old solution index by the new one
			pieces[1] = str(isol)
# get modified first line
			line = ' '.join(pieces) + "\n"
# put modified first line to its place
			lines[0] = line
# make modified solution string
			onesol = "".join(lines)
# write the best solutions first
			rankedsolfh.write(onesol)
# one more solution
			isol += 1
# end of .sol file
		rankedsolfh.write("0\n")
# close file with ranked solutions
		rankedsolfh.close()
# flag that indicates if ranking was achieved
		self.rankedsolutions = True
				
if __name__ == "__main__":
	import eraser
# get PyLSD default values
	default = defaults.Defaults()
# count input file arguments. Should be 0 or 1
	argc = len(sys.argv) - 1
# determine PyLSD input file name and rootname for PyLSD-generated files
	if argc:
# get new rootname from command line argument
		arg = sys.argv[1]
# the argument may have a . inside
		pos = arg.rfind(".")
# the rootname has no . inside 
		rootname = arg[0:pos] if pos != -1 else arg
# input file name from rootname
		fn = arg if os.path.exists(arg) else rootname + default.lsd_extension
	else:
# use default rootname
		rootname = default.rootname
		fn = None
# make a modified defaults object with the current rootname to which an underscore (_) is appended
	default.rootname = rootname + "_"
# get an Eraser object to remove all input and output LSD files that start with rootname
	eras = eraser.Eraser(default)
# erase the input/solution LSD files that start with rootname
	eras.erase()
# make a new LSD object from defaults and input file name (None if standard input)
	lsd = LSD(default, fn)
# look if VMF and VSA expansion yielded at least one LSD data set
	if lsd.notEmpty:
# create LSD input files
		lsd.generate_files()
# run the created input files
		lsd.run_files()
# make a single result file from the ones created by the lsd input file instances
		lsd.collect_results()
# print the final number of solutions
		print ("Overall: %d solution%s in %s" % (lsd.nsol, "s" if lsd.nsol > 1 else "", lsd.solpath))
# post-processing ranks the solutions, creates structure diagrams and displays them
		lsd.postproc()	
# incoherent input file
	else:
# no LSD input file were generated from the PyLSD input
		print ("No output file produced. Check VMF and VSA related commands.")

# TODO
# test test test
# Fix the CNTD(LSD)/PIEC(VSAResolver) problem

# 2 atoms with different chemical shifts are not equivalent for lsd.py
# lsd.c makes equivalence classes of C atoms if CCLA is 1 (defaults is 0)
# lsd.c with CCLA 1 should never consider two X atoms as equivalent if they have
# different chemical shifts or if their attached H (if any) have the same chemical shift.

# Interprete LSD atom list setting commands and implement the "conditional property" (COPR)
# command that is formally identical to the PROP command but with an initial list parameter.
# if the list is not empty then the property command is printed for LSD
# (without the initial list argument) else it is printed as a command.
# Thus a C atom at 110 ppm can have only carbon neighbors if sp2 and
# exactly 2 oxygen neighbors if sp3.
