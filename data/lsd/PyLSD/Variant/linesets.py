"""The LineSets class handles the sets of LSD input lines
in order to resolve the ambiguities in the molecular formula
"""
import status
import Line
import VSAresolver
import re
import sumMF
import shifts
from functools import reduce

class LineSets:
	"""Expansion of molecular formula ambiguities"""
	def __init__(self, text):
		"""transforms text, a list of character lines (each ended by \n),
		into a list of (list of Line objects). Each list corresponds to
		a particular expansion of the molecular formula, when expansion is required."""

# LineSetArray is the list of (list of Line objects) to be produced
		self.LineSetArray = []
# Error warning messages are collected in the errormsg string. For debugging.
		self.errormsg = ""
# Lines0 is the list of the Line objects that correspond to the original text lines.
		self.Lines0 = [Line.decode(line) for line in text]
# list of indexes of MULT lines, useful if there are variants on the number of heavy atoms
# and if chemical shift prediction is possible
		self.imultlines = findallindexeswithtype(self.Lines0, Line.MULTLine)
# list of indexes of SHIX lines for the generation of files with chemical shift values
		self.ishixlines = findallindexeswithtype(self.Lines0, Line.SHIXLine)
# prepare the generation of chemical shift files, if needed
		self.cs = shifts.shifts(self)
# Tries to know whether the molecular formula requires ambiguity expansion
		needsMFexpansion = self.needsMFexpansion()
		print ("needsMFexpansion:", needsMFexpansion)
		if needsMFexpansion:
# expands
			self.MFexpand()
		else:
			print ("No ambiguity in molecular formula")
# the orginal Line objects constitue the single set of required Line objects
			self.LineSetArray.append(self.Lines0)
# print error messages, if any
		if len(self.errormsg):
			print ("*********** Error(s)")
			print (self.errormsg)
# print the obtained list(s) of Line lists.
		print()
		print ("LSD data sets with disambiguated molecular formula, VSA still possible")
		for Lines in self.LineSetArray:
			print()
			for line in Lines:
				print (line.text)
	
	def run(self):
		"""transforms the lists in LineSetArray into fully disambiguated datasets for LSD in alllsddata"""
# alllsddata is a list of texts, each represents a separate input file for LSD
		self.alllsddata = []
# get a list of Line objects
		for Lines in self.LineSetArray:
			print()
			print ("*"*25, " BEGIN LSD_problem_with_VSA")
			for line in Lines:
				print (line.text)
			print ("*"*25, "  END  LSD_problem_with_VSA")
# transforms the current list of Line objects into textual content(s) of LSD input files
			self.procLineSet(Lines)
		print ()
		print ("Fully disambiguated LSD data sets")
		for lsddata in self.alllsddata:
			print ()
			print (lsddata)
	
	def procLineSet(self, Lines):
		"""transforms Lines, a list of Line objects, into textual content(s) of LSD input files"""
# each line is transformed into a list of parameters
		pb = [line.for_resolver() for line in Lines]
# remove empty lists in pb. Only pertinent VSAResolver information is kept
		pb = [x for x in pb if x]
		print ("pb: ", pb)
# resolves VSAs (one possibility for variable status atom -> variable possibilities with fixed status atoms)
		r = VSAresolver.VSAresolver(pb)
		print (r)
# kept track of Line objects that correspond to MULT commands
		multlines = [line for line in Lines if line.updatable()]
# stbyat is the list of sets of MULT lines
		stbyat = r.stbyat
# go on?
		if len(stbyat) == 0:
			return
# scan through input file indexes and corresponding sets of MULT lines
		for sts in stbyat:
# get MULT lines
			newtextlines = r.linearray(sts)
# scan through MULT lines and corresponding Line objects
			for newcontent, multline in zip(newtextlines, multlines):
# replace old content by new content in MULT lines
				multline.update(newcontent)
# new LSD data set by concatenation of textual lines
			lsddata = reduce(lambda x, y: x+y, [line.from_resolver() for line in Lines], "")  #0901 faire une boucle
# append it to the list
			self.alllsddata.append(lsddata)

	def needsMFexpansion(self):
		"""return True if the FORM command in Lines, the original Line object list from the original text input,
		contains variable parts"""
		Lines = self.Lines0
# look in Lines for the list of indexes of FORM lines
		iformlines = findallindexeswithtype(Lines, Line.FORMLine)
# number of FORM lines in original text, should be 1
		numformlines = len(iformlines)
		if numformlines != 1:
			self.errormsg += "Only one FORM command allowed. Presently: %d\n" % numformlines
# !!!!!!!!!!!!!!! the case iformlines would be empty is not considered !!!!!!!!!!!!!!!
		iformline = iformlines[-1]
# formline is the reference of the last FORM command of the data set
		self.iformline = iformline
# get the active FORM Line object and get information out of it
		formline = Lines[iformline]
		print
		print ("formline:", formline)
# remove the enclosing double quotes around the formula textual parameter
		formula = formline.formula.strip("\"")
		print ("formula:", formula)
# get formula pieces
		pieces = formula.split()
# get number of formula pieces
		numpieces = len(pieces)
		if numpieces % 2:
# the number of formula pieces must be even
			self.errormsg += "Odd number of fields in FORM parameter: %d\n" % numpieces
# get number of elements (odd-indexed in pieces)
		numbers = [pieces[i] for i in range(1, numpieces, 2)]
# get corresponding element symbols (even-indexed in pieces)
		elements = [pieces[i] for i in range(0, numpieces, 2)]
# formula is a dict in which the keys are the element symbols
# and the values are the list of the possible atom numbers
		self.formula = {}
# by default, no formula ambiguity is pre-supposed
		retval = False
# loop over element symbols and associated atom numbers
		for elt, num in zip(elements, numbers):
# put the textual number information into pieces
			m = re.match(r"^(\d+)(-(\d+))?$", num)
# an int (group 1) or and int followed by a dash and int (group 3) must be there
			first = int(m.group(1))
			last = m.group(3)
# make an interval of number values, singleton if no int-int pattern
			interval = range(first, int(last)+1) if last else [first]
# detection of bad intervals
			if last and not interval:
				self.errormsg += "Bad range for atom number in FORM: %d-%d\n" % (first, int(last))
# expansion is needed if at least one element has an interval with more than one value
			retval = retval or len(interval) > 1
# stores formula attribute of the current object
			self.formula[elt] = interval
		print ("formula as dict:", self.formula)
# perform some text analysis in order to extract information
# that is needed by all subsequent processing steps and for error detection.
# looking the indexes of the MOMA lines in the original text
		imomalines = findallindexeswithtype(Lines, Line.MOMALine)
# count MOMA lines
		nummomalines = len(imomalines)
		if nummomalines > 1:
			self.errormsg += "No more than one FORM command allowed. Presently: %d\n" % nummomalines
# get index of the active MOMA line. -1 if none.
		imomaline = imomalines[-1] if nummomalines != 0 else -1
# store position of the MOMA line as current object attribute
		self.imomaline = imomaline 
# momaline is the reference to the last MOMA command of the data set. None if none.
		momaline = Lines[imomaline] if imomaline >= 0 else None
		if momaline and not retval:
# A molecular mass was found and there is no molecular formula ambiguity
# check if these two informations are compatible.
# Get flatformula as a dict that gives a number of atoms for each element instead of a range
			flatformula = dict([[elt, nums[0]] for elt, nums in self.formula.items()])
# calculate mass from the FORM line
			massfromform = calcmass(flatformula)
			print ("Mass from formula:", massfromform)
# get minimum and maximum mass from the MOMA line
			massfrommoma = momaline.mass
# look if FORM and MOMA fit together
			goodmass = massfrommoma[0] <= massfromform <= massfrommoma[1]
			if not goodmass:
# bad fit
				self.errormsg += "Mismatch between mass from FORM (%d) and mass range (%d-%d) from MOMA" \
					% (massfromform, massfrommoma[0], massfrommoma[1])
		print ("retval:", retval)
# return True if molecular formula expansion is needed
		return retval

	def MFexpand(self):
# get molecular mass range from MOMA line (mandatory if MF expansion in necessary)
		if self.imomaline < 0:
			self.errormsg += "A MOMA line is required for MF expansion\n"
			return
# get position of MOMA line from current LineSets object
		imomaline = self.imomaline
# get pointer to original Line objects
		Lines = self.Lines0
# get pointer to MOMA line
		momaline = Lines[imomaline]
# get molecular mass constraint
		massfrommoma = momaline.mass
		print ("mass from MOMA line:", massfrommoma)
# looking for an ELEC line, in order to constrain the parity of the molecular valence
		ieleclines = findallindexeswithtype(Lines, Line.ELECLine)
# number of ELEC lines
		numeleclines = len(ieleclines)
		if numeleclines > 1:
			self.errormsg += "No more than one ELEC command allowed. Presently: %d\n" % numeleclines
# index of the last ELEC line in the orginal text, None if none
		ielecline = ieleclines[-1] if numeleclines != 0 else None
# reference to active ELEC line, None if none
		elecline = Lines[ielecline] if ielecline else None
# molecular electric charge, as a list of possible values. [0] if no ELEC line
		elec = elecline.elec if elecline else [0]
		print ("electric charge:", elec)
# make the list of the possible electric charge parity (remainder of its division by 2)
		parities = sorted(list(set([e%2 for e in elec])))
# parity matters if only one parity value is possible
		paritymatters = len(parities) == 1
		print ("paritymatters:", paritymatters)
# elecparity is the parity of the molecular electric charge, only meaningful if paritymatters
		elecparity = parities[0]
		if paritymatters:
			print ("electric charge parity:", elecparity)
# get formula back from current object
		formula = self.formula
		print ("formula:", formula)
# get list of the possible number of H atoms
		numH = formula["H"] if "H" in formula else []
		print ("numH:", numH)
# the information about H is kept in numH, remove H as key in formula
		if numH:
			del formula["H"]
# list of elements for which the exact number of atoms is unknown (not exactly known)
		uElts = [elt for elt in list(formula.keys()) if len(formula[elt]) > 1]
# is H the only element for which the number of atoms in not known?
		if len(uElts) == 0:
# at this point there is an element in an unkown number because MFexpand is run only if expansion is needed.
# It must be H because the list of non-H elements in known numner is empty.
# Expansion of the number of H only
			print ("Ambiguity only on the number of H atoms:", numH)
# MF associate an element symbol with the corresponding unambiguous number of atoms
			MF = dict([[elt, nums[0]] for elt, nums in formula.items()])
# the number of H is contrained if paritymatters
			if paritymatters:
# MFelec_parity is the valence parity to which the electric charge parity is added.
# only the parity of the number of H atoms is not taken into account yet
				MFelec_parity = VSAresolver.MFnumH_parity(MF, elecparity)
# keep from the numH list only the values that will make an even parity
				numH = [h for h in numH if h%2 == MFelec_parity]
				print ("Ambiguity only on the number of H atoms (considering ELEC):", numH)
# get molecular mass without H atoms
			mass_noH = calcmass(MF)
# keep from the numH list only the values that are compatible with the constraint of the MOMA command
			numH = [h for h in numH if massfrommoma[0] <= mass_noH + h <= massfrommoma[1]]
			print ("Ambiguity only on the number of H atoms (considering ELEC and MOMA):", numH)
# all the values in numH are valid number of H atoms
			for h in numH:
# completes the molecular formula with the current number of H aoms
				MF["H"] = h
				print()
				print ("MF:", MF)
# create a textual replacement for the actual FORM line
				newFORMtext = makeFORMtext(MF)
				print ("newFORMtext", newFORMtext)
# make a Line object from the current textual FORM line
				newline = Line.decode(newFORMtext)
# replace the old FORMLine object by the new one
				Lines[self.iformline] = newline
# register the current text as a new list of Lines objects
				self.LineSetArray.append(list(Lines))
		else:
# there are non-H elements for which the number of atoms is unknown
# kElts is the list of the elements for which the number of atoms is known
			kElts = [elt for elt in list(formula.keys()) if len(formula[elt]) == 1]
			print ("elements in known number:", kElts)
			print ("elements in unknown number:", uElts)
# known part of the molecular formula (no H)
			kMF = dict([[elt, formula[elt][0]] for elt in kElts])
			print ("known part of the MF:", kMF)
# unknown part of the molecular formula (no H)
			uMF = dict([[elt, formula[elt]] for elt in uElts])
			print ("unknown part of the MF:", uMF)
# minimal number for elements in unknown number from FORM line
# this number must be greater or equal to the number of MULT lines for the same element
			umin = dict([[elt, uMF[elt][0]] for elt in uElts])
			print ("minimal number for elements in unknown number:", umin)
# list of indexes of MULT lines
			imultlines = self.imultlines
# list of MULT Line objects
			multlines = [Lines[i] for i in imultlines]
# number of MULT lines for the elements in unknown number
			num_mult_lines = dict([[elt, len([line for line in multlines if line.initial_status[1] == elt])] for elt in uElts])
			print ("number of MULT lines for the elements in unknown number", num_mult_lines)
# elements with atoms in unknown number for which the number of MULT lines is too high
			uElt_error = [elt for elt in uElts if umin[elt] < num_mult_lines[elt]]
			if len(uElt_error):
# at least one elements has a FORM/MULT mismatch
				print ("Too many MULT lines according to the FORM line for: " + repr(uElt_error) + "\n")
				self.errormsg += "Too many MULT lines according to the FORM line for: " + repr(uElt_error) + "\n"
				return
# collect implicit default status for the elements in unknown number
			udefstat = dict([[elt, status.defstat[elt]] for elt in uElts])
			print ("implicit default status for the elements in unknown number:", udefstat)
# indexes of DEMU (DEfault MUltiplicity) lines
			idemulines = findallindexeswithtype(Lines, Line.DEMULine)
# list of DEMU lines
			demulines = [Lines[i] for i in idemulines]
# start collecting the default status from the DEMU lines
			explicitdefstats = {}
# scan through DEMU lines
			for line in demulines:
# current element
				elt = line.elt
# corresponding current default status
				defstat = line.defstat
# has a status been already explicity defined from the current element?
				if elt in list(explicitdefstats.keys()):
# append the current status to the already known one
					explicitdefstats[elt].append(line.defstat)
				else:
# first DEMU line for this element
					explicitdefstats[elt] = [defstat]
			print ("explicit defaut status (from DEMU line):", explicitdefstats)
# scan through elements and default status from DEMU lines
			for elt, defstats in explicitdefstats.items():
# number of status for the current element
				n = len(defstats)
				if n > 1:
# multiply defined default status
					print ("More than one default status for %s: %d\n" % (elt, n))
					self.errormsg += "More than one default status for %s: %d\n" % (elt, n)
				if elt in udefstat:
# if the element has already a default default status, remplace it by the one from DEMU
					udefstat[elt] = defstats[-1]
				else:
# useless DEMU line
					print ("Unused default status (DEMU line) for %s\n" % elt)
					self.errormsg += "Unused default status (DEMU line) for %s\n" % elt
			print ("default status for the elements in unknown number:", udefstat)
# default status as [elt, hyb, mult, charge] is converted to textual format
			self.formattedUdefstat = dict([[elt, formatStatforLSD(stat)] for elt, stat in udefstat.items()])
			print ("formatted default status for the elements in unknown number:", self.formattedUdefstat)
# prepare to get first possible atom id, get position for the insertion of new MULT commands in Lines0
# make list of lists of atom ids
			ids = [line.initial_status[0] for line in multlines]
# list of atom ids without duplicates
			ids_nodup = list(set(ids))
# check if there is at least an id is duplicated
			if len(ids) != len(ids_nodup):
				print ("At least one atom id is duplicated\n")
				self.errormsg += "At least one atom id is duplicated\n"
# id duplication is a fatal error
				return
# get first free atom id for the atoms in FORM that have no MULT counterpart
# and that have to be automatically added with their default status
			self.first_free_id = max(ids) + 1
			print ("first free atom id:", self.first_free_id)
# MULT lines are generally grouped together. insertion point
# for eventual auto-generated MULT lines is immediately after the last explicit MULT line.
			self.insert_before_line = imultlines[-1] + 1
			print ("insert new MULT lines before line (start at 0):", self.insert_before_line)
# if valence parity does not matter ther is no point to group the elements with an odd velence first
			if paritymatters:
# calculate the sum of valence parity for the known part of the molecular formula, including the electric charge
				kMFelec_parity = VSAresolver.MFnumH_parity(kMF, elecparity)
				print ("total valence parity of the elements in known number, including total charge:", kMFelec_parity)
# elements in unknown number with odd valence parity
				uElts_odd = [elt for elt in uElts if status.elt_valence_parity[elt] == 1]
				print ("elements in unknown number with odd valence parity:", uElts_odd)
# number of elements with atoms in unknown number and with odd valence parity
				num_parity_odd = len(uElts_odd)
				print ("number of elements with atoms in unknown number and with odd valence parity:", num_parity_odd)
# elements in unknown number with even valence parity
				uElts_even = [elt for elt in uElts if status.elt_valence_parity[elt] == 0]
				print ("elements in unknown number with even valence parity:", uElts_even)
# ordered elements for which atoms are in unknown number, elements with odd valence first,
# for which valence parity constitutes a constraint
				uElts_ordered = uElts_odd + uElts_even
			else:
# no need for element ordering according to valence parity
				uElts_ordered = uElts
			print ("valence-parity-ordered elements in unknown number:", uElts_ordered)
# number of atoms by valence-parity-ordered elements
			unumA = [formula[elt] for elt in uElts_ordered]
			print ("number of atoms by valence-parity-ordered elements", unumA)
# atomic mass by valence-parity-ordered elements
			umassA = [status.masses[elt] for elt in uElts_ordered]
			print ("atomic mass by valence-parity-ordered elements", umassA)
# total mass of the atoms whose elements are in known number
			kmass = calcmass(kMF)
			print ("total mass of the atoms whose elements are in known number:", kmass)
# range of mass for the atoms whose elements are in unknown number
			umass = [massfrommoma[0]-kmass, massfrommoma[1]-kmass]
			print ("range of mass for the atoms whose elements are in unknown number:", umass)
# get the lists	of atom number lists for the elements for which the atoms are in unknown number. number of H is prepended, if any H	
			coeffss = sumMF.sumMF(paritymatters, num_parity_odd, kMFelec_parity, numH, umass, unumA, umassA)
#			coeffss = [[3, 1, 2]]
			for coeffs in coeffss:
# exploit the current list of atom numbers.
# make a copy of the original list of Line objects
				newLines = list(Lines)
				print
# the current molecular formula at least contains the atoms whose number is known
				formulaexact = dict(kMF)
				if numH:
# the exact current number of H is the first value in coeffs. get it and remove it from coeffs
					numHexact = coeffs.pop(0)
					print ("exact number of H:", numHexact)
# update current molecular formula with the number of H atoms, if any
					formulaexact["H"] = numHexact
# scan through elements in unknown number and the currently found numbers
				for i, elt in enumerate(uElts_ordered):
# update current molecular formula
					formulaexact[elt] = coeffs[i]
				print ("formulaexact:", formulaexact)
# get textual version of the exact current molecular formula
				newFORMtext = makeFORMtext(formulaexact)
				print ("newFORMtext", newFORMtext)
# make a new Line object from the new, unambiguous (I mean, exact) molecular formula
				newline = Line.decode(newFORMtext)
# replace the copy of the original FORM Line object by the new one
				newLines[self.iformline] = newline
# prepare the generation of the MULT lines for the atoms that do not have a MULT line yet
				for i, elt in enumerate(uElts_ordered):
# number of MULT lines to be generated, difference between the required number and
# the already existing number. This number is >= 0 because coeffs[i] >= umin[elt] >= num_mult_lines[elt]
					numlines = coeffs[i] - num_mult_lines[elt]
					print ("generate %d default MULT line(s) for %s" % (numlines, elt))
# generate the missing MULT lines and insert them at the right place in the copy of the original list of Lines
					self.generateMULTlines(newLines, elt, numlines)
# this new (list of Line) is appended to the list of (List of Line)
				self.LineSetArray.append(newLines)

	def generateMULTlines(self, newLines, elt, numlines):
		"""insert in newLines, a list of Line objects, numlines new MULT lines for element elt"""
# loops numlines times
		for i in range(numlines):
# get default status for elt
			status = self.formattedUdefstat[elt]
# where (before which line) to insert, index in newLines
			pos = self.insert_before_line
# get new atom index
			id = self.first_free_id
# make the textual version of the new MULT line
			newMULTtext = "MULT " + str(id) + status + " ; auto-generated \n"
# make a new Line object from the textual version
			newline = Line.decode(newMULTtext)
# insert the new line at the place
			newLines.insert(pos, newline)
# update the insertion position for the next MULT line
			self.insert_before_line += 1
# find the next available atom id
			self.first_free_id += 1

	def getSolutionCounter(self, defsolncounter):
# array of original Lines
		Lines = self.Lines0
# looking the indexes of the COUF lines in the original text
		icouflines = findallindexeswithtype(Lines, Line.COUFLine)
# count COUF lines
		numcouflines = len(icouflines)
		if numcouflines > 1:
			self.errormsg += "No more than one COUF command allowed. Presently: %d\n" % numcouflines
# get index of the active COUF line. -1 if none.
		icoufline = icouflines[-1] if numcouflines != 0 else -1
		return defsolncounter if icoufline < 0 else Lines[icoufline].couf
		
	def getStopFile(self, defstopfile):
# array of original Lines
		Lines = self.Lines0
# looking the indexes of the STOF lines in the original text
		istoflines = findallindexeswithtype(Lines, Line.STOFLine)
# count STOF lines
		numstoflines = len(istoflines)
		if numstoflines > 1:
			self.errormsg += "No more than one STOF command allowed. Presently: %d\n" % numstoflines
# get index of the active STOF line. -1 if none.
		istofline = istoflines[-1] if numstoflines != 0 else -1
		return defstopfile if istofline < 0 else Lines[istofline].stof
		
def calcmass(formula):
	"""calculate the mass that corresponds to molecular formula (as dict)"""
# get masses dictionary from status module
	masses = status.masses
# initial mass
	m = 0
# scan through element, corresponding number of atoms pairs
	for elt, num in formula.items():
# update mass with the mass of the current atoms of the current element
		m += num * masses[elt]
# return molecular mass
	return m

def findallindexeswithpredicate(l, f):
	"""return the list of indexes in a list l that corresponds only to the list elements for which predicate f is True"""
	return [i for i, x in enumerate(l) if f(x)]

def findallindexeswithtype(l, t):
	"""return the list of indexes in a list l that have t as type"""
	return findallindexeswithpredicate(l, lambda x: isinstance(x, t))

# def findindexoflastwithtype(l, t):
	# indexes = findallindexeswithtype(l, t)
	# return indexes[-1] if indexes else -1

def makeFORMtext(MF):
	"""transform a molecular formula (as dict) into a textual FORM line"""
# start FORM line with a FORM command, a blank, and an opening double quote
	ch = "FORM \""
# list of textual elements that constitue the arguments of the FORM line
	l = []
# scan through the canonical order of elements C, H, N, O and the others by alphabetical order
	for elt in status.canonical_elt_order:
# do something for the elemnts in the molecular formula
		if elt in MF:
# number of atoms for the current element
			num = MF[elt]
# do something only if this number is not zero
			if num:
# append the element symbol to the list of textual elements
				l.append(elt)
# append the string that corresponds to the number of atoms
				l.append(str(num))
# append the blank separated textual elements to the already formed string
	ch += ' '.join(l)
# end the result with the closing double quote and a new line character
	ch += "\"\n"
	return ch

def formatListforLSD(l, opt=0):
	"""get a textual version of a list in the format of LSD alternatives"""
# empty list
	if len(l) == 0:
# return empty string
		return ""
# list of one value which is equal to zero in an optional context (electric charge)
	if len(l) == 1 and l[0] == 0 and opt:
# return empty string
		return ""
# singleton
	if len(l) == 1:
# return the textual equivalent
		return str(l[0])
# list with more than one element, return a left parenthesis, followed by a blank-separated list
# of textual value equivalents, followed by a right parenthesis
	return "(" + " ".join([str(i) for i in l]) + ")"

def formatStatforLSD(stat):
	"""transforms a [elt, hyb, mult, elec] status quadruplet into a textual version"""
# make a list of textual version of the status parts. starts with the element symbol
	pieces = [stat[0]]
# go on with hybridization state
	pieces.append(formatListforLSD(stat[1]))
# go on with multiplicity
	pieces.append(formatListforLSD(stat[2]))
# go on with optional ("" if 0) electric charge
	pieces.append(formatListforLSD(stat[3], "opt"))
# remove trailing empty string (for 0 as charge)
	if pieces[-1] == "":
		pieces.pop()
# return a blank followed by blank seprated textual status items
	return " " + " ".join(pieces)

if __name__ == "__main__":
	ch = """; Nitromethane and methyl nitrate

CNTD 0
FORM "C 1 H 3 N 1 O 2-3"
DEMU O O (2 3) (0 1) (-1 0)
; DEMU N N (1 2 3) (0 1 2 3) (0 1)
MOMA 1-1000
PIEC 1
ELEC 0
MAXP 1
MAXN 1
MULT 1 C 3 3
MULT 2 N (1 2 3) (0 1 2 3) (0 1)
MULT 3 O (2 3) (0 1) (-1 0)
MULT 4 O (2 3) (0 1) (-1 0)
FULL L1
; LIST L2 3 4
ELEM L2 O
DIFF L1 L2 L3
PROP L2 0 L3
EXIT
"""
	print
	text1 = ch.splitlines(True)
	print ("text1:", text1)
	ls = LineSets(text1)
#	ls.run()
