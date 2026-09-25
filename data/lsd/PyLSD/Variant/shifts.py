"""Production of pilot-files for X chemical shift predictors, from a linesets object"""
import Line
import status
import re
import os

class shifts:
	"""Production of pilot-files for chemical shift predictors, from a linesets object"""
	def __init__(self, ls):
		"""Collect atom numbering and chemical shift from input file"""
		Lines = ls.Lines0
# get MULT line objects
		multlines = [Lines[i] for i in ls.imultlines]
# get X (X is any element but H) chemical shift information
		shixlines = [Lines[i] for i in ls.ishixlines]
# sort the MULT lines according to their LSD numbering
		multlines.sort(key = lambda l: l.initial_status[0])
# associate an LSD atom number (in text format) with the corresponding SD atom number (start at 1)
# prediction will be achieved from an SDF file in which the LSD numbers are forgotten
		lsd2sd = {}
		for i, line in enumerate(multlines):
# current SD number
			sdindex = i + 1
# current LSD number
			lsdindex = line.initial_status[0]
# current element, remove unusual valence indication, if present
			elt = simple_element(line.initial_status[1])
# associates an LSD number (as text) to the correspond (SD number, element) pair
			lsd2sd[str(lsdindex)] = (sdindex, elt)
		print ("lsd2sd:", lsd2sd)
# associates an element with the (SD number, LSD number, chemical shift) triplet
		self.elt2cs = {}
# scan through SHIX lines
		for line in shixlines:
# get LSD index
			lsdindex = line.id
# get chemical shift value
			shiftvalue = line.delta
# get SD number and element for this LSD number and chemical shift
			status = lsd2sd[str(lsdindex)]
# get SD index from lsd2sd
			sdindex = status[0]
# get element index from lsd2sd
			elt = status[1]
# build a new (SD number, LSD number, chemical shift) triplet for the current element
			cs = (sdindex, lsdindex, shiftvalue)
# is there already an atom of such an element with a chemical shift?
			if not elt in self.elt2cs:
# if not, start a new list of (SD number, LSD number, chemical shift) triplets
				self.elt2cs[elt] = []
# associate the element with the new (SD number, LSD number, chemical shift) triplet
			self.elt2cs[elt].append(cs)
		print ("elt2cs", self.elt2cs)
	
	def write_shift_files(self, defaults):
		"""write pilot files for the known predictors"""
# get the list of X elements for which a prediction tool exists
		predicted = defaults.predictors.keys()
# start a dict of pilot-file names with element as key
		self.filenames = {}
# scan through predictable elements
		for elt in predicted:
# are there atoms of known experimental chemical shift whose value is predictible?
			if elt in self.elt2cs:
# if yes, get the list of (SD number, LSD number, chemical shift) triplets for the current predictible element
				cslist = self.elt2cs[elt]
# find a file name for the current pilot file
				filename = os.path.join(defaults.predictorspath, defaults.rootname + elt + ".txt")
				print ("writing to: ", filename)
# remember the name of current new pilot-file, as concerning the current predictible element
				self.filenames[elt] = filename
# open the current pilot file
				fh = open(filename, "w")
# scan through atoms of the current element that have a defined experimental chemical shift
				for cs in cslist:
# print SD number, LSD number and chemical shift in pilot-file
					fh.write("%d %s %.2f\n" % cs)
# close file
				fh.close()

def simple_element(valenced_element):
	"""remove unusual valence indication of an element, if any"""
# check for the absence of valence indication
	if valenced_element in status.all_elements:
# nothing to do
		return valenced_element
	else:
# catch trailing integers
		m = re.match(r'([A-Z][a-z]?)([1-6]?)([1-6]?)([1-6]?)', valenced_element)
# get parts of valenced_element
		fields = list(m.groups())
# keep only the element name
		elt = fields[0]
# return bare element name
		return elt

if __name__ == "__main__":
	import defaults
	from linesets import findallindexeswithtype
	class fake_linesets:
		def __init__(self, text):
			self.Lines0 = [Line.decode(line) for line in text]
			self.imultlines = findallindexeswithtype(self.Lines0, Line.MULTLine)
			self.ishixlines = findallindexeswithtype(self.Lines0, Line.SHIXLine)
			self.cs = shifts(self)

	ch = """; test
MULT 1 C 3 2
MULT 2 C 3 2
MULT 3 C 3 2
MULT 5 F 3 0
MULT 6 N 3 2
MULT 7 N5 1 0

SHIX 1 20.1
SHIX 2 25.2
SHIX 5 85
SHIX 6 -43.5
SHIX 7 123.7
"""

	text = ch.splitlines(True)
	fls = fake_linesets(text)
	defs = defaults.Defaults()
	fls.cs.write_shift_files(defs)
