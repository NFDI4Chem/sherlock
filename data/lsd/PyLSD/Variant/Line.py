"""LSD command modelling"""
import re

class Line:
	"""Any text line, without the final \n. Without the comment in the command lines"""
	def __init__(self, text):
		print ("Line")
		"""new Line object, with initialization"""
		self.text = text
		print ("text:", self.text)
	def __repr__(self):
		"""Line printing"""
		return self.text
	def for_resolver(self):
		"""Only FORM PIEC ELEC MAXP MAXN SHIX and MULT lines provide something to the VSAResolver"""
		return []
	def updatable(self):
		"""Only MULT lines can be updated"""
		return False
	def from_resolver(self):
		"""string representation of a line, with the final \n"""
		return self.__repr__() + '\n'

class EmptyLine(Line):
	def __init__(self, text):
		print ("EmptyLine")
		Line.__init__(self, text)
	# def original(self):
		# return self.text

class CommandLine(Line):
	"""Any LSD command, commented or not (no blank line, no comment line)"""
	def __init__(self,  head, sep, tail, comment):
		"""LSD command: mnemonic, separator, arguments, eventual comment"""
		print ("CommandLine")
		print ("head: ", "***%s***" % head)
		print ("sep: ", "***%s***" % sep)
		print ("tail: ", "***%s***" % tail)
		print ("comment: ", "***%s***" % comment)
		Line.__init__(self, head + sep + tail)
		self.head, self.sep, self.tail, self.comment = head, sep, tail, comment
	def __repr__(self):
		"""print by assembling the command itself with the initial comment"""
		return self.text + self.comment
	# def original(self):
		# return self.head + self.sep + self.tail + self.comment

class FORMLine(CommandLine):
	"""FORM command, specific to VSAResolver, gives the molecular formula"""
	def __init__(self, commandargs, formarg):
		"""initialization with standard command parts and formula specific string"""
		print ("FORMLine")
		print ("commandargs: ", commandargs)
		print ("formarg: ", formarg)
		CommandLine.__init__(self, *commandargs)
		self.formula = formarg
	def for_resolver(self):
		"""FORM command list for VSAResolver"""
		return ["FORM", self.formula]
	def from_resolver(self):
		"""commented line returned because lsd does not know about molecular formula"""
		return "; " + repr(self) + "\n"

class PIECLine(CommandLine):
	"""PIEC command, specific to VSAResolver, gives the max number of molecular pieces"""
	def __init__(self, commandargs, piecarg):
		"""initialization with standard command parts and connexivity specific integer (max pieces)"""
		print ("PIECLine")
		print ("commandargs: ", commandargs)
		print ("piecarg: ", piecarg)
		CommandLine.__init__(self, *commandargs)
		self.pieces = piecarg
	def for_resolver(self):
		"""PIEC command list for VSAResolver."""
		return ["PIEC", 1, self.pieces]
	def from_resolver(self):
		"""commented line returned because lsd does not know about connexivity"""
		return "; " + repr(self) + "\n"

class ELECLine(CommandLine):
	"""ELEC command, specific to VSAResolver, gives the total electric charge"""
	def __init__(self, commandargs, elecarg):
		"""initialization with standard command parts and molecular electric charge"""
		print ("ELECLine")
		print ("commandargs: ", commandargs)
		print ("elecarg: ", elecarg)
		CommandLine.__init__(self, *commandargs)
		self.elec = elecarg
	def for_resolver(self):
		"""ELEC command list for VSAResolver."""
		return ["ELEC", self.elec]
	def from_resolver(self):
		"""commented line returned because lsd does not know about molecular electric charge"""
		return "; " + repr(self) + "\n"

class MAXPLine(CommandLine):
	"""MAXP command, specific to VSAResolver,
	gives the maximum number of positively charged atoms.
	The number is -1 if no constraint on this number is asked for. (default)"""
	def __init__(self, commandargs, maxparg):
		"""initialization with standard command parts and maximum number of positively charged atoms"""
		print ("MAXPLine")
		print ("commandargs: ", commandargs)
		print ("maxparg: ", maxparg)
		CommandLine.__init__(self, *commandargs)
		self.maxp = maxparg
	def for_resolver(self):
		"""MAXP command list for VSAResolver."""
		return ["MAXP", self.maxp]
	def from_resolver(self):
		"""commented line returned because lsd does not know about maximum number of positively charged atoms"""
		return "; " + repr(self) + "\n"

class MAXNLine(CommandLine):
	"""MAXN command, specific to VSAResolver,
	gives the maximum number of negatively charged atoms.
	The number is -1 if no constraint on this number is asked for. (default)"""
	def __init__(self, commandargs, maxnarg):
		"""initialization with standard command parts and maximum number of negatively charged atoms"""
		print ("MAXNLine")
		print ("commandargs: ", commandargs)
		print ("maxnarg: ", maxnarg)
		CommandLine.__init__(self, *commandargs)
		self.maxn = maxnarg
	def for_resolver(self):
		"""MAXN command list for VSAResolver."""
		return ["MAXN", self.maxn]
	def from_resolver(self):
		"""commented line returned because lsd does not know about maximum number of negatively charged atoms"""
		return "; " + repr(self) + "\n"

class MOMALine(CommandLine):
	"""MOMA command, specific to linesets,
	gives the molecular mass or a range of molecular mass."""
	def __init__(self, commandargs, massarg):
		"""initialization with standard command parts and mass information"""
		print ("MOMALine")
		print ("commandargs: ", commandargs)
		print ("massarg: ", massarg)
		CommandLine.__init__(self, *commandargs)
		self.mass = massarg
	def from_resolver(self):
		"""commented line returned because lsd does not know about molecular mass"""
		return "; " + repr(self) + "\n"

class DEMULine(CommandLine):
	"""DEMU command, specific to linesets,
	gives the default MULT command."""
	def __init__(self, commandargs, demuarg):
		"""initialization with standard command parts and status information"""
		print ("DEMULine")
		print ("commandargs: ", commandargs)
		print ("demuarg: ", demuarg)
		CommandLine.__init__(self, *commandargs)
		self.elt = demuarg[0]
		self.defstat = demuarg[1:]
		print ("element:", self.elt)
		print ("default status:", self.defstat)
	def from_resolver(self):
		"""commented line returned because lsd does not know about default MULT parameters of atoms"""
		return "; " + repr(self) + "\n"

class SHIXLine(CommandLine):
	"""SHIX command, indicate the chemical shift value of non-H atoms."""
	def __init__(self, commandargs, shixarg):
		"""initialization with standard command parts and chemical shift information"""
		print ("SHIXLine")
		print ("commandargs: ", commandargs)
		print ("shixarg: ", shixarg)
		CommandLine.__init__(self, *commandargs)
		self.id = shixarg[0]
		self.delta = shixarg[1]
		print ("atom id:", self.id)
		print ("chemical shift:", self.delta)
	def for_resolver(self):
		"""SHIX command list for VSAResolver."""
		return ["SHIX", int(self.id), self.delta]

class COUFLine(CommandLine):
	"""COUF command, indicates the path to the file that contains the number of found solutions"""
	def __init__(self, commandargs, coufarg):
		"""initialization with standard command parts and file path information"""
		print ("COUFLine")
		print ("commandargs: ", commandargs)
		print ("coufarg: ", coufarg)
		CommandLine.__init__(self, *commandargs)
		self.couf = coufarg
		print ("path to solution count file:", self.couf)

class STOFLine(CommandLine):
	"""STOF command, indicates the path to the file that stops lsd"""
	def __init__(self, commandargs, stofarg):
		"""initialization with standard command parts and file path information"""
		print ("STOFLine")
		print ("commandargs: ", commandargs)
		print ("stofarg: ", stofarg)
		CommandLine.__init__(self, *commandargs)
		self.stof = stofarg
		print ("path to lsd stop file:", self.stof)

class MULTLine(CommandLine):
	"""MULT command, for VSAs and FSAs"""
	def __init__(self, commandargs, multargs):
		"""initialization with standard command parts and MULT specific (id, elt, hyb, mult)"""
		print ("MULTLine")
		print ("commandargs: ", commandargs)
		print ("multargs: ", multargs)
		CommandLine.__init__(self, *commandargs)
		self.initial_status = multargs
	def for_resolver(self):
		"""MULT command list for VSAResolver."""
		return ["MULT"] + self.initial_status
	def updatable(self):
		"""the only command to be updatable for the replacement of VSAs by FSAs, the only one lsd understands so far"""
		return True
	def update(self, newtext):
		"""stores FSA data in VSAs"""
		self.text = newtext

def decode(text):
	"""return a Line object from a line in an LSD input file (with VSAs, possibly)"""
	print ()
	print ("*"*10, "New Line")
# remove trailing \n and enclosing white chars
	tr = text.rstrip("\n")
	ts = tr.strip()
# nothing left. this is an empty line
	if len(ts) == 0:
		print ("Empty line")
		return EmptyLine(tr)
# look for a comment line (only comment inside)
#	m = re.match(r'(\s*;.*$)', ts)
	m = re.match(r'(^;.*$)', ts)
	if m:
# a comment line was found
		print ("Comment line")
		return EmptyLine(tr)
# look for an EXIT line
	m = re.match(r'^EXIT$', ts)
	if m:
# an exit line was found
		print ("EXIT line")
		return EmptyLine(tr)
# look for a command line
	m = re.match(r'^([A-Z 23]{4})(\s+)(.*?)(\s*;.*)?$', ts)
	if m:
# a command line was found
		print ("Command line")
# get mnemonic, separator, arguments and comment
		commandargs = [m.group(i) for i in range(1,5)]
# absent fragment are converted to null strings
		#commandargs = map(lambda x: "" if x==None else x, commandargs)
		commandargs = list(map(lambda x: "" if x==None else x, commandargs))
# head is mnemonic
		head = commandargs[0]
		if head == "FORM":
# a formula command was found
			print ("Formula line")
# remove enclosing double quotes around the string and enclosing blanks
			ch = commandargs[2].strip("\"").strip()
# ch is the string representation of the formula for VSAResolver
			L = FORMLine(commandargs, ch)
			return L
		elif head == "PIEC":
# a pieces command was found
			print ("Pieces line")
# the only argument in the maximum allowed number of pieces
			i = int(commandargs[2])
# i is the integer number of pieces for the VSAResolver
			L = PIECLine(commandargs, i)
			return L
		elif head == "MAXP":
# a maxp command was found
			print ("Maxp line")
# the only argument in the maximum allowed number of positively charged atoms
			i = int(commandargs[2])
# i is the maximal number of positively charged atoms for the VSAResolver
			L = MAXPLine(commandargs, i)
			return L
		elif head == "MAXN":
# a maxn command was found
			print ("Maxn line")
# the only argument in the maximum allowed number of negatively charged atoms
			i = int(commandargs[2])
# i is the maximal number of negatively charged atoms for the VSAResolver
			L = MAXNLine(commandargs, i)
			return L
		elif head == "ELEC":
# a elec command was found
			print ("Elec line")
# the only argument is the list of possible molecular electric charge
			ch = commandargs[2].lstrip("(").rstrip(")").split()
			elec = sorted([int(ec) for ec in ch])
# elec is the list of possible electric charge values for the VSAResolver
			L = ELECLine(commandargs, elec)
			return L
		elif head == "MOMA":
# a moma command was found
			print ("Molecular mass line", commandargs[2])
			m2 = re.match(r'^(\d+)(-(\d+))?$', commandargs[2])
			if m2:
				first = int(m2.group(1))
				last = m2.group(3)
				ma = (first, int(last)) if last else (first, first)
				L = MOMALine(commandargs, ma)
				return L
		elif head == "DEMU":
# a demu command was found
			print ("Default MULT line", commandargs[2])
			m2 = re.match(r'^([A-Z][a-z]?)\s+(.*)$', commandargs[2])
			elt = m2.group(1)
			mu = multsplit(m2.group(2))
			if mu:
# BUG fixed (?): an integer in mu[1] and mu[2] must be turned into a list with a single integer.
				if isinstance(mu[1], int):
					mu[1] = [mu[1]]
				if isinstance(mu[2], int):
					mu[2] = [mu[2]]
				mu.insert(0, elt)
				L = DEMULine(commandargs, mu)
				return L
		elif head == "SHIX":
# a shix command was found
			print ("X Chemical shift line", commandargs[2])
			m2 = re.match(r'^(\d+)\s+([+-]?\d+(\.(\d+)?)?)$', commandargs[2])
			if m2:
				id = m2.group(1)
				delta = float(m2.group(2))
				L = SHIXLine(commandargs, [id, delta])
				return L
		elif head == "COUF":
# acouf command was found
			print ("solution count file line", commandargs[2])
# remove enclosing double quotes around the string and enclosing blanks
			ch = commandargs[2].strip("\"").strip()
			L = COUFLine(commandargs, ch)
			return L
		elif head == "STOF":
# acouf command was found
			print ("lsd stop file line", commandargs[2])
# remove enclosing double quotes around the string and enclosing blanks
			ch = commandargs[2].strip("\"").strip()
			L = STOFLine(commandargs, ch)
			return L
		elif head == "MULT":
# a mult command was found
			print ("Multiplicity line", commandargs[2])
			m2 = re.match(r'^(\d+)\s+(.*)$', commandargs[2])
			id = int(m2.group(1))
			print ("id:", id, "followed by:", m2.group(2))
			mu = multsplit(m2.group(2))
			print ("mu:", mu)
			if mu:
				mu.insert(0, id)
				L = MULTLine(commandargs, mu)
				return L
		else:
# any other LSD command
			print ("VSAResolver unused command line")
			L = CommandLine(*commandargs)
			return L
# should not happen
	return "?????"

def multsplit(ch):
	lre = r'^'
	sepre = r'\s+'
	eltre = r'([A-Z][a-z]?[1-6]{0,3})'
	hybre = r'([1-3]|\((?:[1-3]\s*){2,}\))'
	multre = r'([0-3]|\((?:[0-3]\s*){2,}\))'
	chargere = r'(\s+.*)?'
	rre = r'$'
	bigre = lre + eltre + sepre + hybre + sepre + multre + chargere + rre
	m = re.match(bigre, ch)
	if m:
# get MULT arguments
		print ("*"*20, [m.group(i) for i in range(1,4)])
# atom element is a string
		elt = m.group(1)
		print ("elt:", elt)
# get parts of the hyb information
		hyb = m.group(2).lstrip("(").rstrip(")").split()
		nhyb = len(hyb)
# integer or list of integer for hybridization state
		hyb = int(hyb[0]) if nhyb == 1 else [int(hyb[i]) for i in range(nhyb)]
		print ("hyb:", hyb)
# get parts of the mult information
		mult = m.group(3).lstrip("(").rstrip(")").split()
		nmult = len(mult)
# integer or list of integer for multiplicity
		mult = int(mult[0]) if nmult == 1 else [int(mult[i]) for i in range(nmult)]
		print ("mult:", mult)
		ec = m.group(4)
		if ec:
			ec = sorted([int(e) for e in ec.strip().lstrip("(").rstrip(")").split()])
		else:
			ec = [0]
		return [elt, hyb, mult, ec]
	else:
		return []

if __name__ == "__main__":
	line1 = "\n"
	L1 = decode(line1)
	print ("***")
	print (L1)
	print ("***")
	print ("for status resolver:", L1.for_resolver())
	print ("updatable:", L1.updatable())
	print ("from_resolver:", L1.from_resolver())
	# print "original: |%s|%s|" % (L1.original(), line1)

	print()
	line2 = "; SUBS 0\n"
	L2 = decode(line2)
	print ("***")
	print (L2)
	print ("***")
	print ("for status resolver:", L2.for_resolver())
	print ("updatable:", L2.updatable())
	print ("from_resolver:", L2.from_resolver())
	# print "original: |%s|%s|" % (L2.original(), line2)

	print()
	line3 = "FORM \"C 21 H 22 N 2 O 2\" ; molecular formula\n"
	L3 = decode(line3)
	print ("***")
	print (L3)
	print ("***")
	print ("for status resolver:", L3.for_resolver())
	print ("updatable:", L3.updatable())
	print ("from_resolver:", L3.from_resolver())
	# print "original: |%s|%s|" % (L3.original(), line3)

	print()
	line4 = "PIEC 1\n"
	L4 = decode(line4)
	print ("***")
	print (L4)
	print ("***")
	print ("for status resolver:", L4.for_resolver())
	print ("updatable:", L4.updatable())
	print ("from_resolver:", L4.from_resolver())
	# print "original: |%s|%s|" % (L4.original(), line4)

	print()
	line5 = "MULT 22 N (1 2 3) (0 1 2) (-1 0 1)\n"
	L5 = decode(line5)
	print ("***")
	print (L5)
	print ("***")
	print ("for status resolver:", L5.for_resolver())
	print ("updatable:", L5.updatable())
	print ("from_resolver:", L5.from_resolver())
	# print "original: |%s|%s|" % (L5.original(), line5)

	print()
	line6 = "SSTR S1 N 3 0\n"
	L6 = decode(line6)
	print ("***")
	print (L6)
	print ("***")
	print ("for status resolver:", L6.for_resolver())
	print ("updatable:", L6.updatable())
	print ("from_resolver:", L6.from_resolver())
	# print "original: |%s|%s|" % (L6.original(), line6)

	print ()
	line7 = "ELEC (-1 0 1 2)\n"
	L7 = decode(line7)
	print ("***")
	print (L7)
	print ("***")
	print ("for status resolver:", L7.for_resolver())
	print ("updatable:", L7.updatable())
	print ("from_resolver:", L7.from_resolver())
	# print "original: |%s|%s|" % (L7.original(), line7)

	print()
	line8 = "MAXN 1\n"
	L8 = decode(line8)
	print ("***")
	print (L8)
	print ("***")
	print ("for status resolver:", L8.for_resolver())
	print ("updatable:", L8.updatable())
	print ("from_resolver:", L8.from_resolver())
	# print "original: |%s|%s|" % (L8.original(), line8)

	print()
	line9 = "MOMA 200-300\n"
	L9 = decode(line9)
	print ("***")
	print (L9)
	print ("***")
	print ("for status resolver:", L9.for_resolver())
	print ("updatable:", L9.updatable())
	print ("from_resolver:", L9.from_resolver())
	# print "original: |%s|%s|" % (L9.original(), line9)

	print()
	line10 = "DEMU N N35 (1 2 3) (0 1 2) (-1 0 1)\n"
	L10 = decode(line10)
	print ("***")
	print (L10)
	print ("***")
	print ("for status resolver:", L10.for_resolver())
	print ("updatable:", L10.updatable())
	print ("from_resolver:", L10.from_resolver())
	# print "original: |%s|%s|" % (L10.original(), line10)

	print()
	line11 = "SHIX 5 -37.11\n"
	L11 = decode(line11)
	print ("***")
	print (L11)
	print ("***")
	print ("for status resolver:", L11.for_resolver())
	print ("updatable:", L11.updatable())
	print ("from_resolver:", L11.from_resolver())
	# print "original: |%s|%s|" % (L11.original(), line11)

	print()
	line12 = "COUF \"my_count_file\"\n"
	L12 = decode(line12)
	print ("***")
	print (L12)
	print ("***")
	print ("for status resolver:", L12.for_resolver())
	print ("updatable:", L12.updatable())
	print ("from_resolver:", L12.from_resolver())
	# print "original: |%s|%s|" % (L12.original(), line12)

	print()
	line13 = "STOF \"my_stop_file\"\n"
	L13 = decode(line13)
	print ("***")
	print (L13)
	print ("***")
	print ("stop file:", L13.stof)
	print ("for status resolver:", L13.for_resolver())
	print ("updatable:", L13.updatable())
	print ("from_resolver:", L13.from_resolver())
	# print "original: |%s|%s|" % (L13.original(), line13)

