"""Variable Statum A Resolver related code

Usage:
r = VSAresolver(pb)
for sts in r.stbyat:
	lines = r.linearray(sts)
	... do something with the array of LSD MULT lines
	
see examples at the end of this file for the structure of input data pb
"""
import status
import atom
import sumnn
import equiv
import box
import boxfiller
import sumnH
import sumEc
import expstats
from functools import reduce

class VSAresolver:
	"""Resolves status ambiguities of VSAs (Variable status atoms)"""
	def __init__(self, problem):
		"""Does the job, using information as coded in the example(s), see at bottom of file"""
		self.pb = problem
# opcode decoder
		self.selector = {}
		self.selector["FORM"] = self.FORM
		self.selector["MULT"] = self.MULT
		self.selector["PIEC"] = self.PIEC
		self.selector["ELEC"] = self.ELEC
		self.selector["MAXP"] = self.MAXP
		self.selector["MAXN"] = self.MAXN
		self.selector["SHIX"] = self.SHIX
# empty initial list of atoms 
		self.atoms = []
# empty initial list of shifts 
		self.shifts = []
# initial list of possible molecular electric charges
		self.elec = [0]
# initial required maximum number of positively charged atoms. -1 means no control on this number
		self.reqmaxPC = -1
# initial required maximum number of negatively charged atoms. -1 means no control on this number
		self.reqmaxNC = -1
# run the processing of VSAs
		self.run()
	
	def FORM(self, args):
		"""Decodes the molecular formula information item"""
#		print "FORM", args
		self.MF = {}
		self.numH = 0
		items = args[0].split()
		for i in range(0, len(items), 2):
			elt = items[i]
			num = int(items[i + 1])
			if elt == 'H':
				self.numH = num
			else:
				self.MF[elt] = num
	
	def MULT(self, args):
		"""Decodes the atom definition information items"""
#		print "MULT", args
		if len(args) == 4:
			args.append([0])
		id, elt, hyb, mult, ec = tuple(args)
		newatom = atom.Atom(id, elt, hyb, mult, ec)
		print (newatom)
		self.atoms.append(newatom)
	
	def PIEC(self, args):
		"""Decodes the connexity range information item"""
#		print "PIEC", args
		self.pieces = list(range(args[0], args[1]+1))

	def ELEC(self, args):
		"""Decodes the molecular electric charge information item"""
		print ("ELEC", args)
		self.elec = args[0]

	def MAXP(self, args):
		"""Decodes the required maximum number of positively charged atoms"""
		print ("MAXP", args)
		self.reqmaxPC = args[0]

	def MAXN(self, args):
		"""Decodes the required maximum number of negatively charged atoms"""
		print ("MAXN", args)
		self.reqmaxNC = args[0]

	def SHIX(self, args):
		"""Decodes the atom chemical shifts"""
		print ("SHIX", args)
		self.shifts.append([args[0], args[1]])

	def run(self):
		"""Runs the decoding of the current problem and runs the status generation process"""
# stbyat is the result of the creation of an VSAresolver. Filled in step5()
		self.stbyat = []
# scan input data
		for line in self.pb:
			comm = line[0]
			args = line[1:]
			self.selector[comm](args)
# deal with chemical shift
		for id, cshift in self.shifts:
			alist = [a for a in self.atoms if a.id == id]
			a = alist[0]
			a.chemshift = cshift
# have a look to atoms
		for a in self.atoms:
			a.dump()
# check parity consistency between molecular formula and molecular electric charge
		p = MFnumH_parity(self.MF, self.numH)
# retain in good_elec the molecular electric charge values that fit with the molecular formula
		good_elec = [e for e in self.elec if e % 2 == p]
		print ("valid molecular electric charges:", good_elec)
# is good_elec non-empty?
		if good_elec:
# modify self.elec
			self.elec = good_elec
		else:
# no solution to the problem
			return
# calculate the minium and the maximum number of H according to atom status
# the number of H atoms that is given by the molcular formula must
# be enclosed between these two limits. 
#		self.minH = reduce(lambda x, y: x+y, map(lambda a: a.minH, self.atoms))
#		self.maxH = reduce(lambda x, y: x+y, map(lambda a: a.maxH, self.atoms))

		self.minH = sum([a.minH for a in self.atoms])
		self.maxH = sum([a.maxH for a in self.atoms])
		print ("minH:", self.minH, "maxH:", self.maxH, "numH:", self.numH)
		ok = self.minH <= self.numH <= self.maxH
		print ("number of H atoms in the molecular formula is %s correct" % ("" if ok else "not"))
		if not ok:
			return
# calculate the minimum and the maximum molecular electric charge according to atom status.
# At least one of the declared molecular electric charge must
# be enclosed between these two limits. 
		self.minEc = sum([a.minEc for a in self.atoms])
		self.maxEc = sum([a.maxEc for a in self.atoms])
		print ("minEc:", self.minEc, "maxEc:", self.maxEc, "numEc:", self.elec)
		ok = bool([e for e in self.elec if self.minEc <= e <= self.maxEc])
		print ("molecular electric charge is%s correct" % ("" if ok else " not"))
		if not ok:
			return
# checks if the minimal number of positively charged atoms is compatible
# with the required maximum number of positively charged atoms
		if self.reqmaxPC != -1:
			minnumPC = len([a for a in self.atoms if a.mustbePC])
			print ("minnumPC:", minnumPC, "reqmaxPC", self.reqmaxPC)
			ok = minnumPC <= self.reqmaxPC
			print ("maximimal number of positively charged atoms is%s correct" % ("" if ok else " not"))
			if not ok:
				return
# checks if the minimal number of negatively charged atoms is compatible
# with the required maximum number of negatively charged atoms
		if self.reqmaxNC != -1:
			minnumNC = len([a for a in self.atoms if a.mustbeNC])
			print ("minnumNC:", minnumNC, "reqmaxNC", self.reqmaxNC)
			ok = minnumNC <= self.reqmaxNC
			print ("maximimal number of negatively charged atoms is%s correct" % ("" if ok else " not"))
			if not ok:
				return
# collect information about the non VSA part of the molecule
		self.invariant_part()
# collect information about the VSA part of the molecule
		self.variant_part()
		if self.numvariants != 0:
#			pass
# generate the possible solutions for the number of VSAs
# as a function of the number of neighbors
			self.step1()
# generate arrays of number of neighbors by VSA (in the order of VSA equivalence classes
			self.step2()
# expand each array of number of neighbors by VSA into arrays of status
# takes care of H atoms
			self.step3()
# expand each array of arrays of status into arrays of array of status
# takes care of electric charges
			self.step4()
# expand each array of arrays of status into arrays of array of status
# takes care remaining ambiguities on hybridization/valence
			self.step5()
		else:
			self.stbyat = [[]]

	def invariant_part(self):
		"""Computes the total number of H and of neighborhoods in the fixed part
		and number of atoms according to the number of neighbors"""
# list of fixed status atoms (FSAs)
#		self.invariants = filter(lambda x: x.isInvariant, self.atoms)
		self.invariants = [a for a in self.atoms if a.isInvariant]
#		print "Invariants:", map(lambda x: x.id, self.invariants)
		print ("Invariants:", [a.id for a in self.invariants])
# sn: sum of number of neighbors for FSAs only, sH: number of H on FSAs.
# sPC: number of positively charged atoms, sNC: number of negatively charged atoms
		(sn, sH, sE, sNC, sPC) = (0, 0, 0, 0, 0)
# initial list of the number of FSAs ([0]*5)
		nabynn = [0 for i in range(5)]
# initial number of sp2 FSAs 
		nsp2 = 0
# initial number of sp FSAs 
		nsp = 0
# initial number of terminal atoms, with hybridization as index
		nterm = [0]*4
# initial number of non terminal atoms, with hybridization as index
		nnonterm = [0]*4
# scan through FSAs
		for a in self.invariants:
# stats is always a list, with one element for FSAs
			st = a.stats[0]
# number of neighbors for status st
			nn = status.status_nn_one(st)
# number of H atoms for status st
			nH = status.status_nH_one(st)
# hybridization for status st
			hy = status.status_hyb_one(st)
# electric charge for status st
			e = status.status_Ec_one(st)
# update the sum of number of neighbors for FSAs
			sn += nn
# update the sum of number of H atoms for FSAs
			sH += nH
# update the sum of electric charges for FSAs
			sE += e
# update the number of positively charged FSAs
			sPC += 1 if a.mustbePC else 0
# update the number of negatively charged FSAs
			sNC += 1 if a.mustbeNC else 0
# update the table of number of FSAs by nn
			nabynn[nn] += 1
# update the number of sp2 FSAs
			if hy == 2:
				nsp2 += 1
# update the number of sp FSAs
			if hy == 1:
				nsp += 1
# update the number of terminal / non termainal atoms
			if nn == 1:
				nterm[hy] += 1
			else:
				nnonterm[hy] += 1
# update molecule
		self.invariant_nn = sn
		self.invariant_nH = sH
		self.invariant_elec = sE
		self.invariant_nPC = sPC
		self.invariant_nNC = sNC
		self.invariant_nabynn = nabynn
		self.invariant_nsp2 = nsp2
		self.invariant_nsp = nsp
		self.invariant_nterm = nterm
		self.invariant_nnonterm = nnonterm
# get number of H atoms on VSAs
		self.variantH = self.numH - sH
		print (sn, "fixed neighbors", sH, "fixed H atoms", sE, "fixed charges")
		print (self.variantH, "H atoms on VSAs")
		print ("fixed atoms by nn:", nabynn[1:])
		self.variant_elec = [e - self.invariant_elec for e in self.elec]
		print (self.variant_elec, "electric charges on VSAs")
# get max number of positively/negatively charged VSAs
		self.variant_reqmaxPC = self.reqmaxPC - self.invariant_nPC
		self.variant_reqmaxNC = self.reqmaxNC - self.invariant_nNC

	def variant_part(self):
		"""creates the list of variant atoms and count them"""
# list of VSAs
		self.variants = [a for a in self.atoms if a.isVariant]
# ids of variants
		self.variantids = [a.id for a in self.variants]
		print ("Variants:", self.variantids)
# number of VSAs
		self.numvariants = len(self.variants)
		
	def step1(self):
		"""Generate possible arrays of number of VSAs by nn"""
# Computes the min and max of the number of H and of electric charges
# according to the number of neighbors
# and the maximum number of VSA for each possible number of neighbors

# list of (array of minimum number of H by nn) for all VSAs
		mins = [a.minHbynn for a in self.variants]
# global minimum number of H atoms on VSAs by nn
		self.minminHbynn = reduce(lambda x, y: list(map(min, x, y)), mins)
		print ("minminH:", self.minminHbynn[1:])
# list of (array of maximum number of H by nn) for all VSAs
		maxs = [a.maxHbynn for a in self.variants]
# global maximum number of H atoms on VSAs by nn
		self.maxmaxHbynn = reduce(lambda x, y: list(map(max, x, y)), maxs)
		print ("maxmaxH:", self.maxmaxHbynn[1:])

# list of (array of minimum electric charge by nn) for all VSAs
		mins = [a.minEcbynn for a in self.variants]
# global minimum electric charge on VSAs by nn
		self.minminEcbynn = reduce(lambda x, y: list(map(min, x, y)), mins)
		print ("minminEc:", self.minminEcbynn[1:])
# list of (array of electric charges by nn) for all VSAs
		maxs = [a.maxEcbynn for a in self.variants]
# global maximum electric charge on VSAs by nn
		self.maxmaxEcbynn = reduce(lambda x, y: list(map(max, x, y)), maxs)
		print ("maxmaxEc:", self.maxmaxEcbynn[1:])

# and function for boolean scalars
		andscalfunc = lambda x, y: x and y
# elementwise and function for boolean arrays
		andvecfunc = lambda x, y: list(map(andscalfunc, x, y))
# list of (array of positive charge obligation) for all VSAs
		PCs = [a.mustbePCbynn for a in self.variants]
# obligation for VSAs of being positively charged, by nn
		self.mustbePCbynn = [1 if x else 0 for x in reduce(andvecfunc, PCs)]
		print ("mustbePCbynn:", self.mustbePCbynn[1:])
# list of (array of negative charge obligation) for all VSAs
		PCs = [a.mustbeNCbynn for a in self.variants]
# obligation for VSAs of being positively charged, by nn
		self.mustbeNCbynn = [1 if x else 0 for x in reduce(andvecfunc, PCs)]
		print ("mustbeNCbynn:", self.mustbeNCbynn[1:])

# concatenation des listes de nombre de voisins possibles pour tous les VSAs
		allnn = reduce(lambda x, y: x + y, [a.nn for a in self.variants], [])
# Obfuscated python!! A simpler one-liner is possible with len(reduce()) instead of sum(map()).
# for i in range(5) we get a table . so that .[i] is the number of times i is present in allnn
# This gives the maximum number of atoms with i neighbors.
		self.maxnabynn = list(map(lambda i: sum(list(map(lambda x: 1 if x==i else 0, allnn))) , range(5)))
		print ("maximum number of VSAs by number of neighbors:", self.maxnabynn[1:])

# get possible arrays of number of VSAs by nn according to max and parity constraints
		self.nabynn_crude = sumnn.sumnn(self.maxnabynn, self.numvariants, len(self.atoms), self.invariant_nn % 2)
# validates only those that fit with H information
		self.validateHEc()
		print()
		print ("nabynn_validatedHEc:", self.nabynn_validatedHEc)
# validates only those that fit with ring/pieces information
		self.validateRings()
		print()
		print ("nabynn_validatedHEcandRings:", self.nabynn_validatedHEcandRings)
		print ("Solution size:", len(self.nabynn_validatedHEcandRings))

	def validateHEc(self):
		"""transforms nabynn_crude into nabynn_validatedHEc"""
# number of H atoms on VSAs
		varH = self.variantH
		print()
		print ("Validation of nabynn_crude for variantH:", varH)
# array of possible electric charge on VSAs
		varEc = self.variant_elec
		print ("Validation of nabynn_crude for variant_elec:", varEc)
# number of positively charged atom
		if self.reqmaxPC != -1:
			maxPC = self.variant_reqmaxPC
			print ("Validation of nabynn_crude for required maximum positively charged atoms:", maxPC)
# number of negatively charged atom
		if self.reqmaxNC != -1:
			maxNC = self.variant_reqmaxNC
			print ("Validation of nabynn_crude for required maximum negatively charged atoms:", maxNC)
# initially empty list of results
		self.nabynn_validatedHEc = []
# iterates over crude solutions. A nice reduce() could do the job
		for l in self.nabynn_crude:
# minH is the minimum number of H on VSAs if each one beared its minimum number of H
			minH = prodscal(l, self.minminHbynn)
# minH is the maximum number of H on VSAs if each one beared its maximum number of H
			maxH = prodscal(l, self.maxmaxHbynn)
			print()
			print ("H validation for:", l, "minH:", minH, "maxH:", maxH)
# testing whether OK for H number or not
			if minH <= varH <= maxH:
				print (l, "is OK for the number of H")
# minEc is the minimum electric charge on VSAs if each one beared its minimum electric charge
				minEc = prodscal(l, self.minminEcbynn)
# maxEc is the maximum electric charge on VSAs if each one beared its maximum electric charge
				maxEc = prodscal(l, self.maxmaxEcbynn)
				print ("Ec validation for:", l, "minEc:", minEc, "maxEc:", maxEc)
# testing whether OK for electric charge or not
				if [e for e in varEc if minEc <= e <= maxEc]:
					print (l, "is OK for the molecular electric charge")
# numPC_OK/numNC_OK 
					numPC_OK = True
					if self.reqmaxPC != -1:
						minPC = prodscal(l, self.mustbePCbynn)
						numPC_OK = minPC <= maxPC
						print ("minimum number of positively charged atoms:", minPC, " maximum required:", maxPC)
					numNC_OK = True
					if self.reqmaxNC != -1:
						minNC = prodscal(l, self.mustbeNCbynn)
						numNC_OK = minNC <= maxNC
						print ("minimum number of negatively charged atoms:", minNC, " maximum required:", maxNC)
					if numPC_OK and numNC_OK:
						print (l, "is OK for the number of positively/negatively charged atoms")
# append tested solution to nabynn_validatedHEc
						self.nabynn_validatedHEc.append(l)
					else:
						print ("validation of the number of positively/negatively charged atoms failed")
				else:
					print ("validation of the molecular electric charge failed")
			else:
				print ("validation of the number of H failed")

	def validateRings(self):
		print ()
		print ("Validation of nabynn_validateHEc for invariant_nabynn:", self.invariant_nabynn)
# number of nodes of the molecular graph
		nodes = len(self.atoms)
		print ("nodes:", nodes)
# empty solution list
		self.nabynn_validatedHEcandRings = []
		for l in self.nabynn_validatedHEc:
# calculate the number of atoms by nn for both FSAs and VSAs
			fullnabynn = sumvecvec(self.invariant_nabynn, l)
			print()
			print ("l:", l, "fullnabynn:", fullnabynn)
# molecular number of neighbors
			vertices_x_2 = prodscal(fullnabynn, range(5))
# number of vertices in the molecular graph (all bonds simple)
			vertices = vertices_x_2 / 2
			print ("vertices_x_2:", vertices_x_2, "vertices:", vertices)
# possible number of rings according to vertices, nodes and connexivity (number of pieces)
			rings = sumscalvec(vertices - nodes, self.pieces)
# goodrings is the number of times the connexivity values
# give a positive or zero number of rings
			goodrings = sum(list(map(lambda x: 0 if x < 0 else 1, rings)))
			print ("rings:", rings, "goodrings:", goodrings)
# at least one valid number of rings is required
			if goodrings > 0:
				print (l, "is OK")
# one more valid solution
				self.nabynn_validatedHEcandRings.append(l)
			else:
				print (l, "is not OK")

	def step2(self):
		"""Transforms the lists of number of atoms according to nn in nabynn_validatedHEcandRings
		into lists of nn according to VSAs. VSAs are stored in the order
		the classer object find them and is the one of self.flatatomclasses"""

# nnbyatsOK is the result of step 2
		self.nnbyatsOK = []
# nothing to do if nabynn_validatedHEcandRings is empty
		if len(self.nabynn_validatedHEcandRings) == 0:
			return
# find classes of VSAs according to their initial status
		self.classer = equiv.Classer(self.variants)
# makes the list of atom classes (a class is a list of atoms)
		self.atomclasses = [x.members() for x in self.classer.classes()]
# makes the list of atom id classes. For debugging only.
		self.idclasses = [list(map(lambda a: a.id , c)) for c in self.atomclasses]
		print ("Classes of VSAs:", self.idclasses)
# transforms classes into boxes for VSA repartition
		self.boxes = [box.Box(len(c), c[0].nnFlags, True) for c in self.atomclasses]
		for i,b in enumerate(self.boxes): print ("class ", i, b)
# create a box filler 
		self.filler = boxfiller.BoxFiller(self.boxes)
# each element in list nnbyboxes is a list in which element is
# related to a group of equivalent VSAs. Each group is described
# by the list of number of neighbors. Initially empty.
		self.nnbyboxes = []
# transformation of the arrays of number of atoms by nn
# into array of number of atoms by groups
		for nabynn in self.nabynn_validatedHEcandRings:
			print ("nabynn:", nabynn)
# get one or more lists of lists of nn
			fillings = self.filler.fill(nabynn)
			print ("fillings:", fillings)
# aggregates the new ones to the old ones
			self.nnbyboxes += fillings
			print ("nnbyboxes:", self.nnbyboxes)
		print ("final nnbyboxes:", self.nnbyboxes)
# nothing to do if nnbyboxes is empty
		if len(self.nnbyboxes) == 0:
			return
# builds the list of VSAs, reordered by groups
		self.flatatomclasses = listlist2list(self.atomclasses)
		print ("flatatomclasses:", self.flatatomclasses)
# builds the list of VSA ids, reordered by groups
		self.flatidclasses = listlist2list(self.idclasses)
		print ("flatidclasses:", self.flatidclasses)
# minimum number of H by nn, for each VSA in the order of flatatomclasses
		self.flatminHbynn = list(map(lambda a: a.minHbynn, self.flatatomclasses))
		print ("flatminHbynn:", self.flatminHbynn)
# maximum number of H by nn, for each VSA in the order of flatatomclasses
		self.flatmaxHbynn = list(map(lambda a: a.maxHbynn, self.flatatomclasses))
		print ("flatmaxHbynn:", self.flatmaxHbynn)
# minimum electric charge by nn, for each VSA in the order of flatatomclasses
		self.flatminEcbynn = list(map(lambda a: a.minEcbynn, self.flatatomclasses))
		print ("flatminEcbynn:", self.flatminEcbynn)
# maximum electric charge by nn, for each VSA in the order of flatatomclasses
		self.flatmaxEcbynn = list(map(lambda a: a.maxEcbynn, self.flatatomclasses))
		print ("flatmaxEcbynn:", self.flatmaxEcbynn)
# obligation of being positively charged by nn, for each VSA in the order of flatatomclasses
		if self.reqmaxPC != -1:
			self.flatmustbePCbynn = list(map(lambda a: a.mustbePCbynn, self.flatatomclasses))
			print ("flatmustbePCbynn:", self.flatmustbePCbynn)
# obligation of being negatively charged by nn, for each VSA in the order of flatatomclasses
		if self.reqmaxNC != -1:
			self.flatmustbeNCbynn = list(map(lambda a: a.mustbeNCbynn, self.flatatomclasses))
			print ("flatmustbeNCbynn:", self.flatmustbeNCbynn)
# each element in nnbyats contains a list of nn, in the order of flatatomclasses
		self.nnbyats = list(map(lambda x: listlist2list(x), self.nnbyboxes))
		print()
		print ("nnbyats:", self.nnbyats)
# filter out the lists of nn that do not fit with the constraint of the number of H on VSAs
		self.nnbyatsOK = list(filter(lambda x: self.nnbyatisOK(x), self.nnbyats))
		print ()
		print ("nnbyatsOK:", self.nnbyatsOK)
		print ("Solution size:", len(self.nnbyatsOK))
	
	def nnbyatisOK(self, nnbyat):
		"""checks a list of nn of VSAs for the criterion of number of H and of electric charge"""
# prepare validation on the number of H atoms
		varH = self.variantH
# prepare validation on the number of positively/negatively charged atoms
		if self.reqmaxPC != -1:
			maxPC = self.variant_reqmaxPC
		if self.reqmaxNC != -1:
			maxNC = self.variant_reqmaxNC
# minH is the minimum possible number of H
		minH = sum(list(map(lambda x, y: y[x], nnbyat, self.flatminHbynn)))
# maxH is the maximum possible number of H
		maxH = sum(list(map(lambda x, y: y[x], nnbyat, self.flatmaxHbynn)))
# comparison with the expected number of H on VSAs
		isOK = minH <= varH <= maxH
		print()
		print ("nnbyatisOK(): nnbyat:", nnbyat, "minH:", minH, "maxH:", maxH, 
		"variantH:", varH, "isOK:", isOK)
		if isOK:
			varEc = self.variant_elec
# minEc is the minimum possible electric charge
			minEc = sum(list(map(lambda x, y: y[x], nnbyat, self.flatminEcbynn)))
# maxEc is the maximum possible electric charge
			maxEc = sum(list(map(lambda x, y: y[x], nnbyat, self.flatmaxEcbynn)))
# comparison with the expected molecular electric charge
			isOK = bool([e for e in varEc if minEc <= e <= maxEc])
			print ("nnbyatisOK(): nnbyat:", nnbyat, "minEc:", minEc, "maxEc:", maxEc, 
			"variant_elec:", varEc, "isOK:", isOK)
# constraint on the number of positively/negatively charged atoms
			if isOK:
				numPC_OK = True
				if self.reqmaxPC != -1:
					maxPC = self.variant_reqmaxPC
					minPC = sum(list(map(lambda x, y: y[x], nnbyat, self.flatmustbePCbynn)))
					numPC_OK = minPC <= maxPC
					print ("nnbyatisOK(): nnbyat:", nnbyat, "minPC:", minPC, "maxPC:", maxPC, 
						"numPC_OK:", numPC_OK)
				numNC_OK = True
				if self.reqmaxNC != -1:
					maxNC = self.variant_reqmaxNC
					minNC = sum(list(map(lambda x, y: y[x], nnbyat, self.flatmustbeNCbynn)))
					numNC_OK = minNC <= maxNC
					print ("nnbyatisOK(): nnbyat:", nnbyat, "minNC:", minNC, "maxNC:", maxNC, 
						"numNC_OK:", numNC_OK)
				isOK = numPC_OK and numNC_OK
		return isOK

	def step3(self):
		"""expand each array of number of neighbors by VSA into arrays of arrays of status"""
# stsbyat contains lists of atom status lists by VSA, each in the order of variants.
# stsbyat is the result of step3
		self.stsbyat = []
# nothing to do if nnbyatsOK is empty
		if len(self.nnbyatsOK) == 0:
			return
		for nnbyat in self.nnbyatsOK:
			self.expandnn(nnbyat)
		print()
		print ("End of step 3: status sets for variants: ", self.stsbyat)
		print ("Solution size:", len(self.stsbyat))

	def expandnn(self, nnbyat):
		"""expand an array of number of neighbors by VSA into arrays of status"""
# !!!!!!!!! nnbyat is in the order of flatatomclasses
		print()
		print ("Expanding nn:", nnbyat, "for atoms:", self.flatidclasses)
		for a, nn in zip(self.flatatomclasses, nnbyat):
# each atom gets a new list of status that takes into account the current nn
# the initial list of all possible status is still in a.stats0
			a.stats = a.statbynn[nn]
# possible numbers of H for atom a considering its has nn neighbors
#			a.nHnns = sorted(list(set(status.status_nH(a.stats))))
			a.nHnns = a.nHbynn[nn]
			a.nnexact = nn
		variant_ordered_nnbyat = [a.nnexact for a in self.variants]
		print ("Expanding nn:", variant_ordered_nnbyat, "for atoms:", self.variantids)
# find classes of VSAs taking into account that each VSA has a defined number of neighbors
		self.classer = equiv.Classer(self.variants)
# makes the list of atom classes (a class is a list of atoms)
		self.atomclassesnn = [x.members() for x in self.classer.classes()]
# makes the list of atom id classes. For debugging only.
		self.idclassesnn = [list(map(lambda a: a.id , c)) for c in self.atomclassesnn]
		print ("Classes of VSAs:", self.idclassesnn)
# make the list of atom indexes within classes
		self.indexclassesnn = list(map(lambda x: list(range(len(x))), self.atomclassesnn))
# flat list of atoms, in the order of symmetry classes
		self.flatatomclassesnn = listlist2list(self.atomclassesnn)
# flat list of atom indexes within classes, in the order of flatatomclassesnn
		self.flatindexclassesnn = listlist2list(self.indexclassesnn)
		print ("Index within classes:", self.flatindexclassesnn)
# flat list of lists of possible nH, in the order of flatatomclassesnn
		self.flatnHnn = list(map(lambda a: a.nHnns, self.flatatomclassesnn))
		print ("Possible nH:", self.flatnHnn)
# list of flat lists (in the order of flatatomclassesnn) of possible nH
		self.nHnnbyats = sumnH.sumnH(self.flatnHnn, self.flatindexclassesnn, self.variantH)
		print ("Possiblities for exact total nH:", self.nHnnbyats)
# prepare validation on the number of positively/negatively charged atoms
		if self.reqmaxPC != -1:
			maxPC = self.variant_reqmaxPC
		if self.reqmaxNC != -1:
			maxNC = self.variant_reqmaxNC
# storing back nH in atoms
		for nHnnbyat in self.nHnnbyats:
			for a, nH in zip(self.flatatomclassesnn, nHnnbyat):
				a.nHexact = nH
			self.nHbyat = list(map(lambda a: a.nHexact, self.variants))
			print ("Combine nn:", variant_ordered_nnbyat, "and nH:", self.nHbyat)
			comb = list(map(lambda a, nn, nH: a.statbynnbynH[nn][nH], self.variants, variant_ordered_nnbyat, self.nHbyat))
			print ("Combined status list:" , comb)
			lminEc = list(map(lambda a, nn, nH: a.minEcbynnbynH[nn][nH], self.variants, variant_ordered_nnbyat, self.nHbyat))
			print ("Corresponding minimum electric charge:", lminEc)
			lmaxEc = list(map(lambda a, nn, nH: a.maxEcbynnbynH[nn][nH], self.variants, variant_ordered_nnbyat, self.nHbyat))
			print ("Corresponding maximum electric charge:", lmaxEc)
			minEc = sum(lminEc)
			print ("minimum charge on VSAs:", minEc)
			maxEc = sum(lmaxEc)
			print ("maximum charge on VSAs:", maxEc)
			isOK = [e for e in self.variant_elec if minEc <= e <= maxEc]
			if isOK:
				numPC_OK = True
				if self.reqmaxPC != -1:
					lmustbePC = list(map(lambda a, nn, nH: a.mustbePCbynnbynH[nn][nH], self.variants, variant_ordered_nnbyat, self.nHbyat))
					minPC = sum([1 if x else 0 for x in lmustbePC])
					numPC_OK = minPC <= maxPC
					print ("minPC:", minPC, "maxPC:", maxPC, "numPC_OK:", numPC_OK)
				numNC_OK = True
				if self.reqmaxNC != -1:
					lmustbeNC = list(map(lambda a, nn, nH: a.mustbeNCbynnbynH[nn][nH], self.variants, variant_ordered_nnbyat, self.nHbyat))
					minNC = sum([1 if x else 0 for x in lmustbeNC])
					numNC_OK = minNC <= maxNC
					print ("minPC:", minPC, "maxPC:", maxPC, "numPC_OK:", numPC_OK)
				isOK = numPC_OK and numNC_OK
				if isOK:
# the combination provides a list of status lists that is collected in stsbyat
					self.stsbyat.append(comb)
				else:
					print ("validation of the number of positively/negatively charged atoms failed")
			else:
				print ("validation of the molecular electric charge failed")

	def step4(self):
#		pass
		"""driver for expandec (if needed)"""
# initialize list of status lists for VSAs. ecstsbyat is the result of step 4.
		self.ecstsbyat = []
# nothing to do if stsbyat is empty
		if len(self.stsbyat) == 0:
			return
# scan through the arrays of arrays of status of VSAs from step 3
		for sts in self.stsbyat:
			print()
			print ("Charge expansion for status:", sts)
# each element in sts corresponds to an atom with a given nn and a given nH
			sts0 = [x[0] for x in sts]
			self.nns = status.status_nn(sts0)
			print ("number of neighbors:", self.nns)
			self.nHs = status.status_nH(sts0)
			print ("number of hydrogens:", self.nHs)
			Ecs = [a.EcbynnbynH[nn][nH] for a, nn, nH in zip(self.variants, self.nns, self.nHs)]
			print ("possible charges:", Ecs)
# is there a variant atom with a list of status with more than one element?
			if len(Ecs) != sum([len(x) for x in Ecs]):
# insert status in atoms
				for a, s, nn, nH, Ec in zip(self.variants, sts, self.nns, self.nHs, Ecs):
					a.stats = s
					a.nnexact = nn
					a.nHexact = nH
					a.EcnHnns = Ec
# expand charge lists that are not singletons
				self.expandec(Ecs)
			else:
# Ecs is made of singletons only. No electric charge alternative to expand.
# the molecular electric charge and the number of positively/negatively
# charged atoms have already been approved in step 3
				self.ecstsbyat.append(sts)
				print ("No charge expansion required")
		print()
		print ("End of step 4: status sets for variants: ", self.ecstsbyat)
		print ("Solution size:", len(self.ecstsbyat))

	def expandec(self, Ecs):
		"""expand an array of electric charge arrays"""
# discriminate betweens atoms with already known charge and those with unknown charge
		pairs = zip(Ecs, range(len(Ecs)))
# uIdx are the atom indexes (in the order of variants) for which charge is still unkown (uEcs)
		upairs = [x for x in pairs if len(x[0]) > 1]
		uEcs = [x[0] for x in upairs]
		uIdx = [x[1] for x in upairs]
		print ("unknown electric charges:", uEcs, "at indexes:", uIdx)
# kIdx are the atom indexes (in the order of variants) for which charge is already known (kEcs)
		kpairs = [x for x in pairs if len(x[0]) == 1]
		kEcs = [x[0] for x in kpairs]
		kIdx = [x[1] for x in kpairs]
		print ("known electric charges:", kEcs, "at indexes:", kIdx)
# kEcs0 is like kEcs but with each singleton replaced by its element
		kEcs0 = [x[0] for x in kEcs]
		print ("known electric charges:", kEcs0)
# get the electric charge on atoms with known electric charge
		ksumEc = sum(kEcs0)
		print ("electric charge on atoms with known electric charge:", ksumEc)
# possible total electric charge on atoms with unknown charge
		usumEc = [x - ksumEc for x in self.variant_elec]
		print ("possible total electric charge on atoms with unknown charge:", usumEc)
# get the maximum number of unknown positively/negatively charged atoms
		if self.reqmaxPC != -1:
			knumPC = len(list(filter(lambda x: x > 0, kEcs0)))
			umaxPC = self.variant_reqmaxPC - knumPC
			print ("maximum number of unknown positively charged atoms:", umaxPC)
		else:
			umaxPC = -1
		if self.reqmaxNC != -1:
			knumNC = len(list(filter(lambda x: x < 0, kEcs0)))
			umaxNC = self.variant_reqmaxNC - knumNC
			print ("maximum number of unknown negatively charged atoms:", umaxNC)
		else:
			umaxNC = -1
# atoms of unknown charge
		uatoms = [self.variants[i] for i in uIdx]
		print ("atoms with unknown electric charge:", uatoms)
# atoms of known charge
		katoms = [self.variants[i] for i in kIdx]
		print ("atoms with known electric charge:", katoms)
# setting Ecexact for atoms of known charge
		for a in katoms:
			a.Ecexact = a.EcnHnns[0]

		print ("Buiding equivalence classes of uatoms")
# find equivalence classes of uatoms
		self.classer = equiv.Classer(uatoms)
# makes the list of atom classes (a class is a list of atoms)
		self.atomclassesec = [x.members() for x in self.classer.classes()]
# makes the list of atom id classes. For debugging only.
		self.idclassesec = [list(map(lambda a: a.id , c)) for c in self.atomclassesec]
		print ("Classes of uatoms (Atom Ids):", self.idclassesec)
# make the list of atom indexes within classes
		self.indexclassesec = list(map(lambda x: list(range(len(x))), self.atomclassesec))
# flat list of atoms, in the order of symmetry classes
		self.flatatomclassesec = listlist2list(self.atomclassesec)
# flat list of atom indexes within classes, in the order of flatatomclassesnn
		self.flatindexclassesec = listlist2list(self.indexclassesec)
		print ("Index within classes:", self.flatindexclassesec)
# flat list of lists of possible Ec, in the order of flatatomclassesnn
		self.flatEcnHnn = list(map(lambda a: a.EcnHnns, self.flatatomclassesec))
		print ("Possible Ec:", self.flatEcnHnn)
# list of flat lists (in the order of flatatomclassesec) of possible Ec
		self.EcnHnnbyats = sumEc.sumnEc(self.flatEcnHnn, self.flatindexclassesec, usumEc, umaxPC, umaxNC)
		print ("Possiblities for exact total Ec:", self.EcnHnnbyats)
		for EcnHnnbyat in self.EcnHnnbyats:
			for a, Ec in zip(self.flatatomclassesec, EcnHnnbyat):
				a.Ecexact = Ec
			self.Ecs = list(map(lambda a: a.Ecexact, self.variants))
			print ("Combine nn:", self.nns, "and nH:", self.nHs, "and Ec:", self.Ecs)
			comb = list(map(lambda a, nn, nH, Ec: a.statbynnbynHbyEc[nn][nH][Ec + status.chargeOffset], self.variants, self.nns, self.nHs, self.Ecs))
			print ("Combined status list:" , comb)
			self.ecstsbyat.append(comb)


	def step5(self):
		"""driver for expandstats (if needed)"""
# initialize list of status lists for VSAs. stbyat is the result of step 5.
# Already done in run()
#		self.stbyat = []
# nothing to do if ecstsbyat is empty
		if len(self.ecstsbyat) == 0:
			return
# scan through the arrays of arrays of status of VSAs from step 4
		for sts in self.ecstsbyat:
# is there a variant atom with a list of status with more than one element?
			if len(sts) != sum(list(map(lambda x: len(x), sts))):
# the list of status must be expanded
				self.expandstats(sts)
			else:
# sts is made f singletons only. make the list of status in the order of VSAs
				st = listlist2list(sts)
# validate status for hybridization and connexivity
				if self.finalvalid(st):
# a new array of status was found for the VSAs
					self.stbyat.append(st)
		print()
		print ("End of step 5: status arrays for variants: ", self.stbyat)
		print ("Solution size:", len(self.stbyat))

	def expandstats(self, sts):
		"""expand sts, a flat array of arrays of status (in self.ecstsbyat)
		into an array of flat arrays of status
		"""
# store status in atoms, alternatives are still possible.
		for a, st in zip(self.variants, sts):
			a.stats = st
# discriminate betweens atoms with already known status and those with still unknown status
		pairs = zip(sts, list(range(len(sts))))
# uIdx are the atom indexes (in the order of variants) for which status is still unkown (usts)
		upairs = [x for x in pairs if len(x[0]) > 1]
		usts = [x[0] for x in upairs]
		uIdx = [x[1] for x in upairs]
		print()
		print ("unknown status:", usts, "at indexes:", uIdx)
# kIdx are the atom indexes (in the order of variants) for which status is already known (ksts)
		kpairs = [x for x in pairs if len(x[0]) == 1]
		ksts = [x[0] for x in kpairs]
		kIdx = [x[1] for x in kpairs]
		print ("known status:", ksts, "at indexes:", kIdx)
# atoms of unknown status
		uatoms = [self.variants[i] for i in uIdx]
		print ("atoms with unknown status:", uatoms)
# atoms of known status
		katoms = [self.variants[i] for i in kIdx]
		print ("atoms with known status:", katoms)
# setting statexact for atoms of known status
		for a in katoms:
			a.statexact = a.stats[0]
# find classes of VSAs taking into account that each VSA has a defined array of status
		self.classer = equiv.Classer(uatoms)
# makes the list of atom classes (a class is a list of atoms)
		self.atomclassesstat = [x.members() for x in self.classer.classes()]
# makes the list of atom id classes. For debugging only.
		self.idclassesstat = [list(map(lambda a: a.id , c)) for c in self.atomclassesstat]
		print()
		print ("Classes of VSAs:", self.idclassesstat)
# make the list of atom indexes within classes
		self.indexclassesstat = list(map(lambda x: list(range(len(x))), self.atomclassesstat))
# flat list of atoms, in the order of symmetry classes
		self.flatatomclassesstat = listlist2list(self.atomclassesstat)
# flat list of atom indexes within classes, in the order of flatatomclassesstat
		self.flatindexclassesstat = listlist2list(self.indexclassesstat)
		print ("Index within classes:", self.flatindexclassesstat)
# flat list of lists of possible status, in the order of flatatomclassesstat
		self.flatstats = list(map(lambda a: a.stats, self.flatatomclassesstat))
# list of flat lists (in the order of flatatomclassesstat) of possible status
		self.statbyats = expstats.expstats(self.flatstats, self.flatindexclassesstat)
		print ("Possiblities for status:", self.statbyats)
# storing back stat in atoms
		for statusbyat in self.statbyats:
			for a, st in zip(self.flatatomclassesstat, statusbyat):
				a.statexact = st
			self.statbyat = list(map(lambda a: a.statexact, self.variants))
# validate status for hybridization and connexivity
			if self.finalvalid(self.statbyat):
# a new array of status was found for the VSAs
				self.stbyat.append(self.statbyat)

	def finalvalid(self, st):
		"""Validation of st, an array of status, for hybridization and connexivity"""
		print()
		print ("Validating:", st)
# hybridization, in the order of variants
		hybs = status.status_hybr(st)
		print ("Hybridizations:", hybs)
# total number of sp2
		nsp2 = len(list(filter(lambda h: h == 2, hybs))) + self.invariant_nsp2
# total number of sp
		nsp = len(list(filter(lambda h: h == 1, hybs))) + self.invariant_nsp
		print ("sp2:", nsp2, "sp:", nsp)
# the total number of sp2 atoms is always even
		if nsp2 % 2 != 0:
			print ("Removed by bad number of sp2 atoms:", nsp2)
			return False
# the total number of sp atoms is always even if the number of sp2 atoms is zero
		if nsp2 == 0 and nsp % 2 != 0:
			print ("Removed by bad number of sp atoms:", nsp)
			return False
# number of neighbors, in the order of variants
		nns = list(map(status.status_nn_one, st))
# number of variant terminal atoms, by hybridization state
		nterm = list(map(lambda h: len(list(filter(lambda i: nns[i]==1 and hybs[i]==h, range(len(st))))), range(4)))
# number of variant non terminal atoms, by hybridization state
		nnonterm = list(map(lambda h: len(list(filter(lambda i: nns[i]>1 and hybs[i]==h, range(len(st))))), range(4)))
# total number of terminal atoms, by hybridization state
		nterm = list(map(lambda x, y: x+y, nterm, self.invariant_nterm))
# total number of non terminal atoms, by hybridization state
		nnonterm = list(map(lambda x, y: x+y, nnonterm, self.invariant_nnonterm))
		print ("number of terminal atoms, by hyb:", nterm)
		print ("number of non terminal atoms, by hyb:", nnonterm)
# if the number of sp atoms is zero, each pair of terminal sp2 that exceeds the non terminal sp2
# gives rise to a biatomic small molecule
		if nsp == 0:
			nsp2t = nterm[2]
			nsp2nt = nnonterm[2]
# excess of terminal sp2 atoms
			if nsp2t > nsp2nt:
				nfrag = 1 + (nsp2t - nsp2nt)/2
				if nfrag > max(self.pieces):
					print ("Removed by number of fragments: at least", nfrag)
					return False
		print ("OK")
		return True

	def __repr__(self):
		"""returns the string representation of the VSAresolver result"""
		repr = ""
		for sts in self.stbyat:
			bloc = self.linearray(sts)
			repr += "\n"
			for line in bloc:
				repr += line + "\n"
		return repr

	def linearray(self, sts):
		"""transforms an array of status in the order of variants
		in an array of lsd printable MULT lines, in the order of atoms in pb.
		This function is intended for external use.
		"""
		bloc = []
		for a, st in zip(self.variants, sts):
			a.elt, a.nn, a.mult, a.hyb, a.ec, a.valence, a.definit = status.status_all(st)
			a.lsdline = a.printout()
		for a in self.atoms:
			bloc.append(a.lsdline)
		return bloc

def prodscal(u, v):
	"""scalar product of two vectors"""
	return reduce(lambda xx, yy: xx + yy, list(map(lambda x, y: x * y, u, v)))

def sumvecvec(u, v):
	"""vector sum of two vectors"""
	return list(map(lambda x, y: x + y, u, v))

def sumscalvec(x, v):
	"""adds x to all the elements in v"""
	return sumvecvec(v, [x]*len(v))

def listlist2list(l):
	"""transforms l, a list of lists, into a list by concatenation of the lists"""
	return reduce(lambda x, y: x + y, l, [])

def MFnumH_parity(MF, numH):
	"""calculates the parity of the total number of valences of a molecule,
	ignoring its possible global electric charge"""
	parity = numH
	for elt, nelt in MF.items():
		parity += nelt * status.elt_valence_parity[elt]
	parity = parity % 2
	return parity


if __name__ == "__main__":

# Glycerol carbonate
	pb1 = [
["FORM", "C 4 H 6 O 4"],
["PIEC", 1, 1],
["MULT", 1, "C", 2, 0],
["MULT", 2, "C", 3, 1],
["MULT", 3, "C", 3, 2],
["MULT", 4, "C", 3, 2],
["MULT", 5, "O", [2, 3], [0, 1]],
["MULT", 6, "O", [2, 3], [0, 1]],
["MULT", 7, "O", [2, 3], [0, 1]],
["MULT", 8, "O", [2, 3], [0, 1]]
	]

# C4, 4 carbons atoms only, not the explosive.
	pb2 = [
["FORM", "C 4"],
["PIEC", 1, 1],
["MULT", 1, "C", [1, 2, 3], [0, 1, 2, 3]],
["MULT", 2, "C", [1, 2, 3], [0, 1, 2, 3]],
["MULT", 3, "C", [1, 2, 3], [0, 1, 2, 3]],
["MULT", 4, "C", [1, 2, 3], [0, 1, 2, 3]]
	]

# Asparagine
	pb3 = [
["FORM", "C 4 H 8 N 2 O 3"],
["PIEC", 1, 1],
["MULT", 1, "C", 2, 0],
["MULT", 2, "C", 3, 1],
["MULT", 3, "C", 3, 2],
["MULT", 4, "C", 2, 0],
["MULT", 5, "O", [2, 3], [0, 1]],
["MULT", 6, "O", [2, 3], [0, 1]],
["MULT", 7, "O", [2, 3], [0, 1]],
["MULT", 8, "N", [1, 2, 3], [0, 1, 2]],
["MULT", 9, "N", [1, 2, 3], [0, 1, 2]]
	]

# benzene
	pb4 = [
["FORM", "C 6 H 6"],
["PIEC", 1, 1],
["MULT", 1, "C", [1, 2, 3], [0, 1, 2, 3]],
["MULT", 2, "C", [1, 2, 3], [0, 1, 2, 3]],
["MULT", 3, "C", [1, 2, 3], [0, 1, 2, 3]],
["MULT", 4, "C", [1, 2, 3], [0, 1, 2, 3]],
["MULT", 5, "C", [1, 2, 3], [0, 1, 2, 3]],
["MULT", 6, "C", [1, 2, 3], [0, 1, 2, 3]]
	]

# real benzene
	pb5 = [
["FORM", "C 6 H 6"],
["PIEC", 1, 1],
["MULT", 1, "C", 2, 1],
["MULT", 2, "C", 2, 1],
["MULT", 3, "C", 2, 1],
["MULT", 4, "C", 2, 1],
["MULT", 5, "C", 2, 1],
["MULT", 6, "C", 2, 1]
	]

# nitromethane
	pb6 = [
["FORM", "C 1 H 3 N 1 O 2"],
["PIEC", 1, 1],
["ELEC", [0]],
["MULT", 1, "C", 3, 3],
["MULT", 2, "N", [2, 3], [0, 1, 2], [0, 1]],
["MULT", 3, "O", [2, 3], [0, 1], [-1, 0]],
["MULT", 4, "O", [2, 3], [0, 1], [-1, 0]]
	]

# dinitromethane
	pb7 = [
["FORM", "C 1 H 2 N 2 O 4"],
["PIEC", 1, 1],
["ELEC", [-1, 0, 1]],
["MAXP", 0],
["MAXN", 0],
["MULT", 1, "C", 3, 2],
["MULT", 2, "N35", [1, 2, 3], [0, 1, 2], [0, 1]],
["MULT", 3, "N35", [1, 2, 3], [0, 1, 2], [0, 1]],
["MULT", 4, "O", [2, 3], [0, 1], [-1, 0]],
["MULT", 5, "O", [2, 3], [0, 1], [-1, 0]],
["MULT", 6, "O", [2, 3], [0, 1], [-1, 0]],
["MULT", 7, "O", [2, 3], [0, 1], [-1, 0]],
["SHIX", 1, 75.2]
	]

	r = VSAresolver(pb7)
#	print r
