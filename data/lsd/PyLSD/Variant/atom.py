"""Single atom related code"""

import status

class Atom:
	"""Atom description"""
	def __init__(self, id, elt, hyb, mult, ec=[0]):
		"""Build atom from status and calculate neigborhood related properties"""
		self.id = id
		self.elt = elt
		self.hyb = hyb
		self.mult = mult
		self.ec = ec
		self.chemshift = None
		
		self.Ecexact = 0
		self.EcnHnns = 0
		self.statexact = 0

#		print elt, mult, hyb, charge
# stats: list of all status indexes. no duplicate possible
		self.stats = status.status(elt, mult, hyb, ec)
# stats0: list of all status indexes, backup copy. stats may be overwritten.
#		self.stats0 = self.stats
#		print self.stats
# numstats: number of status indexes
		self.numstats = len(self.stats)
# isVariant: True if the atom has more than a single possible status
		self.isVariant = self.numstats > 1
# isInvariant: True if the atom has a single possible status
		self.isInvariant = not self.isVariant
# an invariant may have variants in its definition. now fixing the definition
		if self.isInvariant:
			self.elt, self.nn, self.mult, self.hyb, self.ec, self.valence, self.definit = status.status_all(self.stats[0])
# lsdline: ready-to-print line for LSD code generation
			self.lsdline = self.printout()
# min and max number of H
			self.minH = self.mult
			self.maxH = self.mult
# min and max electric charge
			self.minEc = self.ec
			self.maxEc = self.ec
# mustbePC : is positively charged
			self.mustbePC = self.ec > 0
# mustbeNC : must be (or is) negatively charged
			self.mustbeNC = self.ec < 0
# singl_nn : True is only a single number of neighbors is possible
#		self.single_nn = len(self.nn) > 1
# multiple_nn : True is more than a single nn is possible
#		self.multiple_nn = not self.single_nn
		
# for VSA only
		if self.isVariant:
# nn_full: list of possible number of neighbors, in the same order as in stats. duplicates possible
			self.nn_full = status.status_nn(self.stats)
# nn: like nn_full, but without duplicates
			self.nn = sorted(list(set(self.nn_full)))
# list of possible nH for a given nn
			self.nHbynn = [[] for nn in range(5)]
# list of possible status for any given nn
			self.statbynn = [[] for nn in range(5)]
# list of possible status for any given nn and any given nH
			self.statbynnbynH = [[] for nn in range(5)]
# list of possible status for any given nn and any given nH and aby given electric charge
			self.statbynnbynHbyEc = [[] for nn in range(5)]
# minimum number of H for any given nn
			self.minHbynn = [0]*5
# maximum number of H for any given nn
			self.maxHbynn = [0]*5
# nH_full: list of possible number of hydrogens, in the same order as in stats. duplicates possible
			self.nH_full = status.status_nH(self.stats)
# nH: list of possible number of hydrogens, no duplicate. should be like mult, but sorted.
			self.nH = sorted(list(set(self.nH_full)))
# min and max number of H
			self.minH = min(self.nH)
			self.maxH = max(self.nH)
# Ec_full: list of possible electric charge, in the same order as in stats. duplicates possible
			self.Ec_full = status.status_Ec(self.stats)
# Ec: list of possible electric charges, no duplicate. should be like charge, but sorted.
			self.Ec = sorted(list(set(self.Ec_full)))
# min and max electric charge
			self.minEc = min(self.Ec)
			self.maxEc = max(self.Ec)
# mustbePC : must be positively charged
			self.mustbePC = self.minEc > 0
# mustbeNC : must be negatively charged
			self.mustbeNC = self.maxEc < 0
# list of possible Ec for a given nn
			self.Ecbynn = [[] for nn in range(5)]
# minimum electric charge for any given nn
			self.minEcbynn = [0]*5
# maximum electric charge for any given nn
			self.maxEcbynn = [0]*5
# must be positively charged, for any given nn
			self.mustbePCbynn = [False]*5
# must be negatively charged, for any given nn
			self.mustbeNCbynn = [False]*5
# list of possible electric charge for any given nn and any given nH
			self.EcbynnbynH = [[] for nn in range(5)]
# minimum electric charge for any given nn and any given nH
			self.minEcbynnbynH = [[] for nn in range(5)]
# maximum electric charge for any given nn and any given nH
			self.maxEcbynnbynH = [[] for nn in range(5)]
# obligation of being positively charged for any given nn and any given nH
			self.mustbePCbynnbynH =  [[] for nn in range(5)]
# obligation of being negatively charged for any given nn and any given nH
			self.mustbeNCbynnbynH =  [[] for nn in range(5)]
# loop through atoms status
			for i in range(self.numstats):
# current status
				st = self.stats[i]
# ... corresponding nn
				nn = self.nn_full[i]
# ... corresponding nH
				nH = self.nH_full[i]
# ... corresponding Ec
				Ec = self.Ec_full[i]
# update nHbynn. duplicates possible
				self.nHbynn[nn].append(nH)
# update Ecbynn. duplicates possible
				self.Ecbynn[nn].append(Ec)
# update statbynn. no duplicate possible
				self.statbynn[nn].append(st)
# initialize statbynnbynH[nn]
#				self.statbynnbynH[nn] = [[] for i in range(4)]

# scan through possible number of neighbors
			for nn in self.nn:
# initialize statbynnbynH[nn]
				self.statbynnbynH[nn] = [[] for i in range(4)]
# initialize statbynnbynHbyEc[nn]
				self.statbynnbynHbyEc[nn] = [[] for i in range(4)]
# initialize EcbynnbynH[nn]
				self.EcbynnbynH[nn] = [[] for i in range(4)]
# initialize minEcbynnbynH[nn]
				self.minEcbynnbynH[nn] = [0] * 4
# initialize maxEcbynnbynH[nn]
				self.maxEcbynnbynH[nn] = [0] * 4
# initialize mustbePCbynnbynH[nn]
				self.mustbePCbynnbynH[nn] = [False] * 4
# initialize mustbeNCbynnbynH[nn]
				self.mustbeNCbynnbynH[nn] = [False] * 4
# eliminate duplicates and sort nHbynn
				l = self.nHbynn[nn]
				l = sorted(list(set(l)))
				self.nHbynn[nn] = l
# if for the current nn there are possible number of H...should always be true
				if l:
# determine the minimum and maximum number of H
					self.minHbynn[nn] = min(l)
					self.maxHbynn[nn] = max(l)
# eliminate duplicates and sort Ecbynn
				l = self.Ecbynn[nn]
				l = sorted(list(set(l)))
				self.Ecbynn[nn] = l
# if for the current nn there are possible electric charges...should always be true
				if l:
# determine the minimum and maximum electric charge
					minl = min(l)
					maxl = max(l)
					self.minEcbynn[nn] = minl
					self.maxEcbynn[nn] = maxl
# determine if the atom must be positively/negatively charged
					self.mustbePCbynn[nn] = minl > 0
					self.mustbeNCbynn[nn] = maxl < 0
# l is now the list of status for the current nn
				l = self.statbynn[nn]
# if l is not empty...should always be true
				if l:
# store a sorted list of status
					l = sorted(l)
					self.statbynn[nn] = l
# makes a set of possible status
					setl = set(l)
# scan through the possible nH of this atom
#					for nH in self.nH:
					for nH in self.nHbynn[nn]:
# get the set of all status that fit with this nH
						lstatnH = status.status_id(status.status_mult, nH)
# determine the intersection of the possible status of this atom when it has the current nn
# with the list of possible status that fit with the current nH.
# Store the corresponding sorted list in statbynnbynH
						lstatnnnH = sorted(list(setl & set(lstatnH)))
						self.statbynnbynH[nn][nH] = lstatnnnH
# initialize statbynnbynHbyEc[nn][nH]
						self.statbynnbynHbyEc[nn][nH] = [[] for i in range(status.numChargeIndexes)]
# determine the list of possible electric charges for the nn and the nH of the current atom
						lEcnnnH = sorted(list(set(status.status_Ec(lstatnnnH))))
						self.EcbynnbynH[nn][nH] = lEcnnnH
# min and max electric charge according to number of neighbors and number of H
						minl = min(lEcnnnH)
						maxl = max(lEcnnnH)
						self.minEcbynnbynH[nn][nH] = minl
						self.maxEcbynnbynH[nn][nH] = maxl
# obligation of being positively/negatively charged according to number of neighbors and number of H
						self.mustbePCbynnbynH[nn][nH] = minl > 0
						self.mustbeNCbynnbynH[nn][nH] = maxl < 0
						setll = set(self.statbynnbynH[nn][nH])
						for Ec in self.EcbynnbynH[nn][nH]:
							lstatEc = status.status_id(status.status_charge, Ec + status.chargeOffset)
							lstatnnnHEc = sorted(list(setll & set(lstatEc)))
							self.statbynnbynHbyEc[nn][nH][Ec + status.chargeOffset] = lstatnnnHEc

# nnFlags indicates whether a given nn value is possible or not for the current atom.
# initialized with False for all nn			
			self.nnFlags = [0] + [False]*4
# scan through possible nn values of this atom
			for nn in self.nn:
				self.nnFlags[nn] = True

	def equiv(self, other):
		"""is the list of possible status of the current atom identical
		to the one of another atom?"""
		st1 = self.stats
		st2 = other.stats
		if len(st1) != len(st2):
# equivalence only if the status lists have the same length...
			return False
		for s1, s2 in zip(st1, st2):
			if s1 != s2:
# ... and the same elements
				return False
# equivalence only if chemical shifts are equal
		if self.chemshift != other.chemshift:
			return False
		return True

	def printout(self):
		"""generate an LSD MULT line for a FSA (could have been a VSA in the past)"""
		id = str(self.id)
		symbol = self.elt
		val = self.valence
# test for ambiguous valence atom sur as N, S or P
		if status.elt_bignval[symbol]:
# symbol includes valency, such as N3 or S6
			symbol += str(self.valence)
# replace explicit symbol (with valency) by implicit symbol when possible
			if symbol in status.elt_usual.keys():
				symbol = status.elt_usual[symbol]
		hyb = str(self.hyb)
		mult = str(self.mult)
# assemble standard MULT command
		line = " ".join(["MULT", id, symbol, hyb, mult])
# append charge if needed
		if self.ec != 0:
			line += " " + str(self.ec)
		return line
	
	def dump(self):
		print()
		print ("Atom:", self.id)
		for k in sorted(self.__dict__.keys()):
			print ("\t", k, ":", self.__dict__[k])
#		for k, v in self.__dict__.iteritems():
#			print "\t", k, ":", v

	def __repr__(self):
		return ("Atom("+str(self.id)+", "+str(self.elt)+", "+str(self.hyb)+", "+str(self.mult)+", "+str(self.ec)+")")

if __name__ == "__main__":
	a1 = Atom(1, 'C', 3, 3)
	print ()
	print (a1)
	print ("stats", a1.stats)
	print ("nn", a1.nn)
	print ("nH", a1.mult)
	print ("minH", a1.minH)
	print ("maxH", a1.maxH)
	print (a1.printout())
	a1.dump()

	
	a2 = Atom(2, 'O', [2, 3], [0, 1])
	print()
	print (a2)
	print ("stats", a2.stats)
	print ("statbynn", a2.statbynn)
	print ("statbynnbynH", a2.statbynnbynH)
	print ("nn_full", a2.nn_full)
	print ("nH_full", a2.nH_full)
	print ("nn", a2.nn)
	print ("nnFlags", a2.nnFlags)
	print ("nH", a2.nH)
	print ("nHbynn", a2.nHbynn)
	print ("minHbynn", a2.minHbynn)
	print ("maxHbynn", a2.maxHbynn)

	a3 = Atom(3, 'C', [1, 2, 3], [0, 1, 2, 3])
	print()
	print (a3)
	print ("stats", a3.stats)
	print ("statbynn", a3.statbynn)
	print ("statbynnbynH", a3.statbynnbynH)
	print ("nn_full", a3.nn_full)
	print ("nH_full", a3.nH_full)
	print ("nn", a3.nn)
	print ("nnFlags", a3.nnFlags)
	print ("nH", a3.nH)
	print ("nHbynn", a3.nHbynn)
	print ("minHbynn", a3.minHbynn)
	print ("maxHbynn", a3.maxHbynn)

#	a4 = Atom(4, 'O', [2, 3], [0, 1], [-1, 0, 1])
	a4 = Atom(4, 'O', 3, 0, [-1, 0, 1])
	print ()
	print (a4)
	print ("stats", a4.stats)
	print ("statbynn", a4.statbynn)
	print ("statbynnbynH", a4.statbynnbynH)
	print ("nn_full", a4.nn_full)
	print ("nH_full", a4.nH_full)
	print ("Ec_full", a4.Ec_full)
	print ("nn", a4.nn)
	print ("nnFlags", a4.nnFlags)
	print ("nH", a4.nH)
	print ("nHbynn", a4.nHbynn)
	print ("minHbynn", a4.minHbynn)
	print ("maxHbynn", a4.maxHbynn)
	print ("Ec", a4.Ec)
#	print "canbePC", a4.canbePC
#	print "canbeNC", a4.canbeNC
	print ("mustbePC", a4.mustbePC)
	print ("mustbeNC", a4.mustbeNC)
	print ("Ecbynn", a4.Ecbynn)
	print ("minEcbynn", a4.minEcbynn)
	print ("maxEcbynn", a4.maxEcbynn)
	print ("mustbePCbynn", a4.mustbePCbynn)
	print ("mustbeNCbynn", a4.mustbeNCbynn)
	print ("EcbynnbynH", a4.EcbynnbynH)
	print ("minEcbynnbynH", a4.minEcbynnbynH)
	print ("maxEcbynnbynH", a4.maxEcbynnbynH)
	
	a4.dump()
	
	a5 = Atom(5, 'N', [1, 2, 3], [0, 1, 2], [0, 1])
	a5.dump()

