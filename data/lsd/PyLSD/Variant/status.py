"""Stores atom status information.

Each status is composed of 5 items:
1. Number of neighbors
2. Number of H atoms
3. Hybridization state n, as in sp^n
4. Electric charge
5. Valence
"""

import re
import statusReader
import defaults
from functools import reduce

# chargeOffset: must be chosen so that charge+chargeOffset is >= 0
# status: status collection, by element symbol
# masses: atomic mass collection, the ones for molecular mass calculation
# elt_usual: usual element symbols from symbol followed by valence
defaults = defaults.Defaults()
r = statusReader.statusReader(defaults.statuslist)
status, masses, elt_usual, chargeOffset, numChargeIndexes = r.read()

# construction of all_elements, a list of existing atomic symbols
# starting with "C", "N", "O" and with the other ones (if any)
# sorted by alphabetical order.
setof_all_elements = set(status.keys())

minimal_list = ["C", "N", "O"]
setof_minimal = set(minimal_list)
setof_extras = setof_all_elements - setof_minimal
list_extras = sorted(list(setof_extras))

# list of all known elements
all_elements = minimal_list + list_extras
print ("All elements:", all_elements)

# canonical element ordering
canonical_elt_order = list(all_elements)
canonical_elt_order.insert(1, "H")
print ("Canonical order of elements:", canonical_elt_order)

# number of known elements
numelt = len(all_elements)

# elt_hash associates each element symbol
# to its index in all_elements
elt_hash = {}
for i in range(numelt):
	elt_hash[all_elements[i]] = i

# base concentrates all the status information,
# each status being given an index value
# that is its index in base.
# Information is taken from status[] and the element name
# is prepended to the status description list.
base = []
for elt in all_elements:
	for uplet in status[elt]:
		uplet.insert(0, elt)
		base.append(uplet)

# number of known status
numstat = len(base)

# each particular piece of status in each status (valence excepted)
# is associated to a corresponding set of status index.
# status_elt_tmp is indexed by element index insted of element symbol.
status_elt_tmp = [ set() for i in range(numelt) ]
status_nei = [ set() for i in range(5) ]
status_mult = [ set() for i in range(4) ]
status_hyb = [ set() for i in range(4) ]
status_charge = [ set() for i in range(numChargeIndexes) ]
status_val = [ set() for i in range(7) ]

elt_val = dict([[elt, set()] for elt in all_elements])
elt_hyb = dict([[elt, set()] for elt in all_elements])
elt_mult = dict([[elt, set()] for elt in all_elements])
elt_charge = dict([[elt, set()] for elt in all_elements])

f = open("base.txt", "w")
f.write("1. status id\n")
f.write("2. element\n")
f.write("3. number of neighbors\n")
f.write("4. number of attached H atoms\n")
f.write("5. hybridization\n")
f.write("6. electric charge\n")
f.write("7. valence\n")
f.write("8. definition\n")
f.write("\n")
for i in range(numstat):
	elt, nei, mult, hyb, charge, val, definit = base[i]
	f.write("%3d: %-2s %d %d %d %2d %d %s\n" % (i, elt, nei, mult, hyb, charge, val, definit))
	status_elt_tmp[elt_hash[elt]].add(i)
	status_nei[nei].add(i)
	status_mult[mult].add(i)
	status_hyb[hyb].add(i)
	status_charge[charge+chargeOffset].add(i)
	status_val[val].add(i)
	elt_val[elt].add(val)
	elt_hyb[elt].add(hyb)
	elt_mult[elt].add(mult)
	elt_charge[elt].add(charge)
f.close()

# default status for all elements
defstat = {}
for elt in all_elements:
	vals = sorted(list(elt_val[elt]))
	symbol = elt
	if len(vals) > 1:
		for v in vals:
			symbol += str(v)
	hybs = sorted(list(elt_hyb[elt]))
	mults = sorted(list(elt_mult[elt]))
	ecs = sorted(list(elt_charge[elt]))
	defstat[elt] = [symbol, hybs, mults, ecs]
	print ("Default status for", elt, "is", defstat[elt])

# indexing by element symbol instead of element index
# for status_elt
status_elt = {}
for i in range(numelt):
	status_elt[all_elements[i]] = status_elt_tmp[i]

# association between element and the number of possible valences being > 1
elt_bignval = {}
for elt, val in elt_val.items():
	elt_bignval[elt] = len(val) > 1

# valence parity for the elements
elt_valence_parity = {}
for elt, val in elt_val.items():
	elt_valence_parity[elt] = (list(set([v%2 for v in list(val)])))[0]
elt_valence_parity["H"] = 1
print ("elt_valence_parity:", elt_valence_parity)

# "inverted" dict of usual element symbols
inv_elt_usual = {}
for unusual, usual in elt_usual.items():
	inv_elt_usual[usual] = unusual

# need or_ operator (set union) for the dealing with variable status
from operator import or_
		
def status_id(l, value):
	"""gets status ids for a status piece, either fixed or variable"""
	if isinstance(value, int) or isinstance(value, str):
		return l[value]
	else:
		stats = list(map(lambda x: l[x], value))
		return reduce(or_, stats, set())

def status(elt, mult, hyb, charge):
	"""gets status ids from the 4 status pieces, each being either fixed or variable"""
	if elt in list(inv_elt_usual.keys()):
# N, P, S, Br, I alone hold for N3, P3, S2, Br1, I1
		elt = inv_elt_usual[elt]
	if elt not in all_elements:
# for N3, N5, N35, P3, P5, P35, S2, S4, S6, S24, S26, S46, S246, Br1, Br3, I1, I3, I5
		m = re.match(r'([A-Z][a-z]?)([1-6]?)([1-6]?)([1-6]?)', elt)
		fields = list(m.groups())
		elt = fields[0]
		vals = [int(val) for val in fields[1:] if val]
		sf = set()
		for val in vals:
			sf = sf | status_id(status_val, val)
	else:
# for any other element symbol
		sf = set(range(numstat))
	s = status_id(status_elt, elt)
	sf = sf & s
#	print s
#	print sf
	s = status_id(status_mult, mult)
	sf = sf & s
#	print s
#	print sf
	s = status_id(status_hyb, hyb)
	sf = sf & s
#	print s
#	print sf
	charge = charge+chargeOffset if isinstance(charge, int) else [c+chargeOffset for c in charge]
	s = status_id(status_charge, charge)
	sf = sf & s
#	print s
#	print sf
	return sorted(list(sf))

def status_nn_one(stat):
	"""returns the number of neighbors(nn) that corresponds to stat, a status index"""
	return base[stat][1]

def status_nn(stats):
	"""returns the list of possible nn according to a list of status indexes"""
	return list(map(status_nn_one, stats))

def status_nH_one(stat):
	"""returns the number of hydrogen atoms(nH) that corresponds to stat, a status index"""
	return base[stat][2]

def status_nH(stats):
	"""returns the list of possible nH according to a list of status indexes"""
	return list(map(status_nH_one, stats))

def status_hyb_one(stat):
	"""returns the hybridization state (hyb) that corresponds to stat, a status index"""
	return base[stat][3]

def status_hybr(stats):
	"""returns the list of possible hybridization states according to a list of status indexes"""
	return list(map(status_hyb_one, stats))

def status_Ec_one(stat):
	"""returns the electric charge (Ec) that corresponds to stat, a status index"""
	return base[stat][4]

def status_Ec(stats):
	"""returns the list of possible electric charges according to a list of status indexes"""
	return list(map(status_Ec_one, stats))

def status_all(stat):
	"""returns element, nn, nH, hybridization state, charge, valence and textual definition from status id"""
	return base[stat]

if __name__ == "__main__":
	print ("base = ", base)
	print ("status_elt = ", status_elt)
	print ("status_nei = ", status_nei)
	print ("status_mult = ", status_mult)
	print ("status_hyb = ", status_hyb)
	print ("status_charge = ", status_charge)
	
	print ("status_id(status_mult, 3) = ", status_id(status_mult, 3))
	print ("status_id(status_mult, [1, 2]) = ", status_id(status_mult, [1, 2]))
	print ("status('C', 3, 3, 0) = ", status('C', 3, 3, 0))
	print ("status('O', [0, 1], [2, 3], 0) = ", status('O', [0, 1], [2, 3], 0))
	print ("status('O', [0, 1], [2, 3], [-1, 0, 1]) = ", status('O', [0, 1], [2, 3], [-1, 0, 1]))
	print ("status('O', 1, [2, 3], 0) = ", status('O', 1, [2, 3], 0))
	print ("status('N', [0, 1, 2], [2, 3], 0) = ", status('N', [0, 1, 2], [2, 3], 0))
	print ("status('N', [0, 1, 2], [1, 2, 3], 0) = ", status('N', [0, 1, 2], [1, 2, 3], 0))
	print ("status('N5', [0, 1, 2], [1, 2, 3], 0) = ", status('N5', [0, 1, 2], [1, 2, 3], 0))
	print ("status('N35', [0, 1, 2], [1, 2, 3], 0) = ", status('N35', [0, 1, 2], [1, 2, 3], 0))
	print ("status_nn(status('O', [0, 1], [2, 3], 0)) = ", status_nn(status('O', [0, 1], [2, 3], 0)))
	print ("status_nH(status('O', [0, 1], [2, 3], 0)) = ", status_nH(status('O', [0, 1], [2, 3], 0)))
	print ("status_Ec(status('O', [0, 1], [2, 3], [-1, 0, 1])) = ", status_Ec(status('O', [0, 1], [2, 3], [-1, 0, 1])))
	
