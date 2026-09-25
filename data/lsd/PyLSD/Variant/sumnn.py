def sumnn(v, finalsum, natoms, nnparity):
	"""Generates lists with numbers of VSA according to their number of neighbors
	
	v: list of maximum number of VSA, for each possible number of neighbors
	finalsum: number of VSA
	natoms: total number of atoms
	nnparity: constrains the sum of number of neighbors of VSA to be even(0) or odd(1).
	In each array of the returned value the sum of the values equals to finalsum
	and each value is less than the corresponding value (same index) in v.
	"""
# takes care of a rare situation, in which the number of atoms
# is less or equal to 4. For 4 atoms, an atom
# cannot have 4 neighbors.
	i = 4
	while i > natoms - 1:
		v[i] = 0
		i -= 1
# remove trailing zeros in v		
	i = 4
	while v[i] == 0:
		i -= 1
	imax = i
	print ("imax:", imax)
# w keeps track of the removal of the places at with v contains 0
	i = 1
	w = list(range(5))
# scan through v
	while i <= imax:
		print ("  i:", i, "v["+str(i)+"]:", v[i])
# if current v element is zero
		if v[i] == 0:
# moves v and w values to the left in order to fill the gap
			for j in range((i+1), (imax+1)):
				v[(j-1)] = v[j]
				w[(j-1)] = w[j]
# one less non-zero value in v
# i is not incremented because the shift of v
# may have introduced another zero in v at index i
			imax -= 1
# no zero
		else:
# one position ahead in v
			i += 1
	print ("imax:", imax)
	print (v[1:imax+1], w[1:imax+1])
# v has been reordered so that indexes of interest are between 1 and imax included
# v does not contain zeros. The values that concern the odd parity nn (1 and/or 3) will now have
# the lowest indexes. Indexes between 1 and nodd (none if nodd is 0) correspond to odd 
# nn values. w is the key that restores the initial order in arrays of number of VSA,
# nn by nn. The initial 0 in arrays is always preserved but has no meaning.
# For example: with v=[0, 0, 4, 6, 8] one gets v=[0, 6, 4, 8] and w=[0, 3, 2, 4] and nodd=1.
# The odd number of neighbors must be considered first because parity violatation
# of the number of molecular neighbors (sum of all number of neighbors of all atoms)
# must be even. For the same reason, this part of the algorithm must be aware
# of the parity of the molecular number of neighbors of the fixed status atoms (nnparity).

# list of odd and even indexes
	wodd = []
	weven = []
# list of odd and even values
	vodd = []
	veven = []
# fill these list according to parity of w elements
	for i in range(1,imax+1):
		print ("i:", i, "w["+str(i)+"]:", w[i])
		if w[i] % 2:
			wodd.append(w[i])
			print ("append", w[i], "to wodd")
			vodd.append(v[i])
			print ("append", v[i], "to vodd")
		else:
			weven.append(w[i])
			print ("append", w[i], "to weven")
			veven.append(v[i])
			print ("append", v[i], "to veven")
# reordering, with a heading 0
	v = [0] + vodd + veven
	w = [0] + wodd + weven
# number of odd number of neighbors
	nodd = len(wodd)
	print ("new v:", v[1:imax+1], "w:", w[1:imax+1], "nodd:", nodd)
	print (v)
# v has now been completely reordered
# A cumul vector is now computed. Its first value
# is the sum of all v values. 
	cumul = [sum(v)]
	for i in range(imax):
# stores the maximum number of atoms to be placed
		cumul.append(cumul[i] - v[i+1])
	print ("cumul:", cumul)
# initialiasing of recursion by the gen() function
# result: the list of all possible number of atoms according to number of neighbors
	biglist = []
	curlist = [0]
	cursum = 0
	curi = 1
	curparity = nnparity
# if molecular nn is odd without any VSA having and odd nn,
# then the list of solutions is empty.
# The number  of atoms must be less or equal to the sum all max nn
	if (False if (nodd == 0 and nnparity != 0) else True) and finalsum <= cumul[0]:
		print ("******** Calling gen() for the first time ********")
		print ("v:", v)
		print ("imax:", imax)
		print ("finalsum:", finalsum)
		print ("cumul:", cumul)
		print ("nodd:", nodd)
		print ("*****************")
# call to recursive gen(). v, imax, finalsum, cumul, nodd are not modified by gen()
		gen(v, imax, finalsum, cumul, nodd, curi, curlist, cursum, curparity, biglist)
		print (biglist)
# reorder vector of number of atoms by nn according to w
	reorderedlist = list(map(lambda x: reorder(x, w), biglist))
	print (reorderedlist)
	return reorderedlist

def gen(v, imax, finalsum, cumul, nodd, curi, curlist, cursum, curparity, biglist):
	"""recursive generation of arrays of number of VSAs by number of neighbors"""
	blanks = " "*(2*curi)
	print (blanks, "***gen() called  with ***")
	print (blanks, "curi:", curi)
	print (blanks, "curlist:", curlist)
	print (blanks, "cursum:", cursum)
	print (blanks, "curparity:", curparity)
	print (blanks, "biglist:", biglist)
	print (blanks, "******")
# missing is the number atoms to be distributed now and later
	missing = finalsum - cursum
	print (blanks, missing, "missing atom"+("s" if missing > 1 else ""))
# are we about to calculate the last value is curlist?
	if curi == imax:
# of course, there are missing missing atoms 
		na = missing
		parityOK = True
# parity of molecular nn is a priori correct unless demonstrated not to be.
# parity issue is not solved yet when only dealing with odd nn atoms
		if curi == nodd:
			newparity = (na + curparity) % 2
			print (blanks, "newparity:", newparity)
			parityOK = newparity == 0
# calculate new parity
		print (blanks, "parityOK:", parityOK)
# parity must be 0 (even)
		if parityOK:
# a new array of number of VSAs by nn is found
			newlist = list(curlist)
# ... and completed with the final value
			newlist.append(na)
			print (blanks, "New list was found:", newlist)
# ... and append to the list of the solutions
			biglist.append(newlist)
# we are not calculating the last value in 
	else:
# cumul[curi] would be the number of missing atoms
# if all missing atoms would have to be distributed later
# the number of atoms that is distributed now: na
# na is therefore greater or equal to min1
		min1 = missing - cumul[curi]
# na is less or equal to missing
		max1 = missing
# na is greater or equal to 0
		min2 = 0
# na is less or equal to v[curi] by definition of v
		max2 = v[curi]
# na is therefore greater or equal to minna
		minna = max(min1, min2)
# na is therefore less or equal to maxna
		maxna = min(max1, max2)
		print (blanks, "na between", minna, "and", maxna)
# Testing all possible values for na, the number of atoms with nn == w[curi]
		for na in range(minna, maxna+1):
			print (blanks, "testing na:", na)
# parity of molecular nn is a priori correct unless demonstrated not to be.
			parityOK = True
			newparity = curparity
# parity is updated only for atoms with odd nn
			if curi <=  nodd:
				newparity = (na + curparity) % 2
				print (blanks, "newparity:", newparity)
# parity can be checked only when all odd nn values have been considered
			if curi == nodd:
# parity must be even
				parityOK = newparity == 0
			print (blanks, "parityOK:", parityOK)
# about to consider the next nn value
			if parityOK:
# next nn value is w[newi]
				newi = curi + 1
# create a new list of number of VSAs by nn
				newlist = list(curlist)
				newlist.append(na)
# new number of placed atoms
				newsum = cursum + na
				print (blanks, "will call gen() with newi:", newi, "newlist:", newlist, "newsum:", newsum, "newparity:", newparity)
# recursive call to gen()
				gen(v, imax, finalsum, cumul, nodd, newi, newlist, newsum, newparity, biglist)

def reorder(x, w):
	"""creates a new array of number of VSAs by nn by rordering x according to w (new2old)"""
	print ("reordering: ", x, "according to:", w)
	l = [0]*5
	for i in range(len(x)):
		l[w[i]] = x[i]
	return l

if __name__ == "__main__":
	sumnn([0, 4, 4, 0, 0], 4, 4, 0)
# up to 4 atoms with nn=1 and 4 atoms with nn=2 for 4 variable status atoms and 4 atoms
# The sum of all numbers of neighbors must be even (0).
