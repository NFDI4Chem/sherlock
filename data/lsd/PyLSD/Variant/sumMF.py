"""Generate variable parts of molecular formula (MF) under molecular mass constraint"""
def sumMF(paritymatters, num_parity_odd, parity, numH, mass, numA, massA):
	"""returns lists of element numbers. Paritity is defined as the remainder of the division by 2
	of the sum of the valence of all atoms. A single electric charge deficit or excess is equivalent
	to a parity variation of one unit. Each atomic element is univocally associated to the parity
	of its valence, even though more than one valence is possible like for N(3,5), P(3,5), S(2,4,6).
	The part of the molecular formula under consideration has num_parity_odd elements with parity equal
	to 1 (odd). The parity (parity) of the known part of the MF is equal to the one of the unknown part.
	numH is a list of possible values for the number of elements H in the MF. mass is a list of two
	values, the min and max of the molecular mass (MM) for the unknown part of the molecule.
	numA is a list of lists of possible number of elements of unknown number (H excluded).
	Each list in numA is sorted in ascending order.
	massA if the list of the masses of the elements of unknown number (H excluded).
	The order in numA and massA corresponds to the same element lists. The first num_parity_odd
	elements in this list (H excluded) have an odd parity."""

# final list of lists of number of elements
	result = []
# number of non-H elements in unknown number
	imax = len(numA)
# should not happen
	if imax == 0:
		return result
# numH is an empty list if the complete MF as no H inside.
# If MF has H, eventual alternatives will be considered first,
# before all non-H elements.
	if numH:
# insert number of H and mass of H in numA and massA at the beginning.
		imax += 1
		numA.insert(0, numH)
		massA.insert(0, 1)
# H is an element with an odd valence number. One more of them
		numodd = num_parity_odd + 1
	else:
# nothing is changed if there is no H
		numodd = num_parity_odd
# calculate the minimum and maximum mass than can be associated to each element
	minA = [n[0] * m for n, m in zip(numA, massA)]
	maxA = [n[-1] * m for n, m in zip(numA, massA)]
# calculate the overall mass that is associated to H and elements in unknown number.
	summin = sum(minA)
	summax = sum(maxA)
# impossible to get numbers of elements if the [summin-summax]
# and [mass[0]-mass[1]] intervals de not have any intersection
	if summax < mass[0] or summin > mass[1]:
		return result
# calculate cumulated min and max of masses, by element
	mincums = [summin - minA[0]]
	maxcums = [summax - maxA[0]]
	for i in range(1, imax):
		mincums.append(mincums[i-1] - minA[i])
		maxcums.append(maxcums[i-1] - maxA[i])
# prepare recursive search of the numbers of atoms for all elements
	curlist = []
	curmass = 0
	curi = 0
	curparity = parity % 2
	print ("imax:", imax)
	print ("paritymatters:", paritymatters)
	print ("numodd:", numodd)
	print ("numA:", numA)
	print ("massA:", massA)
	print ("mincums:", mincums)
	print ("maxcums:", maxcums)
# go!
	recur(imax, paritymatters, numodd, numA, massA, mass, mincums, maxcums, curi, curparity, curmass, curlist, result)
	return result
	
def recur(imax, paritymatters, numodd, numA, massA, mass, mincums, maxcums, curi, curparity, curmass, curlist, result):
	"""search of all list of possible element numbers, by element.
	imax is the number of values i numA, massA, mincums, maxcum; each value corresponds to a chemical element.
	the first numodd chemichal elements have an odd valence
	numA is the list the lists of possible number of chemical elements
	massA is the atomic mass of each chemical element
	mass[0] and mass[1] define the inclusive limits of the molecular mass
	mincums, maxcums are the cumulated min and max mass of the remaining elements (those after the current one)
	curi defines the current element under study; positional parameter in numA, massA, mincums, maxcums
	curparity id the valence parity so far
	curmass is the mass so far
	curlist is the list of the elements for which a number has been postulated so far
	result accumulates the curlists when they are full and valid.
	"""
#	print "curi:", curi, " curparity:", curparity, "curmass:", curmass, "curlist:", curlist
	print ("curmass: %3d" % curmass, "  ", curlist)
#	print "curi:", curi, "not paritymatters:", not paritymatters, "numodd:", numodd, "parity:", curparity
# backtrack occurs if paritymatters and if curi is the one at which parity must be checked (numodd)
# and if parity is 1, because it must be the 0 (sum of all atomic valencies is even, charges included)
	if not paritymatters or numodd != curi or curparity == 0:
#		print "curi:", curi, "parity OK"
# all elements have got a defined number of occurences
		if imax == curi:
			print ("curmass: %3d" % curmass, "  ", curlist, "*"*4)
#			print "curi: ", curi, "*"*4, curlist
# create a deep copy of the list so that it will not be modified later and append the copy to the result list
			result.append(list(curlist))
		else:
# try to find n, the number of elements for the element indexed by curi
			nums = numA[curi]
			m = massA[curi]
			mmin = mincums[curi]
			mmax = maxcums[curi]
			for n in nums:
# add n atoms of the current chemical element to the MF, get new mass
				newmass = curmass + m * n
# get the min and max molecular mass considering that the elements that will be considered later
# are in their minimal or maximal number
				newmmin = mmin + newmass
				newmmax = mmax + newmass
#				print "curi:", curi, "newmass:", newmass, "newmmin:", newmmin, "newmmax:", newmmax
				if newmmax < mass[0]:
# maximal mass that can be reached is less than the imposed minimal final mass. n must be bigger. go on looping
					continue
				if newmmin > mass[1]:
# minimal mass than can be reached is more than the imposed maximal final mass. n is already too big. stop looping 
					break
# at this point the value of n is valid. valence parity will be checked later.
# get new molecular valence parity, if this has a meaning
				newparity = (curparity + n) % 2 if paritymatters and curi < numodd else curparity
#				print "curparity:", curparity, "n:", n, "curi:", curi, "numodd:", numodd, "newparity:", newparity
# ready for next chemical element
				newi = curi + 1
# append n to the current list of n values
				curlist.append(n)
# newlist is the same as curlist
				newlist = curlist
# try to find n for the next element, if any
				recur(imax, paritymatters, numodd, numA, massA, mass, mincums, maxcums, newi, newparity, newmass, newlist, result)
# remove n from curlist, necessary for considering the next n value in nums
				curlist.pop()
				print ("curmass: %3d" % curmass, "  ", curlist)
	else:
#		print "curi:", curi, "parity ignored or not OK"
		pass

if __name__ == "__main__":
# C 1 H 1-5 N 1-2 O 2-3, no charge, molecular mass is 61
# paritymatters because the electric charge is defined (0, no charge)
	paritymatters = True
# H excepted, N is the only (1) chemical element in unknown number with an odd valence.
# The index of N is 0. The index of O is then 1.
	num_parity_odd = 1
# The known part of the formula (C 1, no charge) has an even global valence
	parity = 0
# the H part of the molecular formula
	numH = [1, 2, 3, 4, 5]
# the molecular mass is 61. The mass of the known part is 12. The mass of the unkown part is 61-12=49
	mass = [49, 49]
# the N and O part of the molecular formula
	numA = [[1, 2], [2, 3]]
# the atomic mass of N and O
	massA = [14, 16]
# get all possibilities
	result = sumMF(paritymatters, num_parity_odd, parity, numH, mass, numA, massA)
# print result
	print()
	print (result)

# C 1 H 3 N 1 O 2-3, no charge, no mass constraint
	paritymatters = True
# no element in unknown number with an odd valence
	num_parity_odd = 0
# the known part of the molecule is CN without charge. The sum of valences is odd.
	parity = 1
# 3 H atoms
	numH = [3]
# no mass constraint
	mass = [-100, 100]
# the number of O is the only one to be unknown (2 or 3)
	numA = [[2, 3]]
	massA = [16]
	result = sumMF(paritymatters, num_parity_odd, parity, numH, mass, numA, massA)
	print()
	print (result)
