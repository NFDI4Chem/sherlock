"""Expansion of electric charge alternatives"""
def sumnEc(Eclists, indexes, sumEc, maxPC, maxNC):
	"""return lists of flat lists of electric charge values from Eclists,
	a list of electric charge lists that contain alternatives.
	Indexes refer to index of items in Eclists in their equivalence classes.
	The sum of electric charges must be an element in sumEc.
	There are at most maxPC and maxNC positively and negatively charged atoms.
	Any of these values may be -1 if the constraint is disabled"""
# initial result
	result = []
# number of items in Eclists and indexes
	imax = len(Eclists)
# nothing to do
	if imax == 0:
		return result
# maximal/minimal charge for each atom
	maxEcs = [max(l) for l in Eclists]
	minEcs = [min(l) for l in Eclists]
# maximal/minimal charge for the atoms with undefinite electric charge
	cmax = sum(maxEcs)
	cmin = sum(minEcs)
# at least one value in sumEc must be in the correct range, but who knows?
	ok = bool([e for e in sumEc if cmin <= e <= cmax])
	if not ok:
		return result
# calculate cumulative minimum/maximum charge values
	cummaxs = [0]*imax
	cummins = [0]*imax
	for i in range(imax):
		cmax -= maxEcs[i]
		cmin -= minEcs[i]
		cummaxs[i] = cmax
		cummins[i] = cmin
# prepare recursive search for unambiguous electric charge values
# current index in Eclists and indexes
	curi = 0
# sum of electric charges found so far
	cursum = 0
# number of positively/negatively charged atoms so far
	curPC = 0
	curNC = 0
# list of unambiguous electric charge values so far
	curEcs = []
	print ("Eclists:", Eclists)
	print ("indexes:", indexes)
	print ("sumEc:", sumEc, " maxPC:", maxPC, " maxNC:", maxNC)
	print ("imax:", imax, " cummaxs:", cummaxs, " cummins:", cummins)
# run recursive disambiguation
	recur(Eclists, imax, indexes, sumEc, maxPC, maxNC, cummaxs, cummins, curi, cursum, curPC, curNC, curEcs, result)
	return result

def recur(Eclists, imax, indexes, sumEc, maxPC, maxNC, cummaxs, cummins, curi, cursum, curPC, curNC, curEcs, result):
	"""recursive construction of the list of the lists of unambiguous charge values"""
	if curi == imax:
# all electric charge values found, append their list to the result
		result.append(list(curEcs))
		print ("curi:", curi, " curEcs:", curEcs)
	else:
# current list of electric charge alternatives
		Eclist = Eclists[curi]
# index of the current list of electric charge alternatives in its equivalence class
		index = indexes[curi]
# maximum sum of all remaining electric charges
		cummax = cummaxs[curi]
# minimum sum of all remaining electric charges
		cummin = cummins[curi]
# try an electric charge value Ec from the list of alternatives
		for Ec in Eclist:
			if index != 0:
# this is not the head of its equivalence class, Ec cannot be less that the value that was precedently found
				if Ec < curEcs[curi - 1]:
# try next one
					continue
# accept Ec as valid but test for compatibility with the future choices
			newsum = cursum + Ec
# maximum/minimum electric charge if all the remaining ones are at their maximum/minimum value
			highlimit = newsum + cummax
			lowlimit = newsum + cummin
# test if one of the possible total electric charge values fit with the currently reachable interval
			sumOK = bool([e for e in sumEc if lowlimit <= e <= highlimit])
			print ("curi:", curi, " Ec:", Ec)
# test for the number of electrically charged atoms
			if sumOK:
# positively charged atom count a priori OK
				PC_OK = True
# update the number of positively charged atoms
				newPC = curPC + (1 if Ec > 0 else 0)
# the number of positively charged atoms is limited
				if maxPC != -1:
# check if the limitation is respected
					PC_OK = newPC <= maxPC
				NC_OK = True
				newNC = curNC + (1 if Ec < 0 else 0)
				if maxNC != -1:
					NC_OK = newNC <= maxNC
# check is the constrainsts on the number of positively and negatively charged atoms are respected
				if PC_OK and NC_OK:
# prepare for next recursive function call
					newi = curi + 1
# the currently validated electric charge value is appended to the list of the preceding ones
					curEcs.append(Ec)
# symbolic new list of validated electric charge values
					newEcs = curEcs
# recursive call now
					recur(Eclists, imax, indexes, sumEc, maxPC, maxNC, cummaxs, cummins, newi, newsum, newPC, newNC, newEcs, result)
# put the list of validated electric charge values to its former value.
					curEcs.pop()

if __name__ == "__main__":
	def test(Eclists, indexes, sumEc, maxPC, maxNC):
		print()
		print ("Eclists:", Eclists)
		print ("indexes:", indexes)
		print ("sumEc:", sumEc, " maxPC:", maxPC, " maxNC:", maxNC)
		print ("result:", sumnEc(Eclists, indexes, sumEc, maxPC, maxNC))
	
	Eclists = [[0, 1], [0, 1], [-1, 0], [-1, 0], [-1, 0], [-1, 0]]
	indexes = [0, 1, 0, 1, 2, 3]
	sumEc = [-2, -1, 0, 1, 2]
	maxPC = 1
	maxNC = 1
	test(Eclists, indexes, sumEc, maxPC, maxNC)

