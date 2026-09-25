"""status expansion for the atoms for which alternative status still exist"""
def expstats(sts, indexes):
	"""resolution of status ambiguity"""
# result contains the list of all possibilities for non-ambiguous status lists.
# sts is list of alternative status, as lists of possiblities.
# indexes handles equivalence of status lists. Each item in sts has a corresponding
# index in the list of equivalent status lists
	result = []
# number of items in sts and indexes
	imax = len(sts)
# nothing to do
	if imax == 0:
		return result
# list of status values, initially empty
	curlist = []
# current index in sts and indexes
	curi = 0
#	print "First call of recur() with imax:", imax, "sts:", sts, "indexes:", indexes
# recursively fill curlist and copies its value in results as the end
	recur(imax, sts, indexes, curi, curlist, result)
	return result

def recur(imax, sts, indexes, curi, curlist, result):
	"""recursive status list expansion"""
#	print "Recur called with curi:", curi, "curlist:", curlist, "result:", result
# at the end
	if curi == imax:
#		print "solution found:", curlist
# make a copy of the list of unambiguous status and stores it in the result list
		result.append(list(curlist))
# not finished yet
	else:
# current status list to be expanded
		stlist = sts[curi]
# index of the current status list in its equivalence class
		index = indexes[curi]
# loop over status possibilities
		for st in stlist:
# if the current status list is not the first one in its class
			if index != 0:
# and if the current status value is less than the preceding one
				if st < curlist[curi - 1]:
# the current status values has not to be considered. try next one.
					continue
# prepare for next status list
			newi = curi + 1
# remember that st is a valid status value
			curlist.append(st)
# shallow copy of the list of status found so far
			newlist = curlist
# go one step further in sts and indexes, if possible
			recur(imax, sts, indexes, newi, newlist, result)
#			print "back from recur() at curi:", curi
# put curslist in its previous state so that everything goes fine in the loop over status possibilities
			curlist.pop()

if __name__ == "__main__":
	sts = [[2, 3], [2, 3]] 
	indexes = [0, 1]
	print ("sts:", sts)
	print ("indexes:", indexes)
	print ("result:", expstats(sts, indexes))
	
	print
	sts = [[2, 3], [3, 4]] 
	indexes = [0, 0]
	print ("sts:", sts)
	print ("indexes:", indexes)
	print ("result:", expstats(sts, indexes))

	print
	sts = [[2, 3, 4], [2, 3, 4], [2, 3, 4]] 
	indexes = [0, 1, 2]
	print ("sts:", sts)
	print ("indexes:", indexes)
	print ("result:", expstats(sts, indexes))
