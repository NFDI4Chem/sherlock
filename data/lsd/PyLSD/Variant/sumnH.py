def sumnH(nHlists, indexes, sumH):
	result = []
	imax = len(nHlists)
	if imax == 0:
		return result
	minH = list(map(lambda x: min(x), nHlists))
	maxH = list(map(lambda x: max(x), nHlists))
	summinH = sum(minH)
	summaxH = sum(maxH)
	ok = summinH <= sumH <= summaxH
	if not ok:
		return result
	minHcums = [summinH - minH[0]]
	maxHcums = [summaxH - maxH[0]]
	for i in range(1, imax):
		minHcums.append(minHcums[i-1] - minH[i])
		maxHcums.append(maxHcums[i-1] - maxH[i])
	curlist = []
	cursumH = 0
	curi = 0
#	print "imax", imax
#	print "nHlists", nHlists
#	print "indexes", indexes
#	print "minHcums", minHcums
#	print "maxHcums", maxHcums
#	print "sumH", sumH
	recur(imax, nHlists, indexes, minHcums, maxHcums, sumH, curi, cursumH, curlist, result)
	return result
	
def recur(imax, nHlists, indexes, minHcums, maxHcums, sumH, curi, cursumH, curlist, result):
#	print
#	print "curi", curi
#	print "cursumH", cursumH
#	print "curlist", curlist
#	print "result", result
	if curi == imax:
		result.append(list(curlist))
	else:
		nHlist = nHlists[curi]
		index = indexes[curi]
		minHcum = minHcums[curi]
		maxHcum = maxHcums[curi]
#		print "nHlist:", nHlist, "index:", index, "minHcum:", minHcum, "maxHcum:", maxHcum
		for nH in nHlist:
#			print "Trying nH:", nH, "with curi:", curi
			if index != 0:
				if nH < curlist[curi - 1]:
#					print "try next nH, elimination by symmetry"
					continue
			newsumH = cursumH + nH
			remainH = sumH - newsumH
#			print "newsumH:", newsumH, "remainH:", remainH
			if minHcum > remainH:
#				print "nH is too high. loop finished"
				break
			if remainH > maxHcum:
#				print "nH is too low. try next one."
				continue
#			print "nH is OK"
			curlist.append(nH)
			newlist = curlist
			newi = curi + 1
			recur(imax, nHlists, indexes, minHcums, maxHcums, sumH, newi, newsumH, newlist, result)
#			print "back from recur with curi:", curi
			curlist.pop()



if __name__ == "__main__":
	nHlists = [[0, 1], [0, 1], [0], [0]]
	indexes = [0, 1, 0, 1]
	sumH = 1
	print()
	print ("nHlists:", nHlists)
	print ("indexes:", indexes)
	print ("sumH:", sumH)
	print ("Result:", sumnH(nHlists, indexes, sumH))

	nHlists = [[0, 1], [0, 1], [0, 1], [0, 1]]
	indexes = [0, 1, 2, 3]
	sumH = 1
	print()
	print ("nHlists:", nHlists)
	print ("indexes:", indexes)
	print ("sumH:", sumH)
	print ("Result:", sumnH(nHlists, indexes, sumH))



	nHlists = [[0, 1, 2], [0, 1, 2], [0, 1, 2], [0, 1, 2]]
	indexes = [0, 1, 2, 3]
	sumH = 0
	print()
	print ("nHlists:", nHlists)
	print ("indexes:", indexes)
	print ("sumH:", sumH)
	print ("Result:", sumnH(nHlists, indexes, sumH))

	nHlists = [[0, 1, 2], [0, 1, 2], [0, 1], [0, 1]]
	indexes = [0, 1, 0, 1]
	sumH = 0
	print()
	print ("nHlists:", nHlists)
	print ("indexes:", indexes)
	print ("sumH:", sumH)
	print ("Result:", sumnH(nHlists, indexes, sumH))

	nHlists = [[0, 1], [0, 1], [0, 1], [0, 1]]
	indexes = [0, 1, 2, 3]
	sumH = 0
	print()
	print ("nHlists:", nHlists)
	print ("indexes:", indexes)
	print ("sumH:", sumH)
	print ("Result:", sumnH(nHlists, indexes, sumH))



	nHlists = [[0, 1], [0, 1], [0, 1], [0, 1], [0, 1]]
	indexes = [0, 1, 2, 0, 1]
	sumH = 5
	print()
	print ("nHlists:", nHlists)
	print ("indexes:", indexes)
	print ("sumH:", sumH)
	print ("Result:", sumnH(nHlists, indexes, sumH))

	nHlists = [[0, 1], [0, 1], [0], [1, 2], [0, 1]]
	indexes = [0, 1, 0, 0, 0]
	sumH = 5
	print()
	print ("nHlists:", nHlists)
	print ("indexes:", indexes)
	print ("sumH:", sumH)
	print ("Result:", sumnH(nHlists, indexes, sumH))

	nHlists = [[0, 1], [0], [0], [1, 2], [1, 2]]
	indexes = [0, 0, 1, 0, 1]
	sumH = 5
	print()
	print ("nHlists:", nHlists)
	print ("indexes:", indexes)
	print ("sumH:", sumH)
	print ("Result:", sumnH(nHlists, indexes, sumH))

	nHlists = [[0, 1], [0, 1], [0, 1], [1, 2], [0]]
	indexes = [0, 1, 2, 0, 0]
	sumH = 5
	print()
	print ("nHlists:", nHlists)
	print ("indexes:", indexes)
	print ("sumH:", sumH)
	print ("Result:", sumnH(nHlists, indexes, sumH))

	nHlists = [[0, 1], [0, 1], [0, 1], [1, 2], [1, 2]]
	indexes = [0, 1, 2, 0, 1]
	sumH = 5
	print()
	print ("nHlists:", nHlists)
	print ("indexes:", indexes)
	print ("sumH:", sumH)
	print ("Result:", sumnH(nHlists, indexes, sumH))

