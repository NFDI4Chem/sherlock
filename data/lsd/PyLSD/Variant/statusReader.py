"""collect structured information from a status list file"""
class statusReader:
	"""collect structured information from a status list file"""
	def __init__(self, filename):
		"""remembre the name of the file that contains all the possible atom status"""
		self._filename = filename
	
	def read(self):
		"""read and process lines in status file"""
# status dict has element symbols as keys and a list of status data as value
# status data is a list contains, the number of neighbors, the multiplicity,
# the hybridization state, the electric charge, the valence and a textual description
		status = {}
# masses dict has element symbol and integer mass as value
		masses = {}
# elt_usual dict has valency-aware element symbols as key and the corresponding element key as value
# "N3"->"N" says that 3 is the usual valency of N
		elt_usual = {}
# for the guessing of the minimal charge value
		maxminCharge = 0
		minCharge = maxminCharge
# for the guessing of the maximal atom charge value
		maxCharge = -10
# opening of the status file
		with open(self._filename) as f:
# loop over the lines in the status file
			for line in f:
# remove leading and trailing separation characters (including final newline)
				line = line.strip()
				pos = line.find("#")
# try to find a comment in the line
				if pos >= 0:
# remove comment if any
					line = line[0:pos]
				print ("***%s***" % line)
# split line into fields
				it = line.split()
# get number of fields
				nargs = len(it)
# no field?
				if nargs == 0:
# try next line
					continue
				if nargs <= 3:
# the line is a paragraph header. A paragraph is all about a single element
					elt = it[0]
# already found infrmation obout the current element?
					if elt not in list(status.keys()):
# new element, new status information about the current element,
# ignored if this is not the first paragraph for the current element
						status[elt] = []
# get mass of the current element
# ignored if this is not the first paragraph for the current element
						masses[elt] = int(it[1])
						if nargs == 3:
# the current element has more than one possible valence. get usual valence
							val = it[2]
							elt_usual[elt + val] = elt
				else:
# the current line provides a possible status for the current element
# get all (integers) about the status but the textual description
					st = [int(i) for i in it[0:5]]
					charge = st[3]
# update minimal charge
					if charge < minCharge:
						minCharge = charge
# update maximal charge
					if charge > maxCharge:
						maxCharge = charge
# get textual status definition, possibly over more than one field
					definit = "" if nargs < 7 else ' '.join(it[6:])
# append textual definition to status
					st.append(definit)
# append status to the list of possible status of the current element
					status[elt].append(st)
# chargeOffset is 0 if there is no atom with a negative charge,
# else this is the number that must be added to a charge in order
# to get a non negative number, suitable as list index
		chargeOffset = 0 if minCharge == maxminCharge else -minCharge
# minimal size of a list that would be indexed by charge+chargeOffset values
		numChargeIndexes = maxCharge - minCharge + 1
# mass of H, for completeness
		masses["H"] = 1
# return information for status file
		return status, masses, elt_usual, chargeOffset, numChargeIndexes

	def unicity(self, status):
		"""find elements for which the same number of neighbors, the same multiplicity
		the same electric charge may lead to different hybridization/valence values.
		Could be useful to bypass step5 is VSAResolver
		"""
		retval = []
# scan through known element symbols"
		for elt in list(status.keys()):
# dictionary that associates nn/nH/ec triplets to a list of hyb/val pairs
			d = {}
# scan through status of the current element
			for nn, nH, hyb, ec, val, dfn in status[elt]:
# get key value pair in d
				k = str([nn,nH,ec])
				v = str([hyb,val])
				if k in d:
# the elt/nn/nH/ec combination is already known. append the hyb/val pair to the list of hyb/val pairs
					d[k].append(v)
				else:
# create a new hyb/val pair for the current elt/nn/nH/ec combination
					d[k] = [v]
# scan through the elt/nn/nH/ec combinations
			for k in list(d.keys()):
				v = d[k]
# ambiguity of hyb/val if  elt/nn/nH/ec are fixed?
				if len(v) != 1:
# collect an elt, nn/nH/ec, list of hyb/val triplet
					retval.append([elt, k, v])
#			print elt, d
		return retval
				
if __name__ == "__main__":
	import defaults
	defaults = defaults.Defaults()
	r = statusReader(defaults.statuslist)
	status, masses, elt_usual, chargeOffset, numChargeIndexes = r.read()
	print ("status:\n", status, "\n")
	print ("masses:", masses)
	print ("elt_usual:", elt_usual)
	print ("chargeOffset:", chargeOffset)
	print ("numChargeIndexes:", numChargeIndexes)

	print ("non-unicity situations after disambiguation of elt/nn/nH/ec:", r.unicity(status))
