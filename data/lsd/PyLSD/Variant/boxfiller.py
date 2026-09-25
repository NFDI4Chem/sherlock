"""BoxFiller class definition"""

import box

class BoxFiller:
	"""Filling an array of boxes with the exact number of beads.
	
	Each box has its own size and permissions.
	The number of beads of each kind to be distributed in the boxes is given
	as argument to fill().
	The information about beads is given by an array that starts with 0.
	The next value is the number of beads with id egal to 1,
	the next value is for id 2, and so on...
	The beads are distributed in the boxes with ids in ascending order.
	"""
	def __init__(self, boxlist):
		"""A BoxFiller is created from a list of boxes"""
		self.boxlist = boxlist
		self.numboxes = len(boxlist)
	def getcontents(self):
		"""returns the content of all boxes. New list references are created"""
		return list(map(lambda b: list(b.get()), self.boxlist))
	def fill(self, beads):
		"""Fills the boxes with beads"""
# the list of box contents is initially empty
		self.biglist = []
# organizes bead information as an argument to the recur() function
		self.org(beads)
# call the recursive bead distribution function with the initial argument from org()
		self.recur(self.initarg)
		return self.biglist
	def org(self, beads):
		"""reorganizes bead information into initial recur() argument"""
# initarg is an array of quadruplets, initially empty
		self.initarg = []
# scan bead information. n is the number of beads with id i
		for i, n in enumerate(beads):
# is there at least one bead with bead id i ?
			if n:
				what = i
				howmany = n
# occurence is incremented each time a bead with the same id is distributed
				occurence = 1
# meaningless is occurence is 1. Else, is equal to the box index in which
# a bead with the same id has been distributed. Solves the problem
# of the equivalence of the beads with the same id.
				formerbox = 0
# a new quadruplet is added to initarg
				self.initarg.append([what, howmany, occurence, formerbox])
#		print "Initarg:", self.initarg
	def recur(self, arg):
		"""Recursively distributes beads according to argument arg (quadruplets)"""
#		print "enter recur with arg:", arg
		if len(arg) == 0:
#			print "Got it!"
# argument list is empty. All beads placed. Get box contents
			l = self.getcontents()
#			print "current contents:", l
# append the nex content of boxes to the biglist
			self.biglist.append(l)
#			print "New biglist:", self.biglist
		else:
#			print "Not at the end"
# consider here only the information about the beads with the lowest id
			head = arg[0]
# remember about the beads with other ids, if any
			tail = arg[1:]
#			print "head:", head
#			print "tail:", tail
# get quadruplet members
			what, howmany, occurence, formerbox = head
# If we deal with a bead with a new id (occurence == 1) then the first box
# to consider is the one with index 0. Else, the first box is the same
# as the one that was formerly filled
			startbox = 0 if occurence == 1 else formerbox
# tries all boxes with index starting at startbox, up to the last one
			for ibox in range(startbox, self.numboxes):
# get reference to the current box
				curbox = self.boxlist[ibox]
# The current box must be able to receive a new bead
				if curbox.pushable(what):
# place a new bead with what id
					curbox.push(what)
					print ("***", self.getcontents()) #  
# True if all beads with the current id (what) are distributed.
					if howmany == occurence:
# The argument for the next call to recur() is simply tail
						newarg = tail
# There are is beads with the same id to distribute
					else:
# Same id, same number, one more occurence, current box index -> new head
						newhead = [what, howmany, occurence+1, ibox]
# new arg from new head and old tail
						newarg = [newhead] + tail
# recursive call with new arg
					self.recur(newarg)
# remove the bead that was added before adding it into a next box (if any)
					curbox.pop()
					print ("***", self.getcontents()) #  

if __name__ == "__main__":
#	b0 = box.Box(4, [0, True, True, False, False])
#	b1 = box.Box(2, [0, True, True, True, False])
#	f = BoxFiller([b0, b1])
#	l = f.fill([0, 2, 3, 1, 0])

#	b0 = box.Box(1, [0, True, False, False, False])
#	b1 = box.Box(1, [0, False, True, False, False])
#	f = BoxFiller([b0, b1])
#	l = f.fill([0, 1, 1, 0, 0])

#	b0 = box.Box(4, [0, True, True, False, False])
#	f = BoxFiller([b0])
#	l = f.fill([0, 2, 2, 0, 0])

#	b0 = box.Box(4, [0, True, True, False, False])
#	f = BoxFiller([b0])
#	l = f.fill([0, 4, 0, 0, 0])

# 4 identical oxygens, one or two neighbors
	b0 = box.Box(4, [0, True, True, False, False])
# 2 identical nitrogens, one, two or three neighbors
	b1 = box.Box(2, [0, True, True, True, False])
# 6 VSA, 4 O, 2 N
	f = BoxFiller([b0, b1])
# 6 VSA, 4 O, 2 N, 2 atoms with 1 neighbor, 3 atoms with 2 neighbors, 1 atom with 3 neighbors
	l = f.fill([0, 2, 3, 1, 0])
	print ("Fillings:", l)
# Fillings: [[[1, 1, 2, 2], [2, 3]], [[1, 2, 2, 2], [1, 3]]]
# 2 solutions:
#     2 O with 1 neighbor, 2 O with 2 neighbors, 1 N with 2 neighbors, 1 N with 3 neighbors
#     1 O with 1 neighbor, 3 O with 2 neighbors, 1 N with 1 neighbor, 1 N with 3 neighbors
