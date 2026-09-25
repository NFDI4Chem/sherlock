"""Box class definition"""

class Box:
	"""Box: simulates a group of identical VSAs
	
	A box contains cells. Each cell is either empty or contains a bead.
	A bead has a numberic id painted on. Beads can be added and removed
	as if the box behaved like a stack: last in first out.
	A bead cannot be added if the box is full or if its number
	is not compatible with the box.
	A permission table (array of boolean values) indicates if
	a given id (index in permission table) is compatible
	with the box.
	"""
	def __init__(self, size, permissions, verbose=False):
		"""A box is created with a size, permissions and optional verbosity"""
		self.size = size
		self.permissions = permissions
		self.verbose = verbose
# a box is initially empty
		self.level = 0
		self.content = []
	def pushable(self, v):
		"""Can a bead with id v be added to the box?"""
		return self.level < self.size and self.permissions[v]
	def push(self, v):
		"""bead with id v is added to the box. Beware: No control!!"""
		self.content.append(v)
		self.level += 1
	def pop(self):
		"""remove the last added bead. Beware: No control!!"""
		self.content.pop()
		self.level -= 1
	def get(self):
		"""Returns the reference of the content of the box"""
		return self.content
	def __repr__(self):
		"""Box representation"""
		msg = repr(self.content)
		if self.verbose:
			msg += " level: %d size: %d permissions: %s" % (self.level, self.size, repr(self.permissions[1:]))
		return msg

if __name__ == "__main__":
	b = Box(2, [0, True, True, False, False], False)
	print ("b: ", b)
	print ("3 can be pushed? ", b.pushable(3))
	print ("2 can be pushed? ", b.pushable(2))
	print ("Pushing now 2.")
	b.push(2)
	print ("b: ", b)
	print ("1 can be pushed? ", b.pushable(1))
	print ("Pushing now 1.")
	b.push(1)
	print ("b: ", b)
	print ("1 can be pushed? ", b.pushable(1))
	print ("Poping now.")
	b.pop()
	print ("b: ", b)
	print ("Poping now.")
	b.pop()
	print ("b: ", b)

