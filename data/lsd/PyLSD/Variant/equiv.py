"""Handling of equivalence classes"""

class Class:
	"""A single equivalence class"""
	def __init__(self, first):
		"""Class initialization by means of its first representative"""
		self._first = first
		self._members = [first]
	def first(self):
		"""get first class representative"""
		return self._first
	def members(self):
		"""get list of class members"""
		return self._members
	def add(self, o):
		"""add a new element in class"""
		self._members.append(o)

class Classer:
	"""An equivalent class builder
	
	An object to be classed must have an equiv() method
	that establishes an equivalence relationship between
	itself and any other object of the same nature.
	"""
	def __init__(self, objects):
		"""Initialize a class builder with the list of objects to be classed."""
		self._objects = objects
# list of classes is initially empty
		self._classes = []
# scan object list
		for o in objects:
# the current object has not been found to belong to a known class yet
			found = False
# scan through known classes
			for c in self._classes:
# get first representative of the current class
				f = c.first()
# if the current object is equivalent to a representative of a known class
				if o.equiv(f):
# ... a known class has been found for the object
					found = True
					break
			if found:
# add the current object to its class
				c.add(o)
			else:
# create a new class with the current object
				c = Class(o)
# ... and append this new class to the list of known classes
				self._classes.append(c)
	
	def classes(self):
		"""get list of equivalence classes"""
		return self._classes

if __name__ == "__main__":
	class MyInt:
		"""A class of integers, with the equality of the modulo 10 as equivalence relationship"""
		def __init__(self, v):
			"""create new MyInt object from int value"""
			self._value = v
		def equiv(self, other):
			"""equivalence relationship"""
			return self._value % 10 == other._value % 10
		def value(self):
			"""get int value"""
			return self._value
	
	i = MyInt(17)
	j = MyInt(27)
	print ("Are 17 and 27 equivalent?", i.equiv(j))

# new class with 17 and 27
	z = Class(i)
	z.add(j)
# print values in class
	print ([x.value() for x in z.members()])

# create b, list of MyInt with integers from 0 to 24
	a = range(0, 25)
	b = [MyInt(x) for x in a]
# determines classes from b
	c = Classer(b)
# make list of class contents (as integers)
#	d = [map(lambda y: y.value(), x.members()) for x in c.classes()]
	d = [[y.value() for y in x.members()] for x in c.classes()]
	print (d)
