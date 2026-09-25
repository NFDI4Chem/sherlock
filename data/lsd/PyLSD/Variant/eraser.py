"""LSD input and output file erasing before each run of lsd.py"""
import os
import sys
import re

class Eraser:
	"""File erasing and OS dependent variables. Files are in datafolder"""
	def __init__(self, defaults):
# initial part of all relevent files
		self.rootname = defaults.rootname
# where the action takes place
		self.datafolder = os.path.normpath(defaults.datafolder)
		# localfiles = [files for root, dir, files in os.walk(datafolder) if root==datafolder][0]
# get the name of local files
		localfiles = os.listdir(self.datafolder)
		print ("all files:", localfiles)

# special case for Win systems
		self.isWin = defaults.isWin
		print ("isWin:", self.isWin)
# special LSD input file extension
		self.lsd_extension = defaults.lsd_extension
# assemble regular expression for eraser target files
		expr0 = self.rootname + "[0-9]+"
		expr1 = expr0 + self.lsd_extension
		expr2 = expr0 + r"\.sol"
		expr = "^(" + expr1 + "|" + expr2 + ")$" 
		print ("expr:", expr)
# remember regular expression
		self.matchexpr = expr

# retain only the local file names that matches with expected LSD files
		targets = filter(lambda x: self.isTarget(x), localfiles)
# remember the list of the paths of the files that will be deleted
		self.targets = [os.path.join(self.datafolder, target) for target in targets]
		print ("Eraser targets:", self.targets)

	def isTarget(self, filename):
		"""Predicate for file name matching"""
		m = re.match(self.matchexpr, filename)
		return bool(m)
	
	def erase(self):
		"""does the job"""
		for file in self.targets:
			os.unlink(file)
			# print file, os.stat(file).st_size

if __name__ == "__main__":
	import defaults
	
	default = defaults.Defaults
	eraser = Eraser(default)
	print ("file1.lsd:", eraser.isTarget("file1.lsd"))
	print ("file1.sol:", eraser.isTarget("file1.sol"))
	print ("file.lsd:", eraser.isTarget("file.lsd"))
	print ("toto.sol:", eraser.isTarget("toto.sol"))
	eraser.erase()
