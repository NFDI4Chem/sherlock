"""Management of OS and system specific variables"""
import sys
import os

class Defaults:
	"""container for OS and system specific values"""
	def __init__(self):
		self.statuslist = "statuslist.txt"
		self.datafolder = os.path.normpath("../LSD/Data")
		self.lsdbin = os.path.normpath("../LSD")
		self.rootname = "file"
		self.isWin = sys.platform == "win32"
		self.lsd_extension = r".lsd" if self.isWin else ""
		self.nulldev = "NUL" if self.isWin else "/dev/null"
		self.solncounter = "solncounter"
		self.stopfile = "stoplsd"
		self.viewerpath = '"C:\\Program Files\\Ghostgum\\gsview\\gsview64"' if self.isWin else "evince"
		self.predictorspath = os.path.normpath("../Predict")
		self.predictorsbat = {}
		self.predictorsbat["C"] = "predictLSDC"
		self.java = '"C:\\Program Files\\Java\\jre1.8.0_321\\bin\\java.exe"' if self.isWin else "java"
		jarsC = ["predictorc.jar", "cdk-interfaces.jar", "cdk-io.jar", "cdk-ioformats.jar", "cdk-core.jar", "cdk-data.jar", "cdk-standard.jar", "cdk-valencycheck.jar", "vecmath1.2-1.14.jar", "jgrapht.jar"]
		self.javamemory = 1024 # megabytes memory for java
		jarsep = ";" if self.isWin else ":"
		fulljarsC = [os.path.join(self.predictorspath, x) for x in jarsC]
		fulljarsC.append(self.predictorspath)
		alljarsC = jarsep.join(fulljarsC)
		javamemopt = " -Xmx" + str(self.javamemory) + "m"
		self.predictors = {}
		self.predictors["C"] = self.java + javamemopt + " -classpath " + alljarsC + " PredictorLSDC"
