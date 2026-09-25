/* compile */
/* Unix: javac -classpath predictorc.jar:cdk-interfaces.jar:cdk-io.jar:cdk-core.jar:cdk-data.jar:cdk-standard.jar:cdk-valencycheck.jar:vecmath1.2-1.14.jar PredictorLSDC.java */
/* Windows: "C:\Program Files (x86)\Java\jdk1.7.0_02\bin\javac.exe" -classpath predictorc.jar;cdk-interfaces.jar;cdk-io.jar;cdk-core.jar;cdk-data.jar;cdk-standard.jar;cdk-valencycheck.jar;vecmath1.2-1.14.jar PredictorLSDC.java */

/* use */
/* Unix: java -classpath predictorc.jar:cdk-interfaces.jar:cdk-io.jar:cdk-ioformats.jar:cdk-core.jar:cdk-data.jar:cdk-standard.jar:cdk-valencycheck.jar:vecmath1.2-1.14.jar:jgrapht.jar:. PredictorLSDC rootname */
/* Windows: "C:\Program Files (x86)\Java\jre7\bin\java.exe" -classpath predictorc.jar;cdk-interfaces.jar;cdk-io.jar;cdk-ioformats.jar;cdk-core.jar;cdk-data.jar;cdk-standard.jar;cdk-valencycheck.jar;vecmath1.2-1.14.jar;jgrapht.jar;. PredictorLSDC rootname */

/*
 * PredictorLSDC processes a SD File (rootname0.sdf) and a pilot file (rootnameC.txt)
 * in order to produce a ranked SD file (rootnameD.sdf) and a result file (rootnameR.txt).
 * rootnameD.sdf contains the molecules in rootname0.sdf, but sorted so that the list of the distances between
 * the predicted C chemical shift vector of each molecule in rootnameD.sdf and the experimental chemical shift vector
 * in rootnameC.txt is in increasing order. The distance between two chemical shift vectors is the sum
 * of the absolute values of the chemical shift differences (predicted - experimental).
 * Each line in the pilot file contains 3 values: i) the index of an atom for which the chemical shift
 * has to be calculated from a molecule in rootname0.sdf (index start at 1) ii) the index of the
 * corresponding atom in the original LSD command file; it is equal to the index in the SD file
 * if there is no "hole" in the LSD atom numbering iii) the experimental chemical shift to which the predicted
 * one has to be compared. There is no line for a C atom for which no experimental chemical shift in known.
 * The rootnameR.txt file start with three lines that give the indexes of the SD C atoms (line starts with SD-C),
 * the indexes of the corresponding atoms for LSD (line starts with LSD-C), and the corresponding chemical shifts values
 * (line starts with CS-C) Each subsequent line is related to a single molecule in rootnameD.sdf.
 * The first value is the index of the molecule in rootname0.sdf (index start at 0). The first value is the distance between
 * experimental and predicted chemical shift vectors. The next values correspond to the predicted chemical shift vector
 * for the atoms in the order their index appear in the first line. 999.99 means that prediction failed.
 */

import java.io.*;
import java.util.List;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.Locale;
import java.util.Map;
import java.util.StringTokenizer;
import java.lang.String;
import java.util.Collections;

import org.openscience.cdk.DefaultChemObjectBuilder;
import org.openscience.cdk.Molecule;
import org.openscience.cdk.io.iterator.IteratingMDLReader;
import org.openscience.cdk.io.SDFWriter;
import org.openscience.cdk.interfaces.IAtom;
import org.openscience.cdk.interfaces.IAtomType;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.aromaticity.CDKHueckelAromaticityDetector;
import org.openscience.cdk.atomtype.CDKAtomTypeMatcher;
import org.openscience.cdk.tools.CDKHydrogenAdder;
import org.openscience.cdk.tools.manipulator.AtomTypeManipulator;
import org.openscience.nmrshiftdb.PredictionTool;


public class PredictorLSDC {

// A PairDI object stores a chemical shift distance (dkey) and the corresponding
// molecule index in rootname0.sdf
	public static class PairDI {
		private Double _dkey;
		private Integer _ival;
	
		public PairDI(){
		}
		public PairDI(Double dkey, Integer ival){
			_dkey = dkey;
			_ival = ival;
		}
		public Double getdkey() {
			return _dkey;
		}
		public Integer getival() {
			return _ival;
		}
	}

// add Implicit Hydrogens To Satisfy Valency to molecule mol
	public static void addImplicitHydrogensToSatisfyValency(IAtomContainer mol) throws CDKException {
		CDKAtomTypeMatcher matcher = CDKAtomTypeMatcher.getInstance(mol.getBuilder());           
		CDKHydrogenAdder hAdder = CDKHydrogenAdder.getInstance(mol.getBuilder());
		for (IAtom atom : mol.atoms()) {
			IAtomType type = matcher.findMatchingAtomType(mol, atom);
			if(type!=null){
				AtomTypeManipulator.configure(atom, type);
				hAdder.addImplicitHydrogens(mol, atom);
			}
		}
	}

// isCSkey() returns true if ch is a CS key. A CS key starts with CS and is followed by one or more digits.
// CS keys are used to store predicted chemical shift values in the SD property part of a molecule.
// key SDnn is associated to a five-fields information piece (see below) about the atom indexed nn (SD number).
	public static boolean isCSkey(String ch) {
		int l = ch.length();
		if(l < 3 ||  ch.charAt(0) != 'C' || ch.charAt(1) != 'S') {
			return false;
		}
		for(int i=2; i<l; i++) {
			if(!Character.isDigit(ch.charAt(i))) {
				return false;
			}
		}
		return true;
	}

// residue() calculates the distance between the experimental and predicted
// chemical shift vector in molecule mol. The required values are stored as SD properties
	public static double residue(IAtomContainer mol) {
		double result = 0.0;
		String key;
		String value;
		StringTokenizer st = null;
		// String lsdindexstr;
		String expstr;
		// String minstr;
		String meanstr;
		// String maxstr;
		double expval;
		double thval;
		double diff;
		Map<Object,Object> sdFields = mol.getProperties();
		if(sdFields != null){
			for (Object propKey : sdFields.keySet()) {
				key = (String)propKey;
				// System.out.println("Key: "+key);
				if (isCSkey(key)) {
// use the 5-string value that is associated to the current CS key
					value = (String)sdFields.get(propKey);
					// System.out.println("Value: "+value);
					st = new StringTokenizer(value);
					// lsdindexstr = st.nextToken();
// first item is not used here (LSD atom number)
					st.nextToken();
// second item is the experimental chemical shift value
					expstr = st.nextToken();
					// minstr = st.nextToken();
// third item is not used here. Minimum chemical shift value found by predict()
					st.nextToken();
// fourth item is the mean of the chemical shift values found by predict()
					meanstr = st.nextToken();
					// maxstr = st.nextToken();
// fifth item is not used here. Maximum chemical shift value found by predict()
					st.nextToken();
// get double type values for the experimental and predicted (theoretical) values
					expval = Double.parseDouble(expstr);
					thval = Double.parseDouble(meanstr);
					if(thval < 999.9) {
// A non-error value was found (error is 999.99)
						diff = Math.abs(expval - thval);
					} else {
// a very big difference means that something went wrong
						diff = 99999.99;
					}
					result += diff;
					// System.out.println("diff: "+key+" "+expval+" "+thval+" "+String.format(Locale.US, "%8.2f", diff));
				}
			}
		}
		// System.out.println("result: "+String.format(Locale.US, "%8.2f", result));
		return result;
	}

// writeResultHeader writes the first 3 lines in rootnameR.txt (file handle fhout)
// the rootnameC.txt file has nlines lines, whose first, second and third column respectively
// contain the values in SDIndex, LSDIndex and expShift
	public static void writeResultHeader(PrintWriter fhout, int nlines, 
		ArrayList<String> SDIndex, ArrayList<String> LSDIndex, ArrayList<String> expShift) {
		
		fhout.print("SD-C");
		for(int i=0; i<nlines; i++) {
			fhout.print(" "+SDIndex.get(i));
		}
		fhout.println();
		fhout.print("LSD-C");
		for(int i=0; i<nlines; i++) {
			fhout.print(" "+LSDIndex.get(i));
		}
		fhout.println();
		fhout.print("CS-C");
		for(int i=0; i<nlines; i++) {
			fhout.print(" "+expShift.get(i));
		}
		fhout.println();
	}

// writeResult() writes a line in rootnameR.txt (file handle fhout) for molecule mol
// whose index in rootname0.sdf is index, for which nlines experimental chemical shifts are known,
// and whose corresponding SD atom indexes are in SDIndex (as strings)
	public static void writeResult(PrintWriter fhout, Molecule mol, int index, int nlines, ArrayList<String> SDIndex) {
		Map<Object,Object> sdFields;
		String keyProp;
		String valueProp;
		String[] pieces;

// get all SD properties		
		sdFields = mol.getProperties();
// print molecule index
		fhout.print(index);
// get experimental to predicted distance as string
		valueProp = (String)sdFields.get("PREDICTORC_RESIDUE");
// write distance
		fhout.print(" " + valueProp);
// scan through SD atom numbers
		for(int i=0; i<nlines; i++) {
// build key
			keyProp = "CS" + SDIndex.get(i);
// get LSD index, experiment chemical shift, min, mean and max predicted chemical shift
			valueProp = (String)sdFields.get(keyProp);
// get property pieces
			pieces = valueProp.split("\\s");
// print predicted chemical shift (mean value)
			fhout.print(" " + pieces[3]);
		}
// ready for another molecule
		fhout.println();
	}
	
	public static void main(String[] args) throws Exception {
		if(args.length != 1){
// args: rootname
			java.lang.System.exit(1);
		}
// get handle to pilot file rootnameC.txt
		BufferedReader fhpilot = new BufferedReader(new InputStreamReader(new DataInputStream(new FileInputStream(args[0]+"C.txt"))));
		String strLine;
		String[] pieces;
		ArrayList<String> SDIndex = new ArrayList<String>();
		ArrayList<String> LSDIndex = new ArrayList<String>();
		ArrayList<String> expShift = new ArrayList<String>();
// scan through pilot file
		while ((strLine = fhpilot.readLine()) != null) {
			pieces = strLine.split("\\s");
			SDIndex.add(pieces[0]);
			LSDIndex.add(pieces[1]);
			expShift.add(pieces[2]);
		}
		int nlines = LSDIndex.size();
		// for(int i=0; i<nlines; i++) {
			// System.out.println(LSDIndex.get(i)+" "+SDIndex.get(i)+" "+expShift.get(i));
		// }
// get handle to structure file sootname0.sdf
		IteratingMDLReader sdin = new IteratingMDLReader(new FileReader(new File(args[0]+"0.sdf")), DefaultChemObjectBuilder.getInstance());
		Molecule mol = null;
		double[] result;
		IAtom curAtom;
		int curiSD;
		// int curiLSD;
		// double curdexp;
		String curSD;
		String curLSD;
		String curexp;
		String keyProp;
		String valueProp;
		double error;
		PredictionTool predictor = new PredictionTool();
		ArrayList<Molecule> allMols = new ArrayList<Molecule>();
		List<PairDI> liste = new ArrayList<PairDI>();
		int isol = 0;
		double[] badResult = {999.99, 999.99, 999.99};
// scan through molecules
		while (sdin.hasNext()) {
// get current molecule
			mol = (Molecule)sdin.next();
// make it usable for the predictor
			addImplicitHydrogensToSatisfyValency(mol);
			CDKHueckelAromaticityDetector.detectAromaticity(mol);
// scan through pilot file information in order to supplement the current molecule
// with atomic experimental and predicted chemical shift data
			for(int i=0; i<nlines; i++) {
				curSD = SDIndex.get(i);
				curLSD = LSDIndex.get(i);
				curexp = expShift.get(i);
				curiSD = Integer.parseInt(curSD);
				// curiLSD = Integer.parseInt(curLSD);
				// curdexp = Double.parseDouble(curexp);
// get CDK atom index from SD atom index
				curAtom = mol.getAtom(curiSD-1);
// build key for the current SD atom index
				keyProp = "CS" + curSD;
// start property value with corresponding LSD number and exprimental chemical shift
				valueProp = curLSD + " " + curexp + " ";
				try {
// try to get a predicted chemical shift for the current atom in the current molecule
					result = predictor.predict(mol, curAtom);
				}
				catch (Exception e) {
// assign bad values to chemical shifts if prediction went wrong
					result = badResult;
					// e.printStackTrace();
					// System.out.format(Locale.US, "%3d:%3d%3d%8.2f%8.2f%8.2f%8.2f\n", isol+1, curiLSD, curiSD, curdexp, result[0], result[1], result[2]);
				}
// build the end of the current CS property
				valueProp += String.format(Locale.US, "%.2f %.2f %.2f", result[0], result[1], result[2]);
				// System.out.println(keyProp);
				// System.out.println(valueProp);
// add new property to current molecule
				mol.setProperty(keyProp, valueProp);
			}
// get distance between experimental and predicted chemical shift vectors
			error = residue(mol);
// appends distance as a new molecule SD property
			mol.setProperty("PREDICTORC_RESIDUE", String.format(Locale.US, "%.2f", error));
			// allMols.add((Molecule)mol.clone());
// appends the current molecule to the list of molecules
			allMols.add(mol);
// append the new (distance, index) pair in the list of (distance, index) pairs
// for index reordered according to the ascendent order of distances
			liste.add(new PairDI(error, isol));
// ready for next molecule
			isol++;
			// sdout.write(mol);
		}
// no need for the initial molecule set
		sdin.close();
// get sorted order of molecules in allMols
		Collections.sort(liste, new Comparator<PairDI>(){
			public int compare(PairDI p1, PairDI p2) {
				return p1.getdkey().compareTo(p2.getdkey());
			}
		});
// create rootnameD.sdf. Could be avoided because I do not intend anymore
// to do something with this file, contrarily to what I considered first.
		BufferedWriter fw = new BufferedWriter(new FileWriter(new File(args[0]+"D.sdf")));
		SDFWriter sdout = new SDFWriter(fw);
// create result file rootnameR.txt
		PrintWriter fhout = new PrintWriter(new FileWriter(new File(args[0]+"R.txt")));
// print first 3 lines in result file rootnameR.txt
		writeResultHeader(fhout, nlines, SDIndex, LSDIndex, expShift);
		int index;
		PairDI p = null;
// scans through the re-ordered list of (distance, index) pairs
		for(int i=0; i<liste.size(); i++) {
// get pair
			p = liste.get(i);
// get index of molecule, most likely first
			index = Integer.valueOf(p.getival());
// get molecule
			mol = allMols.get(index);
// remember at which place it initially was
			mol.setProperty("INITIAL_INDEX_FROM_LSD", index);
// write molecule to rootnameD.sdf
			sdout.write(mol);
// write a line about rank, distance and predictions in rootnameR.txt
			writeResult(fhout, mol, index, nlines, SDIndex);
		}
// done
		sdout.close();
 		fhout.close();
	}
}
