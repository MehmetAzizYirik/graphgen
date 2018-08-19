package maygenjunit;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.openscience.cdk.depict.DepictionGenerator;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.graph.GraphUtil;
import org.openscience.cdk.graph.invariant.Canon;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond;
import org.openscience.cdk.interfaces.IBond.Order;
import org.openscience.cdk.tools.manipulator.BondManipulator;

import com.google.common.collect.ListMultimap;
import com.google.common.collect.MultimapBuilder;


public class Functions {
	
	private IAtomContainer atomContainer;
	
	public Functions (IAtomContainer container) {
        this.atomContainer = container;
    }

	
	//Calculates the CDK canon symmetry array representing the symmetry class distribution of the molecule.
	public static long[] canonsym(IAtomContainer mol){
		int[][]	g = GraphUtil.toAdjList(mol);
		long[] sym= Canon.symmetry(mol, g);
		return sym;
	}
	
	//Saturation checker, checking the maximum number of connected bonds of atoms.
	public static boolean satcheck(IAtomContainer mol, int i) throws CloneNotSupportedException, CDKException, IOException{
		if ((mol.getAtom(i).getImplicitHydrogenCount()+ordsum(mol,i))>= (int)valences.get(mol.getAtom(i).getSymbol())){ 
			return false;
		}else{
			return true;
		}
	}
	
	//Check whether the class was saturated or not. SortedSet<Integer> c
	public static boolean classsatcheck(IAtomContainer mol, List<Integer> c)throws CloneNotSupportedException, CDKException, IOException{
		boolean check=true;
		for(int i=0; i<c.size();i++) {
			if(satcheck(mol,c.get(i))) {
				check=false;
			}
		}
		return check;
	}
	
	//Summation of the connected bond orders.
	public static int ordsum(IAtomContainer mol, int i){
		int count=0;
		for(IBond bond: mol.getConnectedBondsList(mol.getAtom(i))){
			count=count+bond.getOrder().numeric();
		}
		return count;
	}
	
	// Calculating open sites of atoms.
	public static int opencounter(IAtomContainer mol, int i)throws CloneNotSupportedException, CDKException, IOException{
		int open = valences.get(mol.getAtom(i).getSymbol()).intValue()- ordsum(mol,i) - mol.getAtom(i).getImplicitHydrogenCount(); 
		return open;
	}
	
	//To get the difference  between the open sites.
	public static int opendif( int i, int j) throws CDKException {
		if(i>j){
			return j;
		}else{
			return i;
		}
	}
	
	//To choose the bond order based on the open differences.
	public static Order ordselect( int i) throws CDKException{
		Order ord=null;
		if(i==1){
			ord= Order.SINGLE;
		}else if(i==2){
			ord= Order.DOUBLE;
		}else if(i>=3){
			ord= Order.TRIPLE;
		}
		return ord;
	}
	
	public static final Comparator<String> DESC_ORDER = new Comparator<String>() {
	    public int compare(String e1, String e2) { 
	        return e2.compareTo(e1);
	    }
	};
	
	//Increase bond order of a given pair of edge vertices.
	public static void bondadder(IAtomContainer mol, int f, int l)throws CloneNotSupportedException, CDKException, IOException {
		IBond add = mol.getBond(mol.getAtom(f), mol.getAtom(l)); //bond is a 1D array with two entries; edge vertices
		if(add == null){ //Bondmanipulator returns nullpointerexception.					
			mol.addBond(f, l, IBond.Order.SINGLE);
			mol.getBond(mol.getAtom(f), mol.getAtom(l)).setID("last");
		}
		else{
			BondManipulator.increaseBondOrder(add); 
			mol.getBond(mol.getAtom(f), mol.getAtom(l)).setID("last");
		}
	}
	
	// Molecule depiction
	public static void depict(IAtomContainer mol, String path) throws CloneNotSupportedException, CDKException, IOException{
		DepictionGenerator depict = new DepictionGenerator();
		depict.withCarbonSymbols().withSize(1000, 1000).withZoom(4).depict(mol).writeTo(path);
	}
	
	//This function removes all the  last added bonds from the atomcontainer.
	public static void lastbondclear(IAtomContainer mol) {
		List<IBond> bonds= new ArrayList<IBond>();
		for(IBond bond:mol.bonds()) {
			if(bond.getID()=="last") {
				bonds.add(bond);
				bond.setID("null");
			}
		}
		for(IBond bond: bonds) {
			mol.removeBond(bond);
		}
	}
	
	//The equivalent classes of molecules are ordered and enumerated in ascending order based on their open values and implicit hydrogens; as described in the paper. 
	public static ListMultimap<String, Integer> ecenum(IAtomContainer acontainer) throws CloneNotSupportedException, CDKException, IOException {
		//SortedSetMultimap<String,Integer> classes = TreeMultimap.create();
		ListMultimap<String, Integer> classes = MultimapBuilder.treeKeys(DESC_ORDER).arrayListValues().build();
		long[] sym=Functions.canonsym(acontainer);
		for(int i=0; i<acontainer.getAtomCount();i++){
			if(Functions.satcheck(acontainer, i)==true){	
				//classes.put(acontainer.getAtom(i).getSymbol()+opencounter(acontainer, i)+Long.valueOf(sym[i]).intValue(), i); //The open sites and the sym values are used for labelling. Less interactions for an atom means lower sym values.
				classes.put(acontainer.getAtom(i).getSymbol()+acontainer.getAtom(i).getImplicitHydrogenCount()+Long.valueOf(sym[i]).intValue(), i);
			}
		}	
		return classes;
	}
	
	public static Set<String> kys= new HashSet<String>();
	public static int soy1=0;
	public static int soy2=0;
	public static List<IAtomContainer> aclist1demo= new ArrayList<IAtomContainer>();
	public static  List<IAtomContainer> classext(IAtomContainer mol, String key,Set<String> keys) throws CloneNotSupportedException, CDKException, IOException { //Util.removehyd(mol, Util.hydcount(mol)); //Out of bounds error returned after removing Hs. How can we consider H atoms for CHs ?
		ListMultimap<String, Integer> ec=ecenum(mol); //Equivalence classes are computed.
		List<Integer> ecv=ec.get(key); //The atom indices of the source equivalence class. That is saturated.
		for(Integer source:ecv) {
			soy1++;//Just for the enumeration of the depictions.
			Iterator<String> iterator = keys.iterator();
			while(iterator.hasNext()){
				soy2++; //Just for the enumeration of the depictions.
				String ky = iterator.next();
				int target=ec.get(ky).get(0);
				if(ky!=key && target!=source) {
					if(satcheck(mol,source) && satcheck(mol,target)){ //Checks whether the source and target indices are saturated.
						Order ord= ordselect(opendif(opencounter(mol, source), opencounter(mol, target))); //open sites are counted. Then difference is calculated to have the maksimum order.
						for(int q=0;q<ord.numeric();q++) { //The order is used as the maximum order.
							bondadder(mol, source, target); //Simply increasing the bond order.
						}
						if(classsatcheck(mol,ecv)) { //Checks whether all the indices in the source e.v. are saturated or not.
							//depict(mol,"C:\\Users\\mehme\\Desktop\\maygen-outputs\\sat"+soy2+".png");
							depict(mol,"..\\maygenjunit\\bin\\sat"+soy1+"."+soy2+".png");
							aclist1demo.add(mol);
							lastbondclear(mol); //The last added bonds are removed for the next extension.
							iterator.remove(); //The target e.c. used for the extension is removed.
							classext(mol,key,keys);
						}
					}
				}
			}
		}
		return aclist1demo;
	}
	public static Map<String, Integer> valences; 
	static {
		//The atom valences from CDK.
		valences = new HashMap<String, Integer>();
			
		valences.put("C", 4);
		valences.put("N", 5);
		valences.put("O", 2);
		valences.put("S", 6);
		valences.put("P", 5);
		valences.put("F", 1);
		valences.put("I", 7);
		valences.put("Cl", 5);
		valences.put("Br", 5);
		valences.put("H", 1);
	}
}
