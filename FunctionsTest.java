package maygenjunit;

import java.util.Arrays;
import java.util.List;
import java.util.Set;

import org.junit.Before;
import org.junit.Test;
import static org.junit.Assert.assertTrue;
import org.openscience.cdk.Atom;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IBond.Order;

import com.google.common.collect.ArrayListMultimap;
import com.google.common.collect.ListMultimap;
import static org.junit.Assert.assertEquals;
import java.io.IOException;
import maygenjunit.Functions;



public class FunctionsTest {
	private IAtomContainer atomContainer;
	
	@Before
	/**
	 * This atomcontainer is built as the example input for all the test cases.
	 * @throws Exception
	 */
    public void setUp() throws Exception {
		IAtomContainer acontainer = new org.openscience.cdk.AtomContainer();
		acontainer.addAtom(new Atom("C"));
		acontainer.addAtom(new Atom("C"));
		acontainer.addAtom(new Atom("C"));
		acontainer.addAtom(new Atom("C"));
		acontainer.addAtom(new Atom("C"));
		acontainer.addAtom(new Atom("C"));

		acontainer.getAtom(0).setImplicitHydrogenCount(3);
		acontainer.getAtom(1).setImplicitHydrogenCount(3);
		acontainer.getAtom(2).setImplicitHydrogenCount(2);
		acontainer.getAtom(3).setImplicitHydrogenCount(2);
		acontainer.getAtom(4).setImplicitHydrogenCount(1);
		acontainer.getAtom(5).setImplicitHydrogenCount(1);
		
		acontainer.addBond(4, 5, Order.SINGLE);
		acontainer.addBond(3, 4, Order.SINGLE);
		
		this.atomContainer=acontainer;
    }
	
	@Test
	/**
	 * Function calculates equivalence classes using canonsym, then classify the atoms.
	 * The equivalence classes are enumerated based on implicit hydrogens and symmetry values
	 * of the atoms. 
	 * 
	 * @throws CloneNotSupportedException
	 * @throws CDKException
	 * @throws IOException
	 */
    public void test_ecenum() throws CloneNotSupportedException, CDKException, IOException {	
		ListMultimap<String, Integer> map= ArrayListMultimap.create();
		map.put("C12", 5);
		map.put("C16", 4);
		map.put("C25", 3);
		map.put("C21", 2);
		map.put("C33", 0);
		map.put("C33", 1);
		//System.out.println(Functions.ecenum(atomContainer));
        assertEquals(map,Functions.ecenum(atomContainer));
    }
	
	@Test
	/**
	 * The canon values are calculated by using CDK Canon class. That calculates the symmetry classes.
	 */
	public void test_canonsym() {
		long[] sym= {3,3,1,5,6,2};
		assertEquals(Arrays.toString(sym), Arrays.toString(Functions.canonsym(atomContainer)));
	}
	
	@Test
	/**
	 * This function is for checking the saturation of an atom in an atomcontainer.
	 * @throws CloneNotSupportedException
	 * @throws CDKException
	 * @throws IOException
	 */
	public void test_satcheck() throws CloneNotSupportedException, CDKException, IOException {
		assertTrue(Functions.satcheck(atomContainer, 0)==true);
	}
	
	@Test
	/**
	 * This function checks whether all the indices of an e.c. are saturated or not.
	 */
	public void test_classsatcheck() throws CloneNotSupportedException, CDKException, IOException {
		List<Integer> list = Arrays.asList(0,1);
		//System.out.println(Functions.classsatcheck(atomContainer, list));
		assertTrue(Functions.classsatcheck(atomContainer, list)==false);
	}
	
	@Test
	/**
	 * The function calculates the open sites of an atom. Here, atom at the 1th index has 1 open site.
	 */
	public void test_opencounter() throws CloneNotSupportedException, CDKException, IOException {
		//System.out.println(Functions.opencounter(atomContainer, 1));
		assertEquals(1,Functions.opencounter(atomContainer, 1));
	}
	
	@Test
	/**
	 * For two integers, the minimum one is returned. This function is used for the detection of maximum bond type
	 * between atom pairs based on their open sites.
	 */
	public void test_opendif() throws CDKException {
		assertEquals(2,Functions.opendif(2, 3));
	}
	
	@Test
	/**
	 * The bond order for a interaction between two atoms is calculated based on opendif value.
	 */
	public void test_ordselect() throws CDKException{
		assertEquals(1,Functions.ordselect(1).numeric().intValue());
	}
	
	@Test
	/**
	 * The function is simply for depiction of molecules.
	 */
	public void test_depict() throws CloneNotSupportedException, CDKException, IOException {
		Functions.depict(atomContainer,"..\\maygenjunit\\bin\\example.png");
	}
	
	@Test
	/**
	 * This function takes a an equivalence class to extend a molecule. For all the other
	 * equivalence classes, only one of their elements are used for the extension as the target
	 * atoms. After extending the molecule, this used target equivalence class is removed from
	 * the equivalence class list. Thus, only one element was used for the extension from all the
	 * equivalence classes.
	 */
	
	public void test_classext() throws CloneNotSupportedException, CDKException, IOException {
		ListMultimap<String, Integer> ecs=Functions.ecenum(atomContainer); //The equivalence class based classification of the atoms.
		Set<String> keys=ecs.keySet(); //The list of the e.c. labels.
		Functions.depict(atomContainer,"..\\maygenjunit\\bin\\init.png"); //Depicting the initial structure.
                Functions.classext(atomContainer,"C33",keys); //Extending the mol from its equivalence class "C33".
	}
}
