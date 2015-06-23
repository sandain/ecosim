
import static org.junit.Assert.assertEquals;

import java.io.File;

import org.junit.Before;
import org.junit.After;
import org.junit.Test;
import org.junit.Ignore;
import org.junit.runner.RunWith;
import org.junit.runners.JUnit4;

import ecosim.NewickTree;
import ecosim.NewickTreeNode;
import ecosim.InvalidNewickException;

public class TestNewickTree {

    @Before
    public void setup () throws InvalidNewickException {
        File treeFile = new File (
            "build/tests/java/assets/TestNewickTree.nwk"
        );
        tree = new NewickTree (treeFile);
    }

    @Test
    public void testIsValid () {
        assertEquals (
            "Invalid tree.",
            true,
            tree.isValid ()
        );
    }

    @Test
    public void testToString () {
        assertEquals (
            "Tree mismatch.",
            testNewickTree,
            tree.toString ()
        );
    }

    @Test
    public void testCompareTo () throws InvalidNewickException {
        NewickTree a = new NewickTree (testNewickTree);
        assertEquals (
            "Tree mismatch.",
            0,
            a.compareTo (tree)
        );
    }

    @Test
    public void testReroot () throws InvalidNewickException {
        NewickTree a = new NewickTree (testNewickTree);
        NewickTree b = new NewickTree (
            "(B:0.1,(A:0.1,((C:0.1,D:0.1):0.2,E:0.8):0.1):0.1):0.0;"
        );
        NewickTree c = new NewickTree (
            "((((A:0.1,B:0.2):0.1,E:0.8):0.2,D:0.1):0.05,C:0.05):0.0;"
        );
        a.reroot ("B");
        assertEquals (
            "Rerooted tree mismatch.",
            0,
            a.compareTo (b)
        );
        a.reroot ("C");
        assertEquals (
            "Rerooted tree mismatch.",
            0,
            a.compareTo (c)
        );
        a.reroot ("E");
        assertEquals (
            "Rerooted tree mismatch.",
            0,
            a.compareTo (tree)
        );
    }

    @Test
    public void testGetRoot () {
        NewickTreeNode root = tree.getRoot ();
        assertEquals (
            "Invalid root node.",
            true,
            root.isRootNode ()
        );
    }

    @Test
    public void testRemoveDescendant () throws InvalidNewickException {
        NewickTree a = new NewickTree (testNewickTree);
        NewickTree b = new NewickTree (
            "(((C:0.1,D:0.1):0.2,B:0.3):0.3,E:0.5):0.0;"
        );
        a.removeDescendant ("A" );
        assertEquals (
            "Remove descendant failed:\n" +
            a.toString () + "\n" + b.toString () + "\n",
            0,
            a.compareTo (b)
        );
    }

    @Test
    public void testRemoveOutgroup () throws InvalidNewickException {
        NewickTree a = new NewickTree (testNewickTree);
        NewickTree b = new NewickTree (
            "((A:0.1,B:0.2):0.1,(C:0.1,D:0.1):0.2):0.3;"
        );
        a.removeDescendant ("E");
        assertEquals (
            "Remove outgroup failed:\n" +
            a.toString () + "\n" + b.toString () + "\n",
            0,
            a.compareTo (b)
        );
    }


    private NewickTree tree;

    //       ┌─ A
    //     ┌─┤
    //     │ └── B
    // ┌───┤
    // │   │  ┌─ C
    // │   └──┤
    // │      └─ D
    // │
    // └───── E
    private final String testNewickTree =
        "(((A:0.10000,B:0.20000):0.10000,(C:0.10000,D:0.10000):0.20000):0.30000,E:0.50000):0.00000;";

}
