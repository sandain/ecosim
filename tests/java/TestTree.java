
import static org.junit.Assert.assertEquals;

import java.io.File;

import org.junit.Before;
import org.junit.After;
import org.junit.Test;
import org.junit.Ignore;
import org.junit.runner.RunWith;
import org.junit.runners.JUnit4;

import ecosim.tree.Tree;
import ecosim.tree.Node;
import ecosim.tree.InvalidTreeException;

public class TestTree {

    @Before
    public void setup () throws InvalidTreeException {
        File treeFile = new File (
            "build/tests/java/assets/TestTree.nwk"
        );
        tree = new Tree (treeFile);
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
            testTree,
            tree.toString ()
        );
    }

    @Test
    public void testCompareTo () throws InvalidTreeException {
        Tree a = new Tree (testTree);
        assertEquals (
            "Tree mismatch.",
            0,
            a.compareTo (tree)
        );
    }

    @Test
    public void testReroot () throws InvalidTreeException {
        Tree a = new Tree (testTree);
        Tree b = new Tree (
            "(B:0.1,(A:0.1,((C:0.1,D:0.1):0.2,E:0.8):0.1):0.1):0.0;"
        );
        Tree c = new Tree (
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
        Node root = tree.getRoot ();
        assertEquals (
            "Invalid root node.",
            true,
            root.isRootNode ()
        );
    }

    @Test
    public void testRemoveDescendant () throws InvalidTreeException {
        Tree a = new Tree (testTree);
        Tree b = new Tree (
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
    public void testRemoveOutgroup () throws InvalidTreeException {
        Tree a = new Tree (testTree);
        Tree b = new Tree (
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


    private Tree tree;

    //       ┌─ A
    //     ┌─┤
    //     │ └── B
    // ┌───┤
    // │   │  ┌─ C
    // │   └──┤
    // │      └─ D
    // │
    // └───── E
    private final String testTree =
        "(((A:0.10000,B:0.20000):0.10000,(C:0.10000,D:0.10000):0.20000" +
        "):0.30000,E:0.50000):0.00000;";

}
