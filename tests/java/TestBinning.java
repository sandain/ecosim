
import static org.junit.Assert.assertEquals;

import java.io.File;
import java.util.ArrayList;

import org.junit.Before;
import org.junit.After;
import org.junit.Test;
import org.junit.Ignore;
import org.junit.runner.RunWith;
import org.junit.runners.JUnit4;

import ecosim.Binning;
import ecosim.BinLevel;
import ecosim.MasterVariables;
import ecosim.tree.Tree;
import ecosim.tree.InvalidTreeException;

public class TestBinning {

    @Before
    public void setup () throws InvalidTreeException {
        File treeFile = new File ("build/tests/java/assets/TestTree.nwk");
        Tree tree = new Tree (treeFile);
        binning = new Binning (tree);
    }

    @Test
    public void testCrit () {
        Double[] expected = {
            0.600d, 0.650d, 0.700d, 0.750d, 0.800d,
            0.810d, 0.820d, 0.830d, 0.840d, 0.850d,
            0.860d, 0.870d, 0.880d, 0.890d, 0.900d,
            0.910d, 0.920d, 0.930d, 0.940d, 0.950d,
            0.955d, 0.960d, 0.965d, 0.970d, 0.975d,
            0.980d, 0.985d, 0.990d, 0.995d, 1.000d
        };
        ArrayList<BinLevel> bins = binning.getBins ();
        for (int i = 0; i < bins.size (); i ++) {
            assertEquals (
                "Unexpected sequence criterion level.",
                bins.get (i).getCrit (), expected[i], MasterVariables.EPSILON
            );
        }
    }

    @Test
    public void testLevel () {
        Integer[] expected = {
            3, 3, 4, 4, 5, 5, 5, 5, 5, 5,
            5, 5, 5, 5, 5, 5, 5, 5, 5, 5,
            5, 5, 5, 5, 5, 5, 5, 5, 5, 5
        };
        ArrayList<BinLevel> bins = binning.getBins ();
        for (int i = 0; i < bins.size (); i ++) {
            assertEquals (
                "Unexpected number of bins.",
                bins.get (i).getLevel (), expected[i]
            );
        }
    }

    private Binning binning;

}
