
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import java.io.File;
import java.util.ArrayList;

import org.junit.Before;
import org.junit.Test;

import ecosim.Binning;
import ecosim.BinLevel;
import ecosim.MainVariables;
import ecosim.tree.Tree;
import ecosim.tree.InvalidTreeException;

/**
 *  Tests for Binning.
 *
 *  TestBinning.nwk is an ultrametric tree of five leaves, A to E, whose
 *  internal nodes have maximum leaf distances of 0.0024 (A, B), 0.0064
 *  (A, B, C), 0.0264 (A, B, C, D), and 0.0864 (all). At a sequence
 *  criterion level crit, a node is split when its maximum leaf distance
 *  exceeds 1 - crit, so the number of bins is one more than the number of
 *  those distances that exceed 1 - crit. None of the distances fall on a
 *  bin level.
 */
public class TestBinning {

    @Before
    public void setup () throws InvalidTreeException {
        tree = new Tree (new File ("build/tests/java/assets/TestBinning.nwk"));
        binning = new Binning (tree);
        binning.run ();
    }

    @Test
    public void testBinLevelsSorted () {
        for (int i = 1; i < Binning.binLevels.length; i ++) {
            assertTrue (
                "Bin levels are not strictly increasing.",
                Binning.binLevels[i] > Binning.binLevels[i - 1]
            );
        }
     }

    @Test
    public void testCrit () {
        ArrayList<BinLevel> bins = binning.getBins ();
        assertEquals (
            "Unexpected number of bin levels.",
            Binning.binLevels.length, bins.size ()
        );
        for (int i = 0; i < bins.size (); i ++) {
            assertEquals (
                "Unexpected sequence criterion level.",
                Binning.binLevels[i], bins.get (i).getCrit (),
                MainVariables.EPSILON
            );
        }
    }

    @Test
    public void testLevel () {
        for (BinLevel bin: binning.getBins ()) {
            assertEquals (
                "Unexpected number of bins at " + bin.getCrit () + ".",
                expectedLevel (bin.getCrit ()), (int)bin.getLevel ()
            );
        }
    }

    @Test
    public void testOutgroup () throws InvalidTreeException {
        // With E as the outgroup, E is no longer counted once the root is
        // split, which happens when the largest distance exceeds 1 - crit.
        Tree rerooted = new Tree (
            new File ("build/tests/java/assets/TestBinning.nwk")
        );
        rerooted.reroot ("E");
        Binning outgroupBinning = new Binning (rerooted);
        outgroupBinning.run ();
        for (BinLevel bin: outgroupBinning.getBins ()) {
            int expected = expectedLevel (bin.getCrit ());
            if (expected > 1) expected --;
            assertEquals (
                "Unexpected number of bins at " + bin.getCrit () + ".",
                expected, (int)bin.getLevel ()
            );
        }
    }

    /**
     *  The expected number of bins for TestBinning.nwk at a sequence
     *  criterion level.
     *
     *  @param crit The sequence criterion level.
     *  @return The expected number of bins.
     */
    private int expectedLevel (Double crit) {
        int level = 1;
        for (double distance: distances) {
            if (distance > 1.0d - crit) level ++;
        }
        return level;
    }

    private final double[] distances = { 0.0024d, 0.0064d, 0.0264d, 0.0864d };

    private Tree tree;
    private Binning binning;

}
