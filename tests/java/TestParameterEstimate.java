
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertTrue;

import org.junit.Test;

import ecosim.Binning;
import ecosim.BinLevel;
import ecosim.ParameterEstimate;
import ecosim.ParameterSet;

/**
 *  Tests for ParameterEstimate.
 *
 *  The binning results are built directly, with a sequence length of 1000,
 *  so that a point x sequence differences apart has a sequence criterion
 *  level of 1 - x / 1000. The number of bins at each point is 2^y, so the
 *  points lie exactly on the lines y = 10 - 0.5 x (sigma 0.5) for x <= 8
 *  and y = 7 - 0.125 x (omega 0.125) for x >= 8. The lines cross at y = 6,
 *  so npop is 2^6 = 64.
 */
public class TestParameterEstimate {

    @Test
    public void testTwoLines () {
        ParameterSet result = estimate (twoLines (
            new int[] { 0, 2, 4, 6, 8, 16, 24, 32, 40, 48, 56 }
        ));
        assertEquals ("Omega mismatch.", 0.125d, result.getOmega (), TOLERANCE);
        assertEquals ("Sigma mismatch.", 0.5d, result.getSigma (), TOLERANCE);
        assertEquals ("Npop mismatch.", 64L, (long)result.getNpop ());
    }

    @Test
    public void testSpacing () {
        // The same two lines sampled with fewer points on the sigma line.
        // The estimate must not depend on where the points fall.
        ParameterSet result = estimate (twoLines (
            new int[] { 0, 4, 8, 16, 24, 32, 40, 48, 56 }
        ));
        assertEquals ("Omega mismatch.", 0.125d, result.getOmega (), TOLERANCE);
        assertEquals ("Sigma mismatch.", 0.5d, result.getSigma (), TOLERANCE);
        assertEquals ("Npop mismatch.", 64L, (long)result.getNpop ());
    }

    @Test
    public void testTooFewPoints () {
        // Two lines need at least four points between the flat ends.
        ParameterSet result = estimate (bins (
            new int[] { 0, 8, 16 }, new int[] { 1024, 64, 32 }
        ));
        assertInvalid (result);
    }

    @Test
    public void testNoValidSplit () {
        // A shallow line at low divergence followed by a steep line: no
        // split has the steeper line at low divergence.
        ParameterSet result = estimate (bins (
            new int[] { 0, 8, 16, 18, 20, 22, 24, 26, 28, 30, 32 },
            new int[] { 1024, 512, 256, 128, 64, 32, 16, 8, 4, 2, 1 }
        ));
        assertInvalid (result);
    }

    /**
     *  Assert that an estimate is not valid: omega and sigma are NaN, and
     *  npop is 1.
     *
     *  @param result The estimate.
     */
    private void assertInvalid (ParameterSet result) {
        assertTrue ("Omega is not NaN.", result.getOmega ().isNaN ());
        assertTrue ("Sigma is not NaN.", result.getSigma ().isNaN ());
        assertEquals ("Npop mismatch.", 1L, (long)result.getNpop ());
    }

    /**
     *  Run the parameter estimate on the binning results, and check that it
     *  ran.
     *
     *  @param binning The binning results.
     *  @return The estimate.
     */
    private ParameterSet estimate (Binning binning) {
        ParameterEstimate estimate = new ParameterEstimate (LENGTH, binning);
        estimate.run ();
        assertTrue ("The parameter estimate did not run.", estimate.hasRun ());
        return estimate.getResult ();
    }

    /**
     *  Binning results on the two lines, at the given numbers of sequence
     *  differences.
     *
     *  @param differences The numbers of sequence differences, in order.
     *  @return The binning results.
     */
    private Binning twoLines (int[] differences) {
        int[] levels = new int[differences.length];
        for (int i = 0; i < differences.length; i ++) {
            int x = differences[i];
            double y = x <= 8 ? 10.0d - 0.5d * x : 7.0d - 0.125d * x;
            levels[i] = (int)Math.round (Math.pow (2.0d, y));
        }
        return bins (differences, levels);
    }

    /**
     *  Binning results with the given numbers of bins at the given numbers
     *  of sequence differences. The bin levels are added in increasing
     *  order of sequence criterion level, as Binning produces them.
     *
     *  @param differences The numbers of sequence differences, in order.
     *  @param levels The number of bins at each number of differences.
     *  @return The binning results.
     */
    private Binning bins (int[] differences, int[] levels) {
        Binning binning = new Binning ();
        for (int i = differences.length - 1; i >= 0; i --) {
            Double crit = 1.0d - differences[i] / (double)LENGTH;
            binning.addBinLevel (new BinLevel (crit, levels[i]));
        }
        return binning;
    }

    private static final int LENGTH = 1000;
    private static final double TOLERANCE = 1.0e-9d;

}
