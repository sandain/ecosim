
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertSame;

import java.util.ArrayList;
import java.util.List;

import org.junit.Test;

import ecosim.BinLevel;
import ecosim.ParameterSet;
import ecosim.Summary;

/**
 *  Tests for Summary.
 */
public class TestSummary {

    @Test
    public void testSettersNotifyListeners () {
        Summary summary = new Summary ();
        final List<Summary> received = new ArrayList<Summary> ();
        summary.addListener (new Summary.Listener () {
            public void summaryChanged (Summary s) {
                received.add (s);
            }
        });
        summary.setNu (5);
        summary.setLength (1000);
        summary.setOutgroup ("outgroup");
        summary.setBins (new ArrayList<BinLevel> ());
        summary.setHillclimbing (new ParameterSet ());
        summary.setConfidenceInterval (new ParameterSet[] {
            new ParameterSet (), new ParameterSet ()
        });
        assertEquals ("Listener calls mismatch.", 6, received.size ());
        for (Summary s: received) {
            assertSame ("Listener was not passed the summary.", summary, s);
        }
        assertEquals ("Nu mismatch.", 5, (int)summary.getNu ());
        assertEquals ("Length mismatch.", 1000, (int)summary.getLength ());
        assertEquals (
            "Outgroup mismatch.", "outgroup", summary.getOutgroup ()
        );
    }

    @Test
    public void testRefreshListeners () {
        Summary summary = new Summary ();
        final int[] calls = new int[2];
        summary.addListener (new Summary.Listener () {
            public void summaryChanged (Summary s) {
                calls[0] ++;
            }
        });
        summary.addListener (new Summary.Listener () {
            public void summaryChanged (Summary s) {
                calls[1] ++;
            }
        });
        summary.refreshListeners ();
        assertEquals ("First listener calls mismatch.", 1, calls[0]);
        assertEquals ("Second listener calls mismatch.", 1, calls[1]);
    }

}
