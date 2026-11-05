
import static org.junit.Assert.assertEquals;

import java.util.ArrayList;
import java.util.List;

import org.junit.Test;

import ecosim.Logger;

/**
 *  Tests for Logger.
 */
public class TestLogger {

    @Test
    public void testAppend () {
        Logger log = new Logger ();
        log.append ("one");
        log.appendln (" two");
        String expected = "one two" + System.getProperty ("line.separator");
        assertEquals ("Log mismatch.", expected, log.toString ());
        assertEquals ("Length mismatch.", expected.length (), log.length ());
    }

    @Test
    public void testListeners () {
        Logger log = new Logger ();
        final List<String> first = new ArrayList<String> ();
        final List<String> second = new ArrayList<String> ();
        log.addListener (new Logger.Listener () {
            public void textAppended (String text) {
                first.add (text);
            }
        });
        log.addListener (new Logger.Listener () {
            public void textAppended (String text) {
                second.add (text);
            }
        });
        log.append ("one");
        log.append ("two");
        List<String> expected = new ArrayList<String> ();
        expected.add ("one");
        expected.add ("two");
        assertEquals ("First listener mismatch.", expected, first);
        assertEquals ("Second listener mismatch.", expected, second);
    }

    @Test
    public void testListenerAddedLater () {
        // A listener only receives text appended after it was added.
        Logger log = new Logger ();
        log.append ("before");
        final List<String> received = new ArrayList<String> ();
        log.addListener (new Logger.Listener () {
            public void textAppended (String text) {
                received.add (text);
            }
        });
        log.append ("after");
        assertEquals ("Listener calls mismatch.", 1, received.size ());
        assertEquals ("Listener text mismatch.", "after", received.get (0));
    }

}
