
import static org.junit.Assert.assertEquals;

import org.junit.Before;
import org.junit.After;
import org.junit.Test;
import org.junit.Ignore;
import org.junit.runner.RunWith;
import org.junit.runners.JUnit4;

import ecosim.Sequence;

public class TestSequence {

    @Before
    public void setup () {
        seq = new Sequence ();
        seq.setIdentifier (testIdentifier);
        seq.setDescription (testDescription);
        seq.setSequence (testSequence);
    }

    @Test
    public void testIdentifier () {
        assertEquals ("Identifier not equal.", seq.getIdentifier (), testIdentifier);
    }

    @Test
    public void testDescription () {
        assertEquals ("Description not equal.", seq.getDescription (), testDescription);
    }

    @Test
    public void testSequence () {
        assertEquals ("Sequence not equal.", seq.getSequence (), testSequence);
    }

    private String testIdentifier = "test_sequence";
    private String testDescription = "Test description";
    private String testSequence = "acgttgca";

    private Sequence seq;

}
