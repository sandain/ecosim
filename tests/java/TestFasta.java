
import static org.junit.Assert.assertEquals;

import java.io.File;

import org.junit.Before;
import org.junit.After;
import org.junit.Test;
import org.junit.Ignore;
import org.junit.runner.RunWith;
import org.junit.runners.JUnit4;

import ecosim.Fasta;
import ecosim.Sequence;
import ecosim.InvalidFastaException;

public class TestFasta {

    @Before
    public void setup () throws InvalidFastaException {
        File fastaFile = new File ("build/tests/java/assets/TestFasta.fa");
        fasta = new Fasta (fastaFile);
    }

    @Test
    public void testOutgroup () throws InvalidFastaException {
        Sequence outgroup = fasta.getOutgroup ();
        assertEquals ("Identifier mismatch.", "test_sequence", outgroup.getIdentifier ());
        assertEquals ("Description mismatch.", "Test description", outgroup.getDescription ());
        assertEquals ("Sequence mismatch.", "acgttgca", outgroup.getSequence ());
    }

    @Test
    public void testNextSequence () throws InvalidFastaException {
        Sequence seq;
        // Check the first sequence.
        seq = fasta.nextSequence ();
        assertEquals ("Identifier mismatch.", "gb|CP000239.1|:2536947-2539214", seq.getIdentifier ());
        // Check the second sequence.
        seq = fasta.nextSequence ();
        assertEquals ("Identifier mismatch.", "gb|CP000240.1|:c28351-26084", seq.getIdentifier ());
    }

    @After
    public void teardown () throws InvalidFastaException {
        fasta.close ();
    }

    private Fasta fasta;

}
