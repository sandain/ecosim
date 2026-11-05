
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertFalse;
import static org.junit.Assert.fail;

import java.io.File;

import org.junit.Test;

import ecosim.Fasta;
import ecosim.InvalidFastaException;

public class TestFasta {

    @Test
    public void testOutgroup () throws InvalidFastaException {
        Fasta fasta = new Fasta (asset ("TestFasta.fa"));
        assertEquals (
            "Outgroup mismatch.", "test_sequence", fasta.getOutgroup ()
        );
    }

    @Test
    public void testLength () throws InvalidFastaException {
        Fasta fasta = new Fasta (asset ("TestFasta.fa"));
        assertEquals ("Length mismatch.", 8, fasta.getLength ());
    }

    @Test
    public void testSize () throws InvalidFastaException {
        Fasta fasta = new Fasta (asset ("TestFasta.fa"));
        assertEquals ("Size mismatch.", 3, fasta.size ());
     }

    @Test
    public void testWindowsLineEndings () throws InvalidFastaException {
        Fasta fasta = new Fasta (asset ("TestFastaCRLF.fa"));
        assertEquals (
            "Outgroup mismatch.", "test_sequence", fasta.getOutgroup ()
        );
        assertEquals ("Length mismatch.", 8, fasta.getLength ());
        assertEquals ("Size mismatch.", 3, fasta.size ());
    }

    @Test (expected = InvalidFastaException.class)
    public void testDuplicateIdentifier () throws InvalidFastaException {
        new Fasta (asset ("TestFastaDuplicate.fa"));
    }

    @Test (expected = InvalidFastaException.class)
    public void testUnaligned () throws InvalidFastaException {
        new Fasta (asset ("TestFastaUnaligned.fa"));
     }

    @Test (expected = InvalidFastaException.class)
    public void testNoHeader () throws InvalidFastaException {
        new Fasta (asset ("TestFastaNoHeader.fa"));
     }

    @Test (expected = InvalidFastaException.class)
    public void testNoSequences () throws InvalidFastaException {
        new Fasta (asset ("TestFastaEmpty.fa"));
    }

    @Test (expected = InvalidFastaException.class)
    public void testNoSequenceData () throws InvalidFastaException {
        new Fasta (asset ("TestFastaNoData.fa"));
    }

    @Test (expected = InvalidFastaException.class)
    public void testNoIdentifier () throws InvalidFastaException {
        new Fasta (asset ("TestFastaNoIdentifier.fa"));
    }

    @Test
    public void testMissingFile () {
        File missing = asset ("TestFastaMissing.fa");
        try {
            new Fasta (missing);
            fail ("A missing file was accepted.");
        }
        catch (InvalidFastaException e) {
            // Expected.
        }
        assertFalse ("A missing file was created.", missing.exists ());
    }

    private File asset (String name) {
        return new File ("build/tests/java/assets/" + name);
    }

}
