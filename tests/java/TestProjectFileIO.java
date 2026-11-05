
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertTrue;

import java.io.File;
import java.io.IOException;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.Locale;

import org.junit.After;
import org.junit.Before;
import org.junit.Test;

import ecosim.Binning;
import ecosim.BinLevel;
import ecosim.MainVariables;
import ecosim.ParameterEstimate;
import ecosim.ParameterSet;
import ecosim.ProjectFileIO;
import ecosim.tree.InvalidTreeException;
import ecosim.tree.Tree;

/**
 *  Tests for ProjectFileIO.
 *
 *  A project with a phylogeny and binning results is saved to a temporary
 *  file, then loaded back.
 */
public class TestProjectFileIO {

    @Before
    public void setup () throws IOException, InvalidTreeException {
        mainVariables = new MainVariables ();
        tree = new Tree (new File ("build/tests/java/assets/TestBinning.nwk"));
        binning = new Binning (tree);
        binning.run ();
        projectFile = File.createTempFile ("TestProjectFileIO", ".xml");
        ProjectFileIO saved = new ProjectFileIO (
            mainVariables, null, NU, LENGTH, OUTGROUP, tree, binning,
            null, null, null, null, null, null
        );
        saved.save (projectFile);
        loaded = new ProjectFileIO (mainVariables, null);
        loaded.load (projectFile);
    }

    @After
    public void cleanup () {
        projectFile.delete ();
    }

    @Test
    public void testPhylogeny () {
        assertEquals ("Nu mismatch.", NU, loaded.getNu ());
        assertEquals ("Length mismatch.", LENGTH, loaded.getLength ());
        assertEquals ("Outgroup mismatch.", OUTGROUP, loaded.getOutgroup ());
        assertNotNull ("Tree not loaded.", loaded.getTree ());
        assertEquals (
            "Tree mismatch.", tree.toString (), loaded.getTree ().toString ()
        );
    }

    @Test
    public void testBinning () {
        assertNotNull ("Binning not loaded.", loaded.getBinning ());
        ArrayList<BinLevel> expected = binning.getBins ();
        ArrayList<BinLevel> bins = loaded.getBinning ().getBins ();
        assertEquals ("Bin levels mismatch.", expected.size (), bins.size ());
        assertTrue ("No bin levels loaded.", bins.size () > 0);
        for (int i = 0; i < bins.size (); i ++) {
            assertEquals (
                "Sequence criterion level mismatch.",
                expected.get (i).getCrit (), bins.get (i).getCrit (),
                MainVariables.EPSILON
            );
            assertEquals (
                "Number of bins mismatch.",
                expected.get (i).getLevel (), bins.get (i).getLevel ()
            );
        }
    }

    @Test
    public void testEncoding () throws IOException, InvalidTreeException {
        // The file is written in UTF-8, whatever the platform's default
        // charset, and a non-ASCII name is loaded back unchanged.
        String name = "\u00c4rz";
        Tree nonAscii = new Tree ("((A:0.1,B:0.1):0.1," + name + ":0.2);");
        ProjectFileIO saved = new ProjectFileIO (
            mainVariables, null, 2, LENGTH, name, nonAscii, null,
            null, null, null, null, null, null
        );
        saved.save (projectFile);
        String text = new String (
            Files.readAllBytes (projectFile.toPath ()), StandardCharsets.UTF_8
        );
        assertTrue (
            "Name not written in UTF-8.",
            text.contains ("<outgroup value=\"" + name + "\"/>")
        );
        ProjectFileIO reloaded = new ProjectFileIO (mainVariables, null);
        reloaded.load (projectFile);
        assertEquals ("Outgroup mismatch.", name, reloaded.getOutgroup ());
    }

    @Test
    public void testLocale () {
        // Numbers are written the same way whatever the default locale, so
        // a project saved in a locale with a decimal comma loads back
        // unchanged.
        ParameterEstimate estimate = new ParameterEstimate (LENGTH, binning);
        estimate.setOmega (-0.125d, 7.0d);
        estimate.setSigma (-0.5d, 10.0d);
        estimate.setResult (new ParameterSet (64L, 0.125d, 0.5d, 0.0d));
        estimate.setHasRun (true);
        Locale defaultLocale = Locale.getDefault ();
        try {
            Locale.setDefault (Locale.GERMANY);
            ProjectFileIO saved = new ProjectFileIO (
                mainVariables, null, NU, LENGTH, OUTGROUP, tree, binning,
                estimate, null, null, null, null, null
            );
            saved.save (projectFile);
        }
        finally {
            Locale.setDefault (defaultLocale);
        }
        ProjectFileIO reloaded = new ProjectFileIO (mainVariables, null);
        reloaded.load (projectFile);
        ParameterEstimate result = reloaded.getEstimate ();
        assertNotNull ("Estimate not loaded.", result);
        assertEquals (
            "Omega slope mismatch.", -0.125d, result.getOmega ()[0], TOLERANCE
        );
        assertEquals (
            "Omega intercept mismatch.", 7.0d, result.getOmega ()[1], TOLERANCE
        );
        assertEquals (
            "Sigma slope mismatch.", -0.5d, result.getSigma ()[0], TOLERANCE
        );
        assertEquals (
            "Sigma intercept mismatch.",
            10.0d, result.getSigma ()[1], TOLERANCE
        );
        assertEquals (
            "Omega mismatch.",
            0.125d, result.getResult ().getOmega (), TOLERANCE
        );
        assertEquals (
            "Sigma mismatch.", 0.5d, result.getResult ().getSigma (), TOLERANCE
        );
    }

    private static final Integer NU = 4;
    private static final Integer LENGTH = 1000;
    private static final String OUTGROUP = "E";
    private static final double TOLERANCE = 1.0e-9d;

    private MainVariables mainVariables;
    private Tree tree;
    private Binning binning;
    private File projectFile;
    private ProjectFileIO loaded;

}
