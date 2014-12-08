/*
 *    Ecotype Simulation models the sequence diversity within a bacterial
 *    clade as the evolutionary result of net ecotype formation and periodic
 *    selection, yielding a certain number of ecotypes.
 *
 *    Copyright (C) 2013-2014  Jason M. Wood, Montana State University
 *
 *    This program is free software: you can redistribute it and/or modify
 *    it under the terms of the GNU General Public License as published by
 *    the Free Software Foundation, either version 3 of the License, or
 *    (at your option) any later version.
 *
 *    This program is distributed in the hope that it will be useful,
 *    but WITHOUT ANY WARRANTY; without even the implied warranty of
 *    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *    GNU General Public License for more details.
 *
 *    You should have received a copy of the GNU General Public License
 *    along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */


package ecosim;

import java.io.File;
import java.util.ArrayList;

/**
 *  The shared methods of the simulation used by SimulationCLI and
 *  SimulationGUI.
 *
 *  @author Jason M. Wood
 *  @copyright GNU General Public License
 */
public class Simulation {

    /**
     *  Simulation constructor.
     *
     *  @param masterVariables The MasterVariables object.
     *  @param fastaFile The fasta formated sequence file.
     *  @param newickFile The newick formated tree file.
     */
    public Simulation (MasterVariables masterVariables, File fastaFile,
        File newickFile) {
        this.masterVariables = masterVariables;
        this.fastaFile = fastaFile;
        this.newickFile = newickFile;
        log = masterVariables.getLog ();
        if (fastaFile != null && fastaFile.exists ()) {
            loadSequenceFile ();
        }
        if (newickFile != null && newickFile.exists ()) {
            loadTreeFile ();
        }
    }

    /**
     *  Load the project from a XML formated file.
     *
     *  @param file The file to load the project from.
     */
    public void loadProjectFile (File file) {
        ProjectFileIO projectFileIO = new ProjectFileIO (masterVariables);
        // Load the project file.
        projectFileIO.load (file);
        // Grab the loaded variables.
        nu = projectFileIO.getNu ();
        length = projectFileIO.getLength ();
        outgroup = projectFileIO.getOutgroup ();
        tree = projectFileIO.getTree ();
        binning = projectFileIO.getBinning ();
        bruteforce = projectFileIO.getBruteforce ();
        hillclimb = projectFileIO.getHillclimb ();
        omegaCI = projectFileIO.getOmegaCI ();
        sigmaCI = projectFileIO.getSigmaCI ();
        npopCI = projectFileIO.getNpopCI ();
        demarcation = projectFileIO.getDemarcation ();
    }

    /**
     *  Save the project to a XML formated file.
     *
     *  @param file The file to save the project to.
     */
    public void saveProjectFile (File file) {
        ProjectFileIO projectFileIO = new ProjectFileIO (
            masterVariables, nu, length, outgroup, tree, binning, bruteforce,
            hillclimb, omegaCI, sigmaCI, npopCI, demarcation
        );
        projectFileIO.save (file);
    }

    /**
     *  Load the fasta formated sequence file.
     */
    protected void loadSequenceFile () {
        // Verify that the fasta file exists.
        if (fastaFile == null || ! fastaFile.exists ()) {
            log.append ("Error loading the Fasta file!\n");
            return;
        }
        log.append (
            "Opening sequence file: " + fastaFile.getName () + "\n"
        );
        try {
            fasta = new Fasta (fastaFile);
            nu = fasta.size ();
            length = fasta.length ();
            outgroup = fasta.getIdentifier (0);
            // Output the number of sequences loaded.
            log.append (String.format (
                "  %d environmental sequences.\n" +
                "  %s is the outgroup.\n\n",
                nu, outgroup
            ));
        }
        catch (InvalidFastaException e) {
            System.out.println ("Error loading sequence file.");
        }
    }

    /**
     *  Load the newick formated tree file.
     */
    protected void loadTreeFile () {
        // Verify that the tree file exists.
        if (newickFile == null || ! newickFile.exists ()) {
            log.append ("Error loading the Newick tree!\n");
            return;
        }
        log.append (
            "Opening tree file: " + newickFile.getName () + "\n"
        );
        try {
            tree = new NewickTree (newickFile);
            tree.reroot (outgroup);
            // Store the tree in file called 'outtree'.
            newickFile = new File (
                masterVariables.getWorkingDirectory () + "outtree"
            );
            tree.save (newickFile);
        }
        catch (InvalidNewickException e) {
            System.out.println ("Error loading tree file.");
        }
    }

    /**
     *  Generate a tree using FastTree.
     */
    protected void generateTree () {
        log.append ("Generating a tree using FastTree...\n");
        // Store the tree in file called 'outtree'.
        newickFile = new File (
            masterVariables.getWorkingDirectory () + "outtree"
        );
        // Reroot the tree using the outgroup.
        Execs execs = masterVariables.getExecs ();
        execs.runFastTree (fastaFile, newickFile);
        try {
            tree = new NewickTree (newickFile);
            tree.reroot (outgroup);
            tree.save (newickFile);
        }
        catch (InvalidNewickException e) {
            System.out.println ("Error loading tree file.");
        }
    }

    /**
     *  Run the binning program.
     */
    protected void runBinning () {
        log.append ("Starting binning...\n");
        binning = new Binning (masterVariables, tree);
        binning.run ();
        // Verify that binning program ran correctly.
        if (! binning.hasRun ()) {
            log.append ("  Error running the binning program!\n");
            return;
        }
        // Output the results from binning.
        log.append ("The result from binning:\n");
        ArrayList<BinLevel> bins = binning.getBins ();
        for (int i = 0; i < bins.size (); i ++) {
            log.append ("  " + bins.get (i).toString () + "\n");
        }
        log.append ("\n");
    }

    /**
     *  Run the bruteforce program.
     */
    protected void runBruteforce () {
        log.append ("Starting bruteforce search...\n");
        bruteforce = new Bruteforce (masterVariables, nu, length, binning);
        while (bruteforce.getNumResults () < masterVariables.NUM_SUCCESSES) {
            bruteforce.run ();
            // Verify that bruteforce ran correctly.
            if (! bruteforce.hasRun ()) {
                log.append ("Error running the bruteforce search program!\n");
                return;
            }
            // Verify that there are enough bruteforce results.
            if (bruteforce.getNumResults () < masterVariables.NUM_SUCCESSES) {
                int criterion = masterVariables.getCriterion ();
                log.append (
                    "  Not enough results at the current criterion (" +
                    masterVariables.getCriterionLabel (criterion) + ")"
                );
                if (criterion > 2) {
                    log.append (", lowering the value.\n");
                    criterion --;
                    masterVariables.setCriterion (criterion);
                }
                else {
                    log.append (", aborting.\n");
                    bruteforce.setHasRun (false);
                    return;
                }
            }
        }
        // Output the best bruteforce result.
        ParameterSet bestBruteforceResult = bruteforce.getBestResult ();
        log.append ("The best result from bruteforce:\n");
        log.append (bestBruteforceResult.toString () + "\n\n");
    }

    /**
     *  Run the hillclimbing program.
     */
    protected void runHillclimbing () {
        log.append ("Starting hillclimbing...\n");
        hillclimb = new Hillclimb (
            masterVariables, nu, length, binning, bruteforce.getBestResult ()
        );
        hillclimb.run ();
        // Verify that hillclimbing ran correctly.
        if (! hillclimb.hasRun ()) {
            log.append ("  Error running the hillclimbing program!\n");
            return;
        }
        // Output the hillclimbing result.
        ParameterSet hClimbResult = hillclimb.getResult ();
        log.append ("The result from hillclimb:\n");
        log.append (hClimbResult.toString () + "\n\n");
    }

    /**
     *  Run the omega confidence interval program.
     */
    protected void runOmegaConfidenceInterval () {
        log.append ("Starting omega confidence interval...\n");
        omegaCI = new OmegaConfidenceInterval (
            masterVariables, nu, length, binning, hillclimb.getResult ()
        );
        omegaCI.run ();
        // Verify that omegaCI ran correctly.
        if (! omegaCI.hasRun ()) {
            log.append (
                "  Error running the omega confidence interval program!\n"
            );
            return;
        }
        // Output the omegaCI result.
        log.append ("The result from omegaCI:\n");
        log.append ("  " + omegaCI.toString () + "\n\n");
    }

    /**
     *  Run the sigma confidence interval program.
     */
    protected void runSigmaConfidenceInterval () {
        log.append ("Starting sigma confidence interval...\n");
        sigmaCI = new SigmaConfidenceInterval (
            masterVariables, nu, length, binning, hillclimb.getResult ()
        );
        sigmaCI.run ();
        // Verify that sigmaCI ran correctly.
        if (! sigmaCI.hasRun ()) {
            log.append (
                "  Error running the sigma confidence interval program!\n"
            );
            return;
        }
        // Output the sigmaCI result.
        log.append ("The result from sigmaCI:\n");
        log.append ("  " + sigmaCI.toString () + "\n\n");
    }

    /**
     *  Run the npop confidence interval program.
     */
    protected void runNpopConfidenceInterval () {
        log.append ("Starting npop confidence interval...\n");
        npopCI = new NpopConfidenceInterval (
            masterVariables, nu, length, binning, hillclimb.getResult ()
        );
        npopCI.run ();
        // Verify that npopCI ran correctly.
        if (! npopCI.hasRun ()) {
            log.append (
                "  Error running the npop confidence interval program!\n"
            );
            return;
        }
        // Output the npopCI result.
        log.append ("The result from npopCI:\n");
        log.append ("  " + npopCI.toString () + "\n\n");
    }

    /**
     *  Run the demarcation program.
     */
    protected void runDemarcation () {
        log.append ("Starting demarcation...\n");
        demarcation = new Demarcation (
            masterVariables, fasta, tree, binning, hillclimb.getResult ()
        );
        demarcation.run ();
        // Verify that demarcation ran correctly.
        if (! demarcation.hasRun ()) {
            log.append ("  Error running the demarcation program!\n");
            return;
        }
        // Output the demarcation result.
        log.append ("The result from demarcation:\n");
        log.append (demarcation.toString () + "\n\n");
    }

    protected MasterVariables masterVariables;
    protected File fastaFile;
    protected File newickFile;
    protected Logger log;
    protected Fasta fasta;
    protected Integer nu;
    protected Integer length;
    protected String outgroup;
    protected NewickTree tree;
    protected Binning binning;
    protected Bruteforce bruteforce;
    protected Hillclimb hillclimb;
    protected OmegaConfidenceInterval omegaCI;
    protected SigmaConfidenceInterval sigmaCI;
    protected NpopConfidenceInterval npopCI;
    protected Demarcation demarcation;
}
