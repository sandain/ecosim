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
            log.appendln ("Error loading the Fasta file!");
            return;
        }
        log.appendln ("Opening sequence file: " + fastaFile.getName ());
        try {
            fasta = new Fasta (fastaFile);
            nu = fasta.size ();
            length = fasta.length ();
            outgroup = fasta.getIdentifier (0);
            // Output the number of sequences loaded.
            log.appendln (String.format ("  %d environmental sequences.", nu));
            log.appendln (String.format ("  %s is the outgroup.", outgroup));
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
            log.appendln ("Error loading the Newick tree!");
            return;
        }
        log.appendln ("Opening tree file: " + newickFile.getName ());
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
        log.appendln ("Generating a tree using FastTree...");
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
        log.appendln ("Starting binning...");
        binning = new Binning (masterVariables, tree);
        // Output the results from binning.
        log.appendln ("The result from binning:");
        ArrayList<BinLevel> bins = binning.getBins ();
        for (int i = 0; i < bins.size (); i ++) {
            log.appendln ("  " + bins.get (i).toString ());
        }
        log.appendln ();
    }

    protected void runParameterEstimate () {
        log.appendln ("Estimating parameters...");
        estimate = new ParameterEstimate (length, binning);
        // Output the parameter estimate.
        log.appendln ("The estimated parameters:");
        log.appendln (estimate.toString ());
        log.appendln ();
    }

    /**
     *  Run the bruteforce program.
     */
    protected void runBruteforce () {
        log.appendln ("Starting bruteforce search...");
        bruteforce = new Bruteforce (masterVariables, nu, length, binning);
        while (bruteforce.getNumResults () < masterVariables.NUM_SUCCESSES) {
            bruteforce.run ();
            // Verify that bruteforce ran correctly.
            if (! bruteforce.hasRun ()) {
                log.appendln ("Error running the bruteforce search program!");
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
                    log.appendln (", lowering the value.");
                    criterion --;
                    masterVariables.setCriterion (criterion);
                }
                else {
                    log.appendln (", aborting.");
                    bruteforce.setHasRun (false);
                    return;
                }
            }
        }
        // Output the best bruteforce result.
        ParameterSet bestBruteforceResult = bruteforce.getBestResult ();
        log.appendln ("The best result from bruteforce:");
        log.appendln (bestBruteforceResult.toString ());
        log.appendln ();
    }

    /**
     *  Run the hillclimbing program.
     */
    protected void runHillclimbing () {
        log.appendln ("Starting hillclimbing...");
        hillclimb = new Hillclimb (
            masterVariables, nu, length, binning, bruteforce.getBestResult ()
        );
        hillclimb.run ();
        // Verify that hillclimbing ran correctly.
        if (! hillclimb.hasRun ()) {
            log.appendln ("  Error running the hillclimbing program!");
            return;
        }
        // Output the hillclimbing result.
        log.appendln ("The result from hillclimb:");
        log.appendln (hillclimb.toString ());
        log.appendln ();
    }

    /**
     *  Run the omega confidence interval program.
     */
    protected void runOmegaConfidenceInterval () {
        log.appendln ("Starting omega confidence interval...");
        omegaCI = new OmegaConfidenceInterval (
            masterVariables, nu, length, binning, hillclimb.getResult ()
        );
        omegaCI.run ();
        // Verify that omegaCI ran correctly.
        if (! omegaCI.hasRun ()) {
            log.appendln (
                "  Error running the omega confidence interval program!"
            );
            return;
        }
        // Output the omegaCI result.
        log.appendln ("The result from omegaCI:");
        log.appendln ("  " + omegaCI.toString ());
        log.appendln ();
    }

    /**
     *  Run the sigma confidence interval program.
     */
    protected void runSigmaConfidenceInterval () {
        log.appendln ("Starting sigma confidence interval...");
        sigmaCI = new SigmaConfidenceInterval (
            masterVariables, nu, length, binning, hillclimb.getResult ()
        );
        sigmaCI.run ();
        // Verify that sigmaCI ran correctly.
        if (! sigmaCI.hasRun ()) {
            log.appendln (
                "  Error running the sigma confidence interval program!"
            );
            return;
        }
        // Output the sigmaCI result.
        log.appendln ("The result from sigmaCI:");
        log.appendln ("  " + sigmaCI.toString ());
        log.appendln ();
    }

    /**
     *  Run the npop confidence interval program.
     */
    protected void runNpopConfidenceInterval () {
        log.appendln ("Starting npop confidence interval...");
        npopCI = new NpopConfidenceInterval (
            masterVariables, nu, length, binning, hillclimb.getResult ()
        );
        npopCI.run ();
        // Verify that npopCI ran correctly.
        if (! npopCI.hasRun ()) {
            log.appendln (
                "  Error running the npop confidence interval program!"
            );
            return;
        }
        // Output the npopCI result.
        log.appendln ("The result from npopCI:");
        log.appendln ("  " + npopCI.toString ());
        log.appendln ();
    }

    /**
     *  Run the demarcation program.
     */
    protected void runDemarcation () {
        log.appendln ("Starting demarcation...");
        demarcation = new Demarcation (
            masterVariables, nu, length, outgroup, tree,
            hillclimb.getResult ()
        );
        demarcation.run ();
        // Verify that demarcation ran correctly.
        if (! demarcation.hasRun ()) {
            log.appendln ("  Error running the demarcation program!");
            return;
        }
        // Output the demarcation result.
        log.appendln ("The result from demarcation:");
        log.appendln (demarcation.toString ());
        log.appendln ();
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
    protected ParameterEstimate estimate;
    protected Bruteforce bruteforce;
    protected Hillclimb hillclimb;
    protected OmegaConfidenceInterval omegaCI;
    protected SigmaConfidenceInterval sigmaCI;
    protected NpopConfidenceInterval npopCI;
    protected Demarcation demarcation;
}
