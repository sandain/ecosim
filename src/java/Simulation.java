/*
 *    Ecotype Simulation models the sequence diversity within a bacterial
 *    clade as the evolutionary result of net ecotype formation and periodic
 *    selection, yielding a certain number of ecotypes.
 *
 *    Copyright (C) 2013-2015  Jason M. Wood, Montana State University
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
     *  @param log The Logger object.
     *  @param masterVariables The MasterVariables object.
     *  @param fastaFile The fasta formated sequence file.
     *  @param newickFile The newick formated tree file.
     */
    public Simulation (Logger log, MasterVariables masterVariables,
        File fastaFile, File newickFile) {
        this.log = log;
        this.masterVariables = masterVariables;
        this.fastaFile = fastaFile;
        this.newickFile = newickFile;
        // Set default demarcation method and precision.
        demarcationMethod = Demarcation.DEMARCATION_METHOD_MONOPHYLY;
        demarcationPrecision = Demarcation.DEMARCATION_PRECISION_FINE_SCALE;
        execs = new Execs (log, masterVariables);
        if (fastaFile != null && fastaFile.exists ()) {
            loadSequenceFile (fastaFile);
        }
        if (newickFile != null && newickFile.exists ()) {
            loadTreeFile (newickFile);
        }
        // None of the programs are currently running.
        running = false;
    }

    /**
     *  Exit the simulation.
     */
    public void exit () {
        masterVariables.exit ();
        System.exit (0);
    }

    /**
     *  Get the method used for demarcation.
     *
     *  @return The demarcation method.
     */
    public int getDemarcationMethod () {
        return demarcationMethod;
    }

    /**
     *  Get the precision used for demarcation.
     *
     *  @return The demarcation precision.
     */
    public int getDemarcationPrecision () {
        return demarcationPrecision;
    }

    /**
     *  Set the method used for demarcation.
     *
     *  @param demarcationMethod The demarcation method.
     */
    public void setDemarcationMethod (int demarcationMethod) {
        this.demarcationMethod = demarcationMethod;
    }

    /**
     *  Set the precision used for demarcation.
     *
     *  @param demarcationPrecision The demarcation precision.
     */
    public void setDemarcationPrecision (int demarcationPrecision) {
        this.demarcationPrecision = demarcationPrecision;
    }

    /**
     *  Load the project from a XML formated file.
     *
     *  @param file The file to load the project from.
     */
    public void loadProjectFile (File file) {
        ProjectFileIO projectFileIO = new ProjectFileIO (
            masterVariables, execs
        );
        // Load the project file.
        projectFileIO.load (file);
        // Grab the loaded variables.
        nu = projectFileIO.getNu ();
        length = projectFileIO.getLength ();
        outgroup = projectFileIO.getOutgroup ();
        tree = projectFileIO.getTree ();
        binning = projectFileIO.getBinning ();
        estimate = projectFileIO.getEstimate ();
        hillclimb = projectFileIO.getHillclimb ();
        npopCI = projectFileIO.getNpopCI ();
        omegaCI = projectFileIO.getOmegaCI ();
        sigmaCI = projectFileIO.getSigmaCI ();
        demarcation = projectFileIO.getDemarcation ();
    }

    /**
     *  Save the project to a XML formated file.
     *
     *  @param file The file to save the project to.
     */
    public void saveProjectFile (File file) {
        ProjectFileIO projectFileIO = new ProjectFileIO (
            masterVariables, execs, nu, length, outgroup, tree, binning,
            estimate, hillclimb, npopCI, omegaCI, sigmaCI, demarcation
        );
        projectFileIO.save (file);
    }

    /**
     *  Load the fasta formated sequence file.
     *
     *  @param file The fasta formated sequence file to load.
     */
    public void loadSequenceFile (File file) {
        // Verify that the sequence file exists.
        if (file == null || ! file.exists ()) {
            log.appendln ("Error loading the Fasta file!");
            return;
        }
        log.appendln ("Opening sequence file...");
        try {
            fasta = new Fasta (file);
            Sequence outgroupSequence = fasta.getOutgroup ();
            length = outgroupSequence.length ();
            outgroup = outgroupSequence.getIdentifier ();
            // Update the summary data.
            summary.setLength (length);
            summary.setOutgroup (outgroup);
            // Output the sequence data.
            log.appendln (String.format (
                "  sequence length: %d.", length
            ));
            log.appendln (String.format (
                "  %s is the outgroup.", outgroup
            ));
        }
        catch (InvalidFastaException e) {
            System.out.println ("Error loading sequence file.");
        }
    }

    /**
     *  Load the newick formated tree file.
     *
     *  @param file The newick formated tree file to load.
     */
    public void loadTreeFile (File file) {
        // Verify that the tree file exists.
        if (file == null || ! file.exists ()) {
            log.appendln ("Error loading the Newick tree!");
            return;
        }
        log.appendln ("Opening tree file...");
        try {
            tree = new NewickTree (file);
            tree.reroot (outgroup);
            // Store the tree in file called 'outtree'.
            tree.save (new File (
                masterVariables.getWorkingDirectory () + "outtree"
            ));
            // Get the number of sequences loaded.
            nu = tree.size ();
            // Update the summary data.
            summary.setNu (nu);
            // Output the number of sequences loaded.
            log.appendln (String.format (
                "  %d environmental sequences.", nu
            ));
        }
        catch (InvalidNewickException e) {
            System.out.println ("Error loading tree file.");
        }
    }

    /**
     *  Generate a tree using FastTree and the given sequence file.
     *
     *  @param file The fasta formated sequence file.
     */
    public void generateTree (File file) {
        log.appendln ("Generating a tree using FastTree...");
        // Store the tree in file called 'fasttree'.
        File newickFile = new File (
            masterVariables.getWorkingDirectory () + "fasttree"
        );
        // Generate a tree using FastTree.
        execs.runFastTree (file, newickFile);
        // Load the tree.
        loadTreeFile (newickFile);
    }

    /**
     *  Check whether the simulation is running.
     *
     *  @return True if the simulation is running, otherwise False.
     */
    public boolean isRunning () {
        return running;
    }

    /**
     *  Run the binning program.
     */
    public void runBinning () {
        // Start running binning.
        running = true;
        log.appendln ("Running binning...");
        binning = new Binning (tree);
        // Output the results from binning.
        log.appendln ("The result from binning:");
        ArrayList<BinLevel> bins = binning.getBins ();
        for (int i = 0; i < bins.size (); i ++) {
            log.appendln ("  " + bins.get (i).toString ());
        }
        log.appendln ();
        // Done running binning.
        running = false;
    }

    /**
     *  Run the parameter estimate program.
     */
    public void runParameterEstimate () {
        // Start running the parameter estimate.
        running = true;
        log.appendln ("Running parameter estimate...");
        estimate = new ParameterEstimate (
            masterVariables, execs, nu, length, binning
        );
        estimate.run ();
        // Output the parameter estimate.
        log.appendln ("The estimated parameters:");
        log.appendln (estimate.toString ());
        log.appendln ();
        // Done running the parameter estimate.
        running = false;
    }

    /**
     *  Run the hillclimbing program.
     */
    public void runHillclimbing () {
        // Start running hillclimbing.
        running = true;
        Integer crit = masterVariables.getCriterion ();
        Double likelihood = 0.0d;
        log.appendln ("Running hillclimbing...");
        log.appendln (
            "Starting with precision: " + masterVariables.getCriterionLabel (crit)
        );
        hillclimb = new Hillclimb (
            masterVariables, execs, nu, length, binning, estimate.getResult ()
        );
        while (likelihood < masterVariables.EPSILON)  {
            // Run hillclimbing using the current criterion.
            hillclimb.run ();
            // Verify that hillclimbing ran correctly.
            if (! hillclimb.hasRun ()) {
                log.appendln ("  Error running the hillclimbing program!");
                return;
            }
            likelihood = hillclimb.getResult ().getLikelihood ();
            // Break out of the loop if likelihood > zero.
            if (likelihood > masterVariables.EPSILON) break;
            log.appendln ("Hillclimbing result has zero likelihood.");
            // Check to see if the precision can be reduced.
            if (crit == 0) {
                log.appendln ("  Error, the precision can not be reduced.");
                break;
            }
            crit --;
            String critLabel = masterVariables.getCriterionLabel (crit);
            log.appendln ("Reducing the precision to " + critLabel + ".");
            masterVariables.setCriterion (crit);
        }
        // Output the hillclimbing result.
        log.appendln ("The result from hillclimb:");
        log.appendln (hillclimb.toString ());
        log.appendln ();
        // Done running hillclimbing.
        running = false;
    }

    /**
     *  Run the npop confidence interval program.
     */
    public void runNpopConfidenceInterval () {
        // Start running the confidence interval.
        running = true;
        log.appendln ("Running npop confidence interval...");
        npopCI = new NpopConfidenceInterval (
            masterVariables, execs, nu, length, binning, hillclimb.getResult ()
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
        // Done running the confidence interval.
        running = false;
    }

    /**
     *  Run the omega confidence interval program.
     */
    public void runOmegaConfidenceInterval () {
        // Start running the confidence interval.
        running = true;
        log.appendln ("Running omega confidence interval...");
        omegaCI = new OmegaConfidenceInterval (
            masterVariables, execs, nu, length, binning, hillclimb.getResult ()
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
        // Done running the confidence interval.
        running = false;
    }

    /**
     *  Run the sigma confidence interval program.
     */
    public void runSigmaConfidenceInterval () {
        // Start running the confidence interval.
        running = true;
        log.appendln ("Running sigma confidence interval...");
        sigmaCI = new SigmaConfidenceInterval (
            masterVariables, execs, nu, length, binning, hillclimb.getResult ()
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
        // Done running the confidence interval.
        running = false;
    }

    /**
     *  Run the omega, sigma, and npop confidence interval programs.
     */
    public void runConfidenceIntervals () {
        runNpopConfidenceInterval ();
        runOmegaConfidenceInterval ();
        runSigmaConfidenceInterval ();
    }

    /**
     *  Run the demarcation program.
     */
    public void runDemarcation () {
        // Start running demarcation.
        running = true;
        log.appendln ("Running demarcation...");
        demarcation = new Demarcation (
            masterVariables, execs, nu, length, outgroup, tree,
            hillclimb.getResult (), demarcationMethod, demarcationPrecision
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
        // Done running demarcation.
        running = false;
    }

    protected Logger log;
    protected MasterVariables masterVariables;
    protected Execs execs;
    protected File fastaFile;
    protected File newickFile;
    protected Fasta fasta;
    protected Integer nu;
    protected Integer length;
    protected String outgroup;
    protected NewickTree tree;
    protected Binning binning;
    protected ParameterEstimate estimate;
    protected Hillclimb hillclimb;
    protected NpopConfidenceInterval npopCI;
    protected OmegaConfidenceInterval omegaCI;
    protected SigmaConfidenceInterval sigmaCI;
    protected Demarcation demarcation;
    protected Integer demarcationMethod;
    protected Integer demarcationPrecision;
    protected boolean running;

}
