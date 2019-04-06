/*
 *    Ecotype Simulation models the sequence diversity within a bacterial
 *    clade as the evolutionary result of net ecotype formation and periodic
 *    selection, yielding a certain number of ecotypes.
 *
 *    Copyright (C) 2013-2019  Jason M. Wood, Montana State University
 *
 *    This program is free software: you can redistribute it and/or modify
 *    it under the terms of the GNU General Public License as published by
 *    the Free Software Foundation, either version 2 of the License, or
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

import ecosim.tree.InvalidTreeException;
import ecosim.tree.SVGPainter;
import ecosim.tree.Tree;

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
     */
    public Simulation (Logger log, MasterVariables masterVariables) {
        this.log = log;
        this.masterVariables = masterVariables;
        // Set default demarcation method.
        demarcationPaintMethod = Demarcation.PAINT_METHOD_DEMARCATED;
        demarcationMethod = Demarcation.DEMARCATION_METHOD_MONOPHYLY;
        execs = new Execs (log, masterVariables);
        summary = new Summary ();
        // None of the programs are currently running.
        running = false;
    }

    /**
     *  Exit the simulation.
     */
    public void exit () {
        // Save any results to the output file if provided.
        File outputFile = masterVariables.getOutputFile ();
        if (outputFile != null) {
            saveProjectFile (outputFile);
        }
        masterVariables.exit ();
        System.exit (0);
    }

    /**
     *  Get the Summary object.
     *
     *  @return The Summary object.
     */
    public Summary getSummary () {
        return summary;
    }

    /**
     *  Get the Tree object.
     *
     *  @return The Tree object.
     */
    public Tree getTree () {
        return tree;
    }

    /**
     *  Get the paint method used for demarcation.
     *
     *  @return The demarcation paint method.
     */
    public int getDemarcationPaintMethod () {
        return demarcationPaintMethod;
    }

    /**
     *  Set the paint method used for demarcation.
     *
     *  @param demarcationPaintMethod The demarcation paint method.
     */
    public void setDemarcationPaintMethod (int demarcationPaintMethod) {
        this.demarcationPaintMethod = demarcationPaintMethod;
        demarcation.setPaintMethod (demarcationPaintMethod);
        summary.refreshObservers ();
    }

    /**
     *  Set the scale used for tree display.
     *
     *  @param scale The scale to use.
     */
    public void setScale (int scale) {
        if (tree != null && tree.isValid ()) {
            tree.setScale (scale);
        }
        if (demarcation != null && demarcation.isValid ()) {
            demarcation.setScale (scale);
        }
        summary.refreshObservers ();
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
        ParameterSet[] confidenceInterval = new ParameterSet[] {
            new ParameterSet (), new ParameterSet ()
        };
        log.append ("Opening: " + file.getPath () + "\n\n");
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
        // Update the summary data and output the values stored in the project
        // file to the log.
        if (nu > 0) {
            log.append ("Phylogeny results:\n");
            log.append (String.format (
                "  %,d environmental sequences.\n" +
                "  %,d sequence length.\n" +
                "  %s is the outgroup.\n\n",
                nu, length, outgroup
            ));
            summary.setNu (nu);
            summary.setLength (length);
            summary.setOutgroup (outgroup);
            summary.setTree (tree);
        }
        if (binning != null) {
            log.append ("Binning result:\n");
            ArrayList<BinLevel> bins = binning.getBins ();
            for (int i = 0; i < bins.size (); i ++) {
                log.append ("  " + bins.get (i).toString () + "\n");
            }
            log.append ("\n");
            summary.setBins (bins);
        }
        if (estimate != null) {
            log.append ("Parameter estimate:\n");
            log.append (estimate.toString () + "\n\n");
            summary.setEstimate (estimate);
        }
        if (hillclimb != null && hillclimb.hasRun ()) {
            log.append ("Hillclimbing result:\n");
            log.append (hillclimb.toString () + "\n\n");
            summary.setHillclimbing (hillclimb.getResult ());
        }
        if (npopCI != null && npopCI.hasRun ()) {
            log.append ("Npop confidence interval result:\n");
            log.append ("  " + npopCI.toString () + "\n\n");
            Long[] npop = npopCI.getResult ();
            confidenceInterval[0].setNpop (npop[0]);
            confidenceInterval[1].setNpop (npop[1]);
        }
        if (omegaCI != null && omegaCI.hasRun ()) {
            log.append ("Omega confidence interval result:\n");
            log.append ("  " + omegaCI.toString () + "\n\n");
            Double[] omega = omegaCI.getResult ();
            confidenceInterval[0].setOmega (omega[0]);
            confidenceInterval[1].setOmega (omega[1]);
        }
        if (sigmaCI != null && sigmaCI.hasRun ()) {
            log.append ("Sigma confidence interval result:\n");
            log.append ("  " + sigmaCI.toString () + "\n\n");
            Double[] sigma = sigmaCI.getResult ();
            confidenceInterval[0].setSigma (sigma[0]);
            confidenceInterval[1].setSigma (sigma[1]);
        }
        summary.setConfidenceInterval (confidenceInterval);
        if (demarcation != null && demarcation.hasRun ()) {
            log.append ("Demarcation result:\n");
            log.append (demarcation.toString () + "\n\n");
            summary.setDemarcation (demarcation);
        }
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
        log.append ("Saving to: " + file.getName () + "\n");
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
                "  sequence length: %,d.", length
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
            tree = new Tree (file);
            tree.makeBinary ();
            tree.reroot (tree.getDescendant (outgroup));
            // Output the tree in Newick and SVG formats if debug is enabled.
            if (masterVariables.getDebug ()) {
                String dir = masterVariables.getWorkingDirectory ();
                tree.toNewick (new File (dir + "outtree.nwk"));
                File svg = new File (dir + "outtree.svg");
                tree.paintTree (new SVGPainter (svg));
            }
            // Get the number of sequences loaded.
            nu = tree.size ();
            // Update the summary data.
            summary.setTree (tree);
            summary.setNu (nu);
            // Output the number of sequences loaded.
            log.appendln (String.format (
                "  %,d environmental sequences.", nu
            ));
            // Output the diversity sampled by the sequences.
            log.appendln (String.format (
                "  %.2f diversity sampled.", tree.getDiversity ()
            ));
        }
        catch (InvalidTreeException e) {
            System.out.println ("Error loading tree file.");
            e.printStackTrace ();
        }
    }

    /**
     *  Generate a tree using FastTree and the given sequence file.
     *
     *  @param file The fasta formated sequence file.
     */
    public File generateTree (File file) {
        log.appendln ("Generating a tree using FastTree...");
        // Store the tree in file called 'fasttree'.
        File newickFile = new File (
            masterVariables.getWorkingDirectory () + "fasttree.nwk"
        );
        // Generate a tree using FastTree.
        execs.runFastTree (file, newickFile);
        return newickFile;
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
     *  Check whether the tree has been loaded.
     *
     *  @return True if the tree has been loaded, otherwise False.
     */
    public boolean treeLoaded () {
        return tree != null && tree.isValid ();
    }

    /**
     *  Check whether hillclimbing has been run.
     *
     *  @return True if the hillclimbing has been run, otherwise False.
     */
    public boolean hillclimbHasRun () {
        return hillclimb != null && hillclimb.hasRun ();
    }

    /**
     *  Run the binning program.
     */
    public void runBinning () {
        // Start running binning.
        running = true;
        log.appendln ("Running binning...");
        binning = new Binning (tree);
        binning.run ();
        ArrayList<BinLevel> bins = binning.getBins ();
        // Update the summary data.
        summary.setBins (bins);
        // Output the results from binning.
        log.appendln ("The result from binning:");
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
        estimate = new ParameterEstimate (nu, length, binning);
        estimate.run ();
        // Update the summary data.
        summary.setEstimate (estimate);
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
            "Starting with precision: " +
            masterVariables.getCriterionLabel (crit)
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
            if (crit == 1) {
                log.appendln ("  Error, the precision can not be reduced.");
                break;
            }
            crit --;
            String critLabel = masterVariables.getCriterionLabel (crit);
            log.appendln ("Reducing the precision to " + critLabel + ".");
            masterVariables.setCriterion (crit);
        }
        // Update the summary data.
        summary.setHillclimbing (hillclimb.getResult ());
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
            masterVariables, execs, nu, length, binning,
            hillclimb.getResult ()
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
            masterVariables, execs, nu, length, binning,
            hillclimb.getResult ()
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
            masterVariables, execs, nu, length, binning,
            hillclimb.getResult ()
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
        // Initialize the confidence interval.
        confidenceInterval = new ParameterSet[] {
            new ParameterSet (), new ParameterSet ()
        };
        // Run the npop confidence interval.
        runNpopConfidenceInterval ();
        Long[] npop = npopCI.getResult ();
        confidenceInterval[0].setNpop (npop[0]);
        confidenceInterval[1].setNpop (npop[1]);
        summary.setConfidenceInterval (confidenceInterval);
        // Run the omega confidence interval.
        runOmegaConfidenceInterval ();
        Double[] omega = omegaCI.getResult ();
        confidenceInterval[0].setOmega (omega[0]);
        confidenceInterval[1].setOmega (omega[1]);
        summary.setConfidenceInterval (confidenceInterval);
        // Run the sigma confidence interval.
        runSigmaConfidenceInterval ();
        Double[] sigma = sigmaCI.getResult ();
        confidenceInterval[0].setSigma (sigma[0]);
        confidenceInterval[1].setSigma (sigma[1]);
        summary.setConfidenceInterval (confidenceInterval);
    }

    /**
     *  Run the demarcation program.
     */
    public void runDemarcation () {
        // Start running demarcation.
        running = true;
        log.appendln ("Running demarcation...");
        try {
            demarcation = new Demarcation (
                masterVariables, execs, nu, length, outgroup, tree,
                hillclimb.getResult (), demarcationMethod
            );
            demarcation.setPaintMethod (demarcationPaintMethod);
            demarcation.run ();
            // Verify that demarcation ran correctly.
            if (! demarcation.hasRun ()) {
                log.appendln ("  Error running the demarcation program!");
                return;
            }
        }
        catch (InvalidTreeException e) {
            log.appendln ("  Error running the demarcation program!");
            return;
        }
        // Output the demarcation in SVG format if debugging is turned on.
        if (masterVariables.getDebug ()) {
            String dir = masterVariables.getWorkingDirectory ();
            File svg = new File (dir + "demarcation.svg");
            demarcation.paintTree (new SVGPainter (svg));
        }
        // Update the summary data.
        summary.setDemarcation (demarcation);
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
    protected Summary summary;
    protected Fasta fasta;
    protected Integer nu;
    protected Integer length;
    protected String outgroup;
    protected Tree tree;
    protected Binning binning;
    protected ParameterEstimate estimate;
    protected Hillclimb hillclimb;
    protected ParameterSet[] confidenceInterval;
    protected NpopConfidenceInterval npopCI;
    protected OmegaConfidenceInterval omegaCI;
    protected SigmaConfidenceInterval sigmaCI;
    protected Demarcation demarcation;
    protected Integer demarcationPaintMethod;
    protected Integer demarcationMethod;
    protected boolean running;

}
