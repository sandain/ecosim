/*
 *    Ecotype Simulation models the sequence diversity within a bacterial clade
 *    as the evolutionary result of net ecotype formation and periodic
 *    selection, yielding a certain number of ecotypes.
 *
 *    Copyright (C) 2013  Jason M. Wood, Montana State University
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
    public Simulation(MasterVariables masterVariables, File fastaFile,
        File newickFile) {
        this.masterVariables = masterVariables;
        this.fastaFile = fastaFile;
        this.newickFile = newickFile;
        log = masterVariables.getLog();
        phylogeny = new Phylogeny(masterVariables);
        if (fastaFile != null && fastaFile.exists()) {
            loadSequenceFile();
        }
        if (newickFile != null && newickFile.exists()) {
            loadTreeFile();
        }
    }

    /**
     *  Save the project to an XML formated file.
     *
     *  @param file The file to save the project to.
     */
    public void saveProjectFile(File file) {
        ProjectFileIO projectFileIO = new ProjectFileIO(
            masterVariables, phylogeny, binning, bruteforce,
            hillclimb, omegaCI, sigmaCI, npopCI, demarcation
        );
        projectFileIO.save(file);
    }

    /**
     *  Load the fasta formated sequence file.
     */
    protected void loadSequenceFile() {
        // Verify that the fasta file exists.
        if (fastaFile == null || ! fastaFile.exists()) {
            log.append("Error loading the Fasta file!\n");
            return;
        }
        log.append(
            "Opening sequence file: " + fastaFile.getName() + "\n"
        );
        phylogeny.loadSequenceFile(fastaFile);
    }

    /**
     *  Load the newick formated tree file.
     */
    protected void loadTreeFile() {
        // Verify that the tree file exists.
        if (newickFile == null || ! newickFile.exists()) {
            log.append("Error loading the Newick tree!\n");
            return;
        }
        log.append(
            "Opening tree file: " + newickFile.getName() + "\n"
        );
        phylogeny.loadTreeFile(newickFile);
    }

    /**
     *  Generate a tree using Phylip's Neibor-Joining or Parsimony methods.
     *
     *  @param method The method to use to generate the tree.
     */
    protected void generateTree(String method) {
        log.append(String.format(
            "Generating a %s tree using Phylip...\n",
            method
        ));
        newickFile = phylogeny.generateTree(method);
    }

    /**
     *  Run the phylogeny program.
     */
    protected void runPhylogeny() {
        log.append(
            "Starting the phylogeny program...\n"
        );
        phylogeny.run();
        // Verify that the phylogeny programs ran correctly.
        if (! phylogeny.hasRun()) {
            log.append("  Error running the phylogeny program!\n");
            return;
        }
        log.append(
            "The results from the phylogeny program:\n"
        );
        // Output the number of sequences loaded.
        log.append(String.format(
            "  %d sequences loaded.\n" +
            "  %s is the outgroup.\n\n",
            phylogeny.getNu(), phylogeny.getOutgroupIdentifier()
        ));
    }

    /**
     *  Run the binning program.
     */
    protected void runBinning() {
        log.append("Starting binning...\n");
        binning = new Binning(masterVariables, phylogeny);
        binning.run();
        // Verify that binning program ran correctly.
        if (! binning.hasRun()) {
            log.append("  Error running the binning program!\n");
            return;
        }
        // Output the results from binning.
        log.append("The result from binning:\n");
        ArrayList<BinLevel> bins = binning.getBinLevels();
        for (int i = 0; i < bins.size(); i ++) {
            log.append("  " + bins.get(i).toString() + "\n");
        }
        log.append("\n");
    }

    /**
     *  Run the bruteforce program.
     */
    protected void runBruteforce() {
        log.append("Starting bruteforce search...\n");
        bruteforce = new Bruteforce(masterVariables, phylogeny, binning);
        while (bruteforce.getNumResults() < masterVariables.NUM_SUCCESSES) {
            bruteforce.run();
            // Verify that bruteforce ran correctly.
            if (! bruteforce.hasRun()) {
                log.append("Error running the bruteforce search program!\n");
                return;
            }
            // Verify that there are enough bruteforce results before continuing.
            if (bruteforce.getNumResults() < masterVariables.NUM_SUCCESSES) {
                int criterion = masterVariables.getCriterion();
                log.append(
                    "  Not enough results at the current criterion (" +
                    masterVariables.getCriterionLabel(criterion) + ")"
                );
                if (criterion > 2) {
                    log.append(", lowering the value.\n");
                    criterion --;
                    masterVariables.setCriterion(criterion);
                }
                else {
                    log.append(", aborting.\n");
                    bruteforce.setHasRun(false);
                    return;
                }
            }
        }
        // Output the best bruteforce result.
        ParameterSet bestBruteforceResult = bruteforce.getBestResult();
        log.append("The best result from bruteforce:\n");
        log.append(bestBruteforceResult.toString() + "\n\n");
    }

    /**
     *  Run the hillclimbing program.
     */
    protected void runHillclimbing() {
        log.append("Starting hillclimbing...\n");
        hillclimb = new Hillclimb(
            masterVariables, phylogeny, binning, bruteforce.getBestResult()
        );
        hillclimb.run();
        // Verify that hillclimbing ran correctly.
        if (! hillclimb.hasRun()) {
            log.append("  Error running the hillclimbing program!\n");
            return;
        }
        // Output the hillclimbing result.
        ParameterSet hClimbResult = hillclimb.getResult();
        log.append("The result from hillclimb:\n");
        log.append(hClimbResult.toString() + "\n\n");
    }

    /**
     *  Run the omega confidence interval program.
     */
    protected void runOmegaConfidenceInterval() {
        log.append("Starting omega confidence interval...\n");
        omegaCI = new OmegaConfidenceInterval(
            masterVariables, phylogeny, binning, hillclimb
        );
        omegaCI.run();
        // Verify that omegaCI ran correctly.
        if (! omegaCI.hasRun()) {
            log.append(
                "  Error running the omega confidence interval program!\n"
            );
            return;
        }
        // Output the omegaCI result.
        log.append("The result from omegaCI:\n");
        log.append("  " + omegaCI.toString() + "\n\n");
    }

    /**
     *  Run the sigma confidence interval program.
     */
    protected void runSigmaConfidenceInterval() {
        log.append("Starting sigma confidence interval...\n");
        sigmaCI = new SigmaConfidenceInterval(
            masterVariables, phylogeny, binning, hillclimb
        );
        sigmaCI.run();
        // Verify that sigmaCI ran correctly.
        if (! sigmaCI.hasRun()) {
            log.append(
                "  Error running the sigma confidence interval program!\n"
            );
            return;
        }
        // Output the sigmaCI result.
        log.append("The result from sigmaCI:\n");
        log.append("  " + sigmaCI.toString() + "\n\n");
    }

    /**
     *  Run the npop confidence interval program.
     */
    protected void runNpopConfidenceInterval() {
        log.append("Starting npop confidence interval...\n");
        npopCI = new NpopConfidenceInterval(
            masterVariables, phylogeny, binning, hillclimb
        );
        npopCI.run();
        // Verify that npopCI ran correctly.
        if (! npopCI.hasRun()) {
            log.append(
                "  Error running the npop confidence interval program!\n"
            );
            return;
        }
        // Output the npopCI result.
        log.append("The result from npopCI:\n");
        log.append("  " + npopCI.toString() + "\n\n");
    }

    /**
     *  Run the demarcation program.
     */
    protected void runDemarcation() {
        log.append("Starting demarcation...\n");
        demarcation = new Demarcation(
            masterVariables, phylogeny, binning, hillclimb
        );
        demarcation.run();
        // Verify that demarcation ran correctly.
        if (! demarcation.hasRun()) {
            log.append("  Error running the demarcation program!\n");
            return;
        }
        // Output the demarcation result.
        log.append("The result from demarcation:\n");
        log.append(demarcation.toString() + "\n\n");
    }

    protected MasterVariables masterVariables;
    protected File fastaFile;
    protected File newickFile;
    protected Logger log;
    protected Phylogeny phylogeny;
    protected Binning binning;
    protected Bruteforce bruteforce;
    protected Hillclimb hillclimb;
    protected OmegaConfidenceInterval omegaCI;
    protected SigmaConfidenceInterval sigmaCI;
    protected NpopConfidenceInterval npopCI;
    protected Demarcation demarcation;
}
