/*
 *    Ecotype Simulation models the sequence diversity within a bacterial clade as
 *    the evolutionary result of net ecotype formation, periodic selection,
 *    and drift, yielding a certain number of ecotypes.
 * 
 *    Copyright (C) 2009  Fred Cohan, Wesleyan University
 *                        Carlo Francisco, Wesleyan University
 *                        Danny Krizanc, Wesleyan University
 *                        Andrew Warner, Wesleyan University
 *                        Jason Wood, Montana State University
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
import javax.swing.JTextArea;

/**
 *  Runs the binning process and also runs full likelihoods.
 *
 *  @author Andrew Warner
 */
public class BinningAndFred implements Runnable {

    /**
     *  Constructor for BinnandAndFred.
     * 
     *  @param inFasta The fasta file to work with.
     *  @param masterVariables The MasterVariables.
     *  @param userChoice The user's choice.
     */
    public BinningAndFred(Fasta inFasta, MasterVariables masterVariables, int userChoice) {
        this.inFasta = inFasta;
        this.masterVariables = masterVariables;
        this.userChoice = userChoice;
        log = masterVariables.getLog();
        narr = masterVariables.getNarrator();
        execs = masterVariables.getExecs();
        workingDirectory = masterVariables.getWorkingDirectory();
        binningFasta = new BinningFasta(masterVariables);
        fredValsReader = new FredValsReader(masterVariables);
        fredMethodInput = new File(workingDirectory + "bruteforceIn.dat");
        fredOutput = new File(workingDirectory + "bruteforceOut.dat");
    }

    /**
     *  Constructor for binningandfred from a SAVED file --->>> NOTE! run()
     *  must NOT be called after running this constructor as inFasta has not
     *  been defined.
     *
     *  @param masterVariables The MasterVariables.
     *  @param sequenceVals The sequence values.
     *  @param bins The bins.
     */
    public BinningAndFred(MasterVariables masterVariables, int[] sequenceVals, ArrayList<String> bins) {
        this.masterVariables = masterVariables;
        this.sequenceVals = sequenceVals;
        this.bins = bins;
        log = masterVariables.getLog();
        narr = masterVariables.getNarrator();
        execs = masterVariables.getExecs();
        workingDirectory = masterVariables.getWorkingDirectory();
        fredMethodInput = new File(workingDirectory + "bruteforceIn.dat");
        fredOutput = new File(workingDirectory + "bruteforceOut.dat");
    }

    /**
     *  Run binning.
     */
    public void run() {
        if (userChoice == -1) {
            masterVariables.setSortPercentage(5);
        }
        else {
            masterVariables.setSortPercentage(userChoice);
        }
        // Run fred and get the best value to run the full likelihood.
        log.append("Running binning program...\n");
        narr.println("Running binning program...");
        ArrayList<FredOutVal> results = getFredOutput(inFasta, fredMethodInput);
        bestFred = findBestHclimb(results);
    }

    /**
     *  Runs only the binning program, useful in demarcations.  In the case
     *  of the original run requiring a value of omega, sigma, npop, etc.
     *  we use the run() function, we runs fred method and hillclimbing,
     *  but this method will give us only the binning data.
     */
    public void runBinningOnly() {
        File outFasta = new File(workingDirectory + "sequencesfasta.txt");
        File numbers = new File(workingDirectory + "numbers.dat");
        // Read in the fasta file, transfer it to sequencesfasta.txt, and write numbers.dat.
        sequenceVals = binningFasta.readFile(inFasta, outFasta, numbers);
        // Run the binning program.
        execs.runBinning(); 
        // Read the new length of the sequence after gaps have been removed.
        File rg = new File(workingDirectory + "fasta.dat");
        sequenceVals[1] = binningFasta.readRGlength(rg);
        File outBin = new File(workingDirectory + "output.dat");
        bins = binningFasta.readBins(outBin, narr);
    }

    /**
     *  Runs the full likelihood, with only one data point involved.
     *
     *  pre: out is the acinas input file, value is the best data point from the initial fred method run, sequenceVals.
     *  post: FredOutVal is the full likelihood value.
     *
     *  @param value The binning data.
     *  @return FredOutVal containing the best value for hillclimbing with it's full likelihood.
     */
    public FredOutVal fullLike(FredOutVal value) {
        // Declare default values
        double[] omegaRange = {value.getOmega(), 10000};
        double[] sigmaRange = {value.getSigma(), 10000};
        if (sigmaRange[0] > 1.0e25)
            sigmaRange[1] = 1.0e31;
        int[] npopRange = {value.getNpop(), 10000};
        int[] xnumics = {0, 0, 0, 0};
        int jwhichxavg = masterVariables.getSortPercentage() + 1;
        int nrep = NREP_FULL_LIKE;
        // Write acinas.dat.
        InputWriter.writeFile(fredMethodInput, bins, omegaRange, sigmaRange,
            npopRange, masterVariables.getDriftRange(), xnumics,
            sequenceVals[0], nrep, sequenceVals[1], jwhichxavg);
        narr.println();
        narr.println("The input for full likelihood: ");
        narr.writeInput(fredMethodInput);
        // Run fred method.
        log.append("Running Brute Force Search...\n");
        execs.runBruteforce();
        // Read in results.
        ArrayList<FredOutVal> results = fredValsReader.readFredOutVals(fredOutput);
        // Return the first and only result.
        return results.get(0);
    }

    /**
     *  Get the sequence values.
     *
     *  @return int array containing the sequence values.
     */
    public int[] getSeqVals() {
        return sequenceVals.clone();
    }

    /**
     *  Get the bins.
     *
     *  @return ArrayList<String> containing the bins.
     */
    public ArrayList<String> getBins() {
        return bins;
    }

    /**
     *  Returns the best FredOutVal.
     * 
     *  @return The best FredOutVal.
     */
    public FredOutVal getValue() {
        return bestFred;
    }

    /**
     *  Set the bins. Be careful with this one!
     *
     *  @param bins The bins.
     */
    public void setBins(ArrayList<String> bins) {
        this.bins = bins;
    }

    /**
     *  Set the sequence values.
     *
     *  @param newVals The sequence values.
     */
    public void setSeqVals(int[] newVals) {
        this.sequenceVals = newVals;
    }

    /**
     *  Choose the best fred value according to the sorting algorithm.
     *
     *  @param data The list of unsorted output fred values from a run of fred method.
     *  @return FredOutVal containing the best Fred value.
     */
    private FredOutVal getBestFred(ArrayList<FredOutVal> data) {
        // Sort the values and pick the best one, storing the sorted values in the debug file.
        File debug = new File(workingDirectory + "fredDebug.dat");
        FredOutVal bestFred;
        bestFred = binningFasta.chooseBest(data, debug);
        return bestFred;
    }

    /**
     *  This method takes in a properly formatted fasta file, runs binning on
     *  it, runs Fred, sorts it, and outputs the unsorted output from the run
     *  of fred method.
     *
     *  pre: inFasta is a properly formatted fasta file.
     *  post: The output results from a run of fred method are returned.
     *
     *  @param inFasta A properly formatted fasta file.
     *  @param fredMethodInput The output from fred method.
     *  @return ArrayList<FredOutVal> containing the output results from a run of fred method.
     */
      private ArrayList<FredOutVal> getFredOutput(Fasta inFasta, File fredMethodInput) {
        File outFasta = new File(workingDirectory + "sequencesfasta.txt");
        File numbers = new File(workingDirectory + "numbers.dat");
        // Output for the fasta reader.
        // Read in the fasta file, transfer it to sequencesfasta.txt, and write numbers.dat.
        sequenceVals = binningFasta.readFile(inFasta, outFasta, numbers);
        // Run the binning program.
        execs.runBinning();
        // Read the new length of the sequence after gaps have been removed.
        File rg = new File(workingDirectory + "fasta.dat");
        sequenceVals[1] = binningFasta.readRGlength(rg);
        narr.println("The number of sequences is: "+sequenceVals[0]);
        narr.println("The length of each sequences after removegaps is: " + sequenceVals[1]);
        File outBin = new File(workingDirectory + "output.dat");
        bins = binningFasta.readBins(outBin, narr);
        int [] npopRange = {1, sequenceVals[0]};
        int nu = sequenceVals[0];
        int lenSequence = sequenceVals[1];
        InputWriter.writeFile(fredMethodInput, bins,
            masterVariables.getOmegaRange(), masterVariables.getSigmaRange(),
            npopRange, masterVariables.getDriftRange(),
            masterVariables.getxnumics(), nu, masterVariables.getnrep(),
            lenSequence, -1);
        narr.println();
        narr.println("The input for Brute Force Search is: ");
        narr.writeInput(fredMethodInput);
        // Run fred method.
        log.append("Running Brute Force Search...\n");
        execs.runBruteforce();
        // Read in results.
        ArrayList<FredOutVal> results = fredValsReader.readFredOutVals(fredOutput);
        return results;
    }

    /**
     *  Print the criterion we are sorting at.
     */
    private void printSorting() {
        int sortPer = masterVariables.getSortPercentage();
        String criterion="";
        switch (sortPer) {
            case 0: criterion = "5x";
                    break;
            case 1: criterion = "2x";
                    break;
            case 2: criterion = "1.5x";
                    break;
            case 3: criterion = "1.25x";
                    break;
            case 4: criterion = "1.10x";
                    break;
            case 5: criterion = "1.05x";
                    break;
        }
        log.append("Now checking at " + criterion + "\n");
        narr.println("Now checking at " + criterion + "\n");
    }

    /**
     *  Finds the best value to run hillclimbing on.  It sorts the output data
     *  by the first column from the right that has a nonzero value, then runs
     *  full likelihood on the best data point from that column, then checks
     *  if it can run hill climbing in a reasonable amount of time (<13 hours)
     *  on that data point, if so, it returns that data point, otherwise it
     *  resorts by the next column to the left and continues.
     *
     *  @param data The output from a run of fred's method.
     *  @return FredOutVal containing the best value for hillclimbing.
     */
    private FredOutVal findBestHclimb(ArrayList<FredOutVal> data) {
        FredOutVal bestFred = getBestFred(data);
        int nrep;
        long startTime = -1;
        long stopTime = -1;
        long runTime = -1;
        double hourTime = -1;
        double predictedTime = -1;
        double [] percentages;
        int sortPer;
        FredOutVal actualLikelihoods = null;
        printSorting();
        log.append("(the criterion that values are sorted by)\n");
        narr.println("(the criterion that values are sorted by)");
        while (true) {
            // Keep running full likelihoods and timing hillclimbing until we find a good likelihood to run it on.
            log.append("\n");
            // Don't run full likelihood again if we already have them.
            if (actualLikelihoods == null) {
                log.append("Running full likelihood with the following result:\n");
                log.append(bestFred.toString() + "\n");
                narr.println();
                narr.println("Running full likelihood with the following result:");
                narr.println(bestFred.toString());
                // Run the full likelihood.
                startTime = System.currentTimeMillis();
                actualLikelihoods = fullLike(bestFred);
                stopTime = System.currentTimeMillis();
                runTime = stopTime - startTime;
                hourTime = (double) runTime / 3600000; // Find the number of hours.
                if (userChoice == -1) {
                    masterVariables.setSortPercentage(5);
                    // Try the least precise criterion until we find one that can be done in 8 hours or
                    // less, unless a criterion was pre-specified.
                }
            }
            percentages = actualLikelihoods.getPercentages();
            sortPer = masterVariables.getSortPercentage();
            while (percentages[sortPer] == 0.0 && sortPer > 0) {
                sortPer --;
                masterVariables.setSortPercentage(sortPer);
            }
            // Figure out the number of reps needed to get 50 successes.
            nrep = (int)(50 / percentages[sortPer]);
            predictedTime = NHCLIMB_REPS * (((double) nrep / NREP_FULL_LIKE) * hourTime);
            log.append("\n");
            log.append("Actual likelihood's for hillclimbing: \n");
            log.append(actualLikelihoods.toString() + "\n");
            narr.println("Actual likelihood's for hillclimbing: ");
            narr.println(actualLikelihoods.toString());
            printSorting();
            log.append("\nTesting hillclimbing...\n");
            narr.println();
            narr.println("Testing hillclimbing...");
            String predicted = String.format("Predicted Run time: %.2g hours", predictedTime);
            log.append(predicted + "\n");
            narr.println(predicted);
            if (userChoice != -1) {
                return actualLikelihoods;
            }
            if (predictedTime < numHours) {
                log.append("Can be run in a reasonable amount of time, continuing...\n");
                narr.println("Can be run in a reasonable amount of time, continuing...");
                return actualLikelihoods;
            }
            if (masterVariables.getSortPercentage() == 0) {
                log.append("Cannot decrease precision further, using last data point as hillclimb seed...\n");
                narr.println("Cannot decrease precision further, using last data point as hillclimb seed...");
                return actualLikelihoods;
            }
            log.append("Hill Climbing would take longer than " + numHours + " hours, decreasing precision...\n");
            narr.println("Hill Climbing would take longer than " + numHours + " hours, decreasing precision...");
            masterVariables.setSortPercentage(sortPer - 1);
        }
    }

    /**
     *  The number of repetitions of one hillclimbing cycle, usually 35-45 but 50 is safe.
     */
    private final int NHCLIMB_REPS = 50;

    /**
     *  The number of reps to do for a full likelihood.
     */
    private final int NREP_FULL_LIKE = 10000;

    /**
     *  The number of hours that we would not like to run hillclimbing longer than.
     */
    private int numHours = 8;

    /**
     *  The data from after the binning program runs.
     */
    private ArrayList<String> bins;

    /**
     *  The Sequence value data, where sequenceVals[0] is the number of sequences
     *  and sequenceVals[1] is the length of the sequence.
     */
    private int[] sequenceVals = new int[2];

    /**
     *  The input file for fred method.
     */
    private File fredMethodInput;

    /**
     *  The output from the fred method program.
     */
    private File fredOutput;

    /**
     *  The name of the input fasta file.
     */
    private Fasta inFasta;

    /**
     *  The JTextArea that is outputted to.
     */
    private JTextArea log;

    /**
     *  The narrative writer class.
     */
    private NarrWriter narr;

    /**
     *  The best fred value.
     */
    private FredOutVal bestFred;

    private int userChoice;
    private MasterVariables masterVariables;
    private Execs execs;
    private String workingDirectory;
    private FredValsReader fredValsReader;
    private BinningFasta binningFasta;
}
