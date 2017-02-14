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

/**
 *  Runs hillclimbing on a given FredOutVal.
 *
 *  @author Andrew Warner
 */
public class HillClimbing implements Runnable {

    /**
     *  Constructor for objects of class HillClimbing.
     *
     *  @param value The value to do hillclimbing on, and the eventual result of hillclimbing.
     *  @param masterVariables The MasterVariables.
     *  @param bins The binning data from the Binning process.
     *  @param sequenceVals The sequence values, where sequenceVals[0] is the number of sequences,
     *  and sequenceVals[1] is the length of each sequence after remove gaps.
     *  @param numSuccesses The ideal number of successes.
     */
    public HillClimbing(FredOutVal value, MasterVariables masterVariables, ArrayList<String> bins, int[] sequenceVals, int numSuccesses) {
        this.value = value;
        this.masterVariables = masterVariables;
        this.bins = bins;
        this.sequenceVals = sequenceVals;
        this.value = value;
        this.numSuccesses = numSuccesses;
        narr = masterVariables.getNarrator();
        execs = masterVariables.getExecs();
        workingDirectory = masterVariables.getWorkingDirectory();
        hillIn = new File(workingDirectory + "hClimbIn.dat");
        hillOut = new File(workingDirectory + "hClimbOut.dat");
    }

    /**
     *  Get the current FredOutVal.
     *
     *  @return FredOutVal containing the current value.
     */
    public FredOutVal getValue() {
        return value;
    }

    /**
     *  Runs hillclimbing on a full likelihood Fred value.
     */
    public void run() {
        hillClimb();
    }

    /**
     *  Runs hillclimbing on a full likelihood Fred value.
     */
    private void hillClimb() {
        // Declare default values.
        double [] omegaRange = {value.getOmega(), 10000};
        double[] sigmaRange = {value.getSigma(), 10000};
        int[] npopRange = {value.getNpop(), 10000};
        int[] xnumics = {0, 0, 0, 0};
        int jwhichxavg = masterVariables.getSortPercentage() + 1;
        double[] percentages = value.getPercentages();
        int nrep = (int)(numSuccesses / percentages[masterVariables.getSortPercentage()]);
        // Write hClimbIn.dat.
        InputWriter.writeFile(hillIn, bins, omegaRange, sigmaRange, npopRange, masterVariables.getDriftRange(), xnumics,
            sequenceVals[0], nrep, sequenceVals[1], jwhichxavg);
        narr.println();
        narr.println("The input for hillclimbing: ");
        narr.writeInput(hillIn);
        // Run hillclimb.exe.
        execs.runHillClimb();
        FredValsReader fredValsReader = new FredValsReader(masterVariables);
        value = fredValsReader.readHillClimbOutput(hillOut);
    }

    private String workingDirectory;

    /**
     *  The input file to write to for hill climbing.
     */
    private File hillIn;

    /**
     *  The output file from hillclimbing.
     */
    private File hillOut;

    /**
     *  The narrative writer for the program.
     */
    private NarrWriter narr;

    /**
     *  The value to do hillclimbing on, and the eventual result of hillclimbing.
     */
    private FredOutVal value;

    /**
     *  The binning data from the Binning process.
     */
    private ArrayList<String> bins;

    /**
     *  The sequence values, where sequenceVals[0] is the number of sequences,
     *  and sequenceVals[1] is the length of each sequence after remove gaps.
     */
    private int[] sequenceVals;

    /**
     *  The ideal number of successes.
     */
    private int numSuccesses;

    private Execs execs;

    private MasterVariables masterVariables;
}
