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

import java.util.ArrayList;
import java.util.StringTokenizer;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import javax.swing.JTextArea;

/**
 *  The ConfidenceInterval class is extended by the different types of confidence intervals.
 *
 *  @author Andrew Warner
 */
public class ConfidenceInterval extends Thread {

    /**
     *  Constructor for objects of class ConfidenceInterval.
     *
     *  @param type Type of confidence interval {drift, npop, omega, sigma}.
     *  @param value The inputted value from hillclimbing that includes theoptimized omega, sigma, and npop.
     *  @param results Results of this confidence interval.
     *  @param masterVariables The MasterVariables.
     *  @param input File containing the values fred method to examine.
     *  @param output Narrative to write to.
     */
    protected ConfidenceInterval(String type, FredOutVal value, BinningAndFred results,
        MasterVariables masterVariables, File input, File output) {
        this.type = type;
        this.value = value;
        this.bins = results.getBins();
        this.masterVariables = masterVariables;
        this.input = input;
        this.output = output;
        narr = masterVariables.getNarrator();
        execs = masterVariables.getExecs();
        sequenceVals = results.getSeqVals();
    }

    /**
     *  Writes the input for the confidence interval.
     *
     *  @param input The input file to write to.
     *  @param range The range of variable for the current confidence interval.
     *  @param xnumics The values for xnumics for the given confidence interval.
     *  @param probThresh The probability threshold.
     */
    protected void writeInput(File input, double[] range, int[] xnumics, double probThresh) {
        this.probThresh = probThresh;
        double [] omegaRange, sigmaRange, driftRange;
        int [] npopRange;
        if (type.equals("omega")) {
            omegaRange = range;
        }
        else {
            omegaRange = new double [2];
            omegaRange[0] = value.getOmega();
            omegaRange[1] = 10000;
        }
        if (type.equals("sigma")) {
            sigmaRange = range;
        }
        else {
            sigmaRange = new double[2];
            sigmaRange[0] = value.getSigma();
            sigmaRange[1] = 10000;
        }
        if (type.equals("drift")) {
            driftRange = range;
        }
        else {
            driftRange = masterVariables.getDriftRange();
        }
        if (type.equals("npop") || type.equals("demarcation")) {
            npopRange = new int[2];
            npopRange[0] = (int)range[0];
            npopRange[1] = (int)range[1];
        }
        else {
            npopRange = new int[2];
            npopRange[0] = value.getNpop();
            npopRange[1] = 10000;
        }
        int jwhichxavg = masterVariables.getSortPercentage() + 1;
        double[] percentages = value.getPercentages();
        int numSuccesses = masterVariables.NUM_CI_SUCCESSES;
        int nrep = (int)(numSuccesses / probThresh);
        // Write the input file.
        InputWriter.writeFile(input, bins, omegaRange, sigmaRange, npopRange, driftRange, xnumics, sequenceVals[0],
            nrep, sequenceVals[1], jwhichxavg, probThresh);
    }

    /**
     *  Runs the confidence interval.
     *
     *  @return FredOutVal containing the confidence interval value.
     */
    protected FredOutVal runCI() {
        FredOutVal exitVal = null;
        if (type.equals("drift")) {
            execs.runDriftCI();
            exitVal = getOutputValue();
        }
        else if (type.equals("npop")) {
            execs.runNpopCI();
            exitVal = getOutputValue();
        }
        else if (type.equals("omega")) {
            execs.runOmegaCI();
            exitVal = getOutputValue();
        }
        else if (type.equals("sigma")) {
            execs.runSigmaCI();
            exitVal = getOutputValue();
        }
        else if (type.equals("demarcation")) {
            execs.runDemarcation();
        }
        return exitVal;
    }

    /**
     *  Reads in the output data from the output file for this particular
     *  confidence interval and returns the last data point recorded (either
     *  the first value out of confidence or the last value run; in the case
     *  of sigma if this value is still within confidence and sigma is 100
     *  then sigma is > 100).
     *
     *  @return FredOutVal containing the result from this run of the confidence interval.
     */
    protected FredOutVal getOutputValue() {
        try {
            BufferedReader input = new BufferedReader(new FileReader(output));
            ArrayList<String> values = new ArrayList<String>();
            String line = input.readLine();
            while (line != null) {
                values.add(line);
                line = input.readLine();
            }
            input.close();
            // These values have to be given temp values but they will be defined in the following for loop.
            double omega = -1;
            double sigma = -1;
            double drift = -1;
            double yvalue = -1;
            int npop = -1;
            // Return the last value if that value is in the confidence interval
            // otherwise return the last value that was in the confidence interval.
            for (int i = 0; i < 2; i ++) {
                StringTokenizer tk = new StringTokenizer(values.get(values.size() - (i + 1)));
                omega = (new Double(tk.nextToken())).doubleValue();
                sigma = (new Double(tk.nextToken())).doubleValue();
                npop = (new Integer(tk.nextToken())).intValue();
                drift = (new Double(tk.nextToken())).doubleValue();
                yvalue = (-1) * (new Double(tk.nextToken())).doubleValue();
                if (yvalue > probThresh) {
                    break;
                }
                if (values.size() < 2) {
                    break;
                }
            }
            int sortPer = masterVariables.getSortPercentage();
            double [] percentagesRes = new double[6];
            percentagesRes[sortPer] = yvalue;
            narr.println();
            narr.println("The output file for the " + type + " confidence interval:");
            narr.writeInput(output); // Write the whole output file to the narrative.
            return new FredOutVal(omega, sigma, npop, drift, percentagesRes, masterVariables);
        }
        catch (IOException e) {
            e.printStackTrace();
            return null;
        }
    }

    /**
     *  Sets the value of the current FredOutVal.
     *
     *  @param newVal The value to set.
     */
    public void setValue(FredOutVal newVal) {
        this.value = newVal;
    }

    /**
     *  Returns the last CIReader.
     *
     *  @return The last CIReader.
     */
    public CIReader getLastReader() {
        return lastReader;
    }

    /**
     *  Sets the upper bound.
     *
     *  @param newVal The new value of upperBound.
     */
    public void setUpperBound(double newVal) {
        this.upperBound = newVal;
    }

    /**
     *  The last reader used for the confidence interval.
     */
    protected CIReader lastReader;

    /**
     *  The value from hillclimbing that we will use to find the npop confidence interval.
     */
    protected FredOutVal value;

    /**
     *  The binning data.
     */
    protected ArrayList<String> bins;

    /**
     *  The narrator class.
     */
    protected NarrWriter narr;

    /**
     *  The sequence values.
     */
    protected int [] sequenceVals;

    /**
     *  The input file for the confidence interval program.
     */
    protected File input;

    /**
     *  The range of the confidence interval; ie, the number to divide
     *  the likelihood by to find out which likelihoods are out of confidence.
     */
    protected double percentRange;

    /**
     *  The type of confidence interval this is.
     */
    protected String type;

    /**
     *  The upper bound of the confidence interval program.
     */
    protected double upperBound;

    /**
     *  The lower bound for all confidence intervals is 1.0e-7 (besides npop).
     */
    protected double lowerBound = 1.0e-7;

    /**
     *  The number of values to do within the interval (ie the value to put in
     *  xnumics under the corresponding variable).
     */
    protected int incUpper;

    /**
     *  The number of increments to run for the lower bound.
     */
    protected int incLower;

    /**
     *  The output file that the confidence interval program will write to.
     */
    protected File output;

    /**
     *  The number of increments per order of magnitude.
     */
    protected int incPerOM;

    /**
     *  The probability threshold; ie, the likelihood value from hillclimbing
     *  divided by either 6.83 (two tailed) or 3.87 (one tailed); any likelihood
     *  less than this value is outside the confidence interval.
     */
    protected double probThresh;

    protected MasterVariables masterVariables;

    private Execs execs;

}
