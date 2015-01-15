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
import javax.swing.JTextArea;
import java.util.Hashtable;

/**
 *  This class holds the master variables that may need to be changed later into
 *  the project, so that the programmer will not have to go through each class
 *  and change them manually.
 *
 *  @author Andrew Warner & Jason Wood
 */
public class MasterVariables {

    /**
     *  The MasterVariables constructor.
     */
    public MasterVariables() {
        options = OptionsIO.load();
        debug = false;
        log = new JTextArea();
        String workingDirectory = options.getDirectoryOption(Options.WORKING_DIRECTORY);
        narrOut = new File(workingDirectory + "narrative.txt");
        narr = new NarrWriter(narrOut);
        noOutgroup = new File(workingDirectory + "noOutgroup.dat");
    }

    /**
     *  Returns the percentage to start sorting at, where 2 starts at 1.5x,
     *  1 starts at 2x, etc (see FredOutVal documentation).
     *
     *  @return The sort percentage.
     */
    public int getSortPercentage() {
        return sortPercentage;
    }

    /**
     *  Set the starting percentage for comparison, the default is 2, ie, values
     *  are sorted based on their 1.5x values, then by their 2x values, then by 5x.
     *  A value of 3 here would mean that output values are sorted instead by
     *  their 1.25x values, then by 1.5x, then by 2x, etc.
     *
     *  pre - newVal is between 2 and 5.
     *  post - The new value is stored in sortPercentage.
     *
     *  @param newVal The new value to be stored in sortPercentage.
     */
    public void setSortPercentage(int newVal) {
        sortPercentage = newVal;
    }

    /**
     *  Get the default range for omega.
     *
     *  @return double array containing the default range of omega.
     */
    public double[] getOmegaRange() {
        return omegaRange;
    }

    /**
     *  Get the default range for sigma.
     *
     *  @return double array containing the default range of sigma.
     */
    public double[] getSigmaRange() {
        return sigmaRange;
    }

    /**
     *  Get the default range for drift.
     *
     *  @return double array containing the Default range of drift.
     */
    public double[] getDriftRange() {
        return driftRange;
    }

    /**
     *  Get the default values of xnumics.
     *
     *  @return int array containing the Default values of xnumics.
     */
    public int[] getxnumics() {
        return xnumics;
    }

    /**
     *  Get the default number of reps.
     *
     *  @return The default number of reps.
     */
    public int getnrep() {
        return nrep;
    }

    /**
     *  Return the default log.
     *
     *  @return The default log.
     */
    public JTextArea getLog() {
        return log;
    }

    /**
     *  Returns the default NarrWriter.
     *
     *  @return The default Narrative.
     */
    public NarrWriter getNarrator() {
        return narr;
    }

    /**
     *  Clears the master NarrWriter, and returns it.
     *
     *  @return The an empty Narrative.
     */
    public NarrWriter clearNarrator() {
        narr.close();
        narr = new NarrWriter(narrOut);
        return narr;
    }

    /**
     *  Returns the Execs object.
     *
     *  @return The default Execs.
     */
    public Execs getExecs() {
        return execs;
    }

    public Options getOptions() {
        return options;
    }

    /**
     *  Returns the binary directory.
     *
     *  @return String containing the binary directory.
     */
    public String getBinaryDirectory() {
        return options.getDirectoryOption(Options.BINARY_DIRECTORY);
    }

    /**
     *  Returns the script directory.
     *
     *  @return String containing the script directory.
     */
    public String getScriptDirectory() {
        return options.getDirectoryOption(Options.SCRIPT_DIRECTORY);
    }

    /**
     *  Returns the help directory.
     *
     *  @return String containing the binary directory.
     */
    public String getHelpDirectory() {
        return options.getDirectoryOption(Options.HELP_DIRECTORY);
    }

    /**
     *  Returns the working directory.
     *
     *  @return String containing the working directory.
     */
    public String getWorkingDirectory() {
        return options.getDirectoryOption(Options.WORKING_DIRECTORY);
    }

    /**
     *  Returns the options file.
     *
     *  @return String containing the options file.
     */
    public String getOptionsFile() {
        return options.OPTIONS_FILE;
    }

    /**
     *  Returns the no outgroup file.
     *
     *  @return The no outgroup file.
     */
    public File getNoOutgroup () {
        return noOutgroup;
    }

    /**
     *  Set the options to those specified.
     *
     *  @param options The options to load.
     */
    public void setOptions (Options options) {
        this.options = options;
    }

    /**
     *  Returns the number of threads to use.
     *
     *  @return The number of threads to use.
     */
    public int getNumberThreads () {
        return numThreads;
    }

    /**
     *  Return the current debug status.
     *
     *  @return The current debug status.
     */
    public boolean getDebug() {
        return debug;
    }

    public void setExecs(Execs execs) {
        this.execs = execs;
    }

    /**
     *  Set the number of threads to use.
     *
     *  @param numThreads The new number of threads to use.
     */
    public void setNumberThreads (int numThreads) {
        this.numThreads = numThreads;
    }

    /**
     *  Set the current debug status.
     *
     *  @param debug The new debug status.
     */
    public void setDebug(boolean debug) {
        this.debug = debug;
    }

    /**
     *  The homogen coefficent.
     */
    public static final double HOMGEN_COEFF = 7.79;

    /**
     *  The confidence interval number for any one-tailed confidence interval.
     *  IE: to see if a value is in confidence, we take our likelihood from
     *  hill climbing and divide it by this number.
     */
    public static final double ONETAIL_CI_NUMBER = 3.87;

    /**
     *  The optimal number of successes that we want.
     */
    public static final int NUM_SUCCESSES = 50;

    /**
     *  The number of ideal successes for a run of the confidence interval.
     */
    public static final int NUM_CI_SUCCESSES = 20;

    /**
     *  The value to start from in the drift confidence interval.
     */
    public static final double DRIFT_CI_START = 1.0e4;

    /**
     *  The number to divide the likelihood by to get the confidence interval range.
     */
    public static final double CI_NUMBER = 6.83;

    /**
     *  The default number of threads is equal to the system maximum.
     */
    private int numThreads = Runtime.getRuntime ().availableProcessors ();

    /**
     *  The Ecotype Simulation version number.
     */
    private String version =
        Package.getPackage ("ecosim").getImplementationVersion ();

    /**
     *  The input file with outgroup removed.
     */
    private File noOutgroup;

    /**
     *  The index of the first value in the percentages array (in FredOutVal objects) to start sorting by.
     */
    private int sortPercentage = 5;

    /**
     *  The default range for omega.
     */
    private double[] omegaRange = {0.001, 100.0};

    /**
     *  The default range for sigma.
     */
    private double[] sigmaRange = {0.001, 100.0};

    /**
     *  The default range for drift.
     */
    private double[] driftRange = {1.0e25, 1.0e26};

    /**
     *  The default values for xnumics.
     */
    private int[] xnumics = {20, 20, 6, 0};

    /**
     *  The default number of reps to do of fred method.
     */
    private int nrep = 20;

    private boolean debug;

    private JTextArea log;
    private NarrWriter narr;
    private Execs execs;
    private File narrOut;
    private Options options;

}
