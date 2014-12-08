/*
 *    Ecotype Simulation models the sequence diversity within a bacterial
 *    clade as the evolutionary result of net ecotype formation and periodic
 *    selection, yielding a certain number of ecotypes.
 *
 *    Copyright (C) 2009       Andrew Warner, Wesleyan University
 *    Copyright (C) 2009-2014  Jason M. Wood, Montana State University
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

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;

/**
 *  This class holds the master variables that may need to be changed later
 *  into the project, so that the programmer will not have to go through each
 *  class and change them manually.
 *
 *  @author Andrew Warner
 *  @author Jason M. Wood
 *  @copyright GNU General Public License
 */
public class MasterVariables {

    /**
     *  The MasterVariables constructor.
     */
    public MasterVariables () {
        debug = false;
        // Create a temporary directory to serve as the working directory.
        try {
            Path tempDirectory = Files.createTempDirectory ("es2-");
            workingDirectory = tempDirectory.toString () +
                System.getProperty ("file.separator");
        }
        catch (IOException e) {
            e.printStackTrace ();
        }
        // Assume the the bin, help, and scripts directories are in the
        // current working directory.
        binaryDirectory = System.getProperty ("user.dir") +
            System.getProperty ("file.separator") + "bin" +
            System.getProperty ("file.separator");

        helpDirectory = System.getProperty ("user.dir") +
            System.getProperty ("file.separator") +"help" +
            System.getProperty ("file.separator");

        scriptDirectory = System.getProperty ("user.dir") +
            System.getProperty ("file.separator") + "scripts" +
            System.getProperty ("file.separator");

        // Initialize various objects.
        execs = new Execs (this);
        log = new Logger ();
    }

    /**
     *  Get the Ecotype Simulation version number.
     *
     *  @return The version number.
     */
    public String getVersion () {
        return version;
    }

    /**
     *  Get the current criterion value.
     *
     *  @return The current criterion value.
     */
    public int getCriterion () {
        return criterion;
    }

    /**
     *  Set the criterion value.
     *
     *  @param criterion The new criterion value.
     */
    public void setCriterion (int criterion) {
        this.criterion = criterion;
    }

    /**
     *  Get the label for the current criterion value.
     *
     *  @return The label for the current criterion value.
     */
    public String getCriterionLabel () {
       return getCriterionLabel (criterion);
    }

    /**
     *  Get the label for a criterion value.
     *
     *  @param crit The criterion value of interest.
     *  @return The level for the criterion value of interest.
     */
    public String getCriterionLabel (int crit) {
       return criterionLabels[crit - 1];
    }

    /**
     *  Return the array of criterion labels.
     *
     *  @return The array of criterion labels.
     */
    public String[] getCriterionLabels () {
        return criterionLabels;
    }

    /**
     *  Get the current PCR error value.
     *
     *  @return The current PCR error value.
     */
    public double getPCRError () {
        return pcrerror;
    }

    /**
     *  Set the PCR error value.
     *
     *  @param pcrerror The new PCR error value.
     */
    public void setPCRError (double pcrerror) {
        this.pcrerror = pcrerror;
    }

    /**
     *  Return the default log.
     *
     *  @return The default log.
     */
    public Logger getLog () {
        return log;
    }

    /**
     *  Returns the Execs object.
     *
     *  @return The default Execs.
     */
    public Execs getExecs () {
        return execs;
    }

    /**
     *  Returns the binary directory.
     *
     *  @return String containing the binary directory.
     */
    public String getBinaryDirectory () {
        return binaryDirectory;
    }

    /**
     *  Returns the script directory.
     *
     *  @return String containing the script directory.
     */
    public String getScriptDirectory () {
        return scriptDirectory;
    }

    /**
     *  Returns the help directory.
     *
     *  @return String containing the binary directory.
     */
    public String getHelpDirectory () {
        return helpDirectory;
    }

    /**
     *  Returns the working directory.
     *
     *  @return String containing the working directory.
     */
    public String getWorkingDirectory () {
        return workingDirectory;
    }

    /**
     *  Returns the location of the NJPlot program.
     *
     *  @return String containing the location of the program.
     */
    public String getProgramNJPlot () {
        return programNJPlot;
    }

    /**
     *  Returns the location of the Phylip DNAPars program.
     *
     *  @return String containing the location of the program.
     */
    public String getProgramDNAPars () {
        return programDNAPars;
    }

    /**
     *  Returns the location of the Phylip DNADist program.
     *
     *  @return String containing the location of the program.
     */
    public String getProgramDNADist () {
        return programDNADist;
    }

    /**
     *  Returns the location of the Phylip Neighbor program.
     *
     *  @return String containing the location of the program.
     */
    public String getProgramNeighbor () {
        return programNeighbor;
    }

    /**
     *  Returns the location of the Phylip Retree program.
     *
     *  @return String containing the location of the program.
     */
    public String getProgramRetree () {
        return programRetree;
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
    public boolean getDebug () {
        return debug;
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
    public void setDebug (boolean debug) {
        this.debug = debug;
    }

    /**
     *  Used for equality tests for floating point values.
     *
     *  ie: "a == 100.0" --> "a >= 100.0 - EPSILON"
     */
    public static final double EPSILON = 1.0e-6;

    /**
     *  The optimal number of successes that we want in bruteforce.
     */
    public static final int NUM_SUCCESSES = 50;

    /**
     *  The default criterion value for when auto is selected.
     *
     *  (1: 500%, 2: 200%, 3: 150%, 4: 125%, 5: 110%, 6: 105%)
     */
    private int criterion = 3;

    /**
     *  The criterion value labels.
     */
    private String[] criterionLabels = new String[] {
        "5.00x", "2.00x", "1.50x", "1.25x", "1.10x", "1.05x"
    };

    /**
     *  The default PCR error.
     *
     *  According to Akmaev and Wang, it is 1.37e-4, or 1/7300
     */
    private double pcrerror = 8.33e-6;

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
     *  The command name of the NJPlot program.
     */
    private String programNJPlot = "njplot";

    /**
     *  The command name of the Phylip DNAPars program.
     */
    private String programDNAPars = "dnapars";

    /**
     *  The command name of the Phylip DNADist program.
     */
    private String programDNADist = "dnadist";

    /**
     *  The command name of the Phylip Neighbor program.
     */
    private String programNeighbor = "neighbor";

    /**
     *  The command name of the Phylip Retree program.
     */
    private String programRetree = "retree";

    /**
     *  The location of the binary directory.
     */
    private String binaryDirectory;

    /**
     *  The location of the help directory.
     */
    private String helpDirectory;

    /**
     *  The location of the script directory.
     */
    private String scriptDirectory;

    /**
     *  The location of the working directory.
     */
    private String workingDirectory;

    /**
     *  Display the debug information if true.
     */
    private boolean debug;

    /**
     *  The Logger object.
     */
    private Logger log;

    /**
     *  The Execs object.
     */
    private Execs execs;

}
