/*
 *    Ecotype Simulation models the sequence diversity within a bacterial
 *    clade as the evolutionary result of net ecotype formation and periodic
 *    selection, yielding a certain number of ecotypes.
 *
 *    Copyright (C) 2009       Andrew Warner, Wesleyan University
 *    Copyright (C) 2009-2015  Jason M. Wood, Montana State University
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
            tempDirectory = Files.createTempDirectory ("es2-");
            workingDirectory = tempDirectory.toString () +
                System.getProperty ("file.separator");
        }
        catch (IOException e) {
            System.err.println (
                "Error creating the temporary directory: " + e
            );
            System.exit (1);
        }
        // Assume the the bin and help directories are in the current working
        // directory.
        binaryDirectory = System.getProperty ("user.dir") +
            System.getProperty ("file.separator") + "bin" +
            System.getProperty ("file.separator");

        helpDirectory = System.getProperty ("user.dir") +
            System.getProperty ("file.separator") +"help" +
            System.getProperty ("file.separator");

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
    public String getCriterionLabel (Integer crit) {
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
     *  Returns the binary directory.
     *
     *  @return String containing the binary directory.
     */
    public String getBinaryDirectory () {
        return binaryDirectory;
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
     *  Returns the number of threads to use.
     *
     *  @return The number of threads to use.
     */
    public Integer getNumberThreads () {
        return numThreads;
    }

    /**
     *  Return the current debug status.
     *
     *  @return The current debug status.
     */
    public Boolean getDebug () {
        return debug;
    }

    /**
     *  Set the number of threads to use.
     *
     *  @param numThreads The new number of threads to use.
     */
    public void setNumberThreads (Integer numThreads) {
        this.numThreads = numThreads;
    }

    /**
     *  Set the current debug status.
     *
     *  @param debug The new debug status.
     */
    public void setDebug (Boolean debug) {
        this.debug = debug;
    }

    /**
     *  Delete the temporary directory and its contents.
     */
    public void exit () {
        try {
            if (! debug) {
                // Delete all of the files in the temporary directory.
                for (Path file: Files.newDirectoryStream (tempDirectory)) {
                    Files.delete (file);
                }
                // Delete the temporary directory.
                Files.delete (tempDirectory);
            }
        }
        catch (IOException e) {
            System.err.println (
                "Error deleting the temporary directory: " + e
            );
        }
    }

    /**
     *  Used for equality tests for floating point values.
     *
     *  ie: "a == 100.0" --> "a >= 100.0 - EPSILON"
     */
    public static final Double EPSILON = 1.0e-6;

    /**
     *  The default criterion value for when auto is selected.
     *
     *  (1: 500%, 2: 200%, 3: 150%, 4: 125%, 5: 110%, 6: 105%)
     */
    private Integer criterion = 3;

    /**
     *  The criterion value labels.
     */
    private String[] criterionLabels = new String[] {
        "5.00x", "2.00x", "1.50x", "1.25x", "1.10x", "1.05x"
    };

    /**
     *  The default number of threads is equal to the system maximum.
     */
    private Integer numThreads = Runtime.getRuntime ().availableProcessors ();

    /**
     *  The Ecotype Simulation version number.
     */
    private String version =
        Package.getPackage ("ecosim").getImplementationVersion ();

    /**
     *  The location of the temporary directory.
     */
    private Path tempDirectory;

    /**
     *  The location of the binary directory.
     */
    private String binaryDirectory;

    /**
     *  The location of the help directory.
     */
    private String helpDirectory;

    /**
     *  The location of the working directory.
     */
    private String workingDirectory;

    /**
     *  Display the debug information if true.
     */
    private Boolean debug;

}
