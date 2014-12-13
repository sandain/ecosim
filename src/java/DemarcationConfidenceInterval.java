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

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.StringTokenizer;


/**
 *  Run the demarcation program.
 *
 *
 *  @author Jason M. Wood
 *  @copyright GNU General Public License
 */
public class DemarcationConfidenceInterval {

    /**
     *  Run the demarcation program.
     *
     *  @param masterVariables The MasterVariables object.
     *  @param nu The number of environmental sequences.
     *  @param sampleNu The number of sequences in the subsample.
     *  @param length The length of the sequences being analyzed.
     *  @param binning The Binning object.
     *  @param hillclimbResult The result from hillclimbing.
     */
    public DemarcationConfidenceInterval (MasterVariables masterVariables,
        int nu, int sampleNu, int length, Binning binning,
        ParameterSet<Double> hillclimbResult) {
        this (masterVariables, nu, sampleNu, length, binning, hillclimbResult, "");
    }

    /**
     *  Run the demarcation program.
     *
     *  @param masterVariables The MasterVariables object.
     *  @param nu The number of environmental sequences.
     *  @param sampleNu The number of sequences in the subsample.
     *  @param length The length of the sequences being analyzed.
     *  @param binning The Binning object.
     *  @param hillclimbResult The result from hillclimbing.
     *  @param suffix The suffix to attach to the end of file names.
     */
    public DemarcationConfidenceInterval (MasterVariables masterVariables,
        int nu, int sampleNu, int length, Binning binning,
        ParameterSet<Double> hillclimbResult, String suffix) {
        this.masterVariables = masterVariables;
        this.nu = nu;
        this.sampleNu = sampleNu;
        this.length = length;
        this.binning = binning;
        this.hillclimbResult = hillclimbResult;
        String workingDirectory = masterVariables.getWorkingDirectory ();
        inputFileName = workingDirectory + "demarcationIn" + suffix + ".dat";
        outputFileName = workingDirectory + "demarcationOut" + suffix +
            ".dat";
        hasRun = false;
    }

    /**
     *  Run the demarcation program.
     */
    public void run () {
        Execs execs = masterVariables.getExecs ();
        int hashCode = this.hashCode ();
        File inputFile = new File (inputFileName);
        File outputFile = new File (outputFileName);
        // Write the input values for the program to the input file.
        writeInputFile (inputFile);
        // Run the demarcation program.
        execs.runDemarcation (inputFile, outputFile);
        // Get the output provided by the demarcation program.
        readOutputFile (outputFile);
        // Set the flag stating that the demarcation program has run.
        if (result > 0) {
            hasRun = true;
        }
    }

    /**
     *  Returns true if the demarcation program has been run,
     *  false otherwise.
     *
     *  @return True if the demarcation program has been run,
     *  false otherwise.
     */
    public boolean hasRun () {
        return hasRun;
    }

    /**
     *  Changes the value of hasRun.
     *
     *  @param hasRun The new value of hasRun.
     */
    public void setHasRun (boolean hasRun) {
        this.hasRun = hasRun;
    }

    /**
     *  Returns the result of the demarcation program.
     *
     *  @return The result.
     */
    public int getResult () {
        return result;
    }

    /**
     *  Returns the likelihood of the result of the demarcation program.
     *
     *  @return The likelihood.
     */
    public double getLikelihood () {
        return likelihood;
    }

    /**
     *  Returns the demarcation result as a String.
     *
     *  @return the demarcation result.
     */
    public String toString () {
        int npop = hillclimbResult.getNpop ();
        return String.format ("%d (%d)", npop, result);
    }
    /**
     *  Private method to write the input file for the demarcation program.
     *
     *  @param inputFile The file to write to.
     */
    private void writeInputFile (File inputFile) {
        ArrayList<BinLevel> bins = binning.getBins ();
        BufferedWriter writer = null;
        try {
            writer = new BufferedWriter (new FileWriter (inputFile));
            writer.write (String.format ("%-20d numcrit\n", bins.size ()));
            // Output the crit levels and the number of bins.
            for (int j = 0; j < bins.size (); j ++) {
                writer.write (String.format (
                    "%-20.6f %-20d\n",
                    bins.get (j).getCrit (),
                    bins.get (j).getLevel ()
                ));
            }
            // Write the omega value.
            writer.write (String.format (
                "%-20.5f omega\n", hillclimbResult.getOmega ()
            ));
            // Write the sigma value.
            writer.write (String.format (
                "%-20.5f sigma\n", hillclimbResult.getSigma ()
            ));
            // Estimate the value of npop by multiplying the npop value found
            // in hillclimbing by the ratio of the sample size to the total
            // number of environmental sequences.
            int npop = hillclimbResult.getNpop () * sampleNu / nu;
            if (npop < 1) {
                npop = 1;
            }
            // Write the npop value.
            writer.write (
                String.format ("%-20d npop\n", npop)
            );
            // Write the step value.
            writer.write (String.format ("%-20d step\n", step));
            // Write the nu value.
            writer.write (
                String.format ("%-20d nu\n", sampleNu)
            );
            // Write the nrep value.
            writer.write (String.format ("%-20d nrep\n", nrep));
            // Create the random number seed; an odd integer less than nine
            // digits long.
            long iii = (long)(100000000 * Math.random ());
            if (iii % 2 == 0) {
                iii ++;
            }
            // Write the random number seed.
            writer.write (
                String.format ("%-20d iii (random number seed)\n", iii)
            );
            // Write the length of the sequences.
            writer.write (
                String.format (
                    "%-20d lengthseq (after deleting gaps, etc.)\n",
                    length
                )
            );
            // Write the whichavg value.
            int whichavg = masterVariables.getCriterion ();
            writer.write (String.format ("%-20d whichavg\n", whichavg));
            // Write the likelihoodsolution value.
            writer.write (
                String.format (
                    "%-20.5f likelihoodsolution\n",
                    hillclimbResult.getValue ()
                )
            );
        }
        catch (IOException e) {
            System.out.println ("Error writing the input file.");
        }
        finally {
            if (writer != null) {
                try {
                    writer.close ();
                }
                catch (IOException e) {
                    System.out.println ("Error closing the input file.");
                }
            }
        }
    }

    /**
     *  Private method to read the output file from the demarcation
     *  program.
     *
     *  @param outputFile The file to read from.
     */
    private void readOutputFile (File outputFile) {
        BufferedReader reader = null;
        try {
            reader = new BufferedReader (new FileReader (outputFile));
            String nextLine = reader.readLine ();
            while (nextLine != null) {
                StringTokenizer st = new StringTokenizer (nextLine);
                // The output contains the most likely npop and the likelihood
                // for that value.
                st.nextToken (); // "npop".
                result = new Integer (st.nextToken ()).intValue ();
                st.nextToken (); // "likelihood".
                likelihood = new Double (st.nextToken ()).doubleValue ();
                nextLine = reader.readLine ();
            }
        }
        catch (IOException e) {
            System.out.println ("Error reading the output file.");
        }
        finally {
            if (reader != null) {
                try {
                    reader.close ();
                }
                catch (IOException e) {
                    System.out.println ("Error closing the output file.");
                }
            }
        }
    }

    private String inputFileName;
    private String outputFileName;

    private MasterVariables masterVariables;
    private int nu;
    private int sampleNu;
    private int length;
    private Binning binning;
    private ParameterSet<Double> hillclimbResult;

    private int nrep = 1000;
    private int step = 3;
    private int result = 0;
    private double likelihood = 0.0;

    private boolean hasRun;

}

