/*
 *    Ecotype Simulation models the sequence diversity within a bacterial
 *    clade as the evolutionary result of net ecotype formation and periodic
 *    selection, yielding a certain number of ecotypes.
 *
 *    Copyright (C) 2013-2018  Jason M. Wood, Montana State University
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

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.StringTokenizer;

/**
 *  Run the npop confidence interval program.
 *
 *
 *  @author Jason M. Wood
 *  @copyright GNU General Public License
 */
public class NpopConfidenceInterval implements Runnable {

    /**
     *  Run the npop confidence interval program.
     *
     *  @param masterVariables The MasterVariables object.
     *  @param execs The Execs object.
     *  @param nu The number of environmental sequences.
     *  @param length The length of the sequences being analyzed.
     *  @param binning The Binning object.
     *  @param hillclimbResult The result from hillclimbing.
     */
    public NpopConfidenceInterval (MasterVariables masterVariables,
        Execs execs, Integer nu, Integer length, Binning binning,
        ParameterSet hillclimbResult) {
        this.masterVariables = masterVariables;
        this.execs = execs;
        this.nu = nu;
        this.length = length;
        this.binning = binning;
        this.hillclimbResult = hillclimbResult;
        String workingDirectory = masterVariables.getWorkingDirectory ();
        inputFileName = workingDirectory + "npopIn.dat";
        outputFileName = workingDirectory + "npopOut.dat";
        hasRun = false;
    }

    /**
     *  Run the npop confidence interval program.
     */
    public void run () {
        File inputFile = new File (inputFileName);
        File outputFile = new File (outputFileName);
        // Write the input values for the program to the npopIn.dat file.
        writeInputFile (inputFile);
        // Run the npopCI program.
        execs.runNpopCI (inputFile, outputFile);
        // Get the output provided by the npopCI program.
        readOutputFile (outputFile);
        // Set the flag stating that the confidence interval program has run.
        if (result[0] > 0L && result[1] > 0L) {
            hasRun = true;
        }
    }

    /**
     *  Returns true if the npop confidence interval has been run, false
     *  otherwise.
     *
     *  @return True if the npop confidence interval has been run, false
     *  otherwise.
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
     *  Returns the result of the npop confidence interval program.
     *
     *  @return The result.
     */
    public Long [] getResult () {
        return result;
    }

    /**
     *  Returns the likelihood of the result of the npop confidence interval
     *  program.
     *
     *  @return The likelihood.
     */
    public Double [] getLikelihood () {
        return likelihood;
    }

    /**
     *  Returns the npop confidence interval as a String.
     *
     *  @return the npop confidence interval.
     */
    public String toString () {
        return String.format (
            "%d to %d (%.4g, %.4g)",
            result[0], result[1], likelihood[0], likelihood[1]
        );
    }

    /**
     *  Change the lower result stored in this object.
     *
     *  @param result The new result to store.
     *  @param likelihood The likelihood of the result.
     */
    public void setLowerResult (Long result, Double likelihood) {
        this.result[0] = result;
        this.likelihood[0] = likelihood;
    }

    /**
     *  Change the upper result stored in this object.
     *
     *  @param result The new result to store.
     *  @param likelihood The likelihood of the result.
     */
    public void setUpperResult (Long result, Double likelihood) {
        this.result[1] = result;
        this.likelihood[1] = likelihood;
    }

    /**
     *  Private method to write the input file for the npop confidence
     *  interval program.
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
            // Write the npop value.
            writer.write (
                String.format ("%-20d npop\n", hillclimbResult.getNpop ())
            );
            // Write the step value.
            writer.write (String.format ("%-20d step\n", step));
            // Write the nu value.
            writer.write (
                String.format ("%-20d nu\n", nu)
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
                    hillclimbResult.getLikelihood ()
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
     *  Private method to read the output file from the npop confidence
     *  interval program.
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
                // The first line contains the upper value of the confidence
                // interval for npop, and the likelihood for that values.
                // The second line contains the lower value.
                String upperLower = st.nextToken (); // "upper" or "lower".
                st.nextToken (); // "bound".
                st.nextToken (); // "npop".
                int index;
                switch (upperLower) {
                    case "lower": index = 0;
                                  break;
                    case "upper": index = 1;
                                  break;
                    default:      System.out.println (
                                      "Unexpected error in input file: " +
                                      outputFile.getName ()
                                  );
                                  return;
                }
                result[index] = new Long (st.nextToken ());
                st.nextToken (); // "likelihood".
                likelihood[index] = new Double (st.nextToken ());
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
    private Execs execs;
    private Integer nu;
    private Integer length;
    private Binning binning;
    private ParameterSet hillclimbResult;

    private Integer nrep = 1000;
    private Integer step = 3;
    private Long [] result = { 0L, 0L };
    private Double [] likelihood = { 0.0d, 0.0d };

    private boolean hasRun;

}

