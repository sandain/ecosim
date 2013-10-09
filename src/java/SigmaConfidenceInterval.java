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

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.StringTokenizer;

/**
 *  Run the sigma confidence interval program.
 *
 *  @author Jason M. Wood
 *  @copyright GNU General Public License
 */
public class SigmaConfidenceInterval {

    /**
     *  Run the sigma confidence interval program.
     *
     *  @param masterVariables The MasterVariables object.
     *  @param phylogeny The Phylogeny object.
     *  @param binning The Binning object.
     *  @param hillclimb The Hillclimb object.
     */
    public SigmaConfidenceInterval(MasterVariables masterVariables,
        Phylogeny phylogeny, Binning binning, Hillclimb hillclimb) {
        this(masterVariables, phylogeny, binning, hillclimb, "");
    }

   /**
     *  Run the sigma confidence interval program.
     *
     *  @param masterVariables The MasterVariables object.
     *  @param phylogeny The Phylogeny object.
     *  @param binning The Binning object.
     *  @param hillclimb The Hillclimb object.
     *  @param suffix The suffix to attach to the end of file names.
     */
    public SigmaConfidenceInterval(MasterVariables masterVariables,
        Phylogeny phylogeny, Binning binning, Hillclimb hillclimb,
        String suffix) {
        this.masterVariables = masterVariables;
        this.phylogeny = phylogeny;
        this.binning = binning;
        this.hillclimb = hillclimb;
        String workingDirectory = masterVariables.getWorkingDirectory();
        inputFileName = workingDirectory + "sigmaIn" + suffix + ".dat";
        outputFileName = workingDirectory + "sigmaOut" + suffix + ".dat";
        hasRun = false;
    }

    /**
     *  Run the sigma confidence interval program.
     */
    public void run() {
        Execs execs = masterVariables.getExecs();
        File inputFile = new File(inputFileName);
        File outputFile = new File(outputFileName);
        // Write the input values for the program to the sigmaIn.dat file.
        writeInputFile(inputFile);
        // Run the sigmaCI program.
        execs.runSigmaCI(inputFile, outputFile);
        // Get the output provided by the sigmaCI program.
        readOutputFile(outputFile);
        // Set the flag stating that the confidence interval program has run.
        if (result[0] > 0.0 && result[1] > 0.0) {
            hasRun = true;
        }
    }

    /**
     *  Returns true if the sigma confidence interval has been run, false
     *  otherwise.
     *
     *  @return True if the sigma confidence interval has been run, false
     *  otherwise.
     */
    public boolean hasRun() {
        return hasRun;
    }

    /**
     *  Changes the value of hasRun.
     *
     *  @param hasRun The new value of hasRun.
     */
    public void setHasRun(boolean hasRun) {
        this.hasRun = hasRun;
    }

    /**
     *  Returns the result of the sigma confidence interval program.
     *
     *  @return The result.
     */
    public double [] getResult() {
        return result;
    }

    /**
     *  Returns the likelihood of the result of the sigma confidence interval
     *  program.
     *
     *  @return The likelihood.
     */
    public double [] getLikelihood() {
        return likelihood;
    }

    /**
     *  Returns the sigma confidence interval as a String.
     *
     *  @return the sigma confidence interval.
     */
    public String toString() {
        String high;
        if (result[1] >= 100.0 - masterVariables.EPSILON) {
            high = ">100";
        }
        else {
            high = String.format("%.5g", result[1]);
        }
        return String.format(
            "%.5g to %s (%.5g, %.5g)",
            result[0], high, likelihood[0], likelihood[1]
        );
    }

    /**
     *  Change the lower result stored in this object.
     *
     *  @param result The new result to store.
     *  @param likelihood The likelihood of the result.
     */
    public void setLowerResult(double result, double likelihood) {
        this.result[0] = result;
        this.likelihood[0] = likelihood;
    }

    /**
     *  Change the upper result stored in this object.
     *
     *  @param result The new result to store.
     *  @param likelihood The likelihood of the result.
     */
    public void setUpperResult(double result, double likelihood) {
        this.result[1] = result;
        this.likelihood[1] = likelihood;
    }

    /**
     *  Private method to write the input file for the sigma confidence
     *  interval program.
     *
     *  @param inputFile The file to write to.
     */
    private void writeInputFile(File inputFile) {
        ArrayList<BinLevel> bins = binning.getBinLevels();
        ParameterSet<Double> hillclimbResult = hillclimb.getResult();
        BufferedWriter writer = null;
        try {
            writer = new BufferedWriter(new FileWriter(inputFile));
            writer.write(String.format("%-20d numcrit\n", bins.size()));
            // Output just the number of bins at each crit level.
            for (int i = 0; i < bins.size(); i ++) {
                writer.write(String.format("%-20d\n", bins.get(i).getLevel()));
            }
            // Output the crit levels and the number of bins.
            for (int j = 0; j < bins.size(); j ++) {
                writer.write(String.format(
                    "%-20.6f %-20d\n",
                    bins.get(j).getCrit(),
                    bins.get(j).getLevel()
                ));
            }
            // Write the omega value.
            writer.write(
                String.format("%-20.5f omega\n", hillclimbResult.getOmega())
            );
            // Write the sigma value.
            writer.write(
                String.format("%-20.5f sigma\n", hillclimbResult.getSigma())
            );
            // Write the npop value.
            writer.write(
                String.format("%-20d npop\n", hillclimbResult.getNpop())
            );
            // Write the step value.
            writer.write(
                String.format("%-20d step\n", masterVariables.getStep())
            );
            // Write the nu value.
            writer.write(String.format("%-20d nu\n", phylogeny.getNu()));
            // Write the nrep value.
            writer.write(
                String.format("%-20d nrep\n", masterVariables.getNrep())
            );
            // Create the random number seed; an odd integer less than nine
            // digits long.
            long iii = (long)(100000000 * Math.random());
            if (iii % 2 == 0) {
                iii ++;
            }
            // Write the random number seed.
            writer.write(
                String.format("%-20d iii (random number seed)\n", iii)
            );
            // Write the length of the sequences.
            writer.write(
                String.format(
                    "%-20d lengthseq (after deleting gaps, etc.)\n",
                    phylogeny.length()
                )
            );
            // Write the whichavg value.
            int whichavg = masterVariables.getCriterion();
            writer.write(String.format("%-20d whichavg\n", whichavg));
            // Write the likelihoodsolution value.
            writer.write(
                String.format(
                    "%-20.5f likelihoodsolution\n",
                    hillclimbResult.getValue()
                )
            );
        }
        catch (IOException e) {
            System.out.println("Error writing the input file for the " +
                               "confidence interval program.");
        }
        finally {
            if (writer != null) {
                try {
                    writer.close();
                }
                catch (IOException e) {
                    System.out.println("Error closing the input file for the " +
                                       "confidence interval program.");
                }
            }
        }
    }

    /**
     *  Private method to read the output file from the sigma confidence
     *  interval program.
     *
     *  @param outputFile The file to read from.
     */
    private void readOutputFile(File outputFile) {
        BufferedReader reader = null;
        try {
            reader = new BufferedReader(new FileReader(outputFile));
            String nextLine = reader.readLine();
            while (nextLine != null) {
                StringTokenizer st = new StringTokenizer(nextLine);
                // The first line contains the upper value of the confidence
                // interval for sigma, and the likelihood for that values.  The
                // second line contains the lower value.
                String upperLower = st.nextToken(); // "upper" or "lower".
                st.nextToken(); // "bound".
                st.nextToken(); // "sigma".
                int index;
                switch (upperLower) {
                    case "lower": index = 0;
                                  break;
                    case "upper": index = 1;
                                  break;
                    default:      System.out.println(
                                      "Unexpected error in input file: " +
                                      outputFile.getName()
                                  );
                                  return;
                }
                result[index] = new Double(st.nextToken()).doubleValue();
                st.nextToken(); // "likelihood".
                likelihood[index] = new Double(st.nextToken()).doubleValue();
                nextLine = reader.readLine();
            }
        }
        catch (IOException e) {
            System.out.println("Error reading the output file from the " +
                               "confidence interval program.");
        }
        finally {
            if (reader != null) {
                try {
                    reader.close();
                }
                catch (IOException e) {
                    System.out.println("Error closing the output file from " +
                                       "the confidence interval program.");
                }
            }
        }
    }

    private String inputFileName;
    private String outputFileName;

    private MasterVariables masterVariables;
    private Phylogeny phylogeny;
    private Binning binning;
    private Hillclimb hillclimb;

    private double [] result = { 0.0, 0.0 };
    private double [] likelihood = { 0.0, 0.0 };

    private boolean hasRun;

}

