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
 *  Object to interact with the hillclimbing program.
 *
 *  @author Jason M. Wood
 *  @copyright GNU General Public License
 */
public class Hillclimb implements Runnable {

    /**
     *  Run the hillclimb program.
     *
     *  @param masterVariables The MasterVariables object.
     *  @param phylogeny The Phylogeny object.
     *  @param binning The Binning object.
     *  @param parameterSet The set of parameters to optimize.
     */
    public Hillclimb(MasterVariables masterVariables, Phylogeny phylogeny,
        Binning binning, ParameterSet parameterSet) {
        this(masterVariables, phylogeny, binning, parameterSet, "");
    }

    /**
     *  Run the hillclimb program.
     *
     *  @param masterVariables The MasterVariables object.
     *  @param phylogeny The Phylogeny object.
     *  @param binning The Binning object.
     *  @param parameterSet The set of parameters to optimize.
     *  @param suffix The suffix to attach to the end of file names.
     */
    public Hillclimb(MasterVariables masterVariables, Phylogeny phylogeny,
        Binning binning, ParameterSet parameterSet, String suffix) {
        this.masterVariables = masterVariables;
        this.phylogeny = phylogeny;
        this.binning = binning;
        this.parameterSet = parameterSet;
        String workingDirectory = masterVariables.getWorkingDirectory();
        inputFileName = workingDirectory + "hillclimbIn" + suffix + ".dat";
        outputFileName = workingDirectory + "hillclimbOut" + suffix + ".dat";
        hasRun = false;
    }

    /**
     *  Run the hillclimb program.
     */
    public void run() {
        Execs execs = masterVariables.getExecs();
        File inputFile = new File(inputFileName);
        File outputFile = new File(outputFileName);
        // Write the input values for the program to the hclimbIn.dat file.
        writeInputFile(inputFile);
        // Run the hillclimb program.
        execs.runHillclimb(inputFile, outputFile);
        // Get the output provided by the hillclimb program.
        readOutputFile(outputFile);
        // Set the flag stating that the hillclimb program has been run.
        if (result.getNpop() > 0) {
            hasRun = true;
        }
    }

    /**
     *  Returns true if hillclimb has been run, false otherwise.
     *
     *  @return True if hillclimb has been run, false otherwise.
     */
    public boolean hasRun() {
        return hasRun;
    }

    /**
     *  Returns the result of the hillclimb program.
     *
     *  @return The result.
     */
    public ParameterSet<Double> getResult() {
        return result;
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
     *  Change the result stored in this object.
     *
     *  @param result The new result to store.
     */
    public void setResult(ParameterSet<Double> result) {
        this.result = result;
    }

    /**
     *  Returns the hillclimbing result as a String.
     *
     *  @return the hillclimbing result.
     */
    public String toString() {
        return result.toString();
    }

    /**
     *  Private method to write the input file for the hillclimb program.
     *
     *  @param inputFile The file to write to.
     */
    private void writeInputFile(File inputFile) {
        ArrayList<BinLevel> bins = binning.getBinLevels();
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
                String.format("%-20.5f omega\n", parameterSet.getOmega())
            );
            // Write the sigma value.
            writer.write(
                String.format("%-20.5f sigma\n", parameterSet.getSigma())
            );
            // Write the npop value.
            writer.write(
                String.format("%-20d npop\n", parameterSet.getNpop())
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

        }
        catch (IOException e) {
            System.out.println(
                "Error writing the input file for the hillclimb program."
            );
        }
        finally {
            if (writer != null) {
                try {
                    writer.close();
                }
                catch (IOException e) {
                    System.out.println(
                        "Error closing the input file for the hillclimb " +
                        "program."
                    );
                }
            }
        }
    }

    /**
     *  Private method to read the output file from the hillclimb program.
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
                // There should only be one line containing omega, sigma, npop,
                // and the likelihood of that result.
                Double omega = new Double(st.nextToken());
                Double sigma = new Double(st.nextToken());
                Integer npop = new Integer(st.nextToken());
                Double likelihood = new Double(st.nextToken());
                result = new ParameterSet<Double>(
                    omega, sigma, npop, likelihood
                );
                nextLine = reader.readLine();
            }
        }
        catch (IOException e) {
            System.out.println(
                "Error reading the output file from the hillclimb program."
            );
        }
        finally {
            if (reader != null) {
                try {
                    reader.close();
                }
                catch (IOException e) {
                    System.out.println(
                        "Error closing the output file from the hillclimb " +
                        "program."
                    );
                }
            }
        }
    }

    private String inputFileName;
    private String outputFileName;

    private MasterVariables masterVariables;
    private Phylogeny phylogeny;
    private Binning binning;
    private ParameterSet parameterSet;
    private ParameterSet<Double> result;

    private boolean hasRun;

}
