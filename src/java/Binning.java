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
 *  Object to interact with the binning programs.
 *
 *  @author Jason M. Wood
 *  @copyright GNU General Public License
 */
public class Binning implements Runnable {

    /**
     *  Object to interact with the binning programs.
     *
     *  @param masterVariables The MasterVariables object.
     *  @param phylogeny The Phylogeny object.
     */
    public Binning(MasterVariables masterVariables, Phylogeny phylogeny) {
        this(masterVariables, phylogeny, "");
    }

    /**
     *  Object to interact with the binning programs.
     *
     *  @param masterVariables The MasterVariables object.
     *  @param phylogeny The Phylogeny object.
     *  @param suffix The suffix to attach to the end of file names.
     */
    public Binning(MasterVariables masterVariables, Phylogeny phylogeny,
        String suffix) {
        this.masterVariables = masterVariables;
        this.phylogeny = phylogeny;
        divergenceMatrix = new DivergenceMatrix(masterVariables, phylogeny);
        bins = new ArrayList<BinLevel>();
        String workingDirectory = masterVariables.getWorkingDirectory();
        inputFileName = workingDirectory + "binningIn" + suffix + ".dat";
        binLevelsFileName = workingDirectory + "binlevels" + suffix + ".dat";
        outputFileName = workingDirectory + "binningOut" + suffix + ".dat";
        hasRun = false;
    }

    /**
     *  Run the binning programs.
     */
    public void run() {
        Execs execs = masterVariables.getExecs();
        File inputFile = new File(inputFileName);
        File binLevelsFile = new File(binLevelsFileName);
        File outputFile = new File (outputFileName);
        // Run the divergence matrix program.
        divergenceMatrix.run();
        if (! divergenceMatrix.hasRun()) {
            return;
        }
        // Output the divergence matrix to be used by the binning program.
        writeInputFile(inputFile);
        // Output the binLevels file to be used by the binning program.
        writeBinLevelsFile(binLevelsFile);
        // Run the binning program.
        execs.runBinningdanny(
            inputFile, binLevelsFile, outputFile
        );
        // Read in the bin levels produced by the binning program.
        readOutputFile(outputFile);
        // Set the flag stating that the binning programs have been run.
        if (bins.size() == binLevels.length) {
            hasRun = true;
        }
    }

    /**
     *  Returns true if binning has been run, false otherwise.
     *
     *  @return True if binning has been run, false otherwise.
     */
    public boolean hasRun() {
        return hasRun;
    }

    /**
     *  Returns an ArrayList<BinLevel> containing the bin levels.
     *
     *  @return The bin levels.
     */
    public ArrayList<BinLevel> getBinLevels() {
        return bins;
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
     *  Change the bin levels stored in this object.
     *
     *  @param bins The new bin levels.
     */
    public void setBinLevels(ArrayList<BinLevel> bins) {
        this.bins = bins;
    }

    /**
     *  Add a bin level to this object.
     *
     *  @param binLevel The new BinLevel to add.
     */
    public void addBinLevel(BinLevel binLevel) {
        bins.add(binLevel);
    }

    /**
     *  Returns the binning result as a String.
     *
     *  @return the binning result.
     */
    public String toString() {
        String str = "";
        for (int i = 0; i < bins.size(); i ++) {
            str += bins.get(i).toString();
            if (i < bins.size() - 1) {
                str += ", ";
            }
        }
        return str;
    }

    /**
     *  Private method to write the input file for the binning program.
     *
     *  @param inputFile The file to write to.
     */
    private void writeInputFile(File inputFile) {
        BufferedWriter writer = null;
        float[][] matrix = divergenceMatrix.getMatrix();
        int nu = phylogeny.getNu();
        try {
            writer = new BufferedWriter(new FileWriter(inputFile));
            writer.write(String.format(
                " %12d %12d\n",
                nu,
                phylogeny.length()
            ));
            for (int i = 0; i < nu; i ++) {
                for (int j = 0; j < nu; j ++) {
                    writer.write(String.format(
                        " %7.4f",
                        matrix[i][j]
                    ));
                }
                writer.write("\n");
            }
        }
        catch (IOException e) {
            System.out.println("Error writing the input file for the " +
                               "binning program.");
        }
        finally {
            if (writer != null) {
                try {
                    writer.close();
                }
                catch (IOException e) {
                    System.out.println("Error closing the input file for the " +
                                       "binning program.");
                }
            }
        }
    }

    /**
     *  Private method to write the binlevels file for the binning program.
     *  
     *  @param binLevelsFile The binlevels.dat file.
     */
    private void writeBinLevelsFile(File binLevelsFile) {
        BufferedWriter writer = null;
        try {
            writer = new BufferedWriter(new FileWriter(binLevelsFile));
            writer.write(String.format(
                "%d\n",
                binLevels.length
            ));
            for (int i = 0; i < binLevels.length; i++) {
                writer.write(String.format(
                    "%5.3f\n",
                    binLevels[i]
                ));
            }
            writer.close();
        }
        catch (IOException e) {
            e.printStackTrace();
        }
    }

    /**
     *  Private method to read the bin levels from the output file.
     *
     *  @param outputFile The file to read from.
     */
    private void readOutputFile(File outputFile) {
        BufferedReader reader = null;
        bins = new ArrayList<BinLevel>();
        try {
            reader = new BufferedReader(new FileReader(outputFile));
            String nextLine = reader.readLine();
            while (nextLine != null) {
                StringTokenizer st = new StringTokenizer(nextLine);
                if (st.countTokens() == 2) {
                    Float crit = new Float(st.nextToken());
                    Integer value = new Integer(st.nextToken());
                    bins.add(new BinLevel(crit, value));
                }
                nextLine = reader.readLine();
            }
        }
        catch (IOException e) {
            System.out.println("Error reading the output file from the " +
                               "binning programs.");
        }
        finally {
            if (reader != null) {
                try {
                    reader.close();
                }
                catch (IOException e) {
                    System.out.println("Error closing the output file from " +
                                       "the binning programs.");
                }
            }
        }
    }

    private String inputFileName;
    private String binLevelsFileName;
    private String outputFileName;

    private MasterVariables masterVariables;
    private Phylogeny phylogeny;
    private DivergenceMatrix divergenceMatrix;

    private ArrayList<BinLevel> bins;

    private boolean hasRun;

    /**
     *  The default bin levels.
     */
    private double[] binLevels = {
        0.800, 0.850, 0.900, 0.950, 0.960, 0.970, 0.980, 0.990, 0.995, 1.000
    };

}
