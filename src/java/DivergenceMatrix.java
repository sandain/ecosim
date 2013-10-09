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
 *  Object to interact with the divergencematrix program.
 *
 *  @author Jason M. Wood
 *  @copyright GNU General Public License
 */
public class DivergenceMatrix implements Runnable {

    /**
     *  Object to interact with the divergencematrix program.
     *
     *  @param masterVariables The MasterVariables object.
     *  @param phylogeny The Phylogeny object.
     */
    public DivergenceMatrix(MasterVariables masterVariables,
        Phylogeny phylogeny) {
        this(masterVariables, phylogeny, "");
    }

    /**
     *  Object to interact with the divergencematrix program.
     *
     *  @param masterVariables The MasterVariables object.
     *  @param phylogeny The Phylogeny object.
     *  @param suffix The suffix to attach to the end of file names.
     */
    public DivergenceMatrix(MasterVariables masterVariables,
        Phylogeny phylogeny, String suffix) {
        this.masterVariables = masterVariables;
        this.phylogeny = phylogeny;

        String workingDirectory = masterVariables.getWorkingDirectory();
        inputFileName = 
            workingDirectory + "divergencematrixIn" + suffix + ".dat";
        outputFileName = 
            workingDirectory + "divergencematrixOut" + suffix + ".dat";
        hasRun = false;
    }

    /**
     *  Run the divergencematrix program.
     */
    public void run() {
        Execs execs = masterVariables.getExecs();
        File inputFile = new File(inputFileName);
        File outputFile = new File(outputFileName);
        // Write the input values for the program to the
        // divergencematrixIn.dat file.
        writeInputFile(inputFile);
        // Run the divergencematrix program.
        execs.runDivergencematrix(inputFile, outputFile);
        // Get the output provided by the divergencematrix program.
        readOutputFile(outputFile);
        // Set the flag stating that the divergencematrix program has been run.
        hasRun = true;
    }

    /**
     *  Returns the divergence matrix.
     *
     *  @return The divergence matrix.
     */
    public float[][] getMatrix() {
        return matrix;
    }

    /**
     *  Changes the divergence matrix.
     *
     *  @param matrix The new divergence matrix.
     */
    public void setMatrix(float[][] matrix) {
        this.matrix = matrix;
    }

    /**
     *  Returns true if divergencematrix has been run, false otherwise.
     *
     *  @return True if divergencematrix has been run, false otherwise.
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
     *  Private method to write the input file for the divergencematrix
     *  program.
     *
     *  @param inputFile The file to write to.
     */
    private void writeInputFile(File inputFile) {
        BufferedWriter writer = null;
        ArrayList<String> seqs = phylogeny.getSequences();
        try {
            writer = new BufferedWriter(new FileWriter(inputFile));
            writer.write(String.format(
                " %12d %12d\n",
                phylogeny.getNu(),
                phylogeny.length()
            ));
            for (int i = 0; i < seqs.size(); i ++) {
                writer.write(seqs.get(i) + "\n");
            }
        }
        catch (IOException e) {
            System.out.println("Error writing the input file for the " +
                               "divergence matrix program.");
        }
        finally {
            if (writer != null) {
                try {
                    writer.close();
                }
                catch (IOException e) {
                    System.out.println("Error closing the input file for the " +
                                       "divergence matrix program.");
                }
            }
        }
    }

    /**
     *  Private method to read the output file from the divergencematrix
     *  program.
     *
     *  @param outputFile The file to read from.
     */
    private void readOutputFile(File outputFile) {
        BufferedReader reader = null;
        int nu = phylogeny.getNu();
        matrix = new float[nu][nu];
        try {
            reader = new BufferedReader(new FileReader(outputFile));
            String nextLine = reader.readLine();
            nextLine = reader.readLine(); // nu, lengthsequence
            int i = 0;
            while (nextLine != null) {
                StringTokenizer st = new StringTokenizer(nextLine);
                for (int j = 0; j < nu; j ++) {
                    String buffer = st.nextToken();
                    matrix[i][j] = new Float(buffer).floatValue();
                }
                i ++;
                nextLine = reader.readLine();
            }
        }
        catch (IOException e) {
            System.out.println("Error reading the output file from the " +
                               "divergence matrix program.");
        }
        finally {
            if (reader != null) {
                try {
                    reader.close();
                }
                catch (IOException e) {
                    System.out.println("Error closing the output file from " +
                                       "the divergence matrix program.");
                }
            }
        }
    }

    private float[][] matrix;

    private MasterVariables masterVariables;
    private Phylogeny phylogeny;

    private String inputFileName;
    private String outputFileName;

    private boolean hasRun;

}
