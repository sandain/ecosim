/*
 *    Ecotype Simulation models the sequence diversity within a bacterial
 *    clade as the evolutionary result of net ecotype formation and periodic
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
     *  @param fasta The Fasta object.
     */
    public Binning (MasterVariables masterVariables, Fasta fasta) {
        this (masterVariables, fasta, "");
    }

    /**
     *  Object to interact with the binning programs.
     *
     *  @param masterVariables The MasterVariables object.
     *  @param fasta The Fasta object.
     *  @param suffix The suffix to attach to the end of file names.
     */
    public Binning (MasterVariables masterVariables, Fasta fasta,
        String suffix) {
        this.masterVariables = masterVariables;
        this.fasta = fasta;
        bins = new ArrayList<BinLevel> ();
        String workingDirectory = masterVariables.getWorkingDirectory ();
        sequencesFileName = workingDirectory + "sequences" + suffix + ".dat";
        numbersFileName = workingDirectory + "numbers" + suffix + ".dat";
        rgFileName = workingDirectory + "removegaps" + suffix + ".dat";
        populationFileName = workingDirectory +
            "population" + suffix + ".dat";
        nameofstrainsFileName = workingDirectory +
            "namesofstrains" + suffix + ".dat";
        pcrerrorFileName = workingDirectory +
            "pcrerror" + suffix + ".dat";
        correctpcrFileName = workingDirectory +
            "correctpcr" + suffix + ".out";
        binningInputFileName = workingDirectory +
            "binningIn" + suffix + ".dat";
        binLevelsFileName = workingDirectory + "binlevels" + suffix + ".dat";
        binningOutputFileName = workingDirectory +
            "binningOut" + suffix + ".dat";
        divergenceMatrixInputFileName = workingDirectory +
            "divergencematrixIn" + suffix + ".dat";
        divergenceMatrixOutputFileName = workingDirectory +
            "divergencematrixOut" + suffix + ".dat";
        hasRun = false;
    }

    /**
     *  Run the binning programs.
     */
    public void run () {
        Execs execs = masterVariables.getExecs ();
        File sequencesFile = new File (sequencesFileName);
        File numbersFile = new File (numbersFileName);
        File rgFile = new File (rgFileName);
        File populationFile = new File (populationFileName);
        File nameofstrainsFile = new File (nameofstrainsFileName);
        File pcrerrorFile = new File (pcrerrorFileName);
        File correctpcrFile = new File (correctpcrFileName);
        File binningInputFile = new File (binningInputFileName);
        File binLevelsFile = new File (binLevelsFileName);
        File binningOutputFile = new File (binningOutputFileName);
        File divergenceMatrixInputFile = new File (
            divergenceMatrixInputFileName
        );
        File divergenceMatrixOutputFile = new File (
            divergenceMatrixOutputFileName
        );
        // Output the sequences.dat and numbers.dat files to be used by the
        // removegaps program.
        fasta.save (sequencesFile);
        writeNumbersFile (numbersFile, fasta.size (), fasta.length ());
        // Run the fasta file through the removegaps program.
        execs.runRemovegaps (sequencesFile, numbersFile, rgFile);
        // Run the readsynec program.
        // XXX readsynec seems to just make sequences lowercase, done in
        // Fasta.
        execs.runReadsynec (rgFile, populationFile, nameofstrainsFile);
        // Output the pcrerror.dat file to be used by the correctpcr program.
        writePCRErrorFile (pcrerrorFile, masterVariables.getPCRError ());
        // Run the correctpcr program.
        execs.runCorrectpcr (populationFile, pcrerrorFile, correctpcrFile);
        // Get the output provided by the correctpcr program.
        readCorrectPCROutputFile (correctpcrFile);
        // Write the input values for the program to the
        // divergencematrixIn.dat file.
        writeDivergenceMatrixInputFile (divergenceMatrixInputFile);
        // Run the divergencematrix program.
        execs.runDivergencematrix (
            divergenceMatrixInputFile, divergenceMatrixOutputFile
        );
        // Get the output provided by the divergencematrix program.
        readDivergenceMatrixOutputFile (divergenceMatrixOutputFile);
        // Output the divergence matrix to be used by the binning program.
        writeBinningInputFile (binningInputFile);
        // Output the binLevels file to be used by the binning program.
        writeBinLevelsFile (binLevelsFile);
        // Run the binning program.
        execs.runBinningdanny (
            binningInputFile, binLevelsFile, binningOutputFile
        );
        // Read in the bin levels produced by the binning program.
        readBinningOutputFile (binningOutputFile);
        // Set the flag stating that the binning programs have been run.
        if (bins.size () == binLevels.length) {
            hasRun = true;
        }
    }

    /**
     *  Returns true if binning has been run, false otherwise.
     *
     *  @return True if binning has been run, false otherwise.
     */
    public boolean hasRun () {
        return hasRun;
    }

    /**
     *  Returns an ArrayList<BinLevel> containing the bin levels.
     *
     *  @return The bin levels.
     */
    public ArrayList<BinLevel> getBinLevels () {
        return bins;
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
     *  Change the bin levels stored in this object.
     *
     *  @param bins The new bin levels.
     */
    public void setBinLevels (ArrayList<BinLevel> bins) {
        this.bins = bins;
    }

    /**
     *  Add a bin level to this object.
     *
     *  @param binLevel The new BinLevel to add.
     */
    public void addBinLevel (BinLevel binLevel) {
        bins.add (binLevel);
    }

    /**
     *  Returns the binning result as a String.
     *
     *  @return the binning result.
     */
    public String toString () {
        String str = "";
        for (int i = 0; i < bins.size (); i ++) {
            str += bins.get (i).toString ();
            if (i < bins.size () - 1) {
                str += ", ";
            }
        }
        return str;
    }

    /**
     *  Write the numbers file.
     *
     *  @param numbers The numbers.dat file.
     *  @param size The number of environmental sequences.
     *  @param length The length of the environmental sequences.
     */
    private void writeNumbersFile (File numbers, int size, int length) {
        BufferedWriter writer = null;
        try {
            writer = new BufferedWriter (new FileWriter (numbers));
            writer.write (String.format ("%d, %d\n", size, length));
        }
        catch (IOException e) {
            System.out.println ("Error writing the numbers.dat file.");
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
     *  Write the pcrerror file.
     *
     *  @param PCRErrorFile The pcrerror.dat file.
     *  @param PCRError The PCR error.
     */
    private void writePCRErrorFile (File PCRErrorFile, double PCRError) {
        // Create the random number seed; an odd less than 9 digits long
        long randValue = (long) (100000000 * Math.random ());
        if (randValue % 2 == 0) {
            randValue ++;
        }
        try {
            BufferedWriter writer = new BufferedWriter (
                new FileWriter (PCRErrorFile)
            );
            writer.write ("" + PCRError);
            writer.newLine ();
            writer.write ("" + randValue);
            writer.newLine ();
            writer.close ();
        }
        catch (IOException e) {
            e.printStackTrace ();
        }
    }

    /**
     *  Private method to write the input file for the binning program.
     *
     *  @param inputFile The file to write to.
     */
    private void writeBinningInputFile (File inputFile) {
        BufferedWriter writer = null;
        int nu = fasta.size ();
        try {
            writer = new BufferedWriter (new FileWriter (inputFile));
            writer.write (String.format (
                " %12d %12d\n",
                nu,
                fasta.length ()
            ));
            for (int i = 0; i < nu; i ++) {
                for (int j = 0; j < nu; j ++) {
                    writer.write (String.format (
                        " %7.4f",
                        matrix[i][j]
                    ));
                }
                writer.write ("\n");
            }
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
     *  Private method to write the binlevels file for the binning program.
     *
     *  @param binLevelsFile The binlevels.dat file.
     */
    private void writeBinLevelsFile (File binLevelsFile) {
        BufferedWriter writer = null;
        try {
            writer = new BufferedWriter (new FileWriter (binLevelsFile));
            writer.write (String.format (
                "%d\n",
                binLevels.length
            ));
            for (int i = 0; i < binLevels.length; i++) {
                writer.write (String.format (
                    "%5.3f\n",
                    binLevels[i]
                ));
            }
            writer.close ();
        }
        catch (IOException e) {
            e.printStackTrace ();
        }
    }

    /**
     *  Private method to write the input file for the divergencematrix
     *  program.
     *
     *  @param inputFile The file to write to.
     */
    private void writeDivergenceMatrixInputFile (File inputFile) {
        BufferedWriter writer = null;
        ArrayList<String> seqs = fasta.getSequences ();
        try {
            writer = new BufferedWriter (new FileWriter (inputFile));
            writer.write (String.format (
                " %12d %12d\n",
                fasta.size (),
                fasta.length ()
            ));
            for (int i = 0; i < seqs.size (); i ++) {
                writer.write (seqs.get (i) + "\n");
            }
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
     *  Private method to read the output file from the correctpcr
     *  program.
     *
     *  @param outputFile The file to read from.
     */
    private void readCorrectPCROutputFile (File outputFile) {
        BufferedReader reader = null;
        try {
            reader = new BufferedReader (new FileReader (outputFile));
            String nextLine = reader.readLine (); // nu, lengthsequence
            int i = 0;
            nextLine = reader.readLine ();
            while (nextLine != null) {
                StringTokenizer st = new StringTokenizer (nextLine);
                // Modify the sequence stored in the fasta file.
                fasta.setSequence (i, st.nextToken ());
                i ++;
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

    /**
     *  Private method to read the bin levels from the output file.
     *
     *  @param outputFile The file to read from.
     */
    private void readBinningOutputFile (File outputFile) {
        BufferedReader reader = null;
        bins = new ArrayList<BinLevel>();
        try {
            reader = new BufferedReader (new FileReader (outputFile));
            String nextLine = reader.readLine ();
            while (nextLine != null) {
                StringTokenizer st = new StringTokenizer (nextLine);
                if (st.countTokens () == 2) {
                    Float crit = new Float (st.nextToken ());
                    Integer value = new Integer (st.nextToken ());
                    bins.add (new BinLevel (crit, value));
                }
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

    /**
     *  Private method to read the output file from the divergencematrix
     *  program.
     *
     *  @param outputFile The file to read from.
     */
    private void readDivergenceMatrixOutputFile (File outputFile) {
        BufferedReader reader = null;
        int nu = fasta.size ();
        matrix = new float[nu][nu];
        try {
            reader = new BufferedReader (new FileReader (outputFile));
            String nextLine = reader.readLine ();
            nextLine = reader.readLine (); // nu, lengthsequence
            int i = 0;
            while (nextLine != null) {
                StringTokenizer st = new StringTokenizer (nextLine);
                for (int j = 0; j < nu; j ++) {
                    String buffer = st.nextToken ();
                    matrix[i][j] = new Float (buffer).floatValue ();
                }
                i ++;
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

    private String sequencesFileName;
    private String numbersFileName;
    private String rgFileName;
    private String populationFileName;
    private String nameofstrainsFileName;
    private String pcrerrorFileName;
    private String correctpcrFileName;
    private String binningInputFileName;
    private String binLevelsFileName;
    private String binningOutputFileName;
    private String divergenceMatrixInputFileName;
    private String divergenceMatrixOutputFileName;


    private MasterVariables masterVariables;
    private Fasta fasta;

    private ArrayList<BinLevel> bins;
    private float[][] matrix;

    private boolean hasRun;

    /**
     *  The default bin levels.
     */
    private double[] binLevels = {
        0.800, 0.850, 0.900, 0.950, 0.960, 0.970, 0.980, 0.990, 0.995, 1.000
    };

}
