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

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintStream;
import javax.swing.JTextArea;

/**
 *  Holds the executable methods for the cohan program.
 *
 *  @author Andrew Warner & Jason Wood
 */
public class Execs {

    /**
     *  Detects the Operating System that this class is running on,
     *  and executes the native fortran applications.
     *
     *  @param masterVariables The MasterVariables.
     */
    public Execs(MasterVariables masterVariables) {
        this.masterVariables = masterVariables;
        options = masterVariables.getOptions();
        // Grab the running environment.
        String osName = System.getProperty("os.name").toLowerCase();
        String osArch = System.getProperty("os.arch").toLowerCase();
        String osVersion = System.getProperty("os.version");
        binaryDirectory = masterVariables.getBinaryDirectory ();
        workingDirectory = masterVariables.getWorkingDirectory ();
        // Setup the rest of the variables.
        log = masterVariables.getLog();
        pcrerror = new File (workingDirectory + "pcrerror.dat");
        // Check which OS we are running on.
        if (osName.contains("windows")) {
            binaryExtension = ".exe";
            // Use the Phylip batch script.
            phylipScript = workingDirectory + "phylip.bat";
            writePhylipScript(phylipBatchScript);
        }
        else if (osName.contains("linux")) {
            if (osArch.contains("i386")) {
                binaryExtension = ".i386";
            }
            else if (osArch.contains("amd64")) {
                binaryExtension = ".amd64";
            }
            else {
                log.append ("Unsupported Linux architecture, contact the developers.");
                log.append (System.getProperty ("line.separator"));
                log.append("Architecture detected: " + osArch);
                log.append (System.getProperty ("line.separator"));
            }
            // Use the Phylip shell script.
            phylipScript = workingDirectory + "phylip.sh";
            writePhylipScript(phylipShellScript);
        }
        else if (osName.contains("mac")) {
            binaryExtension = ".app";
            // Use the Phylip shell script.
            phylipScript = workingDirectory + "phylip.sh";
            writePhylipScript(phylipShellScript);
        }
        else {
            log.append("Unsupported OS, contact the developers.");
            log.append (System.getProperty ("line.separator"));
            log.append("OS detected: " + osName);
            log.append (System.getProperty ("line.separator"));
            log.append("Architecture detected: " + osArch);
            log.append (System.getProperty ("line.separator"));
        }
    }

    /**
     *  Runs the initial Bruteforce method to find a good value with which to do hill climbing.
     *
     *  @return The exit value, -1 if there was an error.
     */
    public int runBruteforce() {
        String[] command = {
            binaryDirectory + "bruteforce" + binaryExtension,
            workingDirectory + "bruteforceIn.dat",
            workingDirectory + "bruteforceOut.dat",
            workingDirectory + "time.dat",
            Integer.toString (masterVariables.getNumberThreads ()),
            Boolean.toString(masterVariables.getDebug())
        };
        PrintStream errorStream = null;
        PrintStream outputStream = null;
        // Catch program output if debugging is enabled.
        if (masterVariables.getDebug ()) {
            errorStream = System.err;
            outputStream = System.out;
        }
        return runApplication (
            command,
            errorStream,
            "ERROR (Brute Force Search)>",
            outputStream,
            "Brute Force Search>",
            true
        );
    }

    /**
     *  Runs the binning process.
     *
     *  Pre: sequencesfasta and numbers.dat are formatted correctly.  Also
     *  the programs removegaps.exe, readsynec.exe, etc must all be in the
     *  root folder.
     *
     *  Post: output.dat contains the binning data
     *
     *  @return The exit value.
     */
    public int runBinning() {
        int exitVal = -1;
        checkBinLevels();
        makeRandPCRError();
        exitVal = runRemovegaps();
        exitVal += runReadsynec();
        exitVal += runCorrectpcr();
        exitVal += runDivergencematrix();
        exitVal += runBinningdanny();
        return exitVal;
    }

    /**
     *  Runs the hillclimb program.
     *
     *  @return The exit value.
     */
    public int runHillClimb() {
        String[] command = {
            binaryDirectory + "hillclimb" + binaryExtension,
            workingDirectory + "hClimbIn.dat",
            workingDirectory + "hClimbOut.dat",
            Integer.toString (masterVariables.getNumberThreads ()),
            Boolean.toString(masterVariables.getDebug())
        };
        PrintStream errorStream = null;
        PrintStream outputStream = null;
        // Catch program output if debugging is enabled.
        if (masterVariables.getDebug ()) {
            errorStream = System.err;
            outputStream = System.out;
        }
        return runApplication (
            command,
            errorStream,
            "ERROR (Hill Climb)>",
            outputStream,
            "Hill Climb>",
            true
        );
    }

    /**
     *  Opens the tree using NJPlot so that the user can see it to start analysis.
     *
     *  @param treeFile The tree file to open in NJPlot.
     *  @return The exit value.
     */
    public int openTree(String treeFile) {
        String[] command = {
            options.getHelperProgram(Options.NJPLOT),
            treeFile
        };
        return runApplication (
            command,
            System.err,
            "ERROR (NJPlot)>",
            System.out,
            "NJPlot>",
            false
        );
    }

    /**
     *  Runs the dnapars application
     *
     *  @return The exit value.
     */
    public int runDNAPars() {
        return runPhylip (
            Options.DNAPARS,
            "V" + System.getProperty ("line.separator") +
            "1" + System.getProperty ("line.separator") +
            "Y" + System.getProperty ("line.separator")
        );
    }

    /**
     *  Runs the dnadist application.
     *
     *  @return The exit value.
     */
    public int runDNADist() {
        return runPhylip (
            Options.DNADIST,
            "Y" + System.getProperty ("line.separator")
        );
    }

    /**
     *  Runs the Neighbor application.
     *
     *  @return The exit value.
     */
    public int runNJ() {
        return runPhylip (
            Options.NEIGHBOR,
            "Y" + System.getProperty ("line.separator")
        );
    }

    /**
     *  Runs the Drift Confidence Interval application.
     *
     *  @return The exit value.
     */
    public int runDriftCI() {
        String[] command = {
            binaryDirectory + "driftCI" + binaryExtension,
            workingDirectory + "driftIn.dat",
            workingDirectory + "driftOut.dat",
            Integer.toString (masterVariables.getNumberThreads ()),
            Boolean.toString(masterVariables.getDebug())
        };
        PrintStream errorStream = null;
        PrintStream outputStream = null;
        // Catch program output if debugging is enabled.
        if (masterVariables.getDebug ()) {
            errorStream = System.err;
            outputStream = System.out;
        }
        return runApplication (
            command,
            errorStream,
            "ERROR (Drift CI)>",
            outputStream,
            "Drift CI>",
            true
        );
    }

    /**
     *  Runs the Npop Confidence Interval application.
     *
     *  @return The exit value.
     */
    public int runNpopCI() {
        String[] command = {
            binaryDirectory + "npopCI" + binaryExtension,
            workingDirectory + "npopIn.dat",
            workingDirectory + "npopOut.dat",
            Integer.toString (masterVariables.getNumberThreads ()),
            Boolean.toString(masterVariables.getDebug())
        };
        PrintStream errorStream = null;
        PrintStream outputStream = null;
        // Catch program output if debugging is enabled.
        if (masterVariables.getDebug ()) {
            errorStream = System.err;
            outputStream = System.out;
        }
        return runApplication (
            command,
            errorStream,
            "ERROR (Npop CI)>",
            outputStream,
            "Npop CI>",
            true
        );
    }

    /**
     *  Runs the Demarcations Confidence Interval application.
     *
     *  @return The exit value.
     */
    public int runDemarcation() {
        String[] command = {
            binaryDirectory + "demarcation" + binaryExtension,
            workingDirectory + "demarcationIn.dat",
            workingDirectory + "demarcationOut.dat",
            Integer.toString (masterVariables.getNumberThreads ()),
            Boolean.toString(masterVariables.getDebug())
        };
        PrintStream errorStream = null;
        PrintStream outputStream = null;
        // Catch program output if debugging is enabled.
        if (masterVariables.getDebug ()) {
            errorStream = System.err;
            outputStream = System.out;
        }
        return runApplication (
            command,
            errorStream,
            "ERROR (Demarcation)>",
            outputStream,
            "Demarcation>",
            true
        );
    }

    /**
     *  Runs the Omega Confidence Interval application.
     *
     *  @return The exit value.
     */
    public int runOmegaCI() {
        String[] command = {
            binaryDirectory + "omegaCI" + binaryExtension,
            workingDirectory + "omegaIn.dat",
            workingDirectory + "omegaOut.dat",
            Integer.toString (masterVariables.getNumberThreads ()),
            Boolean.toString(masterVariables.getDebug())
        };
        PrintStream errorStream = null;
        PrintStream outputStream = null;
        // Catch program output if debugging is enabled.
        if (masterVariables.getDebug ()) {
            errorStream = System.err;
            outputStream = System.out;
        }
        return runApplication (
            command,
            errorStream,
            "ERROR (Omega CI)>",
            outputStream,
            "Omega CI>",
            true
        );
    }

    /**
     *  Runs the Sigma Confidence Interval application.
     *
     *  @return The exit value.
     */
    public int runSigmaCI() {
        String[] command = {
            binaryDirectory + "sigmaCI" + binaryExtension,
            workingDirectory + "sigmaIn.dat",
            workingDirectory + "sigmaOut.dat",
            Integer.toString (masterVariables.getNumberThreads ()),
            Boolean.toString(masterVariables.getDebug())
        };
        PrintStream errorStream = null;
        PrintStream outputStream = null;
        // Catch program output if debugging is enabled.
        if (masterVariables.getDebug ()) {
            errorStream = System.err;
            outputStream = System.out;
        }
        return runApplication (
            command,
            errorStream,
            "ERROR (Sigma CI)>",
            outputStream,
            "Sigma CI>",
            true
        );
    }

    /**
     *  Change the PCRError value in the pcrerror.dat file.
     *
     *  @param PCRError The new value to be assigned to pcrerror.
     */
    public void changePCRError(double PCRError) {
        // Create the random number seed; an odd less than 9 digits long
        long randValue = (long)(100000000 * Math.random());
        if (randValue % 2 == 0) {
            randValue ++;
        }
        try {
            BufferedWriter output = new BufferedWriter(new FileWriter(this.pcrerror));
            output.write("" + PCRError);
            output.newLine();
            output.write("" + randValue);
            output.newLine();
            output.close();
        }
        catch (IOException e) {
            e.printStackTrace();
        }
    }

    /**
     *  Run the Removegaps program.
     *
     *  @return The exit value.
     */
    private int runRemovegaps() {
        String[] command = {
            binaryDirectory + "removegaps" + binaryExtension,
            workingDirectory + "sequencesfasta.txt",
            workingDirectory + "numbers.dat",
            workingDirectory + "fasta.dat",
            Boolean.toString(masterVariables.getDebug())
        };
        PrintStream errorStream = null;
        PrintStream outputStream = null;
        // Catch program error output if debugging is enabled.
        if (masterVariables.getDebug ()) {
            errorStream = System.err;
            outputStream = System.out;
        }
        return runApplication (
            command,
            errorStream,
            "ERROR (Remove Gaps)>",
            outputStream,
            "Remove Gaps>",
            true
        );
    }

    /**
     *  Run the Readsynec program.
     *
     *  @return The exit value.
     */
    private int runReadsynec() {
        String[] command = {
            binaryDirectory + "readsynec" + binaryExtension,
            workingDirectory + "fasta.dat",
            workingDirectory + "population.dat",
            workingDirectory + "namesofstrains.dat",
            Boolean.toString(masterVariables.getDebug())
        };
        PrintStream errorStream = null;
        PrintStream outputStream = null;
        // Catch program error output if debugging is enabled.
        if (masterVariables.getDebug ()) {
            errorStream = System.err;
            outputStream = System.out;
        }
        return runApplication (
            command,
            errorStream,
            "ERROR (Read Synec)>",
            outputStream,
            "Read Synec>",
            true
        );
    }

    /**
     *  Run the Correctpcr program.
     *
     *  @return The exit value.
     */
    private int runCorrectpcr() {
        String[] command = {
            binaryDirectory + "correctpcr" + binaryExtension,
            workingDirectory + "population.dat",
            workingDirectory + "pcrerror.dat",
            workingDirectory + "correctpcr.out",
            Boolean.toString(masterVariables.getDebug())
        };
        PrintStream errorStream = null;
        PrintStream outputStream = null;
        // Catch program error output if debugging is enabled.
        if (masterVariables.getDebug ()) {
            errorStream = System.err;
            outputStream = System.out;
        }
        return runApplication (
            command,
            errorStream,
            "ERROR (Correct PCR)>",
            outputStream,
            "Correct PCR>",
            true
        );
    }

    /**
     *  Run the Divergencematrix program.
     *
     *  @return The exit value.
     */
    private int runDivergencematrix() {
        String[] command = {
            binaryDirectory + "divergencematrix" + binaryExtension,
            workingDirectory  + "correctpcr.out",
            workingDirectory + "identitymatrix.dat",
            Boolean.toString(masterVariables.getDebug())
        };
        PrintStream errorStream = null;
        PrintStream outputStream = null;
        // Catch program error output if debugging is enabled.
        if (masterVariables.getDebug ()) {
            errorStream = System.err;
            outputStream = System.out;
        }
        return runApplication (
            command,
            errorStream,
            "ERROR (Divergence Matrix)>",
            outputStream,
            "Divergence Matrix>",
            true
        );
    }

    /**
     *  Run the Binningdanny program.
     *
     *  @return The exit value.
     */
    private int runBinningdanny() {
        String[] command = {
            binaryDirectory + "binningdanny" + binaryExtension,
            workingDirectory + "identitymatrix.dat",
            workingDirectory + "binlevels.dat",
            workingDirectory + "output.dat",
            Boolean.toString(masterVariables.getDebug())
        };
        PrintStream errorStream = null;
        PrintStream outputStream = null;
        // Catch program error output if debugging is enabled.
        if (masterVariables.getDebug ()) {
            errorStream = System.err;
            outputStream = System.out;
        }
        return runApplication (
            command,
            errorStream,
            "ERROR (Binning)>",
            outputStream,
            "Binning>",
            true
        );
    }

    /**
     *  Run a Phylip program.
     *
     *  @param program The Phylip program to run.
     *  @param arguments The arguments to the Phylip program.
     *  @return The exit value.
     */
    private int runPhylip(String program, String arguments) {
        File inputFile = new File (workingDirectory + "input");
        BufferedWriter writer = null;
        try {
            writer = new BufferedWriter(new FileWriter(inputFile));
            writer.write(arguments);
        }
        catch (IOException e) {
            e.printStackTrace();
        }
        finally {
            if (writer != null) {
                try {
                    writer.close();
                }
                catch (IOException e) {
                    System.out.println(
                        "Error closing the input file for Phylip."
                    );
                }
            }
        }
        String[] command = {
            phylipScript,
            options.getHelperProgram(program),
            workingDirectory,
            Boolean.toString(masterVariables.getDebug())
        };
        PrintStream errorStream = null;
        PrintStream outputStream = null;
        // Catch program error output if debugging is enabled.
        if (masterVariables.getDebug ()) {
            errorStream = System.err;
            outputStream = System.out;
        }
        return runApplication (
            command,
            errorStream,
            "ERROR (Phylip " + program + ")>",
            outputStream,
            "Phylip " + program + ">",
            true
        );
    }


    /**
     *  Runs the provided application with the provided args.
     *  If the wait boolean is set, waits for the application to finish.
     *
     *  @param command A String array containing the path and filename of the
     *  application, and any arguments.
     *  @param errorStream The IO Stream to print error messages to.
     *  @param errorMessage The title for error messages.
     *  @param outputStream The IO Stream to print standard messages to.
     *  @param outputMessage The title for the output messages.
     *  @param wait Set to TRUE to wait for application to exit.
     *  @return The exit value.
     */
    private int runApplication (
        String[] command, PrintStream errorStream, String errorMessage,
        PrintStream outputStream, String outputMessage, boolean wait
    ) {
        int exitVal = -1;
        try {
            ProcessBuilder pb = new ProcessBuilder (command);
            Process p = pb.start ();
            StreamGobbler errorGobbler = null;
            StreamGobbler outputGobbler = null;
            // Display debugging output if needed.
            if (masterVariables.getDebug ()) {
                System.out.print ("Execute:");
                for (int i = 0; i < command.length; i ++) {
                    System.out.print (" " + command[i]);
                }
                System.out.print ("\n");
            }
            if (errorStream != null) {
                // Grab error messages.
                errorGobbler = new StreamGobbler (
                    p.getErrorStream (),
                    errorStream,
                    errorMessage
                );
                errorGobbler.start ();
            }
            if (outputStream != null) {
                // Grab output messages.
                outputGobbler = new StreamGobbler (
                    p.getInputStream (),
                    outputStream,
                    outputMessage
                );
                outputGobbler.start ();
            }
            // Wait for application to finish if needed.
            if (wait) {
                exitVal = p.waitFor ();
                // Also wait for the StreamGobbler threads to finish.
                if (errorGobbler != null) errorGobbler.join ();
                if (outputGobbler != null) outputGobbler.join ();
            }
            else {
                exitVal = 0;
            }
        }
        catch (IOException e) {
            e.printStackTrace ();
        }
        catch (InterruptedException e) {
            e.printStackTrace ();
        }
        return exitVal;
    }

   /**
    *  Check that the binlevels.dat file exists, if not create a default one.
    */
    private void checkBinLevels() {
        File binLevels = new File(workingDirectory + "binlevels.dat");
        if (! binLevels.exists()) {
            try {
                BufferedWriter output = new BufferedWriter(new FileWriter(workingDirectory + "binlevels.dat"));
                output.write(Integer.toString(defaultBinLevels.length));
                output.newLine();
                for (int i = 0; i < defaultBinLevels.length; i++) {
                    output.write(Double.toString(defaultBinLevels[i]));
                    output.newLine();
                }
                output.close();
            }
            catch (IOException e) {
                e.printStackTrace();
            }
        }
    }

    /**
     *  Generates a random number for the pcrerror file each time binning is run.
     */
    private void makeRandPCRError() {
        double PCRError = defaultPCRError;
        String PCRErrorIn = "";
        if (pcrerror.exists()) {
            try {
                BufferedReader input = new BufferedReader(new FileReader(this.pcrerror));
                PCRErrorIn = input.readLine();
                input.close();
            }
            catch (IOException e) {
                e.printStackTrace();
            }
            try {
                PCRError = (new Double(PCRErrorIn)).doubleValue();
            }
            catch (NumberFormatException e) {
                PCRError = defaultPCRError;
            }
        }
        changePCRError(PCRError);
    }

    /**
     *  Write the script to run Phylip.
     *
     *  @param script The Phylip script to write.
     */
    private void writePhylipScript(String script) {
        File scriptFile = new File(phylipScript);
        BufferedWriter writer = null;
        try {
            writer = new BufferedWriter(new FileWriter(scriptFile));
            writer.write(script);
        }
        catch (IOException e) {
            System.out.println("Error writing the script for Phylip.");
        }
        finally {
            if (writer != null) {
                try {
                    writer.close();
                }
                catch (IOException e) {
                    System.out.println("Error closing the script for Phylip");
                }
            }
        }
        scriptFile.setExecutable(true);
    }

    private MasterVariables masterVariables;
    private Options options;
    private JTextArea log;
    private File pcrerror;
    private String phylipScript;
    private String binaryExtension;
    private String binaryDirectory;
    private String workingDirectory;
    private double defaultPCRError = 8.33e-6;
    private double defaultBinLevels[] = {
        0.650, 0.700, 0.750, 0.800, 0.850, 0.9000, 0.9500, 0.9600, 0.9700, 0.9800, 0.9900, 0.9950, 1.0000
    };
    private static final String phylipBatchScript =
        "@echo off" +
        System.getProperty ("line.separator") +
        "cd /D \"%2\"" +
        System.getProperty ("line.separator") +
        "if exist outfile del outfile" +
        System.getProperty ("line.separator") +
        "if exist outtree del outtree" +
        System.getProperty ("line.separator") +
        "if \"%3\"==\"true\" (" +
        System.getProperty ("line.separator") +
        "  type input | \"%1\"" +
        System.getProperty ("line.separator") +
        ") else (" +
        System.getProperty ("line.separator") +
        "  type input | \"%1\" > screenout" +
        System.getProperty ("line.separator") +
        "  del screenout" +
        System.getProperty ("line.separator") +
        ")" +
        System.getProperty ("line.separator");

    private static final String phylipShellScript =
        "#!/bin/sh" +
        System.getProperty ("line.separator") +
        "cd \"$2\"" +
        System.getProperty ("line.separator") +
        "rm -f outfile outtree" +
        System.getProperty ("line.separator") +
        "COMMAND=\"$1\"" +
        System.getProperty ("line.separator") +
        "if [ -f /etc/debian_version ]; then" +
        System.getProperty ("line.separator") +
        "  COMMAND=\"phylip $1\"" +
        System.getProperty ("line.separator") +
        "fi" +
        System.getProperty ("line.separator") +
        "if [ \"$3\" = \"true\" ]; then" +
        System.getProperty ("line.separator") +
        "  cat input | $COMMAND" +
        System.getProperty ("line.separator") +
        "else" +
        System.getProperty ("line.separator") +
        "  cat input | $COMMAND > /dev/null" +
        System.getProperty ("line.separator") +
        "fi" +
        System.getProperty ("line.separator");

}
