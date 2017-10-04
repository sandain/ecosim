/*
 *    Ecotype Simulation models the sequence diversity within a bacterial
 *    clade as the evolutionary result of net ecotype formation and periodic
 *    selection, yielding a certain number of ecotypes.
 *
 *    Copyright (C) 2009       Andrew Warner, Wesleyan University
 *    Copyright (C) 2009-2016  Jason M. Wood, Montana State University
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
import java.io.PrintStream;
import java.io.IOException;

/**
 *  Holds the executable methods for the cohan program.
 *
 *  @author Andrew Warner
 *  @author Jason M. Wood
 *  @copyright GNU General Public License
 */
public class Execs {

    /**
     *  Detects the Operating System that this class is running on,
     *  and executes the native fortran applications.
     *
     *  @param log The Logger object.
     *  @param masterVariables The MasterVariables.
     */
    public Execs (Logger log, MasterVariables masterVariables) {
        this.log = log;
        this.masterVariables = masterVariables;
        // Grab the running environment.
        String osName = System.getProperty ("os.name").toLowerCase ();
        String osArch = System.getProperty ("os.arch").toLowerCase ();
        String osVersion = System.getProperty ("os.version");
        // Setup the rest of the variables.
        binaryDirectory = masterVariables.getBinaryDirectory ();
        // Check which OS we are running on.
        if (osName.contains ("windows")) {
            binaryExtension = ".exe";
        }
        else if (osName.contains ("linux")) {
            binaryExtension = "";
        }
        else if (osName.contains ("mac")) {
            binaryExtension = ".app";
        }
        else {
            log.append ("Unsupported OS, contact the developers.\n");
            log.append ("OS detected: " + osName + "\n");
            log.append ("Architecture detected: " + osArch + "\n");
        }
    }

    /**
     *  Runs the hillclimb program.
     *
     *  @param input The hillclimb input file.
     *  @param output The hillclimb output file.
     *  @return The exit value.
     */
    public int runHillclimb (File input, File output) {
        String[] command = {
            binaryDirectory + "hillclimb" + binaryExtension,
            input.getAbsolutePath (),
            output.getAbsolutePath (),
            Integer.toString (masterVariables.getNumberThreads ()),
            Boolean.toString (masterVariables.getDebug ())
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
     *  Runs the Npop Confidence Interval application.
     *
     *  @param input The npopCI input file.
     *  @param output The npopCI output file.
     *  @return The exit value.
     */
    public int runNpopCI (File input, File output) {
        String[] command = {
            binaryDirectory + "npopCI" + binaryExtension,
            input.getAbsolutePath (),
            output.getAbsolutePath (),
            Integer.toString (masterVariables.getNumberThreads ()),
            Boolean.toString (masterVariables.getDebug ())
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
     *  Runs the Demarcation application.
     *
     *  @param input The demarcation input file.
     *  @param output The demarcation output file.
     *  @return The exit value.
     */
    public int runDemarcation (File input, File output) {
        String[] command = {
            binaryDirectory + "demarcation" + binaryExtension,
            input.getAbsolutePath (),
            output.getAbsolutePath (),
            Integer.toString (masterVariables.getNumberThreads ()),
            Boolean.toString (masterVariables.getDebug ())
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
     *  @param input The omegaCI input file.
     *  @param output The omegaCI output file.
     *  @return The exit value.
     */
    public int runOmegaCI (File input, File output) {
        String[] command = {
            binaryDirectory + "omegaCI" + binaryExtension,
            input.getAbsolutePath (),
            output.getAbsolutePath (),
            Integer.toString (masterVariables.getNumberThreads ()),
            Boolean.toString (masterVariables.getDebug ())
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
     *  @param input The sigmaCI input file.
     *  @param output The sigmaCI output file.
     *  @return The exit value.
     */
    public int runSigmaCI (File input, File output) {
        String[] command = {
            binaryDirectory + "sigmaCI" + binaryExtension,
            input.getAbsolutePath (),
            output.getAbsolutePath (),
            Integer.toString (masterVariables.getNumberThreads ()),
            Boolean.toString (masterVariables.getDebug ())
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
     *  Runs FastTree on the Fasta formated input file to generate a Newick
     *  formated output file.
     *
     *  @param input The Fasta formated input file.
     *  @param output The Newick formated output file.
     *  @return The exit value.
     */
    public int runFastTree (File input, File output) {
        String[] command = {
            binaryDirectory + "fasttree" + binaryExtension,
            "-nt",
            input.getAbsolutePath (),
        };
        PrintStream errorStream = null;
        PrintStream outputStream = null;
        // Catch program error output if debugging is enabled.
        if (masterVariables.getDebug ()) {
            errorStream = System.err;
        }
        // Catch program output in the output file.
        try {
            outputStream = new PrintStream (output);
        }
        catch (IOException e) {
            e.printStackTrace ();
        }
        return runApplication (
            command,
            errorStream,
            "FastTree>",
            outputStream,
            "",
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

    private MasterVariables masterVariables;
    private Logger log;
    private String binaryExtension;
    private String binaryDirectory;

}
