/*
 *    Ecotype Simulation models the sequence diversity within a bacterial
 *    clade as the evolutionary result of net ecotype formation and periodic
 *    selection, yielding a certain number of ecotypes.
 *
 *    Copyright (C) 2009       Andrew Warner, Wesleyan University
 *    Copyright (C) 2009-2013  Jason M. Wood, Montana State University
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
     *  @param masterVariables The MasterVariables.
     */
    public Execs (MasterVariables masterVariables) {
        this.masterVariables = masterVariables;
        // Grab the running environment.
        String osName = System.getProperty ("os.name").toLowerCase ();
        String osArch = System.getProperty ("os.arch").toLowerCase ();
        String osVersion = System.getProperty ("os.version");
        // Setup the rest of the variables.
        binaryDirectory = masterVariables.getBinaryDirectory ();
        scriptDirectory = masterVariables.getScriptDirectory ();
        workingDirectory = masterVariables.getWorkingDirectory ();
        log = masterVariables.getLog ();
        // Check which OS we are running on.
        if (osName.contains ("windows")) {
            binaryExtension = ".exe";
            // Use the Phylip batch script.
            phylipScript = workingDirectory + "phylip.bat";
            writePhylipScript (phylipBatchScript);
        }
        else if (osName.contains ("linux")) {
            if (osArch.contains ("i386")) {
                binaryExtension = ".i386";
            }
            else if (osArch.contains ("amd64")) {
                binaryExtension = ".amd64";
            }
            else {
                log.append (
                    "Unsupported architecture, contact the developers.\n" +
                    "Architecture detected: " + osArch + "\n"
                );
            }
            // Use the Phylip shell script.
            phylipScript = workingDirectory + "phylip.sh";
            writePhylipScript (phylipShellScript);
        }
        else if (osName.contains ("mac")) {
            binaryExtension = ".app";
            // Use the Phylip shell script.
            phylipScript = workingDirectory + "phylip.sh";
            writePhylipScript (phylipShellScript);
        }
        else {
            log.append ("Unsupported OS, contact the developers.\n");
            log.append ("OS detected: " + osName + "\n");
            log.append ("Architecture detected: " + osArch + "\n");
        }
    }

    /**
     *  Runs the initial Fred method to find a good value with which to do
     *  hillclimbing.
     *
     *  @param input The bruteforce input file.
     *  @param output The bruteforce output file.
     *  @return The exit value, -1 if there was an error.
     */
    public int runBruteforce (File input, File output) {
        String[] command = {
            binaryDirectory + "bruteforce" + binaryExtension,
            input.getAbsolutePath (),
            output.getAbsolutePath (),
            Integer.toString (masterVariables.getNumberThreads ()),
            Boolean.toString (masterVariables.getDebug ())
        };
        return runApplication ("Brute Force Search", command, true);
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
        return runApplication ("Hill Climb", command, true);
    }

    /**
     *  Opens the tree using NJPlot so that the user can see it to start
     *  analysis.
     *
     *  @param treeFile The tree file to open in NJPlot.
     *  @return The exit value.
     */
    public int openTree (File treeFile) {
        String[] command = {
            masterVariables.getProgramNJPlot (),
            treeFile.getPath ()
        };
        return runApplication ("NJPlot", command, false);
    }

    /**
     *  Runs the dnapars application
     *
     *  @return The exit value.
     */
    public int runDNAPars () {
        return runPhylip (masterVariables.getProgramDNAPars (), "V\n1\nY\n");
    }

    /**
     *  Runs the dnadist application.
     *
     *  @return The exit value.
     */
    public int runDNADist () {
        return runPhylip (masterVariables.getProgramDNADist (), "Y\n");
    }

    /**
     *  Runs the Neighbor application.
     *
     *  @return The exit value.
     */
    public int runNJ () {
        return runPhylip (masterVariables.getProgramNeighbor (), "Y\n");
    }

    /**
     *  Runs the Retree application.
     *
     *  @param nu The number of environmental sequences including the
     *      outgroup.
     *  @return The exit value.
     */
    public int runRetree (int nu) {
        String arguments = String.format (
            "Y\nO\n%d\nW\nR\nQ\n",
            nu
        );
        return runPhylip (masterVariables.getProgramRetree (), arguments);
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
        return runApplication ("Npop CI", command, true);
    }

    /**
     *  Runs the Demarcation Confidence Interval application.
     *
     *  @param input The demarcationCI input file.
     *  @param output The demarcationCI output file.
     *  @return The exit value.
     */
    public int runDemarcationCI (File input, File output) {
        String[] command = {
            binaryDirectory + "demarcationCI" + binaryExtension,
            input.getAbsolutePath (),
            output.getAbsolutePath (),
            Integer.toString (masterVariables.getNumberThreads ()),
            Boolean.toString (masterVariables.getDebug ())
        };
        return runApplication ("Demarcation CI", command, true);
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
        return runApplication ("Omega CI", command, true);
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
        return runApplication ("Sigma CI", command, true);
    }

    /**
     *  Run the Removegaps program.
     *
     *  @param fasta The fasta formated input file.
     *  @param numbers The numbers file
     *  @param removegapsOut The removegaps output file.
     *  @return The exit value.
     */
    public int runRemovegaps (File fasta, File numbers, File removegapsOut) {
        String[] command = {
            binaryDirectory + "removegaps" + binaryExtension,
            fasta.getAbsolutePath (),
            numbers.getAbsolutePath (),
            removegapsOut.getAbsolutePath (),
            Boolean.toString (masterVariables.getDebug ())
        };
        return runApplication ("Remove Gaps", command, true);
    }

    /**
     *  Run the Readsynec program.
     *
     *  @param removegapsOut The removegaps file.
     *  @param population The population file.
     *  @param nameofstrains The nameofstrains file.
     *  @return The exit value.
     */
    public int runReadsynec (File removegapsOut, File population,
        File nameofstrains) {
        String[] command = {
            binaryDirectory + "readsynec" + binaryExtension,
            removegapsOut.getAbsolutePath (),
            population.getAbsolutePath (),
            nameofstrains.getAbsolutePath (),
            Boolean.toString (masterVariables.getDebug ())
        };
        return runApplication ("Read Synec", command, true);
    }

    /**
     *  Run the Correctpcr program.
     *
     *  @param population The population file.
     *  @param pcrerror The pcrerror file.
     *  @param correctpcr The correctpcr output file.
     *  @return The exit value.
     */
    public int runCorrectpcr (File population, File pcrerror,
        File correctpcr) {
        String[] command = {
            binaryDirectory + "correctpcr" + binaryExtension,
            population.getAbsolutePath (),
            pcrerror.getAbsolutePath (),
            correctpcr.getAbsolutePath (),
            Boolean.toString (masterVariables.getDebug ())
        };
        return runApplication ("Correct PCR", command, true);
    }

    /**
     *  Run the Divergencematrix program.
     *
     *  @param correctpcr The correctpcr file.
     *  @param divergencematrix The divergence matrix output file.
     *  @return The exit value.
     */
    public int runDivergencematrix (File correctpcr, File divergencematrix) {
        String[] command = {
            binaryDirectory + "divergencematrix" + binaryExtension,
            correctpcr.getAbsolutePath (),
            divergencematrix.getAbsolutePath (),
            Boolean.toString (masterVariables.getDebug ())
        };
        return runApplication ("Divergence Matrix", command, true);
    }

    /**
     *  Run the Binningdanny program.
     *
     *  @param divergencematrix The divergence matrix file.
     *  @param binlevels The binlevels file.
     *  @param binningdanny The binningdanny file.
     *  @return The exit value.
     */
    public int runBinningdanny (File divergencematrix, File binlevels,
        File binningdanny) {
        String[] command = {
            binaryDirectory + "binningdanny" + binaryExtension,
            divergencematrix.getAbsolutePath (),
            binlevels.getAbsolutePath (),
            binningdanny.getAbsolutePath (),
            Boolean.toString (masterVariables.getDebug ())
        };
        return runApplication ("Binning", command, true);
    }

    /**
     *  Run a Phylip program.
     *
     *  @param program The Phylip program to run.
     *  @param arguments The arguments to the Phylip program.
     *  @return The exit value.
     */
    private int runPhylip (String program, String arguments) {
        File inputFile = new File (workingDirectory + "input");
        BufferedWriter writer = null;
        try {
            writer = new BufferedWriter (new FileWriter (inputFile));
            writer.write (arguments);
        }
        catch (IOException e) {
            e.printStackTrace ();
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
        String[] command = {
            phylipScript,
            program,
            workingDirectory,
            Boolean.toString (masterVariables.getDebug ())
        };
        return runApplication (program, command, true);
    }

    /**
     *  Runs the provided application with the provided args.
     *  If the wait boolean is set, waits for the application to finish.
     *
     *  @param name The name of the application to run.
     *  @param command A String array containing the path and filename of the
     *  application, and any arguments.
     *  @param wait Set to TRUE to wait for application to exit.
     *  @return The exit value.
     */
    private int runApplication (String name, String[] command, boolean wait) {
        int exitVal = -1;
        try {
            Process p = new ProcessBuilder (command).start ();
            StreamGobbler errorGobbler = null;
            StreamGobbler outputGobbler = null;
            // Display debugging output if needed.
            if (masterVariables.getDebug ()) {
                System.out.print ("Execute:");
                for (int i = 0; i < command.length; i ++) {
                    System.out.print (" " + command[i]);
                }
                System.out.print ("\n");
                // Grab error messages.
                errorGobbler = new StreamGobbler (
                    p.getErrorStream (),
                    "ERROR (" + name + ")"
                );
                errorGobbler.start ();
                // Grab output messages.
                outputGobbler = new StreamGobbler (
                    p.getInputStream (), name
                );
                outputGobbler.start ();
            }
            // Wait for application to finish if needed.
            if (wait) {
                exitVal = p.waitFor ();
                // Also wait for the StreamGobbler threads to finish.
                if (masterVariables.getDebug ()) {
                    errorGobbler.join ();
                    outputGobbler.join ();
                }
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
     *  Write the script to run Phylip.
     *
     *  @param script The Phylip script to write.
     */
    private void writePhylipScript (String script) {
        File scriptFile = new File (phylipScript);
        BufferedWriter writer = null;
        try {
            writer = new BufferedWriter (new FileWriter (scriptFile));
            writer.write (script);
        }
        catch (IOException e) {
            System.out.println ("Error writing the script for Phylip.");
        }
        finally {
            if (writer != null) {
                try {
                    writer.close();
                }
                catch (IOException e) {
                    System.out.println (
                        "Error closing the script for Phylip"
                    );
                }
            }
        }
        scriptFile.setExecutable(true);
    }

    private MasterVariables masterVariables;
    private Logger log;
    private String phylipScript;
    private String binaryExtension;
    private String binaryDirectory;
    private String scriptDirectory;
    private String workingDirectory;

    private static final String phylipBatchScript =
        "@echo off\n" +
        "cd \"%2\"\n" +
        "if exist outfile del outfile\n" +
        "if exist outtree del outtree\n" +
        "if \"%3\"==\"true\" (\n" +
        "  type input | \"%1\"\n" +
        ") else (\n" +
        "  type input | \"%1\" > screenout\n" +
        "  del screenout\n" +
        ")\n";

    private static final String phylipShellScript =
        "#!/bin/sh\n" +
        "cd \"$2\"\n" +
        "rm -f outfile outtree\n" +
        "COMMAND=\"$1\"\n" +
        "if [ -f /etc/debian_version ]; then\n" +
        "  COMMAND=\"phylip $1\"\n" +
        "fi\n" +
        "if [ \"$3\" = \"true\" ]; then\n" +
        "  cat input | $COMMAND\n" +
        "else\n" +
        "  cat input | $COMMAND > /dev/null\n" +
        "fi\n";

}
