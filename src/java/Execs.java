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
        binaryDirectory = options.getDirectoryOption(Options.BINARY_DIRECTORY);
        scriptDirectory = options.getDirectoryOption(Options.SCRIPT_DIRECTORY);
        workingDirectory = options.getDirectoryOption(Options.WORKING_DIRECTORY);
        // Setup the rest of the variables.
        log = masterVariables.getLog();
        pcrerror = new File (workingDirectory + "pcrerror.dat");
        // Check which OS we are running on.
        if (osName.contains("windows")) {
            binaryExtension = ".exe";
            scriptExtension = ".bat";
        }
        else if (osName.contains("linux")) {
            scriptExtension = ".sh";
            if (osArch.contains("i386")) {
                binaryExtension = ".i386";
            }
            else if (osArch.contains("amd64")) {
                binaryExtension = ".amd64";
            }
            else {
                log.append("Unsupported Linux architecture, contact the developers.\n");
                log.append("Architecture detected: " + osArch + "\n");
            }
        }
        else if (osName.contains("mac")) {
            binaryExtension = ".app";
            scriptExtension = ".sh";
        }
        else {
            log.append("Unsupported OS, contact the developers.\n");
            log.append("OS detected: " + osName + "\n");
            log.append("Architecture detected: " + osArch + "\n");
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
        return runApplication("Brute Force", command, true);
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
        return runApplication("Hill Climb", command, true);
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
        return runApplication("NJPlot", command, false);
    }

    /**
     *  Runs the dnapars application
     *
     *  @return The exit value.
     */
    public int runDNAPars() {
        return runPhylip(Options.DNAPARS);
    }

    /**
     *  Runs the dnadist application.
     *
     *  @return The exit value.
     */
    public int runDNADist() {
        return runPhylip(Options.DNADIST);
    }

    /**
     *  Runs the Neighbor application.
     *
     *  @return The exit value.
     */
    public int runNJ() {
        return runPhylip(Options.NEIGHBOR);
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
        return runApplication("Drift CI", command, true);
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
        return runApplication("Npop CI", command, true);
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
        return runApplication("Demarcation", command, true);
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
        return runApplication("Omega CI", command, true);
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
        return runApplication("Sigma CI", command, true);
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
        return runApplication("Remove Gaps", command, true);
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
        return runApplication("Read Synec", command, true);
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
        return runApplication("Correct PCR", command, true);
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
        return runApplication("Divergence Matrix", command, true);
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
        return runApplication("Binning", command, true);
    }

    /**
     *  Run a Phylip program.
     *
     *  @param program The Phylip program to run.
     *  @return The exit value.
     */
    private int runPhylip(String program) {
        String[] command = {
            scriptDirectory + "phylip" + scriptExtension,
            options.getHelperProgram(program),
            workingDirectory
        };
        return runApplication("Phylip", command, true);
    }

    /**
     *  Runs the provided application with the provided args.
     *  If the wait boolean is set, waits for the application to finish.
     *
     *  @param name The name  of the application.
     *  @param command A String array containing the path and filename of the application, and any arguments.
     *  @param wait Set to TRUE to wait for application to exit.
     *  @return The exit value.
     */
    private int runApplication(String name, String[] command, boolean wait) {
        int exitVal = -1;
        try {
            Process p = new ProcessBuilder(command).start();
            StreamGobbler errorGobbler = null;
            StreamGobbler outputGobbler = null;
            // Display debugging output if needed.
            if (masterVariables.getDebug()) {
                log.append("Execute: " + command[0] + "\n");
                // Grab error messages.
                errorGobbler = new StreamGobbler(p.getErrorStream(), "ERROR (" + name + ")");
                errorGobbler.start();
                // Grab output messages.
                outputGobbler = new StreamGobbler(p.getInputStream(), name);
                outputGobbler.start();
            }
            // Wait for application to finish if needed.
            if (wait) {
                exitVal = p.waitFor();
                // Also wait for the StreamGobbler threads to finish.
                if (masterVariables.getDebug()) {
                    errorGobbler.join();
                    outputGobbler.join();
                }
            }
            else {
                exitVal = 0;
            }
        }
        catch (IOException e) {
            e.printStackTrace();
        }
        catch (InterruptedException e) {
            e.printStackTrace();
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

    private MasterVariables masterVariables;
    private Options options;
    private JTextArea log;
    private File pcrerror;
    private String binaryExtension;
    private String scriptExtension;
    private String binaryDirectory;
    private String scriptDirectory;
    private String workingDirectory;
    private double defaultPCRError = 8.33e-6;
    private double defaultBinLevels[] = {
        0.650, 0.700, 0.750, 0.800, 0.850, 0.9000, 0.9500, 0.9600, 0.9700, 0.9800, 0.9900, 0.9950, 1.0000
    };
}
