/*
 *    Ecotype Simulation models the sequence diversity within a bacterial clade
 *    as the evolutionary result of net ecotype formation and periodic
 *    selection, yielding a certain number of ecotypes.
 *
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

import java.io.File;

/**
 * @mainpage Ecotype %Simulation
 *
 * Ecotype %Simulation models the sequence diversity within a microbial clade
 * as the evolutionary result of net ecotype formation (omega) and periodic
 * selection (sigma), yielding a number of putative ecotypes (npop).  Ecotype
 * %Simulation also demarcates sequences belonging to each putative ecotype.
 *
 * @par
 * We define ecotype here as a phylogenetic group of closely related organisms
 * that are ecologically interchangeable, in that members of the same ecotype
 * share genetic adaptations to a particular set of habitats, resources, and
 * conditions, while simultaneously remaining ecologically distinct from each
 * other.
 *
 * @tableofcontents
 *
 * @section usage Usage
 * @code
 * java -jar EcoSim.jar [OPTIONS] /path/to/sequence_fasta /path/to/sequence_tree output_xml
 *
 * Sequences should be aligned and in a fasta formated file, with the outgroup
 * listed first.
 *
 * The newick formated tree file must include the same leaf node names and
 * number as the sequences in the fasta file.
 *
 * Output is saved in XML format.
 *
 * Options:
 *   -d, --debug        : Display debugging output.
 *   -h, --help         : Display helpful information.
 *   -n, --nogui        : Hide the default GUI.  Implies --runall.
 *   -r, --runall       : Run everything, including demarcation.
 *   -t=n, --threads=n  : Set the number of threads (n) to start, default
 *                        to system maximum.
 *   -v, --version      : Display the version number.
 * @endcode
 *
 * @section fortran Fortran
 * @subsection fortran-programs Programs
 * @li @b ::readsynec - Verifies sequence data input.
 * @li @b ::removegaps - Removes gaps in sequence data.
 * @li @b ::correctpcr - Corrects the sequence data for PCR error.
 * @li @b ::divergencematrix - Calculates the divergence matrix.
 * @li @b ::binningdanny - Bins the sequence data using the divergence matrix.
 * @li @b ::bruteforce - Brute force search of the omega/sigma/npop space.
 * @li @b ::hillclimb - Performs a hillclimbing on the bruteforce result.
 * @li @b ::omegaci - Calculates the confidence interval of the omega value.
 * @li @b ::sigmaci - Calculates the confidence interval of the sigma value.
 * @li @b ::npopci - Calculates the confidence interval of the npop value.
 * @li @b ::demarcationci - Calculates the lower bound of the npop CI.
 *
 * @subsection fortran-libraries Libraries
 * @li @b ::darray - Stores dynamic arrays.
 * @li @b ::methods - Common methods of the simulation.
 * @li @b ::simplexmethod - The Nelder-Mead Simplex Method.
 * @li @b ::ziggurat - The Ziggurat Random Number Generator.
 *
 * @section java Java
 * @li @b EcotypeSimulation - The main object of the Ecotype %Simulation.
 * @li @b BinLevel - Stores the bin levels for the Binning object.
 * @li @b Binning - Object to interact with the binning programs.
 * @li @b Bruteforce - Object to interact with the bruteforce search program.
 * @li @b Demarcation - Demarcates ecotypes based on the hillclimbing values
 *        and the phylogeny of the sequences.
 * @li @b DemarcationConfidenceInterval - Run the demarcation confidence
 *        interval program.
 * @li @b DivergenceMatrix - Object to interact with the divergencematrix
 *        program.
 * @li @b Execs - Holds the executable methods for the cohan program.
 * @li @b Fasta - Handles the input and output of fasta formatted text files.
 * @li @b FileChooser - Defines a custom file chooser.
 * @li @b Heapsorter - Runs the heapsort on a given set of data.
 * @li @b HelpAboutWindow - Displays a Help/About GUI window.
 * @li @b Hillclimb - Object to interact with the hillclimbing program.
 * @li @b InvalidNewickException - Report a malformed Newick tree.
 * @li @b Likelihood - Stores likelihood values.
 * @li @b Logger - Display text to the user.
 * @li @b MasterVariables - Common variables used through the program.
 * @li @b NewickTree - Interact with newick phylogenetic trees.
 * @li @b NewickTreeNode - A node of a newick phylogenetic tree.
 * @li @b NpopConfidenceInterval - Run the npop confidence interval program.
 * @li @b OmegaConfidenceInterval - Run the omega confidence interval program.
 * @li @b ParameterSet - An object to store the parameter values.
 * @li @b Phylogeny - Object to interact with the phylogeny programs.
 * @li @b ProjectFileIO - Perform IO operations for the XML project file.
 * @li @b SigmaConfidenceInterval - Run the sigma confidence interval program.
 * @li @b Simulation - The shared methods of the simulation.
 * @li @b SimulationCLI - The Ecotype %Simulation CLI.
 * @li @b SimulationGUI - The Ecotype %Simulation GUI.
 * @li @b StreamGobbler - Captures output from the Fortran programs.
 */

/**
 *  The main class for instantiating one of the Ecotype %Simulation
 *  interfaces provided by SimulationCLI and SimulationGUI.
 *
 *  @author Jason M. Wood
 *  @copyright GNU General Public License
 */
public class EcotypeSimulation implements Runnable {

    /**
     *  EcotypeSimulation constructor.
     *
     *  @param args The command line arguments.
     */
    public EcotypeSimulation(String[] args) {
        // Initialize the master variables.
        masterVariables = new MasterVariables();
        execs = masterVariables.getExecs();
        log = masterVariables.getLog();
        // Initialize variables.
        noGUI = false;
        runAll = false;
        // Check for command line arguments.
        checkArguments(args);
    }

    /**
     *  Run Ecotype Simulation.
     */
    public void run() {
        // Start one of the Ecotype Simulation interfaces.
        Simulation simulation;
        if (noGUI) {
            // Start the command line interface (CLI).
            simulation = new SimulationCLI(
                masterVariables, fastaFile, newickFile
            );
        }
        else {
            // Start the graphical user interface (GUI).
            simulation = new SimulationGUI(
                masterVariables, fastaFile, newickFile
            );
        }
        // Run the phylogeny program if the fasta and newick files are loaded.
        if (fastaFile != null && fastaFile.exists() && 
            newickFile != null && newickFile.exists()) {
            simulation.runPhylogeny();
        }
        // Launch NJPlot to view the tree.
        if (! noGUI && newickFile != null) {
            log.append("Displaying the tree with NJplot.\n\n");
            execs.openTree(newickFile);
        }
        // If the runAll flag has been set, run the simulation.
        if (runAll) {
            simulation.runBinning();
            simulation.runBruteforce();
            simulation.runHillclimbing();
            simulation.runNpopConfidenceInterval();
            simulation.runOmegaConfidenceInterval();
            simulation.runSigmaConfidenceInterval();
            simulation.runDemarcation();
            if (outputFile != null) {
                simulation.saveProjectFile(outputFile);
            }
        }
    }

    /**
     *  Start an instance of EcotypeSimulation.
     *
     *  @param args The command line arguments.
     */
    public static void main(final String[] args) {
        EcotypeSimulation es = new EcotypeSimulation(args);
        new Thread(es).start();
    }

    /**
     *  Check for command line arguments.
     *
     *  @param args The command line arguments.
     */
    private void checkArguments (String[] args) {
        boolean error = false;
        for (int i = 0; i < args.length; i ++) {
            String[] keyValue = args[i].split("=");
            String key = keyValue[0];
            String value = "";
            // Look for key=value pairs in the arguments.
            if (keyValue.length > 1) {
                value = keyValue[1];
            }
            // Catch a space after the equal sign.
            else if (args.length > i + 1 && args[i].endsWith("=")) {
                value = args[i + 1];
                i ++;
            }
            // Catch a space before the equal sign.
            else if (args.length > i + 1 && args[i + 1].startsWith("=") && ! args[i + 1].equals("=")) {
                value = args[i + 1].substring(1);
                i ++;
            }
            // Catch spaces surrounding the equal sign.
            else if (args.length > i + 2 && args[i + 1].equals("=")) {
                value = args[i + 2];
                i += 2;
            }
            // Check for various arguments.
            switch(key) {
                case "-d":
                case "--debug":
                    masterVariables.setDebug(true);
                    break;
                case "-h":
                case "--help":
                    System.out.println(String.format(
                        "Ecotype Simulation %s\n\n%s",
                        masterVariables.getVersion(), usage
                    ));
                    noGUI = true;
                    break;
                case "-n":
                case "--nogui":
                    noGUI = true;
                    runAll = true;
                    break;
                case "-r":
                case "--runall":
                    runAll = true;
                    break;
                case "-t":
                case "--threads":
                    if (value.length() > 0) {
                        int procs = Runtime.getRuntime().availableProcessors();
                        int num = 0;
                        try {
                            num = new Integer(value).intValue();
                        }
                        catch (NumberFormatException e) {
                            System.out.println(String.format(
                                "Syntax error: Expected a number.\n%s\n%s",
                                e, usage
                            ));
                            System.exit(1);
                        }
                        if (num >= 1 && num <= procs) {
                            masterVariables.setNumberThreads(num);
                        }
                        else {
                            System.out.println(String.format(
                                "Syntax error: Invalid number of threads " + 
                                "specified: %d of %d possible.\n",
                                num, procs
                            ));
                            System.exit(1);
                        }
                    }
                    else {
                        // Number of threads not provided, print an error.
                        System.out.print(String.format(
                            "Syntax error: Number of threads missing.\n%s",
                            usage
                        ));
                        System.exit(1);
                    }
                    break;
                case "-v":
                case "--version":
                    log.append(String.format(
                        "Ecotype Simulation %s\n", masterVariables.getVersion()
                    ));
                    noGUI = true;
                    break;
                default:
                    // Look for file name arguments.
                    File arg = new File(key);
                    // Expect the fasta file first.
                    if (fastaFile == null && arg.isFile()) {
                        fastaFile = arg;
                        break;
                    }
                    else if (fastaFile == null) {
                        // File not found, print an error.
                        System.out.print(String.format(
                            "Syntax error: File %s not found.\n%s",
                            key, usage
                        ));
                        System.exit(1);
                    }
                    // Expect the newick file second.
                    if (newickFile == null && arg.isFile()) {
                        newickFile = arg;
                        break;
                    }
                    else if (newickFile == null) {
                        // File not found, print an error.
                        System.out.print(String.format(
                            "Syntax error: File %s not found.\n%s",
                            key, usage
                        ));
                        System.exit(1);
                    }
                    // Expect the save file last.
                    if (fastaFile != null && newickFile != null &&
                        outputFile == null) {
                        outputFile = arg;
                        break;
                    }
                    if (key.length() > 0) {
                        // Argument unknown, print an error.
                        System.out.print(String.format(
                            "Syntax error: Unknown argument %s supplied.\n%s",
                            key, usage
                        ));
                        System.exit(1);
                    }
            }
        }
    }

    private boolean noGUI;
    private boolean runAll;
    private MasterVariables masterVariables;
    private Execs execs;
    private Logger log;
    private File fastaFile;
    private File newickFile;
    private File outputFile;

    private String usage =
        "Usage:\n" +
        "  java -jar EcoSim.jar [OPTIONS] /path/to/sequence_fasta " +
        "/path/to/sequence_tree output_xml\n" +
        "\n" +
        "  Sequences should be aligned and in a fasta formated file, with\n" +
        "  the outgroup listed first.\n" +
        "\n" +
        "  The newick formated tree file must include the same leaf node\n" +
        "  names and number as the sequences in the fasta file.\n" +
        "\n" +
        "  Output is saved in XML format.\n" +
        "\n" +
        "  Options:\n" +
        "    -d, --debug        : Display debugging output.\n" +
        "    -h, --help         : Display helpful information.\n" +
        "    -n, --nogui        : Hide the default GUI.  Implies --runall.\n" +
        "    -r, --runall       : Run everything, including demarcation.\n" +
        "    -t=n, --threads=n  : Set the number of threads (n) to start, default" + 
                                " to system maximum.\n" +
        "    -v, --version      : Display the version number.\n";

}
