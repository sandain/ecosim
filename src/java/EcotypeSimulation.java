/*
 *    Ecotype Simulation models the sequence diversity within a bacterial
 *    clade as the evolutionary result of net ecotype formation and periodic
 *    selection, yielding a certain number of ecotypes.
 *
 *    Copyright (C) 2009-2015  Jason M. Wood, Montana State University
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

import ecosim.gui.MainWindow;

import java.io.File;
import java.util.Observable;
import java.util.Observer;

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
 *
 *     java -jar EcoSim.jar [OPTIONS] sequence_fasta sequence_tree output_xml
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
 *
 *     -d, --debug        : Display debugging output.
 *     -h, --help         : Display helpful information.
 *     -n, --nogui        : Hide the default GUI.  Implies --runall.
 *     -r, --runall       : Run everything, including demarcation.
 *     -t=n, --threads=n  : Set the number of threads (n) to start, default
 *                          to system maximum.
 *     -v, --version      : Display the version number.
 *
 * @section fortran Fortran
 * @subsection fortran-programs Programs
 * @li @b ::fredmethod - Test the omega/sigma/npop space.
 * @li @b ::hillclimb - Performs a hillclimbing on the bruteforce result.
 * @li @b ::omegaci - Calculates the confidence interval of the omega value.
 * @li @b ::sigmaci - Calculates the confidence interval of the sigma value.
 * @li @b ::npopci - Calculates the confidence interval of the npop value.
 * @li @b ::demarcation - Calculates npop values used to demarcate ecotypes.
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
 * @li @b Binning - Object to run the binning algorithm.
 * @li @b Demarcation - Demarcates ecotypes based on the hillclimbing values
 *        and the phylogeny of the sequences using the ::demarcation program.
 * @li @b Execs - Holds the executable methods for the various programs.
 * @li @b Fasta - Handles the input and output of fasta formatted text files.
 * @li @b FredMethod - Object to interact with the ::fredmethod program.
 * @li @b Heapsorter - Runs the heapsort on a given set of data.
 * @li @b Hillclimb - Object to interact with the ::hillclimb program.
 * @li @b InvalidFastaException - Report a malformed Fasta file.
 * @li @b Logger - Display text to the user.
 * @li @b MasterVariables - Common variables used through the program.
 * @li @b NpopConfidenceInterval - Run the ::npopci program.
 * @li @b OmegaConfidenceInterval - Run the ::omegaci program.
 * @li @b ParameterEstimate - An object to estimate the parameter values.
 * @li @b ParameterSet - An object to store the parameter values.
 * @li @b ProjectFileIO - Perform IO operations for the XML project file.
 * @li @b SigmaConfidenceInterval - Run the ::sigmaci program.
 * @li @b Simulation - The shared methods of the simulation.
 * @li @b StreamGobbler - Captures output from the Fortran programs.
 * @li @b Summary - An object to hold summary data.
 * @li @b gui.ButtonPane - Defines the main button panel for the GUI.
 * @li @b gui.FileChooser - Defines a custom file chooser for the GUI.
 * @li @b gui.HelpAboutWindow - Displays a Help/About GUI window.
 * @li @b gui.LoggerPane - Defines a custom panel to display the Logger.
 * @li @b gui.MainWindow - Defines the main GUI window.
 * @li @b gui.MenuBar - Defines a custom menu bar for the GUI.
 * @li @b gui.OptionsPane - Defines a custom panel to display the options.
 * @li @b gui.SummaryPane - Defines a custom panel to display the summary.
 * @li @b tree.InvalidTreeException - Report a malformed tree.
 * @li @b tree.NewickReader - Read a Newick formatted tree.
 * @li @b tree.NewickWriter - Write the tree in Newick format.
 * @li @b tree.Node - A node of a phylogenetic tree.
 * @li @b tree.SVGWriter - Write the tree in SVG format.
 * @li @b tree.Tree - Interact with phylogenetic trees.
 */

/**
 *  The main class for instantiating Ecotype %Simulation.
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
    public EcotypeSimulation (String[] args) {
        // Initialize the variables.
        log = new Logger ();
        masterVariables = new MasterVariables ();
        noGUI = false;
        runAll = false;
        // Check for command line arguments.
        checkArguments (args);
        // Create the Simulation.
        simulation = new Simulation (
            log, masterVariables, fastaFile, newickFile
        );
    }

    /**
     *  Run Ecotype Simulation.
     */
    public void run () {
        // Display the program name and version.
        System.out.print (String.format (
            "Ecotype Simulation %s\n\n", masterVariables.getVersion ()
        ));
        // Startup the CLI or the GUI.
        setupInterface ();
        // Generate the tree if the fasta file was provided with a tree.
        if (fastaFile != null && fastaFile.exists () &&
            (newickFile == null || ! newickFile.exists ())
        ) {
            simulation.generateTree (fastaFile);
        }
        // Run the binning and parameter estimate programs.
        if (newickFile != null && newickFile.exists ()) {
            simulation.runBinning ();
            simulation.runParameterEstimate ();
        }
        // If the runAll flag has been set, run the simulation.
        if (runAll) {
            simulation.runHillclimbing ();
            simulation.runConfidenceIntervals ();
            simulation.runDemarcation ();
            if (outputFile != null) {
                simulation.saveProjectFile (outputFile);
            }
        }
        // Exit the simulation if the GUI wasn't started.
        if (noGUI) simulation.exit ();
    }

    /**
     *  Start an instance of EcotypeSimulation.
     *
     *  @param args The command line arguments.
     */
    public static void main (final String[] args) {
        EcotypeSimulation es = new EcotypeSimulation (args);
        Thread t = new Thread (es);
        t.start ();
    }

    /**
     *  Check for command line arguments.
     *
     *  @param args The command line arguments.
     */
    private void checkArguments (String[] args) {
        boolean error = false;
        for (int i = 0; i < args.length; i ++) {
            String[] keyValue = args[i].split ("=");
            String key = keyValue[0];
            String value = "";
            // Look for key=value pairs in the arguments.
            if (keyValue.length > 1) {
                value = keyValue[1];
            }
            // Catch a space after the equal sign.
            else if (args.length > i + 1 && args[i].endsWith ("=")) {
                value = args[i + 1];
                i ++;
            }
            // Catch a space before the equal sign.
            else if (args.length > i + 1 && args[i + 1].startsWith ("=") &&
                ! args[i + 1].equals ("=")) {
                value = args[i + 1].substring (1);
                i ++;
            }
            // Catch spaces surrounding the equal sign.
            else if (args.length > i + 2 && args[i + 1].equals ("=")) {
                value = args[i + 2];
                i += 2;
            }
            // Check for various arguments.
            switch (key) {
                case "-d":
                case "--debug":
                    masterVariables.setDebug (true);
                    break;
                case "-h":
                case "--help":
                    System.out.println (String.format (
                        "Ecotype Simulation %s\n\n%s",
                        masterVariables.getVersion (), usage
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
                    if (value.length () > 0) {
                        Runtime r = Runtime.getRuntime ();
                        int procs = r.availableProcessors ();
                        int num = 0;
                        try {
                            num = new Integer (value).intValue ();
                        }
                        catch (NumberFormatException e) {
                            System.out.println (String.format (
                                "Syntax error: Expected a number.\n%s\n%s",
                                e, usage
                            ));
                            System.exit (1);
                        }
                        if (num >= 1 && num <= procs) {
                            masterVariables.setNumberThreads (num);
                        }
                        else {
                            System.out.println (String.format (
                                "Syntax error: Invalid number of threads " +
                                "specified: %d of %d possible.\n",
                                num, procs
                            ));
                            System.exit (1);
                        }
                    }
                    else {
                        // Number of threads not provided, print an error.
                        System.out.println (String.format (
                            "Syntax error: Number of threads missing.\n%s",
                            usage
                        ));
                        System.exit (1);
                    }
                    break;
                case "-v":
                case "--version":
                    log.append (String.format (
                        "Ecotype Simulation %s\n",
                        masterVariables.getVersion ()
                    ));
                    noGUI = true;
                    break;
                default:
                    // Look for unrecognized options.
                    if (key.startsWith ("-")) {
                        System.out.println (String.format (
                            "Syntax error: Option %s not recognized.\n%s",
                            key, usage
                        ));
                        System.exit (1);
                    }
                    // Look for file name arguments.
                    File arg = new File (key);
                    // Expect the fasta file first.
                    if (fastaFile == null && arg.isFile ()) {
                        fastaFile = arg;
                        break;
                    }
                    else if (fastaFile == null) {
                        // File not found, print an error.
                        System.out.println (String.format (
                            "Syntax error: File %s not found.\n%s",
                            key, usage
                        ));
                        System.exit (1);
                    }
                    // Expect the newick file second.
                    if (newickFile == null && arg.isFile ()) {
                        newickFile = arg;
                        break;
                    }
                    else if (newickFile == null) {
                        // File not found, print an error.
                        System.out.println (String.format (
                            "Syntax error: File %s not found.\n%s",
                            key, usage
                        ));
                        System.exit (1);
                    }
                    // Expect the save file last.
                    if (fastaFile != null && newickFile != null &&
                        outputFile == null) {
                        outputFile = arg;
                        break;
                    }
                    if (key.length () > 0) {
                        // Argument unknown, print an error.
                        System.out.println (String.format (
                            "Syntax error: Unknown argument %s supplied.\n%s",
                            key, usage
                        ));
                        System.exit (1);
                    }
            }
        }
    }

    /**
     *  Startup either the command line interface (CLI) or the graphical
     *  user interface (GUI) depending on the noGUI flag.
     */
    private void setupInterface () {
        if (noGUI) {
            // The CLI only requires adding an observer to the log.  When
            // the log changes, the new text is printed to the terminal.
            log.addObserver (new Observer () {
                public void update (Observable o, Object str) {
                    System.out.print ((String) str);
                }
            });
        }
        else {
            // Start the main window of the GUI.
            MainWindow gui = new MainWindow (
                log, masterVariables, simulation
            );
            gui.run ();
        }
    }

    private boolean noGUI;
    private boolean runAll;
    private Logger log;
    private MasterVariables masterVariables;
    private Simulation simulation;
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
        "    -n, --nogui        : Hide the default GUI.  Implies" +
                                " --runall.\n" +
        "    -r, --runall       : Run everything, including demarcation.\n" +
        "    -t=n, --threads=n  : Set the number of threads (n) to start," +
                                " default to system maximum.\n" +
        "    -v, --version      : Display the version number.\n";

}
