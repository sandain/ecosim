/*
 *    Ecotype Simulation models the sequence diversity within a bacterial
 *    clade as the evolutionary result of net ecotype formation and periodic
 *    selection, yielding a certain number of ecotypes.
 *
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

import ecosim.gui.MainWindow;

import java.io.File;
import java.util.Observable;
import java.util.Observer;
import javax.swing.SwingUtilities;

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
 *     java -jar ecosim.jar [OPTIONS]
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
 *     -i, --input=[file]     : A XML formated save file for input.
 *     -o, --output=[file]    : A XML formated save for for output.
 *     -s, --sequences=[file] : A Fasta formated file for input.
 *     -p, --phylogeny=[file] : A Newick formatted file for input.
 *     -d, --debug            : Display debugging output.
 *     -h, --help             : Display helpful information.
 *     -n, --nogui            : Hide the default GUI.  Implies --runall.
 *     -r, --runall           : Run everything, including demarcation.
 *     -t, --threads=[n]      : Set the number of threads (n) to start,
 *                              default to system maximum.
 *     -v, --version          : Display the version number.
 *
 * @section fortran Fortran
 * @subsection fortran-programs Programs
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
 * @li @b api.Painter - Defines a custom method to paint on a surface.
 * @li @b gui.ButtonPane - Defines the main button panel for the GUI.
 * @li @b gui.FileChooser - Defines a custom file chooser for the GUI.
 * @li @b gui.HelpAboutWindow - Displays a Help/About GUI window.
 * @li @b gui.LoggerPane - Defines a custom panel to display the Logger.
 * @li @b gui.MainWindow - Defines the main GUI window.
 * @li @b gui.MenuBar - Defines a custom menu bar for the GUI.
 * @li @b gui.OptionsPane - Defines a custom panel to display the options.
 * @li @b gui.SummaryPane - Defines a custom panel to display the summary.
 * @li @b gui.TiledPainter - Defines a custom tile-based painter for the GUI.
 * @li @b tree.InvalidTreeException - Report a malformed tree.
 * @li @b tree.NewickReader - Read a Newick formatted tree.
 * @li @b tree.NewickWriter - Write the tree in Newick format.
 * @li @b tree.Node - A node of a phylogenetic tree.
 * @li @b tree.SVGPainter - A text-based painter to save a tree in SVG format.
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
        this.args = args;
        // Initialize the variables.
        log = new Logger ();
        masterVariables = new MasterVariables ();
        simulation = new Simulation (log, masterVariables);
        noGUI = false;
        runAll = false;
    }

    /**
     *  Run Ecotype Simulation.
     */
    public void run () {
        // Display the program name and version.
        System.out.print (String.format (
            "Ecotype Simulation %s\n\n", masterVariables.getVersion ()
        ));
        // Check for command line arguments.
        checkArguments (args);
        // Startup the CLI or the GUI.
        setupInterface ();
        // Load the input file, or the sequence and phylogeny files.
        if (inputFile != null) {
            if (! inputFile.exists ()) {
                System.err.println ("Error, input file does not exist.");
                System.exit (1);
            }
            // Load the input file.
            simulation.loadProjectFile (inputFile);
        }
        else if (fastaFile != null) {
            if (! fastaFile.exists ()) {
                System.err.println ("Error, sequence file does not exist.");
                System.exit (1);
            }
            // Load the sequence file.
            simulation.loadSequenceFile (fastaFile);
            // Generate the phylogeny file if not provided.
            if (newickFile == null || ! newickFile.exists ()) {
                newickFile = simulation.generateTree (fastaFile);
            }
            // Load the phylogeny file.
            simulation.loadTreeFile (newickFile);
            // Run the binning and parameter estimate programs.
            simulation.runBinning ();
            simulation.runParameterEstimate ();
        }
        // Run the simulation if requested.
        if (runAll && simulation.treeLoaded ()) {
            simulation.runHillclimbing ();
            simulation.runConfidenceIntervals ();
            simulation.runDemarcation ();
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
                    System.out.println (usage);
                    System.exit (0);
                    break;
                case "-i":
                case "--input":
                    if (value.length () > 0) {
                        inputFile = new File (value);
                    }
                    break;
                case "-n":
                case "--nogui":
                    noGUI = true;
                    runAll = true;
                    break;
                case "-o":
                case "--output":
                    if (value.length () > 0) {
                        masterVariables.setOutputFile (new File (value));
                    }
                    break;
                case "-p":
                case "--phylogeny":
                    if (value.length () > 0) {
                        newickFile = new File (value);
                    }
                    break;

                case "-r":
                case "--runall":
                    runAll = true;
                    break;
                case "-s":
                case "--sequences":
                    if (value.length () > 0) {
                        fastaFile = new File (value);
                    }
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
                    // The version has already been printed, just exit.
                    System.exit (0);
                    break;
                default:
                    // Look for unrecognized options.
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
            SwingUtilities.invokeLater (gui);
        }
    }

    private String[] args;
    private boolean noGUI;
    private boolean runAll;
    private Logger log;
    private MasterVariables masterVariables;
    private Simulation simulation;
    private File fastaFile;
    private File newickFile;
    private File inputFile;

    private String usage =
        "Usage:\n" +
        "  java -jar ecosim.jar [OPTIONS]\n" +
        "\n" +
        "  Sequences should be aligned and in a Fasta formated file, with\n" +
        "  the outgroup listed first.\n" +
        "\n" +
        "  The Newick formated tree file must include the same leaf node\n" +
        "  names and number as the sequences in the Fasta file.\n" +
        "\n" +
        "  Output is saved in XML format.\n" +
        "\n" +
        "  Options:\n" +
        "    -i, --input=[file]     : A XML formated save file for input.\n" +
        "    -o, --output=[file]    : A XML formated save for for output.\n" +
        "    -s, --sequences=[file] : A Fasta formated file for input.\n" +
        "    -p, --phylogeny=[file] : A Newick formatted file for input.\n" +
        "    -d, --debug            : Display debugging output.\n" +
        "    -h, --help             : Display helpful information.\n" +
        "    -n, --nogui            : Hide the default GUI.  Implies" +
                                    " --runall.\n" +
        "    -r, --runall           : Run everything, including" +
                                    " demarcation.\n" +
        "    -t, --threads=[n]      : Set the number of threads (n) to" +
                                    " start, default to system maximum.\n" +
        "    -v, --version          : Display the version number.\n";

}
