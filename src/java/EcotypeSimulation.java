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

import java.awt.EventQueue;
import java.awt.BorderLayout;
import java.awt.Dimension;
import java.awt.GridLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.WindowEvent;
import java.awt.event.WindowAdapter;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.StringTokenizer;
import javax.swing.Box;
import javax.swing.JButton;
import javax.swing.JComboBox;
import javax.swing.DefaultComboBoxModel;
import javax.swing.JComponent;
import javax.swing.JDialog;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JMenu;
import javax.swing.JMenuBar;
import javax.swing.JMenuItem;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JScrollPane;
import javax.swing.JTextArea;

/**
  *  EcotypeSimulation
  *
  */
public class EcotypeSimulation extends JFrame {

    /**
     *  EcotypeSimulation constructor.
     */
    public EcotypeSimulation() {
        masterVariables = new MasterVariables();
        helpAbout = new HelpAboutWindow(masterVariables);
        log = masterVariables.getLog();
        narr = masterVariables.getNarrator();
        workingDirectory = masterVariables.getWorkingDirectory();
        execs = new Execs(masterVariables);
        masterVariables.setExecs(execs);
        initComponents();
    }

    /**
     *  EcotypeSimulation constructor.
     *
     *  @param args The command line arguments.
     */
    public EcotypeSimulation(String[] args) {
        masterVariables = new MasterVariables();
        helpAbout = new HelpAboutWindow(masterVariables);
        log = masterVariables.getLog();
        narr = masterVariables.getNarrator();
        workingDirectory = masterVariables.getWorkingDirectory();
        checkArguments(args);
        execs = new Execs(masterVariables);
        masterVariables.setExecs(execs);
        initComponents();
    }

    /**
     *  Start an instance of EcotypeSimulation.
     *
     *  @param args The command line arguments.
     */
    public static void main(final String[] args) {
        Runnable r = new Runnable() {
            public void run() {
                new EcotypeSimulation(args);
            }
        };
        EventQueue.invokeLater(r);
    }

    /**
     *  Runs the drift confidence interval.
     *
     *  @return FredOutVal that is null if no better value is found, but if during the confidence interval a set of
     *  parameters with a better likelihood is found (via) hillclimbing it will return that value.
     */
    public FredOutVal runDriftConfidenceInterval() {
       log.append ("Starting drift confidence interval...");
       log.append (System.getProperty ("line.separator"));
       narr.println("Starting drift confidence interval...");
       ArrayList<String> bins = values.getBins();
       int[] sequenceVals = values.getSeqVals();
       File driftIn = new File(workingDirectory + "driftIn.dat");
       File driftOut = new File(workingDirectory + "driftOut.dat");
       DriftConfidence drift = new DriftConfidence(hClimbResult, values, masterVariables, driftIn, driftOut);
       // Run the lower bound.
       FredOutVal driftRes = drift.lowerBound();
       log.append ("The result from driftCI lower bound:");
       log.append (System.getProperty ("line.separator"));
       log.append ("  " + driftRes.toString());
       log.append (System.getProperty ("line.separator"));
       narr.println("The result from driftCI:");
       narr.println("  " + driftRes.toString());
       driftConfidenceInterval = 1 / driftRes.getDrift();
       return null;
    }

    /**
     *  Runs the omega confidence interval.
     *
     *  @return FredOutVal that is null if no value with a better likelihood was found, otherwise
     *  it returns the value with that better likelihood.
     */
    public FredOutVal runOmegaConfidenceInterval() {
       log.append ("Starting omega confidence interval...");
       log.append (System.getProperty ("line.separator"));
       narr.println("Starting omega confidence interval...");
       ArrayList<String> bins = values.getBins();
       int[] sequenceVals = values.getSeqVals();
       File omegaIn = new File(workingDirectory + "omegaIn.dat");
       File omegaOut = new File(workingDirectory + "omegaOut.dat");
       OmegaConfidence omega = new OmegaConfidence(hClimbResult, values, masterVariables, omegaIn, omegaOut);
       // Run the lower bound.
       FredOutVal omegaRes = omega.lowerBound();
       log.append ("The result from omegaCI lower bound:");
       log.append (System.getProperty ("line.separator"));
       log.append ("  " + omegaRes.toString());
       log.append (System.getProperty ("line.separator"));
       narr.println("The result from omegaCI:");
       narr.println("  " + omegaRes.toString());
       // Run the upper bound.
       FredOutVal omegaUpRes = omega.upperBound();
       log.append ("The result from omegaCI upper bound:");
       log.append (System.getProperty ("line.separator"));
       log.append ("  " + omegaUpRes.toString());
       log.append (System.getProperty ("line.separator"));
       narr.println("The result from omegaCI upper bound:");
       narr.println("  " + omegaUpRes.toString());
       omegaConfidenceInterval = new double[2];
       omegaConfidenceInterval[0] = omegaRes.getOmega();
       omegaConfidenceInterval[1] = omegaUpRes.getOmega();
       return null;
    }

    /**
     *  This method is called from within the constructor to initialize the form.
     */
    private void initComponents() {
        runToHClimb = new JButton();
        npopConf = new JButton();
        sigmaConf = new JButton();
        omegaConf = new JButton();
        driftConf = new JButton();
        runAll = new JButton();
        critSelector = new JComboBox();
        criterion = new JLabel();
        runAllDemarcs = new JButton();
        openFastaFile = new JMenuItem();
        openSaveFile = new JMenuItem();
        saveFile = new JMenuItem();
        changePCRError = new JMenuItem();
        exitProgram = new JMenuItem();
        demarcations = new JMenuItem();
        demarcsAuto = new JMenuItem();
        userGuideWindow = new JMenuItem();
        aboutWindow = new JMenuItem();
        optionsWindow = new JMenuItem();
        // Start the HelpAboutWindow thread.
        helpAbout = new HelpAboutWindow(masterVariables);
        new Thread(helpAbout).start();
        /**
         *  Need this part to get the narrator to close!
         */
        addWindowListener(new WindowAdapter() {
            public void windowActivated(WindowEvent we) {
                log.repaint();
            }

            public void windowClosing(WindowEvent we) {
                printResults();
                if (narr != null)
                    narr.close();
                System.exit(0);
            }
        });
        setTitle("Ecotype Simulation");
        setJMenuBar(makeMenuBar());
        setLayout(new BorderLayout(15, 15));
        setMinimumSize(new Dimension(640, 400));
        setPreferredSize(new Dimension(640, 700));
        add(makeButtonPane(), BorderLayout.NORTH);
        add(makeLogPane(), BorderLayout.CENTER);
        pack();
        setVisible(true);
    }

    private boolean checkArguments (String[] args) {
        boolean error = false;
        String errorMessage = "";
        errorMessage += "Syntax error: Unknown argument supplied.";
        errorMessage += System.getProperty ("line.separator");
        errorMessage += "Usage: java -jar EcoSim.jar [OPTION]...";
        errorMessage += System.getProperty ("line.separator");
        errorMessage += "  Options";
        errorMessage += System.getProperty ("line.separator");
        errorMessage += "    -d, --debug   : Displays debugging output.";
        errorMessage += System.getProperty ("line.separator");
        errorMessage += "    -a, --about   : Displays the about window.";
        errorMessage += System.getProperty ("line.separator");
        errorMessage += "    -h, --help    : Displays the help window.";
        errorMessage += System.getProperty ("line.separator");
        errorMessage += "    -l, --license : Displays the license window.";
        errorMessage += System.getProperty ("line.separator");
        errorMessage += "    -v, --version : Displays the version window.";
        errorMessage += System.getProperty ("line.separator");
        for (String arg: args) {
            if (arg.equals("-d") || arg.equals("--debug")) {
                masterVariables.setDebug(true);
            }
            else if (arg.equals("-a") || arg.equals("--about")) {
                helpAboutActionPerformed(HelpAboutWindow.about);
            }
            else if (arg.equals("-h") || arg.equals("--help")) {
                helpAboutActionPerformed(HelpAboutWindow.userGuide);
            }
            else if (arg.equals("-l") || arg.equals("--license")) {
                helpAboutActionPerformed(HelpAboutWindow.license);
            }
            else if (arg.equals("-v") || arg.equals("--version")) {
                helpAboutActionPerformed(HelpAboutWindow.version);
            }
            else {
                error = true;
                System.out.print(errorMessage);
            }
        }
        return error;
    }

    private JMenuBar makeMenuBar() {
        JMenuBar menuBar = new JMenuBar();
        JMenu fileMenu = new JMenu();
        JMenu otherProgMenu = new JMenu();
        JMenu helpAboutMenu = new JMenu();
        fileMenu.setText("File");
        openFastaFile.setText("Open Fasta File");
        openFastaFile.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent evt) {
                openFastaFileActionPerformed(evt);
            }
        });
        openSaveFile.setText("Open Saved File");
        openSaveFile.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent evt) {
                openSaveFileActionPerformed(evt);
            }
        });
        saveFile.setText("Save Current Project");
        saveFile.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent evt) {
                saveFileActionPerformed(evt);
            }
        });
        changePCRError.setText("Modify PCR Error");
        changePCRError.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent evt) {
                changePCRErrorActionPerformed(evt);
            }
        });
        optionsWindow.setText("Options");
        optionsWindow.addActionListener( new ActionListener() {
            public void actionPerformed(ActionEvent evt) {
                optionsWindowActionPerformed(evt);
            }
        });
        exitProgram.setText("Exit");
        exitProgram.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent evt) {
                exitProgramActionPerformed(evt);
            }
        });
        // Add Menu Items to the File Menu.
        fileMenu.add(openFastaFile);
        fileMenu.add(openSaveFile);
        fileMenu.add(saveFile);
        fileMenu.add(changePCRError);
        fileMenu.add(optionsWindow);
        fileMenu.add(exitProgram);
        // Create Other Programs Menu.
        otherProgMenu.setText("Other Programs");
        demarcations.setText("Demarcations");
        demarcations.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent evt) {
                demarcationsActionPerformed(evt);
            }
        });
        demarcsAuto.setText("Demarcations (Auto)");
        demarcsAuto.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent evt) {
                demarcsAutoActionPerformed(evt);
            }
        });
        // Create Help About Menu.
        helpAboutMenu.setText("Help");
        userGuideWindow.setText(HelpAboutWindow.userGuide);
        userGuideWindow.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent evt) {
                helpAboutActionPerformed(evt.getActionCommand());
            }
        });
        aboutWindow.setText(HelpAboutWindow.about);
        aboutWindow.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent evt) {
                helpAboutActionPerformed(evt.getActionCommand());
            }
        });
        helpAboutMenu.add(userGuideWindow);
        helpAboutMenu.add(aboutWindow);
        // Add Menu Items to the Other Programs Menu.
        otherProgMenu.add(demarcations);
        otherProgMenu.add(demarcsAuto);
        // Add Menus to the MenuBar.
        menuBar.add(fileMenu);
        menuBar.add(otherProgMenu);
        menuBar.add(helpAboutMenu);
        return menuBar;
    }

    private JComponent makeButtonPane() {
        JPanel pane = new JPanel();
        JPanel buttons = new JPanel();
        // Create a border around the button area.
        pane.setLayout(new BorderLayout());
        pane.add(Box.createRigidArea(new Dimension(0, 15)), BorderLayout.NORTH);
        pane.add(Box.createRigidArea(new Dimension(0, 15)), BorderLayout.SOUTH);
        pane.add(Box.createRigidArea(new Dimension(15, 0)), BorderLayout.EAST);
        pane.add(Box.createRigidArea(new Dimension(15, 0)), BorderLayout.WEST);
        // Add the buttons.
        buttons.setLayout(new GridLayout(1, 2, 15, 15));
        buttons.add(makeSimulationButtonPane());
        buttons.add(makeConfidenceIntervalButtonPane());
        pane.add(buttons, BorderLayout.CENTER);
        return pane;
    }

    private JComponent makeSimulationButtonPane() {
        JPanel pane = new JPanel();
        pane.setLayout(new GridLayout(4,1,15,15));
        // Create Buttons.
        runToHClimb.setText("Run through HillClimbing");
        runToHClimb.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent evt) {
                runToHClimbActionPerformed(evt);
            }
        });
        runAll.setText("Run Everything (no demarcations)");
        runAll.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent evt) {
                runAllActionPerformed(evt);
            }
        });
        runAllDemarcs.setText("Run Everything (with demarcations)");
        runAllDemarcs.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent evt) {
                runAllDemarcsActionPerformed(evt);
            }
        });
        // Create Criterion Selector.
        JPanel crit = new JPanel();
        critSelector.setModel(new DefaultComboBoxModel(
          new String[] { "auto", "5x", "2x", "1.5x", "1.25x", "1.10x", "1.05x" })
        );
        JLabel critLabel = new JLabel("Select Criterion: ", JLabel.TRAILING);
        criterion.setText("Select Criterion");
        crit.add(critLabel);
        crit.add(critSelector);
        // Add everything to the pane.
        pane.add(runToHClimb);
        pane.add(runAll);
        pane.add(runAllDemarcs);
        pane.add(crit);
        return pane;
    }

    private JComponent makeConfidenceIntervalButtonPane() {
        JPanel pane = new JPanel();
        pane.setLayout(new GridLayout(4,1,15,15));
        // Create Buttons.
        npopConf.setText("Run Npop Confidence Interval");
        npopConf.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent evt) {
                npopConfActionPerformed(evt);
            }
        });
        sigmaConf.setText("Run Sigma Confidence Interval");
        sigmaConf.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent evt) {
                sigmaConfActionPerformed(evt);
            }
        });
        omegaConf.setText("Run Omega Confidence Interval");
        omegaConf.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent evt) {
                omegaConfActionPerformed(evt);
            }
        });
        driftConf.setText("Run Drift Confidence Interval");
        driftConf.addActionListener(new ActionListener() {
            public void actionPerformed(ActionEvent evt) {
                driftConfActionPerformed(evt);
            }
        });
        // Add the Buttons to the Pane.
        pane.add(npopConf);
        pane.add(sigmaConf);
        pane.add(omegaConf);
        pane.add(driftConf);
        return pane;
    }

    private JComponent makeLogPane() {
        JScrollPane pane = new JScrollPane();
        log.setColumns(20);
        log.setEditable(false);
        log.setRows(5);
        log.setDoubleBuffered(true);
        pane.setViewportView(log);
        return pane;
    }

    private void changePCRErrorActionPerformed(ActionEvent evt) {
        if (evt.getSource() == changePCRError) {
            Thread thread = new Thread() {
                public void run() {
                    double value = -1;
                    while (true) {
                        String userInput = JOptionPane.showInputDialog("Please input a value for pcrerror:");
                        if (userInput == null)
                            break;
                        try {
                            value = (new Double(userInput)).doubleValue();
                        }
                        catch (NumberFormatException e) {
                            continue;
                        }
                        break;
                    }
                    execs.changePCRError(value);
                }
            };
            thread.start();
        }
    }

    private void exitProgramActionPerformed(ActionEvent evt) {
        if (evt.getSource() == exitProgram) {
            printResults();
            if (narr != null)
                narr.close();
            System.exit(0);
        }
    }

    private void demarcationsActionPerformed(ActionEvent evt) {
        if (evt.getSource() == demarcations) {
            if (readyForCI) {
                Thread thread = new Thread() {
                    public void run() {
                        new Demarcations(hClimbResult, values, masterVariables).setVisible(true);
                    }
                };
                thread.start();
            }
            else {
                log.append (
                    "Please run through hillclimbing or load a saved " +
                    "file before running demarcations."
                );
                log.append (System.getProperty ("line.separator"));
            }
        }
    }

    private void demarcsAutoActionPerformed(ActionEvent evt) {
        if (evt.getSource() == demarcsAuto) {
            if (readyForCI) {
                FileChooser fc = new FileChooser("fasta");
                int returnVal = fc.showOpenDialog(this);
                if (returnVal == FileChooser.APPROVE_OPTION) {
                    File inputFile = fc.getSelectedFile();
                    inputFasta = new Fasta(inputFile);
                    noOutgroupFasta = removeOutgroup(inputFasta);
                    log.append ("Opening: " + inputFile.getName());
                    log.append (System.getProperty ("line.separator"));
                    // Test that the file extension is correct and that the file can be read.
                    if (! inputFasta.isValid()) {
                        log.append ("That is not a valid fasta file.");
                        log.append (System.getProperty ("line.separator"));
                    }
                    else {
                        if (chooseTree() != -1) {
                             Thread thread = new Thread() {
                                public void run() {
                                    new DemarcationsAuto(hClimbResult, values, inputFasta, tree, masterVariables).setVisible(true);
                                }
                            };
                            thread.start();
                        }
                    }
                }
            }
            else {
                log.append (
                    "Please run through hillclimbing or load a saved file " +
                    "before running demarcations."
                );
                log.append (System.getProperty ("line.separator"));
            }
        }
    }

    private int chooseTree() {
        Object[] options = {"Phylip", "Newick"};
        int useNewick = JOptionPane.showOptionDialog(EcotypeSimulation.this,
                        "Generate Neighbor Joining/Parsimony tree from Fasta file using Phylip "
                        + "or provide own Newick tree?",
                        "Tree Generation",
                        JOptionPane.YES_NO_OPTION,
                        JOptionPane.QUESTION_MESSAGE,
                        null, options, options[0]);
        if (useNewick == 1) {
            FileChooser fc = new FileChooser("newick");
            int returnVal = fc.showOpenDialog(this);
            if (returnVal == FileChooser.APPROVE_OPTION) {
                tree = fc.getSelectedFile();
                log.append ("Opening: " + tree.getName ());
                log.append (System.getProperty ("line.separator"));
            }
            else {
                useNewick = -1;
            }
        }
        else {
            tree = null;
        }
        return useNewick;
    }

    /**
     *  Opens the Options Window.
     */
    private void optionsWindowActionPerformed(ActionEvent evt) {
        if (evt.getSource() == optionsWindow) {
            Thread thread = new Thread() {
                public void run() {
                    OptionsWindow options = new OptionsWindow(masterVariables);
                }
            };
            thread.start();
        }
    }

    /**
     *  Opens the Help/About Window.
     *
     *  @param item
     */
    private void helpAboutActionPerformed(String item) {
        if (item == HelpAboutWindow.userGuide) {
            helpAbout.setVisible(HelpAboutWindow.userGuide);
        }
        else if (item == HelpAboutWindow.about) {
            helpAbout.setVisible(HelpAboutWindow.about);
        }
    }

    /**
     *  Save the file appropriately.
     */
    private void saveFileActionPerformed(ActionEvent evt) {
        if (evt.getSource() == saveFile) {
            FileChooser fc = new FileChooser("cpm");
            int returnVal = fc.showSaveDialog(this);
            if (returnVal == FileChooser.APPROVE_OPTION) {
                File userFile = fc.getSelectedFile();
                String path = userFile.getPath();
                //if it doesn't already have the cpm extension, add it
                if (! ((path.substring(path.length() - 4, path.length())).equals(".cpm"))) {
                    userFile = new File(userFile.getPath()+".cpm");
                }
                log.append ("Saving to: " + userFile.getName ());
                log.append (System.getProperty ("line.separator"));
                narr.println("Saving to: " + userFile.getName());
                if (!saveFile(userFile)) {
                    log.append ("Please run up to hillclimbing before saving to a file.");
                    log.append (System.getProperty ("line.separator"));
                }
            }
        }
    }

    private void openSaveFileActionPerformed(ActionEvent evt) {
        if (evt.getSource() == openSaveFile) {
            FileChooser fc = new FileChooser("cpm");
            int returnVal = fc.showOpenDialog(this);
            if (returnVal == FileChooser.APPROVE_OPTION) {
                File savedFile = fc.getSelectedFile();
                log.append ("Opening: " + savedFile.getName());
                log.append (System.getProperty ("line.separator"));
                String name = savedFile.getName();
                // Test that the file extension is correct and that the file can be read.
                if (! savedFile.canRead() || ! ((name.substring(name.length() - 4, name.length()).equals(".cpm"))) ||
                    ! recoverSavedData(savedFile)) {
                    log.append (
                        "That is not a valid saved file, please choose" +
                        " a file previously saved in this program."
                    );
                    log.append (System.getProperty ("line.separator"));
                }
            }
        }
    }

    private void openFastaFileActionPerformed(ActionEvent evt) {
        if (evt.getSource() == openFastaFile) {
            FileChooser fc = new FileChooser("fasta");
            int returnVal = fc.showOpenDialog(this);
            if (returnVal == FileChooser.APPROVE_OPTION) {
                File inputFile = fc.getSelectedFile();
                inputFasta = new Fasta(inputFile);
                noOutgroupFasta = removeOutgroup(inputFasta);
                log.append ("Opening: " + inputFile.getName ());
                log.append (System.getProperty ("line.separator"));
                if (! inputFasta.isValid()) {
                    log.append (
                        "That is not a valid fasta file, please choose a " +
                        "properly formatted fasta file."
                    );
                    log.append (System.getProperty ("line.separator"));
                    return;
                }
                narr.println("Opening: " + inputFile.getName());
            }
        }
    }

    private void runAllActionPerformed(ActionEvent evt) {
        if (evt.getSource() == runAll) {
            if (! inputFasta.isValid()) {
                log.append ("Please open a valid fasta file first.");
                log.append (System.getProperty ("line.separator"));
                return;
            }
            Thread thread = new Thread() {
                public void run() {
                    if (running) {
                        log.append ("Already running...");
                        log.append (System.getProperty ("line.separator"));
                    }
                    else {
                        running = true;
                        int userChoice = userSortPercentage();
                        values = new BinningAndFred(noOutgroupFasta, masterVariables, userChoice);
                        values.run();
                        hClimbResult = hillClimbing(values.getValue());
                        runNpopConfidenceInterval();
                        runSigmaConfidenceInterval();
                        runOmegaConfidenceInterval();
                        runDriftConfidenceInterval();
                        readyForCI = true;
                        printResults();
                        running = false;
                    }
                }
            };
            thread.start();
        }
    }

    private void runToHClimbActionPerformed(ActionEvent evt) {
        if (evt.getSource() == runToHClimb) {
            if (! inputFasta.isValid()) {
                log.append ("Please open a valid fasta file first.");
                log.append (System.getProperty ("line.separator"));
                return;
            }
            Thread thread = new Thread() {
                public void run() {
                    if (running) {
                        log.append ("Already running...");
                        log.append (System.getProperty ("line.separator"));
                    }
                    else {
                        running = true;
                        int userChoice = userSortPercentage();
                        values = new BinningAndFred(noOutgroupFasta, masterVariables, userChoice);
                        values.run();
                        hClimbResult =  hillClimbing(values.getValue());
                        readyForCI = true;
                        running = false;
                    }
                }
            };
            thread.start();
        }
    }

    private void sigmaConfActionPerformed(ActionEvent evt) {
        if (evt.getSource() == sigmaConf) {
            if (readyForCI) {
                Thread thread = new Thread() {
                    public void run() {
                        if (running) {
                            log.append ("Already running...");
                            log.append (System.getProperty ("line.separator"));
                        }
                        else {
                            running = true;
                            FredOutVal result = runSigmaConfidenceInterval();
                            if (result != null) {
                                if (userApproval(result)) {
                                    readyForCI = false;
                                    hClimbResult = hillClimbing(result);
                                    readyForCI = true;
                                }
                            }
                            running = false;
                        }
                    }
                };
                thread.start();
            }
            else {
                log.append (
                    "You must run binning and hillclimbing from a " +
                    "fasta file before you can run a confidence interval."
                );
                log.append (System.getProperty ("line.separator"));
            }
        }
    }

    private void omegaConfActionPerformed(ActionEvent evt) {
        if (evt.getSource() == omegaConf) {
            if (readyForCI) {
                Thread thread = new Thread() {
                    public void run() {
                        if (running) {
                            log.append ("Already running...");
                            log.append (System.getProperty ("line.separator"));
                        }
                        else {
                            running = true;
                            FredOutVal result = runOmegaConfidenceInterval();
                            if (result != null) {
                                if (userApproval(result)) {
                                    readyForCI = false;
                                    hClimbResult = hillClimbing(result);
                                    readyForCI = true;
                                }
                            }
                            running = false;
                        }
                    }
                };
                thread.start();
            }
            else {
                log.append (
                    "You must run binning and hillclimbing from a " +
                    "fasta file before you can run a confidence interval."
                );
                log.append (System.getProperty ("line.separator"));
            }
        }
    }

    private void npopConfActionPerformed(ActionEvent evt) {
        if (evt.getSource() == npopConf) {
            if (readyForCI) {
                Thread thread = new Thread() {
                    public void run() {
                        if (running) {
                            log.append ("Already running...");
                            log.append (System.getProperty ("line.separator"));
                        }
                        else {
                            running = true;
                            FredOutVal result = runNpopConfidenceInterval();
                            if (result != null) {
                                if (userApproval(result)) {
                                    readyForCI = false;
                                    hClimbResult = hillClimbing(result);
                                    readyForCI = true;
                                }
                            }
                            running = false;
                        }
                    }
                };
                thread.start();
            }
            else {
                log.append (
                    "You must run binning and hillclimbing from a " +
                    "fasta file before you can run a confidence interval."
                );
                log.append (System.getProperty ("line.separator"));
            }
        }
    }

    private void driftConfActionPerformed(ActionEvent evt) {
        if (evt.getSource() == driftConf) {
            if (readyForCI) {
                Thread thread = new Thread() {
                    public void run() {
                        if (running) {
                            log.append ("Already running...");
                            log.append (System.getProperty ("line.separator"));
                        }
                        else {
                            running = true;
                            FredOutVal result = runDriftConfidenceInterval();
                            if (result != null) {
                                if (userApproval(result)) {
                                    readyForCI = false;
                                    hClimbResult = hillClimbing(result);
                                    readyForCI = true;
                                }
                            }
                            running = false;
                        }
                    }
                };
                thread.start();
            }
            else {
                log.append (
                    "You must run binning and hillclimbing from a " +
                    "fasta file before you can run a confidence interval."
                );
                log.append (System.getProperty ("line.separator"));
            }
        }
    }

    private void runAllDemarcsActionPerformed(ActionEvent evt) {
        if (evt.getSource() == runAllDemarcs) {
            if (! inputFasta.isValid()) {
                log.append ("Please open a valid fasta file first.");
                log.append (System.getProperty ("line.separator"));
                return;
            }
            Thread thread = new Thread() {
                public void run() {
                    if (running) {
                        log.append ("Already running...");
                        log.append (System.getProperty ("line.separator"));
                    }
                    else {
                        running = true;
                        if (chooseTree() != -1) {
                            int userChoice = userSortPercentage();
                            values = new BinningAndFred(noOutgroupFasta, masterVariables, userChoice);
                            values.run();
                            hClimbResult =  hillClimbing(values.getValue());
                            runNpopConfidenceInterval();
                            runSigmaConfidenceInterval();
                            runOmegaConfidenceInterval();
                            runDriftConfidenceInterval();
                            readyForCI = true;
                            printResults();
                            new DemarcationsAuto(hClimbResult, values, inputFasta, tree, masterVariables).setVisible(true);
                        }
                        running = false;
                    }
                }
            };
            thread.start();
        }
    }

    /**
     *  Remove the outgroup (the first gene) from the input file.
     *  Note: all input files must have the first gene be an outgroup! this may
     *  be generated by whoever is using the program; if there is no outgroup,
     *  the first sequence in the file will be removed.
     *
     *  @param inFasta The FASTA formated input file.
     *  @return The FASTA file with the outgroup removed.
     */
    private Fasta removeOutgroup(Fasta inFasta) {
        Fasta outFasta = new Fasta(inFasta);
        outFasta.remove(0);
        outFasta.save(masterVariables.getNoOutgroup());
        return outFasta;
    }

    /**
     *  Runs the sigma confidence interval.
     *
     *  @return null if no value with a better likelihood than the one found
     *  in hillclimbing, otherwise it returns that better value.
     */
    private FredOutVal runSigmaConfidenceInterval() {
       log.append ("Starting sigma confidence interval: ");
       log.append (System.getProperty ("line.separator"));
       narr.println("Starting sigma confidence interval: ");
       ArrayList<String> bins = values.getBins();
       int[] sequenceVals = values.getSeqVals();
       // Run the new sigma confidence interval from the estimated interval.
       File sigmaIn = new File(workingDirectory + "sigmaIn.dat");
       File sigmaOut = new File(workingDirectory + "sigmaOut.dat");
       SigmaConfidence sigma = new SigmaConfidence(hClimbResult, values, masterVariables, sigmaIn, sigmaOut);
       // Run the lower bound.
       FredOutVal sigmaRes = sigma.lowerBound();
       log.append ("The result from the lower bound of sigmaCI:");
       log.append (System.getProperty ("line.separator"));
       log.append ("  " + sigmaRes.toString ());
       log.append (System.getProperty ("line.separator"));
       narr.println("The result from the lower bound of sigmaCI:");
       narr.println("  " + sigmaRes.toString());
       // Run the upper bound.
       FredOutVal sigmaUpRes = sigma.upperBound();
       // If sigma is greater than 100, the sigma confidence interval class will have returned null for the answer to the upper bound.
       if (sigmaUpRes == null) {
           log.append ("The upper bound of sigma is greater than 100");
           log.append (System.getProperty ("line.separator"));
           narr.println("The upper bound of sigma is greater than 100");
           sigmaConfidenceInterval = new double[2];
           sigmaConfidenceInterval[0] = sigmaRes.getSigma();
           sigmaConfidenceInterval[1] = 100.0;
           return null;
        }
       log.append ("The result from sigmaCI upper bound:");
       log.append (System.getProperty ("line.separator"));
       log.append ("  " + sigmaUpRes.toString ());
       log.append (System.getProperty ("line.separator"));
       narr.println("The result from sigmaCI upper bound:");
       narr.println("  " + sigmaUpRes.toString());
       sigmaConfidenceInterval = new double[2];
       sigmaConfidenceInterval[0] = sigmaRes.getSigma();
       sigmaConfidenceInterval[1] = sigmaUpRes.getSigma();
       return null;
    }

    /**
     *  Runs the npop confidence interval for the given hillclimbing value.
     *
     *  @return FredOutVal with a better likelihood if there is one found
     *  during the confidence interval run, otherwise null
     */
    private FredOutVal runNpopConfidenceInterval() {
       log.append ("Starting npop confidence interval...");
       log.append (System.getProperty ("line.separator"));
       narr.println("Starting npop confidence interval...");
       ArrayList<String> bins = values.getBins();
       int[] sequenceVals = values.getSeqVals();
       File npopIn = new File(workingDirectory + "npopIn.dat");
       File npopOut = new File(workingDirectory + "npopOut.dat");
       NpopConfidence npop = new NpopConfidence(hClimbResult, values, masterVariables, npopIn, npopOut);
       // Run the lower bound.
       FredOutVal npopRes = npop.lowerBound();
       log.append ("The result from npopCI lower bound:");
       log.append (System.getProperty ("line.separator"));
       log.append ("  " + npopRes.toString ());
       log.append (System.getProperty ("line.separator"));
       narr.println("The result from npopCI:");
       narr.println("  " + npopRes.toString());
       // Run the upper bound.
       FredOutVal npopUpRes = npop.upperBound();
       log.append ("The result from npopCI upper bound:");
       log.append (System.getProperty ("line.separator"));
       log.append ("  " + npopUpRes.toString ());
       log.append (System.getProperty ("line.separator"));
       narr.println("The result from npopCI upper bound:");
       narr.println("  " + npopUpRes.toString());
       npopConfidenceInterval = new int[2];
       npopConfidenceInterval[0] = npopRes.getNpop();
       npopConfidenceInterval[1] = npopUpRes.getNpop();
       return null;
    }

    /**
     *  Gets the user's approval to re-run hillclimbing on a new value before continuing.
     *
     *  @param newValue The new value to run hillclimbing on.
     *  @return True if the user selectes "YES" otherwise false.
     */
    private boolean userApproval(FredOutVal newValue) {
        String prompt = "";
        prompt += "The program has found a value from the confidence interval that";
        prompt += System.getProperty ("line.separator");
        prompt += "has a better likelihood than the value from hill climbing.";
        prompt += System.getProperty ("line.separator");
        prompt += "The Value is:" + System.getProperty ("line.separator");
        prompt += "omega: " + newValue.getOmega ();
        prompt += System.getProperty ("line.separator");
        prompt += "sigma: " + newValue.getSigma ();
        prompt += System.getProperty ("line.separator");
        prompt += "npop: " + newValue.getNpop ();
        prompt += System.getProperty ("line.separator");
        prompt += "drift: " + newValue.getDrift ();
        prompt += System.getProperty ("line.separator");
        prompt += "Re-do hillclimbing?";
        JOptionPane useNewVal = new JOptionPane(prompt,JOptionPane.QUESTION_MESSAGE, JOptionPane.YES_NO_OPTION);
        JDialog dialog = useNewVal.createDialog(EcotypeSimulation.this, "Better parameters found");
        dialog.setDefaultCloseOperation(JDialog.DO_NOTHING_ON_CLOSE);
        dialog.setVisible(true);
        Integer selectedValue = (Integer)useNewVal.getValue();
        if (selectedValue.intValue() == JOptionPane.YES_OPTION) {
            log.append ("User elected to re-do hillclimbing on the value:");
            log.append (System.getProperty ("line.separator"));
            log.append ("  " + newValue.toString());
            log.append (System.getProperty ("line.separator"));
            narr.println("User elected to re-do hillclimbing on the value:");
            narr.println("  " + newValue.toString());
            return true;
        }
        else {
            log.append ("User elected to continue with the current value from hillclimbing");
            log.append (System.getProperty ("line.separator"));
            narr.println("User elected to continue with the current value from hillclimbing");
            return false;
        }
    }

    /**
     *  Runs hillclimbing.
     *
     *  @param value The value to run hillclimbing on.
     *  @return The optimized set of parameters from hillclimbing.
     */
    private FredOutVal hillClimbing(FredOutVal value) {
        FredOutVal bestFred = value;
        ArrayList<String> bins = values.getBins();
        int[] sequenceVals = values.getSeqVals();
        HillClimbing hillOne = new HillClimbing(bestFred, masterVariables, bins, sequenceVals, masterVariables.NUM_SUCCESSES);
        // Run and time hillclimbing.
        long startTime, stopTime, runTime;
        double hourTime;
        log.append ("Running hillclimbing...");
        log.append (System.getProperty ("line.separator"));
        narr.println("Running hillclimbing...");
        startTime = System.currentTimeMillis();
        hillOne.run();
        stopTime = System.currentTimeMillis();
        runTime = stopTime - startTime;
        hourTime = (double)runTime / 3600000;
        log.append (String.format (
            "Actual runtime was: %.2g hours.", hourTime
        ));
        log.append (System.getProperty ("line.separator"));
        FredOutVal hClimbResult = hillOne.getValue();
        String omegaResult = String.format("omega: %.3g", hClimbResult.getOmega());
        String sigmaResult = String.format("sigma: %.3g", hClimbResult.getSigma());
        String npopResult = String.format("npop: %1d", hClimbResult.getNpop());
        log.append (System.getProperty ("line.separator"));
        log.append ("The values from hill climbing:");
        log.append (System.getProperty ("line.separator"));
        log.append ("  " + omegaResult);
        log.append (System.getProperty ("line.separator"));
        log.append ("  " + sigmaResult);
        log.append (System.getProperty ("line.separator"));
        log.append ("  " + npopResult);
        log.append (System.getProperty ("line.separator"));
        narr.println();
        narr.println("The values from hill climbing: ");
        narr.println(omegaResult);
        narr.println(sigmaResult);
        narr.println(npopResult);
        hClimbResult = values.fullLike(hClimbResult);
        log.append ("The full likelihood from hill climbing:");
        log.append (System.getProperty ("line.separator"));
        log.append ("  " + hClimbResult.toString ());
        log.append (System.getProperty ("line.separator"));
        narr.println("The full likelihood from hill climbing:");
        narr.println("  " + hClimbResult.toString());
        return hClimbResult;
    }

    /**
     *  Print out the final results from hillclimbing and confidence intervals.
     */
    private void printResults() {
        if (hClimbResult != null) {
            log.append ("Final results");
            log.append (System.getProperty ("line.separator"));
            log.append ("The result from hillclimbing was:");
            log.append (System.getProperty ("line.separator"));
            log.append (hClimbResult.toString ());
            log.append (System.getProperty ("line.separator"));
            narr.println("Final results");
            narr.println("The result from hillclimbing was:");
            narr.println(hClimbResult.toString());
        }
        if (omegaConfidenceInterval != null) {
            String omega = String.format("%.3g", hClimbResult.getOmega());
            String omegaHigh;
            String omegaLow;
            if (omegaConfidenceInterval[1] == 100.0) {
                omegaHigh = ">100";
            }
            else {
                omegaHigh = String.format("%.3g", omegaConfidenceInterval[1]);
            }
            if (omegaConfidenceInterval[0] < 2e-7) {
                omegaLow = "<2e-7";
            }
            else {
                omegaLow = String.format("%.3g", omegaConfidenceInterval[0]);
            }
            log.append ("omega: " + omega + " (" + omegaLow + " to " + omegaHigh + ")");
            log.append (System.getProperty ("line.separator"));
            narr.println("omega: " + omega + " (" + omegaLow + " to " + omegaHigh + ")");
        }
        if (sigmaConfidenceInterval != null) {
            String sigma = String.format("%.3g", hClimbResult.getSigma());
            String sigmaLow;
            String sigmaHigh;
            if (sigmaConfidenceInterval[1] == 100.0) {
                sigmaHigh = ">100";
            }
            else {
                sigmaHigh = String.format("%.3g", sigmaConfidenceInterval[1]);
            }
            if (sigmaConfidenceInterval[0] < 2e-7) {
                sigmaLow = "<1e-7";
            }
            else {
                sigmaLow = String.format("%.3g", sigmaConfidenceInterval[0]);
            }
            log.append ("sigma: " + sigma + " (" + sigmaLow + " to " + sigmaHigh + ")");
            log.append (System.getProperty ("line.separator"));
            narr.println("sigma: " + sigma + " (" + sigmaLow + " to " + sigmaHigh + ")");
        }
        if (npopConfidenceInterval != null) {
            log.append ("npop: " + hClimbResult.getNpop() + " (" + npopConfidenceInterval[0] + " to " + npopConfidenceInterval[1] + ")");
            log.append (System.getProperty ("line.separator"));
            narr.println("npop: " + hClimbResult.getNpop() + " (" + npopConfidenceInterval[0] + " to " + npopConfidenceInterval[1] + ")");
        }
        if (driftConfidenceInterval != 0) {
            String driftHigh;
            if ((1 / driftConfidenceInterval) < 2e-7) {
                driftHigh = ">1e7";
            }
            else {
                driftHigh = String.format("%.3g", driftConfidenceInterval);
            }
            log.append ("drift: 0 (0, " + driftHigh + ")");
            log.append (System.getProperty ("line.separator"));
            narr.println("drift: 0 (0, " + driftHigh + ")");
        }
    }

    /**
     *  Returns a user selected criterion for sorting if selected.
     *
     *  @return -1 if auto is selected, otherwise it returns the sorting
     *  percentage that the user has selected.
     */
    private int userSortPercentage() {
        String userChoice = (String)critSelector.getSelectedItem();
        if (userChoice.equals("5x"))
            return 0;
        if (userChoice.equals("2x"))
            return 1;
        if (userChoice.equals("1.5x"))
            return 2;
        if (userChoice.equals("1.25x"))
            return 3;
        if (userChoice.equals("1.10x"))
            return 4;
        if (userChoice.equals("1.05x"))
            return 5;
        else return -1;
    }

    /**
     *  Saves the progress of the current program run to an output file
     *  Note: users should NOT save files in the main directory in order
     *  to avoid overwriting important program files
     *
     *  @param output The name of the output file to save to.
     *  @return True if there is data to save, false if not.
     */
    private boolean saveFile(File output) {
        if (values == null) {
            return false;
        }
        try {
            BufferedWriter saveOut = new BufferedWriter(new FileWriter(output));
            // Write the log.
            saveOut.write("CohanLabProg save file" + System.getProperty ("line.separator"));
            String allText = log.getText();
            saveOut.write(allText + System.getProperty ("line.separator"));
            // Write the narrator.
            saveOut.write("narr" + System.getProperty ("line.separator"));
            String allNarrText = narr.getText();
            saveOut.write(allNarrText + System.getProperty ("line.separator"));
            // Write the bins.
            saveOut.write("bins" + System.getProperty ("line.separator"));
            ArrayList<String> bins = values.getBins();
            for (int i = 0; i < bins.size(); i ++) {
                saveOut.write(bins.get(i) + System.getProperty ("line.separator"));
            }
            // Write the sequenceVals.
            saveOut.write("sequenceVals" + System.getProperty ("line.separator"));
            int[] sequenceVals = values.getSeqVals();
            saveOut.write(sequenceVals[0] + " " + sequenceVals[1] + System.getProperty ("line.separator"));
            // Write the sorting percentage.
            int sortPer = masterVariables.getSortPercentage();
            saveOut.write("sortPer" + System.getProperty ("line.separator"));
            saveOut.write(sortPer+System.getProperty ("line.separator"));
            // Write the hClimbResult.
            saveOut.write("hClimbResult" + System.getProperty ("line.separator"));
            double[] percentages = hClimbResult.getPercentages();
            String save = "";
            save += hClimbResult.getOmega() + " " + hClimbResult.getSigma() + " " + hClimbResult.getNpop();
            save += " " + hClimbResult.getDrift() + " " + percentages[0] + " " + percentages[1] + " ";
            save += percentages[2] + " " + percentages[3] + " " + percentages[4] + " " + percentages[5] + System.getProperty ("line.separator");
            saveOut.write(save);
            if (omegaConfidenceInterval != null) {
                saveOut.write("omega ci" + System.getProperty ("line.separator"));
                saveOut.write(omegaConfidenceInterval[0] + " " + omegaConfidenceInterval[1] + System.getProperty ("line.separator"));
            }
            if (sigmaConfidenceInterval != null) {
                saveOut.write("sigma ci" + System.getProperty ("line.separator"));
                saveOut.write(sigmaConfidenceInterval[0] + " " + sigmaConfidenceInterval[1] + System.getProperty ("line.separator"));
            }
            if (npopConfidenceInterval != null) {
                saveOut.write("npop ci" + System.getProperty ("line.separator"));
                saveOut.write(npopConfidenceInterval[0] + " " + npopConfidenceInterval[1] + System.getProperty ("line.separator"));
            }
            if (driftConfidenceInterval != 0) {
                saveOut.write("drift ci" + System.getProperty ("line.separator"));
                saveOut.write(driftConfidenceInterval + System.getProperty ("line.separator"));
            }
            saveOut.close();
        }
        catch (IOException e) {
            e.printStackTrace();
        }
        return true;
    }

    /**
     *  Load data from a previously saved file.
     *
     *  @param savedFile A previously saved file from this program.
     *  @return boolean True if successfully recoverd, false if not.
     */
    private boolean recoverSavedData(File savedFile) {
        try {
            BufferedReader savedIn = new BufferedReader(new FileReader(savedFile));
            String line = savedIn.readLine();
            if (! line.equals("CohanLabProg save file")) {
                return false;
            }
            log.setText("");
            line = savedIn.readLine();
            while (! (line.equals("narr"))) {
                log.append (line);
                log.append (System.getProperty ("line.separator"));
                line = savedIn.readLine();
            }
            // Delete the old narrative file and rewrite the new one.
            narr = masterVariables.clearNarrator();
            line = savedIn.readLine();
            while (! (line.equals("bins"))) {
                narr.println(line);
                line = savedIn.readLine();
            }
            // Read the bins in.
            line = savedIn.readLine();
            ArrayList<String> bins = new ArrayList<String>();
            while (! (line.equals("sequenceVals"))) {
                bins.add(line);
                line = savedIn.readLine();
            }
            // Read in the sequencevals.
            line = savedIn.readLine();
            StringTokenizer tk = new StringTokenizer(line);
            int [] sequenceVals = new int[2];
            sequenceVals[0] = (new Integer(tk.nextToken())).intValue();
            sequenceVals[1] = (new Integer(tk.nextToken())).intValue();
            // Read the sort percentage.
            line = savedIn.readLine();
            line = savedIn.readLine();
            tk = new StringTokenizer(line);
            int sortPer = (new Integer(tk.nextToken())).intValue();
            masterVariables.setSortPercentage(sortPer);
            // Read the hill climbing result.
            line = savedIn.readLine();
            line = savedIn.readLine();
            tk = new StringTokenizer(line);
            double omega = (new Double(tk.nextToken())).doubleValue();
            double sigma = (new Double(tk.nextToken())).doubleValue();
            int npop = (new Integer(tk.nextToken())).intValue();
            double drift = (new Double(tk.nextToken())).doubleValue();
            double[] percentages = new double[6];
            for (int i = 0; i < percentages.length; i ++) {
                percentages[i] = (new Double(tk.nextToken())).doubleValue();
            }
            hClimbResult = new FredOutVal(omega, sigma, npop, drift, percentages, masterVariables);
            // Create the new values object.
            values = new BinningAndFred(masterVariables, sequenceVals, bins);
            // Set ready for ci.
            readyForCI = true;
            // Read in potential confidence intervals.
            line = savedIn.readLine();
            while (line != null) {
                String ciType = line;
                line = savedIn.readLine();
                tk = new StringTokenizer(line);
                if (ciType.equals("omega ci")) {
                    omegaConfidenceInterval = new double[2];
                    omegaConfidenceInterval[0] = (new Double(tk.nextToken())).doubleValue();
                    omegaConfidenceInterval[1] = (new Double(tk.nextToken())).doubleValue();
                }
                else if (ciType.equals("sigma ci")) {
                    sigmaConfidenceInterval = new double[2];
                    sigmaConfidenceInterval[0] = (new Double(tk.nextToken())).doubleValue();
                    sigmaConfidenceInterval[1] = (new Double(tk.nextToken())).doubleValue();
                }
                else if (ciType.equals("npop ci")) {
                    npopConfidenceInterval = new int[2];
                    npopConfidenceInterval[0] = (new Integer(tk.nextToken())).intValue();
                    npopConfidenceInterval[1] = (new Integer(tk.nextToken())).intValue();
                }
                else if (ciType.equals("drift ci")) {
                    driftConfidenceInterval = (new Double(tk.nextToken())).doubleValue();
                }
                line = savedIn.readLine();
            }
            savedIn.close();
        }
        catch (IOException e) {
            e.printStackTrace();
        }
        return true;
    }

    private HelpAboutWindow helpAbout;
    private MasterVariables masterVariables;
    private String workingDirectory;
    private Execs execs;
    private NarrWriter narr;
    private File tree;
    private double[] omegaConfidenceInterval;
    private double[] sigmaConfidenceInterval;
    private double driftConfidenceInterval;
    private int[] npopConfidenceInterval;
    private BinningAndFred values;
    private FredOutVal hClimbResult;
    private Fasta inputFasta;
    private Fasta noOutgroupFasta;
    private boolean readyForCI = false;
    private boolean running = false;
    private JMenuItem userGuideWindow;
    private JMenuItem aboutWindow;
    private JMenuItem optionsWindow;
    private JMenuItem changePCRError;
    private JComboBox critSelector;
    private JLabel criterion;
    private JMenuItem demarcations;
    private JButton driftConf;
    private JMenuItem exitProgram;
    private JTextArea log;
    private JButton npopConf;
    private JButton omegaConf;
    private JMenuItem openFastaFile;
    private JMenuItem openSaveFile;
    private JButton runAll;
    private JButton runAllDemarcs;
    private JButton runToHClimb;
    private JMenuItem saveFile;
    private JButton sigmaConf;
    private JMenuItem demarcsAuto;
}
