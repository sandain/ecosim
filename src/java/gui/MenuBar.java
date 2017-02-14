/*
 *    Ecotype Simulation models the sequence diversity within a bacterial
 *    clade as the evolutionary result of net ecotype formation and periodic
 *    selection, yielding a certain number of ecotypes.
 *
 *    Copyright (C) 2015-2016  Jason M. Wood, Montana State University
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


package ecosim.gui;

import ecosim.Logger;
import ecosim.Simulation;

import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.File;

import javax.swing.JMenu;
import javax.swing.JMenuBar;
import javax.swing.JMenuItem;
import javax.swing.JOptionPane;

/**
 *  Create a JMenuBar to display the menu.
 *
 *  @author Jason M. Wood
 *  @copyright GNU General Public License
 */
public class MenuBar extends JMenuBar {

    /**
     *  Create a JMenuBar to display the menu.
     *
     *  @param log The logger.
     *  @param simulation The simulation.
     *  @param helpAbout The Help/About window.
     */
    public MenuBar (
        Logger log, Simulation simulation, HelpAboutWindow helpAbout
    ) {
        this.log = log;
        this.simulation = simulation;
        JMenuItem openSequenceFile = new JMenuItem ();
        openSequenceFile.setText ("Load Sequence File");
        openSequenceFile.addActionListener (new ActionListener () {
            public void actionPerformed (ActionEvent evt) {
                openSequenceFileActionPerformed ();
            }
        });
        JMenuItem loadProjectFile = new JMenuItem ();
        loadProjectFile.setText ("Load Project File");
        loadProjectFile.addActionListener (new ActionListener () {
            public void actionPerformed (ActionEvent evt) {
                loadProjectFileActionPerformed ();
            }
        });
        JMenuItem saveProjectFile = new JMenuItem ();
        saveProjectFile.setText ("Save Project File");
        saveProjectFile.addActionListener (new ActionListener () {
            public void actionPerformed (ActionEvent evt) {
                saveProjectFileActionPerformed ();
            }
        });
        JMenuItem exitProgram = new JMenuItem ();
        exitProgram.setText ("Exit");
        exitProgram.addActionListener (new ActionListener () {
            public void actionPerformed (ActionEvent evt) {
                simulation.exit ();
            }
        });
        // Add Menu Items to the File Menu.
        JMenu fileMenu = new JMenu ();
        fileMenu.setText ("File");
        fileMenu.add (openSequenceFile);
        fileMenu.addSeparator ();
        fileMenu.add (loadProjectFile);
        fileMenu.add (saveProjectFile);
        fileMenu.addSeparator ();
        fileMenu.addSeparator ();
        fileMenu.add (exitProgram);
        // Create Help About Menu.
        JMenu helpAboutMenu = new JMenu ();
        helpAboutMenu.setText ("Help");
        JMenuItem aboutWindow = new JMenuItem ();
        aboutWindow.setText (HelpAboutWindow.ABOUT);
        aboutWindow.addActionListener (new ActionListener () {
            public void actionPerformed (ActionEvent evt) {
                helpAbout.setVisible (HelpAboutWindow.ABOUT);
            }
        });
        JMenuItem userGuideWindow = new JMenuItem ();
        userGuideWindow.setText (HelpAboutWindow.USER_GUIDE);
        userGuideWindow.addActionListener (new ActionListener () {
            public void actionPerformed (ActionEvent evt) {
                helpAbout.setVisible (HelpAboutWindow.USER_GUIDE);
            }
        });
        JMenuItem licenseWindow = new JMenuItem ();
        licenseWindow.setText (HelpAboutWindow.LICENSE);
        licenseWindow.addActionListener (new ActionListener () {
            public void actionPerformed (ActionEvent evt) {
                helpAbout.setVisible (HelpAboutWindow.LICENSE);
            }
        });
        helpAboutMenu.add (aboutWindow);
        helpAboutMenu.add (userGuideWindow);
        helpAboutMenu.add (licenseWindow);
        // Add the menus to the MenuBar.
        add (fileMenu);
        add (helpAboutMenu);
    }

   /**
    *  The user has asked to load a sequence file.
    */
    private void openSequenceFileActionPerformed () {
        // Make sure the simulation isn't already running.
        if (simulation.isRunning ()) {
            log.append ("The simulation is currently running...\n");
            return;
        }
        // Open the fasta file chooser dialog.
        FileChooser fc = new FileChooser ("fasta");
        int returnVal = fc.showOpenDialog (this);
        File file = fc.getSelectedFile ();
        if (returnVal != FileChooser.APPROVE_OPTION) return;
        // Ask the user if they want to provide or generate a tree.
        String[] options = { "Generate", "Newick" };
        int type = 0; // Default to using the Parsimony method.
        type = JOptionPane.showOptionDialog (
            this,
            "Generate a tree with FastTree or use a Newick formatted file?",
            "Tree Type",
            JOptionPane.YES_NO_OPTION,
            JOptionPane.QUESTION_MESSAGE,
            null,
            options,
            options[type]
        );
        final File treeFile;
        switch (options[type]) {
            case "Generate":
                // Generate a tree with FastTree.
                treeFile = simulation.generateTree (file);
                break;
            case "Newick":
                // Open the newick file chooser dialog.
                fc = new FileChooser ("newick");
                returnVal = fc.showOpenDialog (this);
                if (returnVal != FileChooser.APPROVE_OPTION) return;
                treeFile = fc.getSelectedFile ();
                break;
            default:
                treeFile = null;
                break;
        }
        Thread t = new Thread (
            new Runnable () {
                public void run () {
                    if (treeFile == null) return;
                    // Load the sequence and tree files.
                    simulation.loadSequenceFile (file);
                    simulation.loadTreeFile (treeFile);
                    // Run binning and estimate the parameters.
                    simulation.runBinning ();
                    simulation.runParameterEstimate ();
                }
            }
        );
        t.start ();
    }

    /**
     *  The user has asked to load a previously saved project file.
     */
    private void loadProjectFileActionPerformed () {
        FileChooser fc = new FileChooser ("xml");
        int returnVal = fc.showOpenDialog (this);
        if (returnVal == FileChooser.APPROVE_OPTION) {
            File userFile = fc.getSelectedFile ();
            // Test that the file can be read.
            if (! userFile.canRead ()) {
                log.append ("Unable to read from file!\n");
                return;
            }
            // Load the project file.
            simulation.loadProjectFile (userFile);
        }
    }

    /**
     *  The user has asked to save a project file.
     */
    private void saveProjectFileActionPerformed () {
        FileChooser fc = new FileChooser ("xml");
        int returnVal = fc.showSaveDialog (this);
        if (returnVal == FileChooser.APPROVE_OPTION) {
            File userFile = fc.getSelectedFile ();
            String path = userFile.getPath ();
            // if it doesn't already have the xml extension, add it.
            String ext = path.substring (path.length () - 4, path.length ());
            if (! ext.equals (".xml")) {
                userFile = new File (userFile.getPath () + ".xml");
            }
            simulation.saveProjectFile (userFile);
        }
    }

    private Logger log;
    private Simulation simulation;

}
