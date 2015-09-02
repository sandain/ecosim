/*
 *    Ecotype Simulation models the sequence diversity within a bacterial
 *    clade as the evolutionary result of net ecotype formation and periodic
 *    selection, yielding a certain number of ecotypes.
 *
 *    Copyright (C) 2015  Jason M. Wood, Montana State University
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


package ecosim.gui;

import ecosim.Logger;
import ecosim.MasterVariables;
import ecosim.Simulation;

import java.awt.BorderLayout;
import java.awt.Dimension;
import java.awt.GridLayout;
import java.awt.event.WindowEvent;
import java.awt.event.WindowAdapter;
import javax.swing.BorderFactory;
import javax.swing.JPanel;
import javax.swing.JFrame;
import javax.swing.JTabbedPane;

/**
 *  Create a JFrame to display the main window.
 *
 *  @author Jason M. Wood
 *  @copyright GNU General Public License
 */
public class MainWindow extends JFrame {

    public MainWindow (
        Logger log, MasterVariables masterVariables, Simulation simulation
    ) {
        this.log = log;
        this.masterVariables = masterVariables;
        this.simulation = simulation;
        // Start a thread for the Help/About window.
        HelpAboutWindow helpAbout = new HelpAboutWindow (masterVariables);
        new Thread (helpAbout).start ();
        // Setup the main window.
        setTitle ("Ecotype Simulation");
        setJMenuBar (new MenuBar (log, simulation, helpAbout));
        setLayout (new BorderLayout (15, 15));
        setMinimumSize (new Dimension (525, 475));
        setPreferredSize (new Dimension (600, 600));
        add (makeTopPane (), BorderLayout.NORTH);
        add (makeBottomPane (), BorderLayout.CENTER);
        pack ();
        // Watch for window activation and closing events.
        addWindowListener (new WindowAdapter () {
            public void windowClosing (WindowEvent we) {
                simulation.exit ();
            }
        });
    }

    /**
     *  Make the top pane.
     */
    private JPanel makeTopPane () {
        JPanel pane = new JPanel ();
        // Create a border around the top pane.
        pane.setBorder (BorderFactory.createEmptyBorder (15, 15, 15, 15));
        // Add the left and right panes.
        pane.setLayout (new GridLayout (1, 2, 15, 15));
        pane.add (new ButtonPane (log, simulation));
        pane.add (new OptionsPane (log, simulation));
        return pane;
    }

    /**
     *  Make the bottom pane.
     */
    private JTabbedPane makeBottomPane () {
        JTabbedPane pane = new JTabbedPane ();
        pane.addTab ("Summary", new SummaryPane (simulation.getSummary ()));
        pane.addTab ("Log", new LoggerPane (log));
        return pane;
    }

    private Logger log;
    private MasterVariables masterVariables;
    private Simulation simulation;

}
