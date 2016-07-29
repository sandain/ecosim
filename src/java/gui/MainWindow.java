/*
 *    Ecotype Simulation models the sequence diversity within a bacterial
 *    clade as the evolutionary result of net ecotype formation and periodic
 *    selection, yielding a certain number of ecotypes.
 *
 *    Copyright (C) 2015-2016  Jason M. Wood, Montana State University
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

import ecosim.Demarcation;
import ecosim.Logger;
import ecosim.MasterVariables;
import ecosim.Simulation;
import ecosim.Summary;
import ecosim.tree.Tree;
import ecosim.tree.SVGPainter;

import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Dimension;
import java.awt.GridLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.awt.event.WindowEvent;
import java.awt.event.WindowAdapter;
import java.io.File;
import java.util.Observable;
import java.util.Observer;
import javax.swing.BorderFactory;
import javax.swing.JFrame;
import javax.swing.JPanel;
import javax.swing.JPopupMenu;
import javax.swing.JMenuItem;
import javax.swing.JScrollPane;
import javax.swing.JTabbedPane;
import javax.swing.JViewport;

/**
 *  Create a JFrame to display the main window.
 *
 *  @author Jason M. Wood
 *  @copyright GNU General Public License
 */
public class MainWindow extends JFrame implements Runnable {

    public MainWindow (
        Logger log, MasterVariables masterVariables, Simulation simulation
    ) {
        this.log = log;
        this.masterVariables = masterVariables;
        this.simulation = simulation;
    }

    /**
     *  Run the main GUI window of the program.
     */
    public void run () {
        // Start a thread for the Help/About window.
        HelpAboutWindow helpAbout = new HelpAboutWindow (masterVariables);
        helpAbout.run ();
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
        setVisible (true);
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
        Summary summary = simulation.getSummary ();
        JTabbedPane pane = new JTabbedPane ();
        pane.addTab ("Summary", new SummaryPane (summary));
        pane.addTab ("Log", new LoggerPane (log));
        // Watch for changes to the Summary object.
        summary.addObserver (new Observer () {
            private boolean treeDisplayed = false;
            private boolean demarcationDisplayed = false;
            private JScrollPane treeScroll = new JScrollPane ();
            private JScrollPane demarcationScroll = new JScrollPane ();
            public void update (Observable o, Object obj) {
                Summary s = (Summary)obj;
                // Update the tree scroll pane.
                Tree t = s.getTree ();
                if (t == null || ! t.isValid ()) return;
                if (! treeDisplayed) {
                    JViewport viewport = treeScroll.getViewport ();
                    viewport.setOpaque (true);
                    viewport.setBackground (Color.WHITE);
                    pane.addTab ("Phylogeny", treeScroll);
                    treeDisplayed = true;
                }
                // Create a TiledPainter to paint the tree with.
                TiledPainter treePainter = new TiledPainter ();
                // Paint the tree.
                t.paintTree (treePainter);
                treeScroll.setViewportView (treePainter);
                // Add a mouse listener to the treePainter.
                treePainter.addMouseListener (new MouseAdapter () {
                    @Override
                    public void mousePressed (MouseEvent e) {
                      if (e.isPopupTrigger ()) {
                            JPopupMenu popup = new JPopupMenu ();
                            JMenuItem menu = new JMenuItem ("Save as SVG");
                            menu.addActionListener (new ActionListener () {
                                public void actionPerformed (ActionEvent evt) {
                                    saveSVGActionPerformed (t);
                                }
                            });
                            popup.add (menu);
                            popup.show (
                                e.getComponent (), e.getX (), e.getY ()
                            );
                        }
                    }
                });
                // Update the demarcation scroll pane.
                Demarcation d = s.getDemarcation ();
                if (d == null || ! d.isValid ()) return;
                if (! demarcationDisplayed) {
                    JViewport viewport = demarcationScroll.getViewport ();
                    viewport.setOpaque (true);
                    viewport.setBackground (Color.WHITE);
                    pane.addTab ("Demarcation", demarcationScroll);
                    demarcationDisplayed = true;
                }
                // Create a TiledPainter to paint the tree with.
                TiledPainter demarcationPainter = new TiledPainter ();
                // Paint the tree.
                d.paintTree (demarcationPainter);
                demarcationScroll.setViewportView (demarcationPainter);
                // Add a mouse listener to the demarcationPainter.
                demarcationPainter.addMouseListener (new MouseAdapter () {
                    @Override
                    public void mousePressed (MouseEvent e) {
                      if (e.isPopupTrigger ()) {
                            JPopupMenu popup = new JPopupMenu ();
                            JMenuItem menu = new JMenuItem ("Save as SVG");
                            menu.addActionListener (new ActionListener () {
                                public void actionPerformed (ActionEvent evt) {
                                    saveSVGActionPerformed (d);
                                }
                            });
                            popup.add (menu);
                            popup.show (
                                e.getComponent (), e.getX (), e.getY ()
                            );
                        }
                    }
                });
                // Repaint the pane.
                pane.repaint ();
            }
        });
        summary.refreshObservers ();
        return pane;
    }

    /**
     *  Save the tree as an SVG file.
     *
     *  @param tree The tree to save as an SVG file.
     */
    private void saveSVGActionPerformed (Tree tree) {
        FileChooser fc = new FileChooser ("svg");
        int returnVal = fc.showSaveDialog (this);
        File file = fc.getSelectedFile ();
        if (returnVal != FileChooser.APPROVE_OPTION) return;
        String fname = file.getAbsolutePath ();
        if (! fname.endsWith (".svg")) file = new File (fname + ".svg");
        tree.paintTree (new SVGPainter (file));
    }

    private Logger log;
    private MasterVariables masterVariables;
    private Simulation simulation;

}
