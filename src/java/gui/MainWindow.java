/*
 *    Ecotype Simulation models the sequence diversity within a bacterial
 *    clade as the evolutionary result of net ecotype formation and periodic
 *    selection, yielding a certain number of ecotypes.
 *
 *    Copyright (C) 2015-2019  Jason M. Wood, Montana State University
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
import java.awt.Image;
import java.awt.Toolkit;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.MouseAdapter;
import java.awt.event.MouseEvent;
import java.awt.event.WindowEvent;
import java.awt.event.WindowAdapter;
import java.io.File;
import java.util.Arrays;
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
        // Create a list of icon images.
        setIconImages (Arrays.asList (images));
        // Setup the main window.
        setTitle ("Ecotype Simulation");
        setJMenuBar (
            new MenuBar (log, masterVariables, simulation, helpAbout)
        );
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
                // Create a pop up menu to offer a Save as SVG prompt.
                JPopupMenu treePopup = createSaveAsSVGPopupMenu (t);
                treePainter.setComponentPopupMenu (treePopup);
                // Add a mouse listener to the treePainter.
                MouseAdapter tma = createPopupMouseAdapter (treePopup);
                treePainter.addMouseListener (tma);
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
                // Create a pop up menu to offer a Save as SVG prompt.
                JPopupMenu demarcationPopup = createSaveAsSVGPopupMenu (d);
                demarcationPainter.setComponentPopupMenu (demarcationPopup);
                // Add a mouse listener to the demarcationPainter.
                MouseAdapter dma = createPopupMouseAdapter (demarcationPopup);
                demarcationPainter.addMouseListener (dma);
                // Repaint the pane.
                pane.repaint ();
            }
        });
        summary.refreshObservers ();
        return pane;
    }

    /**
     *  Create a JPopupMenu to offer a Save as SVG prompt for Tree objects.
     *
     *  @param tree The Tree to offer to save as SVG.
     */
    private JPopupMenu createSaveAsSVGPopupMenu (Tree tree) {
        JPopupMenu popup = new JPopupMenu ();
        JMenuItem menu = new JMenuItem ("Save as SVG");
        menu.addActionListener (new ActionListener () {
            @Override
            public void actionPerformed (ActionEvent evt) {
                saveSVGActionPerformed (tree);
            }
        });
        popup.add (menu);
        return popup;
    }

    /**
     *  Create a MouseAdapter for a popup menu.
     *
     *  @param popup The popup menu.
     */
    private MouseAdapter createPopupMouseAdapter (JPopupMenu popup) {
        MouseAdapter adapter = new MouseAdapter () {
            @Override
            public void mousePressed (MouseEvent e) {
                if (e.isPopupTrigger ()) {
                    popup.show (e.getComponent (), e.getX (), e.getY ());
                }
            }
        };
        return adapter;
    }

    /**
     *  Save the tree as an SVG file.
     *
     *  @param tree The tree to save as an SVG file.
     */
    private void saveSVGActionPerformed (Tree tree) {
        FileChooser fc = new FileChooser (
            masterVariables.getCurrentDirectory (), "svg"
        );
        int returnVal = fc.showSaveDialog (this);
        File file = fc.getSelectedFile ();
        if (returnVal != FileChooser.APPROVE_OPTION) return;
        masterVariables.setCurrentDirectory (file.getParent ());
        String fname = file.getAbsolutePath ();
        if (! fname.endsWith (".svg")) file = new File (fname + ".svg");
        tree.paintTree (new SVGPainter (file));
    }

    private Logger log;
    private MasterVariables masterVariables;
    private Simulation simulation;

    private Image[] images = new Image[] {
        Toolkit.getDefaultToolkit().getImage(getClass().getResource("/icons/ES.128x128.png")),
        Toolkit.getDefaultToolkit().getImage(getClass().getResource("/icons/ES.64x64.png")),
        Toolkit.getDefaultToolkit().getImage(getClass().getResource("/icons/ES.32x32.png"))
    };

}
