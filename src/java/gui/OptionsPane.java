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

import ecosim.Demarcation;
import ecosim.Logger;
import ecosim.Simulation;

import java.awt.BorderLayout;
import java.awt.GridLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import javax.swing.BorderFactory;
import javax.swing.ButtonGroup;
import javax.swing.JLabel;
import javax.swing.JLayeredPane;
import javax.swing.JPanel;
import javax.swing.JRadioButton;

/**
 *  Create a JPanel to display the options.
 *
 *  @author Jason M. Wood
 *  @copyright GNU General Public License
 */
public class OptionsPane extends JPanel {

    /**
     *  Create a JPanel to display the options.
     *
     *  @param log The logger.
     *  @param simulation The simulation.
     */
    public OptionsPane (Logger log, Simulation simulation) {
        this.log = log;
        this.simulation = simulation;
        setBorder (BorderFactory.createTitledBorder (
            "Demarcation Options"
        ));
        setLayout (new GridLayout (2,1,15,15));
        // Create the demarcation display option buttions.
        JRadioButton bars = new JRadioButton ("Bars", true);
        bars.addActionListener (new ActionListener () {
            public void actionPerformed (ActionEvent evt) {
                setBarsDemarcationPaintMethodActionPerformed ();
            }
        });
        JRadioButton triangles = new JRadioButton ("Triangles", false);
        triangles.addActionListener (new ActionListener () {
            public void actionPerformed (ActionEvent evt) {
                setTrianglesDemarcationPaintMethodActionPerformed ();
            }
        });
        ButtonGroup demarcationPaintMethod = new ButtonGroup ();
        demarcationPaintMethod.add (bars);
        demarcationPaintMethod.add (triangles);
        // Create the demarcation paint method options pane.
        JLayeredPane paintMethodSelector = new JLayeredPane ();
        paintMethodSelector.setLayout (new GridLayout (2,1,0,0));
        paintMethodSelector.add (bars);
        paintMethodSelector.add (triangles);
        JLayeredPane paintMethodOption = new JLayeredPane ();
        paintMethodOption.setLayout (new BorderLayout ());
        paintMethodOption.add (new JLabel ("Display Option:"), BorderLayout.NORTH);
        paintMethodOption.add (paintMethodSelector, BorderLayout.WEST);
        // Add the paint method option pane to the top right pane.
        add (paintMethodOption);
    }

    /**
     *  The user has asked to change the paint method to bars.
     */
    private void setBarsDemarcationPaintMethodActionPerformed () {
        int bars = Demarcation.PAINT_METHOD_DEMARCATED;
        if (simulation.getDemarcationPaintMethod () == bars) return;
        simulation.setDemarcationPaintMethod (bars);
        log.appendln ("Demarcation paint method changed to: bars.");
    }

    /**
     *  The user has asked to change the paint method to triangles.
     */
    private void setTrianglesDemarcationPaintMethodActionPerformed () {
        int triangles = Demarcation.PAINT_METHOD_COLLAPSED;
        if (simulation.getDemarcationPaintMethod () == triangles) return;
        simulation.setDemarcationPaintMethod (triangles);
        log.appendln ("Demarcation paint method changed to: triangles.");
    }

    private Logger log;
    private Simulation simulation;

}
