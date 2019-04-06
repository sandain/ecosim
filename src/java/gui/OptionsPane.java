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
import javax.swing.JSlider;
import javax.swing.event.ChangeEvent;
import javax.swing.event.ChangeListener;

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
            "Display Options"
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
        JLabel displayLabel = new JLabel ("Demarcation Display:");
        JLayeredPane paintMethodSelector = new JLayeredPane ();
        paintMethodSelector.setLayout (new GridLayout (2,1,0,0));
        paintMethodSelector.add (bars);
        paintMethodSelector.add (triangles);
        JLayeredPane paintMethodOption = new JLayeredPane ();
        paintMethodOption.setLayout (new BorderLayout ());
        paintMethodOption.add (displayLabel, BorderLayout.NORTH);
        paintMethodOption.add (paintMethodSelector, BorderLayout.WEST);
        // Create the scale slider bar.
        JLabel scaleLabel = new JLabel ("Scale:");
        JSlider scaleSlider = new JSlider (
            JSlider.HORIZONTAL, 1000, 9000, 5000
        );
        scaleSlider.addChangeListener (new ChangeListener () {
            public void stateChanged (ChangeEvent evt) {
                JSlider scale = (JSlider)evt.getSource ();
                setScaleActionPerformed ((int)scale.getValue ());
            }
        });
        scaleSlider.setMajorTickSpacing (1000);
        scaleSlider.setMinorTickSpacing (500);
        scaleSlider.setPaintTicks (true);
        JLayeredPane scaleOption = new JLayeredPane ();
        scaleOption.setLayout (new BorderLayout ());
        scaleOption.add (scaleLabel, BorderLayout.NORTH);
        scaleOption.add (scaleSlider, BorderLayout.CENTER);
        // Add the paint method option pane to the top right pane.
        add (paintMethodOption);
        add (scaleOption);
    }

    /**
     *  The user has asked to change the paint method to bars.
     */
    private void setBarsDemarcationPaintMethodActionPerformed () {
        int bars = Demarcation.PAINT_METHOD_DEMARCATED;
        if (simulation.getDemarcationPaintMethod () == bars) return;
        simulation.setDemarcationPaintMethod (bars);
    }

    /**
     *  The user has asked to change the paint method to triangles.
     */
    private void setTrianglesDemarcationPaintMethodActionPerformed () {
        int triangles = Demarcation.PAINT_METHOD_COLLAPSED;
        if (simulation.getDemarcationPaintMethod () == triangles) return;
        simulation.setDemarcationPaintMethod (triangles);
    }

    /**
     *  The user has asked to change the scale.
     */
    private void setScaleActionPerformed (int scale) {
        simulation.setScale (scale);
    }

    private Logger log;
    private Simulation simulation;

}
