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
        // Create the demarcation precision buttons.
        JRadioButton coarseScale = new JRadioButton ("Coarse-scale", false);
        coarseScale.addActionListener (new ActionListener () {
            public void actionPerformed (ActionEvent evt) {
                setCoarseScaleDemarcationActionPerformed ();
            }
        });
        JRadioButton fineScale = new JRadioButton ("Fine-scale", true);
        fineScale.addActionListener (new ActionListener () {
            public void actionPerformed (ActionEvent evt) {
                setFineScaleDemarcationActionPerformed ();
            }
        });
        ButtonGroup demarcationPrecision = new ButtonGroup ();
        demarcationPrecision.add (coarseScale);
        demarcationPrecision.add (fineScale);
        // Create the demarcation method buttons.
        JRadioButton monophyletic = new JRadioButton (
            "Require monophyletic clades", true
        );
        monophyletic.addActionListener (new ActionListener () {
            public void actionPerformed (ActionEvent evt) {
                setMonophyleticDemarcationActionPerformed ();
            }
        });
        JRadioButton polyphyletic = new JRadioButton (
            "Allow polyphyletic clades", false
        );
        polyphyletic.addActionListener (new ActionListener () {
            public void actionPerformed (ActionEvent evt) {
                setPolyphyleticDemarcationActionPerformed ();
            }
        });
        ButtonGroup demarcationMethod = new ButtonGroup ();
        demarcationMethod.add (monophyletic);
        demarcationMethod.add (polyphyletic);
        // Create the demarcation precision pane.
        JLayeredPane precisionSelector = new JLayeredPane ();
        precisionSelector.setLayout (new GridLayout (2,1,0,0));
        precisionSelector.add (fineScale);
        precisionSelector.add (coarseScale);
        JLayeredPane precision = new JLayeredPane ();
        precision.setLayout (new BorderLayout ());
        precision.add (new JLabel ("Precision:"), BorderLayout.NORTH);
        precision.add (precisionSelector, BorderLayout.WEST);
        // Create the demarcation method pane.
        JLayeredPane methodSelector = new JLayeredPane ();
        methodSelector.setLayout (new GridLayout (2,1,0,0));
        methodSelector.add (monophyletic);
        methodSelector.add (polyphyletic);
        JLayeredPane method = new JLayeredPane ();
        method.setLayout (new BorderLayout ());
        method.add (new JLabel ("Method:"), BorderLayout.NORTH);
        method.add (methodSelector, BorderLayout.WEST);
        // Add the precision and method panes to the top right pane.
        add (precision);
        add (method);
    }

    /**
     *  The user has asked to change the precision to coarse-scale.
     */
    private void setCoarseScaleDemarcationActionPerformed () {
        int precision = Demarcation.DEMARCATION_PRECISION_COARSE_SCALE;
        if (simulation.getDemarcationPrecision () == precision) return;
        simulation.setDemarcationPrecision (precision);
        log.appendln ("Demarcation precision changed to: coarse-scale.");
    }

    /**
     *  The user has asked to change the precision to fine-scale.
     */
    private void setFineScaleDemarcationActionPerformed () {
        int precision = Demarcation.DEMARCATION_PRECISION_FINE_SCALE;
        if (simulation.getDemarcationPrecision () == precision) return;
        simulation.setDemarcationPrecision (precision);
        log.appendln ("Demarcation precision changed to: fine-scale.");
    }

    /**
     *  The user has asked to change the method to monophyletic.
     */
    private void setMonophyleticDemarcationActionPerformed () {
        int method = Demarcation.DEMARCATION_METHOD_MONOPHYLY;
        if (simulation.getDemarcationMethod () == method) return;
        simulation.setDemarcationMethod (method);
        log.appendln ("Demarcation method changed to: monophyletic.");
    }

    /**
     *  The user has asked to change the method to polyphyletic.
     */
    private void setPolyphyleticDemarcationActionPerformed () {
        int method = Demarcation.DEMARCATION_METHOD_POLYPHYLY;
        if (simulation.getDemarcationMethod () == method) return;
        simulation.setDemarcationMethod (method);
        log.appendln ("Demarcation method changed to: polyphyletic.");
    }

    private Logger log;
    private Simulation simulation;

}
