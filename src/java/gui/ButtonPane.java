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

import ecosim.Logger;
import ecosim.Simulation;

import java.awt.GridLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import javax.swing.JButton;
import javax.swing.JPanel;

/**
 *  Create a JPanel to display the buttons.
 *
 *  @author Jason M. Wood
 *  @copyright GNU General Public License
 */
public class ButtonPane extends JPanel {

    /**
     *  Create a JPanel to display the buttons.
     *
     *  @param log The logger.
     *  @param simulation The simulation.
     */
    public ButtonPane (Logger log, Simulation simulation) {
        this.log = log;
        this.simulation = simulation;
        setLayout (new GridLayout (4,1,15,15));
        // Create and add the buttons to the button pane.
        JButton runHillclimbing = new JButton ();
        runHillclimbing.setText ("Run Hillclimbing");
        runHillclimbing.addActionListener (new ActionListener () {
            public void actionPerformed (ActionEvent evt) {
                runHillclimbingActionPerformed ();
            }
        });
        add (runHillclimbing);
        JButton runCI = new JButton ();
        runCI.setText ("Run Confidence Intervals");
        runCI.addActionListener (new ActionListener () {
            public void actionPerformed (ActionEvent evt) {
                runConfidenceIntervalsActionPerformed ();
            }
        });
        add (runCI);
        JButton runDemarcation = new JButton ();
        runDemarcation.setText ("Run Demarcation");
        runDemarcation.addActionListener (new ActionListener () {
            public void actionPerformed (ActionEvent evt) {
                runDemarcationActionPerformed ();
            }
        });
        add (runDemarcation);
        JButton runAll = new JButton ();
        runAll.setText ("Run Everything");
        runAll.addActionListener (new ActionListener () {
            public void actionPerformed (ActionEvent evt) {
                runAllActionPerformed ();
            }
        });
        add (runAll);
    }

    /**
     *  The user has asked to run the hillclimbing program.
     */
    private void runHillclimbingActionPerformed () {
        if (! simulation.treeLoaded ()) {
            log.append ("Please load a valid Newick formatted phylogeny first.\n");
            return;
        }
        Thread thread = new Thread () {
            public void run () {
                if (simulation.isRunning ()) {
                    log.append ("Already running...\n");
                    return;
                }
                simulation.runHillclimbing ();
            }
        };
        thread.start ();
    }

    /**
     *  The user has asked to run the confidence interval programs.
     */
    private void runConfidenceIntervalsActionPerformed () {
        Thread thread = new Thread () {
            public void run () {
                if (simulation.isRunning ()) {
                    log.append ("Already running...\n");
                    return;
                }
                if (! simulation.hillclimbHasRun ()) {
                    log.append ("Please run through hillclimbing first.\n");
                    return;
                }
                simulation.runConfidenceIntervals ();
            }
        };
        thread.start ();
    }

    /**
     *  The user has asked to run the demarcation program.
     */
    private void runDemarcationActionPerformed () {
        Thread thread = new Thread () {
            public void run () {
                if (simulation.isRunning ()) {
                    log.append ("Already running...\n");
                    return;
                }
                if (! simulation.hillclimbHasRun ()) {
                    log.append ("Please run through hillclimbing first.\n");
                    return;
                }
                simulation.runDemarcation ();
            }
        };
        thread.start ();
    }

    /**
     *  The user has asked to run all of the programs.
     */
    private void runAllActionPerformed () {
        Thread thread = new Thread () {
            public void run () {
                if (simulation.isRunning ()) {
                    log.append ("Already running...\n");
                    return;
                }
                simulation.runHillclimbing ();
                simulation.runConfidenceIntervals ();
                simulation.runDemarcation ();
            }
        };
        thread.start ();
    }

    private Logger log;
    private Simulation simulation;

}
