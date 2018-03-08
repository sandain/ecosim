/*
 *    Ecotype Simulation models the sequence diversity within a bacterial
 *    clade as the evolutionary result of net ecotype formation and periodic
 *    selection, yielding a certain number of ecotypes.
 *
 *    Copyright (C) 2015-2018  Jason M. Wood, Montana State University
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

import java.awt.Font;
import java.util.Observable;
import java.util.Observer;
import javax.swing.JScrollPane;
import javax.swing.JTextArea;

/**
 *  Create a JScrollPane to display the log.
 *
 *  @author Jason M. Wood
 *  @copyright GNU General Public License
 */
public class LoggerPane extends JScrollPane {

    /**
     *  Create a JScrollPane to display the log.
     *
     *  @param log The Logger to display.
     */
    public LoggerPane (Logger log) {
        // Setup the log text area.
        final JTextArea textArea = new JTextArea ();
        textArea.setColumns (20);
        textArea.setRows (5);
        textArea.setFont (new Font ("monospaced", Font.PLAIN, 12));
        textArea.setEditable (false);
        textArea.setDoubleBuffered (true);
        textArea.append (log.toString ());
        // Update the log text area when the log changes.
        log.addObserver (new Observer () {
            public void update (Observable o, Object str) {
                textArea.append ((String) str);
                // Auto update the caret position.
                textArea.setCaretPosition (log.length ());
                // Repaint the log area.
                textArea.repaint ();
            }
        });
        // Add the log text area to the pane.
        setViewportView (textArea);
    }

    private Logger log;

}
