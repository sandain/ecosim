/*
 *    Ecotype Simulation models the sequence diversity within a bacterial
 *    clade as the evolutionary result of net ecotype formation and periodic
 *    selection, yielding a certain number of ecotypes.
 *
 *    Copyright (C) 2013-2026  Jason M. Wood
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


package ecosim;

import java.util.List;
import java.util.concurrent.CopyOnWriteArrayList;

/**
 *  The logger class is used to display text to the user.  Listeners are
 *  passed the text each time text is appended to the log.
 *
 *  @author Jason M. Wood
 *  @copyright GNU General Public License
 */
public class Logger {

    /**
     *  Receives the text appended to a Logger.
     */
    public interface Listener {

        /**
         *  Called each time text is appended to the log.
         *
         *  @param text The text appended to the log.
         */
        void textAppended (String text);

    }

    /**
     *  Logger constructor.
     */
    public Logger () {
        log = "";
        listeners = new CopyOnWriteArrayList<Listener> ();
    }

    /**
     *  Add a listener, to be passed the text each time text is appended to
     *  the log.  Listeners are called in the order they were added.
     *
     *  @param listener The listener to add.
     */
    public void addListener (Listener listener) {
        listeners.add (listener);
    }

    /**
     *  Append text to the log.
     *
     *  @param str The text to append to the log.
     */
    public void append (String str) {
        // Append the text to the log.
        log += str;
        // Pass the text to the listeners.
        for (Listener listener: listeners) {
            listener.textAppended (str);
        }
    }

    /**
     *  Append text to the log.
     *
     *  @param str The text to append to the log.
     */
    public void appendln (String str) {
        append (str + System.getProperty ("line.separator"));
    }

    /**
     *  Append a new line to the log.
     */
    public void appendln () {
        append (System.getProperty ("line.separator"));
    }


    /**
     *  Return the length in number of characters in the log.
     *
     *  @return the length of the log.
     */
    public int length () {
        return log.length ();
    }

    /**
     *  Return the text currently in the log.
     */
    public String toString () {
        return log;
    }

    private String log;
    private List<Listener> listeners;

}
