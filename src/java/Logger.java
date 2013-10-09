/*
 *    Ecotype Simulation models the sequence diversity within a bacterial clade
 *    as the evolutionary result of net ecotype formation and periodic
 *    selection, yielding a certain number of ecotypes.
 *
 *    Copyright (C) 2013  Jason M. Wood, Montana State University
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


package ecosim;

import java.util.Observable;

/**
 *  The logger class is used to display text to the user.  This class utilizes
 *  the observer pattern to update listeners when the log has changed.
 *
 *  @author Jason M. Wood
 *  @copyright GNU General Public License
 */
public class Logger extends Observable {

    /**
     *  Logger constructor.
     */
    public Logger() {
        log = "";
    }

    /**
     *  Append text to the log.
     *
     *  @param str The text to append to the log.
     */
   public void append(String str) {
        // Append the text to the log.
        log += str;
        // Mark this log as changed.
        setChanged();
        // Notify observers of the change.
        notifyObservers(str);
    }

    /**
     *  Return the length in number of characters in the log.
     *
     *  @return the length of the log.
     */
    public int length() {
        return log.length();
    }

    /**
     *  Return the text currently in the log.
     */
    public String toString() {
        return log;
    }

    private String log;

}
