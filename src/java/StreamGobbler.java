/*
 *    Ecotype Simulation models the sequence diversity within a bacterial clade
 *    as the evolutionary result of net ecotype formation and periodic
 *    selection, yielding a certain number of ecotypes.
 *
 *    Copyright (C) 2009       Andrew Warner, Wesleyan University
 *    Copyright (C) 2009-2013  Jason M. Wood, Montana State University
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

import java.io.BufferedReader;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.IOException;

/**
 *  This code comes from an article on java world.  I had a lot of trouble
 *  getting windows processes to run but this article explained it.
 *  http://www.javaworld.com/javaworld/jw-12-2000/jw-1229-traps.html?page=4
 *
 *  This class basically takes in an input stream and reads and outputs it.
 *
 *  @author Andrew Warner
 *  @author Jason M. Wood
 *  @copyright GNU General Public License
 */
class StreamGobbler extends Thread {

    /**
     *  Constructor for StreamGobbler.
     *
     *  @param is The input stream to use.
     *  @param type The type of stream.
     */
    public StreamGobbler(InputStream is, String type) {
        this.is = is;
        this.type = type;
    }

    /**
     *  Start this StreamGobbler thread.
     */
    public void run() {
        try {
            InputStreamReader isr = new InputStreamReader(is);
            BufferedReader br = new BufferedReader(isr);
            String line = null;
            while ((line = br.readLine()) != null) {
                System.out.println(type + ">" + line);
            }
         }
         catch (IOException ioe) {
            ioe.printStackTrace();
         }
    }

    private InputStream is;
    private String type;
}
