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


package ecosim.tree;

import java.io.BufferedWriter;
import java.io.IOException;
import java.io.Writer;

/**
 *  Converts a Node based tree into Newick format.
 *
 *  @author Jason M. Wood
 *  @copyright GNU General Public License
 */

public class NewickWriter extends BufferedWriter {

    public NewickWriter (Writer writer) {
        super (writer);
    }

    /**
     *  Write the tree data in Newick format.
     *
     *  @param tree The tree to write.
     */
    public void write (Tree tree) throws IOException {
        write (tree.toString () + System.getProperty ("line.separator"));
    }

}
