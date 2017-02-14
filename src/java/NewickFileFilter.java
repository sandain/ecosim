/*
 *    Ecotype Simulation models the sequence diversity within a bacterial clade as
 *    the evolutionary result of net ecotype formation, periodic selection,
 *    and drift, yielding a certain number of ecotypes.
 * 
 *    Copyright (C) 2009  Fred Cohan, Wesleyan University
 *                        Carlo Francisco, Wesleyan University
 *                        Danny Krizanc, Wesleyan University
 *                        Andrew Warner, Wesleyan University
 *                        Jason Wood, Montana State University
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

import java.io.File;
import javax.swing.filechooser.FileFilter;

/**
 *  Class NewickFileFilter defines a custom FileFilter for Newick Tree files.
 *
 *  @author Jason Wood
 */
public class NewickFileFilter extends FileFilter {

    /**
     *  Accept all directories and files ending in common Newick Tree file extensions.
     *
     *  @param f File to check for the correct ending.
     *  @return True if accepted, false otherwise.
     */
    public boolean accept(File f) {
        String extension = null;
        boolean accepted = false;
        int i = f.getName().lastIndexOf('.');
        if (i > 0 && i < f.getName().length() - 1) {
            extension = f.getName().substring(i + 1).toLowerCase();
        }
        if (f.isDirectory()) {
            accepted = true;
        }
        else if (extension != null && (extension.equals("tree") || extension.equals("newick"))) {
            accepted = true;
        }
        return accepted;
    }

   /**
    *  Set the description to be displayed for this FileFilter.
    *
    *  @return String containing the description for this FileFilter.
    */
    public String getDescription() {
        return "Newick Tree Files";
    }
}
