/*
 *    Ecotype Simulation models the sequence diversity within a bacterial
 *    clade as the evolutionary result of net ecotype formation and periodic
 *    selection, yielding a certain number of ecotypes.
 *
 *    Copyright (C) 2009-2019  Jason M. Wood, Montana State University
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

import java.io.File;
import java.util.Arrays;
import java.util.ArrayList;
import javax.swing.JFileChooser;
import javax.swing.filechooser.FileFilter;

/**
 *  Defines a custom file chooser.
 *
 *  @author Jason M. Wood
 *  @copyright GNU General Public License
 */
public class FileChooser extends JFileChooser {

    /**
     *  Create a default FileChooser.
     */
    public FileChooser () {
        this (System.getProperty ("user.home"), "default");
    }

    /**
     *  Create a FileChooser for the given type.
     *
     *  @param currentDirectory The current directory.
     *  @param type The type of FileChooser to create.
     */
    public FileChooser (String currentDirectory, String type) {
        setCurrentDirectory (new File (currentDirectory));
        makeFileChooser (type);
    }

    /**
     *  Make the FileChooser for the given type.
     *
     *  @param type The type of FileChooser to create.
     */
    private void makeFileChooser (String type) {
        String[] valid;
        String description;
        int mode;
        switch (type) {
            case "xml": 
                valid = xmlExtensions;
                description = xmlDescription;
                mode = FILES_ONLY;
                break;
            case "svg":
                valid = svgExtensions;
                description = svgDescription;
                mode = FILES_ONLY;
                break;
            case "csv":
                valid = csvExtensions;
                description = csvDescription;
                mode = FILES_ONLY;
                break;
            case "fasta":
                valid = fastaExtensions;
                description = fastaDescription;
                mode = FILES_ONLY;
                break;
            case "newick":
                valid = newickExtensions;
                description = newickDescription;
                mode = FILES_ONLY;
                break;
            case "directory":
                valid = defaultExtensions;
                description = directoryDescription;
                mode = DIRECTORIES_ONLY;
                break;
            default:
                valid = defaultExtensions;
                description = defaultDescription;
                mode = FILES_ONLY;
        }
        FileFilter ff = new FileFilter () {
            public boolean accept (File f) {
                // Return true if the file is a directory.
                if (f.isDirectory ()) return true;
                // Find the extension of the file.
                String extension = null;
                int i = f.getName ().lastIndexOf ('.');
                if (i > 0 && i < f.getName ().length () - 1) {
                    extension = f.getName ().substring (i + 1).toLowerCase ();
                }
                // Return false if the extension was not found.
                if (extension == null) return false;
                // Transform the valid options array into an ArrayList.
                ArrayList<String> validArray = new ArrayList<String> (
                    Arrays.asList (valid)
                );
                // Return true if the valid options array contains a wildcard.
                if (validArray.contains ("*")) return true;
                // Return true if the array of valid options contains the
                // extension, otherwise return false.
                return validArray.contains (extension);
            }
            public String getDescription () {
                return description;
            }
        };
        addChoosableFileFilter (ff);
        setFileFilter (ff);
        setFileSelectionMode (mode);
    }

    private static final String[] defaultExtensions = {
        "*"
    };

    private static final String[] xmlExtensions = {
        "xml"
    };

    private static final String[] svgExtensions = {
        "svg"
    };

    private static final String[] csvExtensions = {
        "csv"
    };

    private static final String[] fastaExtensions = {
        "fa", "faa", "fas", "fasta", "fna", "fsa", "fst", "fta", "mpfa", "txt"
    };

    private static final String[] newickExtensions = {
        "nwk", "newick", "phy", "txt"
    };

    private static final String defaultDescription = "Default";
    private static final String directoryDescription = "Directories";
    private static final String xmlDescription = "XML Project Files";
    private static final String svgDescription = "SVG Image Files";
    private static final String csvDescription = "CSV Save Files";
    private static final String fastaDescription = "FASTA Files";
    private static final String newickDescription = "Newick Tree Files";

}
