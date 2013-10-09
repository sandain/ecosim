/*
 *    Ecotype Simulation models the sequence diversity within a bacterial clade
 *    as the evolutionary result of net ecotype formation and periodic
 *    selection, yielding a certain number of ecotypes.
 *
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

import java.io.File;
import java.util.Arrays;
import java.util.ArrayList;

/**
 *  Defines a custom file chooser.
 *
 *
 *  @author Jason M. Wood
 *  @copyright GNU General Public License
 */
public class FileChooser extends javax.swing.JFileChooser {

    public FileChooser () {
        initialize("default");
    }

    public FileChooser (String type) {
        initialize(type);
    }

    public FileChooser (File directory) {
        super(directory);
        initialize("default");
    }

    public FileChooser (String type, File directory) {
        super(directory);
        initialize(type);
    }

    private void initialize (String type) {
        int mode;
        FileFilter ff;
        if (type.equals("xml")) {
            valid = new ArrayList<String>(Arrays.asList(xmlExtensions));
            ff = new FileFilter(xmlDescription);
            mode = FILES_ONLY;
        }
        else if (type.equals("csv")) {
            valid = new ArrayList<String>(Arrays.asList(csvExtensions));
            ff = new FileFilter(csvDescription);
            mode = FILES_ONLY;
        }
        else if (type.equals("fasta")) {
            valid = new ArrayList<String>(Arrays.asList(fastaExtensions));
            ff = new FileFilter(fastaDescription);
            mode = FILES_ONLY;
        }
        else if (type.equals("newick")) {
            valid = new ArrayList<String>(Arrays.asList(newickExtensions));
            ff = new FileFilter(newickDescription);
            mode = FILES_ONLY;
        }
        else if (type.equals("directory")) {
            valid = new ArrayList<String>();
            ff = new FileFilter("Directories");
            mode = DIRECTORIES_ONLY;
        }
        else {
            valid = new ArrayList<String>();
            ff = new FileFilter("Default");
            mode = FILES_ONLY;
        }
        if (valid.size() > 0) {
          addChoosableFileFilter(ff);
          setFileFilter(ff);
        }
        setFileSelectionMode(mode);
    }

    /**
     *  Defines a custom file filter.
     */
    private class FileFilter extends javax.swing.filechooser.FileFilter {

        private FileFilter (String description) {
            this.description = description;
        }

        /**
         *  Accept a File based on the type of FileChooser and the extension of
         *  the file.  Also accepts directories.
         *
         *  @param f File to check.
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
            else if (extension != null) {
                accepted = valid.contains(extension);
            }
            return accepted;
        }

        public String getDescription() {
            return description;
        }
        private String description;
    }

    private ArrayList<String> valid;

    private String[] xmlExtensions = { "xml" };
    private String[] csvExtensions = { "csv" };
    private String[] fastaExtensions = {
        "fa", "faa", "fas", "fasta", "fna", "fsa", "fst", "fta", "mpfa", "txt"
    };
    private String[] newickExtensions = {
        "nwk", "newick", "phy", "txt"
    };

    private String xmlDescription = "XML Project Files";
    private String csvDescription = "CSV Save Files";
    private String fastaDescription = "FASTA Files";
    private String newickDescription = "Newick Tree Files";

}
