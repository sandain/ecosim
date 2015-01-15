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
import java.util.HashMap;
import java.util.Set;

public class Options {


    public Options() {
        initialize();
    }

    public String[] getTabs() {
        return tabs;
    }

    public String[] getDirectoryOptionList() {
        return directoryOptionList;        
    }

    public String[] getHelperProgramList() {
        return helperProgramList;
    }

    /**
     *  Return the specified directory option.
     *
     *  @param option One of BINARY_DIRECTORY, SCRIPT_DIRECTORY, HELP_DIRECTORY, WORKING_DIRECTORY
     */
    public String getDirectoryOption(String option) {
        String directory = directoryOptions.get(option);
        if (! directory.endsWith("/") && ! directory.endsWith("\\")) {
            directory = directory + System.getProperty("file.separator");
        }
        checkDirectory(directory);
        return directory;
    }

    /**
     *  Return the specified helper program.
     *
     *  @param program One of NJPLOT, DNAPARS, DNADIST, NEIGHBOR
     */
    public String getHelperProgram(String program) {
        return helperPrograms.get(program);
    }

    public String[] getKeys(String tab) {
        String[] keys = { "" };
        if (tab.equals(DIRECTORY_OPTIONS)) {
            keys = directoryOptionList;
        }
        else if (tab.equals(HELPER_PROGRAMS)) {
            keys = helperProgramList;
        }
        return keys;
    }

    public String getValue(String tab, String label) {
        String value = "";
        if (tab.equals(DIRECTORY_OPTIONS)) {
            value = directoryOptions.get(label);
        }
        else if (tab.equals(HELPER_PROGRAMS)) {
            value = helperPrograms.get(label);
        }
        return value;
    }

    public void setDirectoryOption(String label, String value) {
        directoryOptions.put(label, value); 
    }

    public void setHelperProgram(String label, String value) {
        helperPrograms.put(label, value);
    }

    public static String OPTIONS_FILE = "options.xml";

    public static String DIRECTORY_OPTIONS = "Directory Options";
    public static String HELPER_PROGRAMS = "Helper Programs";

    // Directory Option Labels.
    public static String BINARY_DIRECTORY = "Binary Directory";
    public static String SCRIPT_DIRECTORY = "Script Directory";
    public static String HELP_DIRECTORY = "Help Directory";
    public static String WORKING_DIRECTORY = "Working Directory";

    // Helper Program Labels
    public static String NJPLOT = "NJplot";
    public static String DNAPARS = "Phylip DNApars";
    public static String DNADIST = "Phylip DNAdist";
    public static String NEIGHBOR = "Phylip Neighbor";

    /**
     *  Initialize the data.
     */
    private void initialize() {
        directoryOptions = new HashMap<String, String>(defaultDirectoryOptions);
        helperPrograms = new HashMap<String, String>(defaultHelperPrograms);
    }

    /**
     *  Check to make sure that the given directory structure exists.  If it
     *  doesn't, attempt to create it.
     *
     *  @param directory The directory structure to check.
     *  @return True if the directory already exists, or was successfully created.  False otherwise.
     */
    private boolean checkDirectory(String directory) {
        File f = new File(directory);
        boolean exists = false;
        if (! f.exists()) {
            try {
                exists = f.mkdirs();
            }
            catch (SecurityException e) {
                e.printStackTrace();
            }
        }
        else {
            exists = true;
        }
        return exists;
    }

    private String[] tabs = {
        DIRECTORY_OPTIONS,
        HELPER_PROGRAMS,
    };

    private String[] directoryOptionList = {
        BINARY_DIRECTORY,
        SCRIPT_DIRECTORY,
        HELP_DIRECTORY,
        WORKING_DIRECTORY
    };

    private String[] helperProgramList = {
        NJPLOT,
        DNAPARS,
        DNADIST,
        NEIGHBOR
    };

    private HashMap<String, String> directoryOptions;
    private HashMap<String, String> helperPrograms;

    private static HashMap<String, String> defaultDirectoryOptions;
    static {
        defaultDirectoryOptions = new HashMap<String, String>();
        defaultDirectoryOptions.put (
            BINARY_DIRECTORY,
            System.getProperty("user.dir") + System.getProperty("file.separator") + "bin" + System.getProperty("file.separator")
        );
        defaultDirectoryOptions.put (
            SCRIPT_DIRECTORY,
            System.getProperty("user.dir") + System.getProperty("file.separator") + "scripts" + System.getProperty("file.separator")
        );
        defaultDirectoryOptions.put (
            HELP_DIRECTORY,
            System.getProperty("user.dir") + System.getProperty("file.separator") + "help" + System.getProperty("file.separator")
        );
        defaultDirectoryOptions.put (
            WORKING_DIRECTORY,
            System.getProperty("user.dir") + System.getProperty("file.separator") + "working" + System.getProperty("file.separator")
        );
    }

    private static HashMap<String, String> defaultHelperPrograms;
    static {
        defaultHelperPrograms = new HashMap<String, String>();
        defaultHelperPrograms.put (
            NJPLOT,
            "njplot"
        );
        defaultHelperPrograms.put (
            DNAPARS,
            "dnapars"
        );
        defaultHelperPrograms.put (
            DNADIST,
            "dnadist"
        );
        defaultHelperPrograms.put (
            NEIGHBOR,
            "neighbor"
        );
    }

}
