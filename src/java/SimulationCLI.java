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

import java.io.File;
import java.util.Observable;
import java.util.Observer;

/**
 *  The Ecotype %Simulation Command Line Interface (CLI).
 *
 *  @author Jason M. Wood
 *  @copyright GNU General Public License
 */
public class SimulationCLI extends Simulation {

    /**
     *  SimulationCLI constructor.  Used for the command line interface
     *  of Ecotype Simulation.
     *
     *  @param masterVariables The MasterVariables object.
     *  @param fastaFile The fasta formated sequence file.
     *  @param newickFile The newick formated tree file.
     */
    public SimulationCLI(MasterVariables masterVariables,
        File fastaFile, File newickFile) {
        super(masterVariables, fastaFile, newickFile);
        log = masterVariables.getLog();
        // Display the program name and version.
        System.out.print(String.format(
            "Ecotype Simulation %s\n\n", masterVariables.getVersion()
        ));
        // When the log changes, print the text to the terminal.
        log.addObserver(new Observer() {
            public void update(Observable o, Object str) {
                System.out.print((String)str);
            }
        });
    }

    private Logger log;

}
