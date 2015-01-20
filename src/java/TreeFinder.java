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

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.StringTokenizer;

/**
 *  Creates a Newick tree using Phylip's parsimony or neighbor joining methods, and displays the resulting tree with NJPlot.
 * 
 *  @author Andrew Warner
 */
public class TreeFinder {

    /**
     *  Creates a Newick tree using Phylip's parsimony or neighbor joining methods, and displays the resulting tree with NJPlot.
     * 
     *  @param inFastaNRG The input fasta File with all the sequences, the gaps in the sequences of this file are not removed.
     *  @param type The String containing the type of tree, either "pars" for parsimony or "nj" for neighbor joining.
     *  @param masterVariables The MasterVariables object.
     */
    public TreeFinder(Fasta inFastaNRG, String type, MasterVariables masterVariables) {
        Execs execs = masterVariables.getExecs();
        String workingDirectory = masterVariables.getWorkingDirectory();
        File infile = new File (workingDirectory + "infile");
        File outfile = new File (workingDirectory + "outfile");
        String outtree = workingDirectory + "outtree";
        File numbers = new File (workingDirectory + "numbers.dat");
        // Create the Phylip infile.
        writeTreeInput(numbers, inFastaNRG, infile);
        // Create the tree using either the parsimony or neighbor joining method.
        if (type.equals("pars")) {
            execs.runDNAPars();
        }
        else {
            execs.runDNADist();
            outfile.renameTo(infile);
            execs.runNJ();
        }
        // Load the Phlyip outtree.
        getTree(inFastaNRG, new File(outtree));
        // Display the tree using NJPlot.
        execs.openTree(outtree);
    }

    /**
     *  Get the values.
     *
     *  @return ArrayList<String> The values.
     */
    public ArrayList<String> getValues() {
        return values;
    }

    /**
     *  Gets the tree from the corresponding newick tree file given and 
     *  adds all sequences in order to the values ArrayList.
     *
     *  @param fasta The Fasta containing the sequence data from the tree.
     *  @param input The File containing the input tree.
     */
    private void getTree(Fasta fasta, File input) {
        try {
            NewickTree tree = new NewickTree(input);
            ArrayList<NewickTreeNode> leaves = tree.getDescendants();
            ArrayList<String> ids = fasta.getIds();
            values = new ArrayList<String>();
            for (int i = 0; i < leaves.size(); i ++) {
                int index = new Integer(leaves.get(i).getName()).intValue();
                leaves.get(i).setName(ids.get(index));
                values.add(ids.get(index));
            }
            tree.save(input);
        }
        catch (InvalidNewickException e) {
            System.err.println(e);
        }
    }

    /**
     *  Writes an infile for the Phylip programs.
     *
     *  @param numbers The File containing the sequence values for the tree.
     *  @param fasta The Fasta containing the sequence data for the tree.
     *  @param infile The File to write to.
     */
    private void writeTreeInput(File numbers, Fasta fasta, File infile) {
        ArrayList<String> ids = fasta.getIds();
        int[] sequenceVals = findSeqVals(numbers);
        try {
            BufferedWriter output = new BufferedWriter(new FileWriter(infile));
            output.write(sequenceVals[0] + "    " + sequenceVals[1]);
            output.newLine();
            for (int i = 0; i < ids.size(); i ++) {
                output.write (String.format (
                    "%-10d    %s", i, fasta.get(i)
                ));
                output.newLine ();
            }
            output.close();
        }
        catch (IOException e) {
            e.printStackTrace();
        }
    }

    /**
     *  Finds the sequence values from numbers.dat, which must be formatted
     *  to have the number of sequences followed by a command followed by the
     *  length of each sequence before remove gaps on the first line.
     *
     *  @param numbers The File containing the numbers to load.
     *  @return int[] An array such that [0] is the number of sequences and [1] is the length of each sequence before removegaps.
     */
    private int[] findSeqVals(File numbers) {
        try {
            int [] sequenceVals = new int[2];
            BufferedReader input = new BufferedReader(new FileReader(numbers));
            String line = input.readLine();
            StringTokenizer tk = new StringTokenizer(line);
            String numSeqs = tk.nextToken();
            numSeqs = numSeqs.substring(0, numSeqs.length() - 1);
            sequenceVals[0] = (new Integer(numSeqs)).intValue();
            sequenceVals[1] = (new Integer(tk.nextToken())).intValue();
            input.close();
            return sequenceVals;
        }
        catch (IOException e) {
            e.printStackTrace();
        }
        return null;
    }

    /**
     * The input file for the DNAPars program
     */
    private File infile;

    /**
     * The values from the tree, in order
     */
    private ArrayList<String> values = null;

}
