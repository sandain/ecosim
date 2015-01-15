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
import java.util.ArrayList;

/**
 *  Selects certain genes from the original input file and bins them, then runs
 *  the npop confidence interval with that binning data.
 *
 *  @author Andrew Warner
 */
public class SelectSeqBins {

    /**
     *  Constructor for objects of class SelectSeqBins.
     *
     *  @param fastaRG The input Fasta with gaps already removed.
     *  @param sequences The list of sequences to be re-binned.
     *  @param masterVariables The master variables.
     */
    public SelectSeqBins(Fasta fastaRG, ArrayList<String> sequences, MasterVariables masterVariables) {
        this.fastaRG = fastaRG;
        this.sequences = sequences;
        this.masterVariables = masterVariables;
        narr = masterVariables.getNarrator();
        execs = masterVariables.getExecs();
        workingDirectory = masterVariables.getWorkingDirectory();
        binningFasta = new BinningFasta(masterVariables);
        binInput = new File(workingDirectory + "sequencesfasta.txt");
        binOutput = new File(workingDirectory + "output.dat");
        numbers = new File(workingDirectory + "numbers.dat");
        writeSequences(binInput);
        execs.runBinning();
        bins = binningFasta.readBins(binOutput, narr);
        //write to the narrative file
        narr.println("The selected sequences to be run: ");
        for (int j = 0; j < sequences.size(); j ++) {
            narr.println(sequences.get(j));
        }
    }

    /**
     *  This constructor just saves the listed sequences to the file "output".
     *
     *  @param fastaRG The input Fasta with gaps already removed.
     *  @param sequences The list of sequences to be re-binned.
     *  @param output The fasta output file to write to.
     *  @param masterVariables The MasterVariables.
     */
    public SelectSeqBins(Fasta fastaRG, ArrayList<String> sequences, File output, MasterVariables masterVariables) {
        this.fastaRG = fastaRG;
        this.sequences = sequences;
        this.masterVariables = masterVariables;
        narr = masterVariables.getNarrator();
        execs = masterVariables.getExecs();
        workingDirectory = masterVariables.getWorkingDirectory();
        writeSequences(output);
    }

    /**
     *  Get the bins from this run.
     *
     *  @return An ArrayList<String> of the bins from this run.
     */
    public ArrayList<String> getBins() {
        return bins;
    }

    /**
     *  Get the sequence values.
     *
     *  @return An int array of sequence values.
     */
    public int[] getSeqVals() {
        return sequenceVals;
    }

    /**
     *  Retrieve the appropriate sequence data from the input file and write
     *  it to the input for the binning program (writing numbers.dat in the process).
     *
     *  @param outFile The output file to write the sequences to.
     */
    private void writeSequences(File outFile) {
        Fasta outFasta = new Fasta(fastaRG, sequences);
        outFasta.save(outFile);
        int[] values = {outFasta.size(), outFasta.length()};
        binningFasta.writeNumbersDat(numbers, values);
        sequenceVals = values;
    }

    /**
     * see constructor documentation for descriptions of these variables.
     */
    private Fasta fastaRG;
    private File binInput;
    private File binOutput;
    private File numbers;
    private ArrayList<String> sequences;
    private NarrWriter narr;
    private ArrayList<String> bins = new ArrayList<String>();
    private int[] sequenceVals = new int[2];
    private Execs execs;
    private String workingDirectory;
    private MasterVariables masterVariables;
    private BinningFasta binningFasta;
}
