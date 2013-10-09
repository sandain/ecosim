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

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.StringTokenizer;
import java.util.HashMap;

/**
 *  Object to interact with the phylogeny programs.
 *
 *  @author Jason M. Wood
 *  @copyright GNU General Public License
 */
public class Phylogeny implements Runnable {

    /**
     *  Object to interact with the phylogeny programs.
     *
     *  @param masterVariables The MasterVariables object.
     */
    public Phylogeny(MasterVariables masterVariables) {
        this(masterVariables, "");
    }

    /**
     *  Object to interact with the phylogeny programs.
     *
     *  @param masterVariables The MasterVariables object.
     *  @param suffix The suffix to attach to the end of file names.
     */
    public Phylogeny(MasterVariables masterVariables, String suffix) {
        this.masterVariables = masterVariables;
        String workingDirectory = masterVariables.getWorkingDirectory();
        sequencesFileName = workingDirectory + "sequences" + suffix + ".dat";
        numbersFileName = workingDirectory + "numbers" + suffix + ".dat";
        rgFileName = workingDirectory + "removegaps" + suffix + ".dat";
        populationFileName = workingDirectory + "population" + suffix + ".dat";
        nameofstrainsFileName = workingDirectory + "namesofstrains" + suffix + ".dat";
        pcrerrorFileName = workingDirectory + "pcrerror" + suffix + ".dat";
        correctpcrFileName = workingDirectory + "correctpcr" + suffix + ".out";

        fasta = new Fasta();
        newickTree = null;

        hasRun = false;
    }

    /**
     *  Run the phylogeny programs.
     */
    public void run() {
        Execs execs = masterVariables.getExecs();
        String workingDirectory = masterVariables.getWorkingDirectory();
        File sequencesFile = new File(sequencesFileName);
        File numbersFile = new File(numbersFileName);
        File rgFile = new File(rgFileName);
        File populationFile = new File(populationFileName);
        File nameofstrainsFile = new File(nameofstrainsFileName);
        File pcrerrorFile = new File(pcrerrorFileName);
        File correctpcrFile = new File(correctpcrFileName);
        // Verify the sequence file exists.
        if (fasta == null || ! fasta.isValid()) {
            return;
        }
        // Verify the newick tree exists.
        if (newickTree == null || ! newickTree.isValid()) {
            return;
        }
        // Output the sequences.dat and numbers.dat files to be used by the
        // removegaps program.
        fasta.save(sequencesFile);
        writeNumbersFile(numbersFile, fasta.size(), fasta.length());
        // Run the fasta file through the removegaps program.
        execs.runRemovegaps(sequencesFile, numbersFile, rgFile);


        // Run the readsynec program.
        // XXX readsynec seems to just make sequences lowercase, done in Fasta.
        execs.runReadsynec(rgFile, populationFile, nameofstrainsFile);


        // Output the pcrerror.dat file to be used by the correctpcr program.
        writePCRErrorFile(pcrerrorFile, masterVariables.getPCRError());
        // Run the correctpcr program.
        execs.runCorrectpcr(populationFile, pcrerrorFile, correctpcrFile);
        // Get the output provided by the correctpcr program.
        readCorrectPCROutputFile(correctpcrFile);


        // Set the flag stating that the phylogeny programs have been run.
        hasRun = verifyPhylogeny();
    }

    /**
     *  Returns the identifier of the outgroup.
     *
     *  @return The identifier of the outgroup.
     */
    public String getOutgroupIdentifier() {
        String identifier = "";
        if (fasta.size() > 0) {
            identifier = fasta.getIdentifier(0);
        }
        return identifier;
    }

    /**
     *  Returns the sequence of the outgroup.
     *
     *  @return The sequence of the outgroup.
     */
    public String getOutgroupSequence() {
        return fasta.getSequence(0);
    }

    /**
     *  Returns the length of the sequences being analyzed.
     *
     *  @return the length of the sequences being analyzed.
     */
    public int length() {
        return fasta.length();
    }

    /**
     *  Returns the number of environmental sequences.
     *
     *  @return The number of environmental sequences.
     */
    public int getNu() {
        return fasta.size();
    }

    /**
     *  Returns the identifier at the provided index.
     *
     *  @param index The index of the identifier needed.
     *  @return The identifier at the provided index.
     */
    public String getIdentifier(int index) {
        return fasta.getIdentifier(index);
    }

    /**
     *  Returns the sequence with the provided identifier.
     *
     *  @param id The identifier of the sequence needed.
     *  @return The sequence with the provided identifier.
     */
    public String getSequence(String id) {
        return fasta.getSequence(id);
    }

    /**
     *  Returns the sequence at the provided index.
     *
     *  @param index The index of the sequence needed.
     *  @return The sequence at the provided index.
     */
    public String getSequence(int index) {
        return fasta.getSequence(index);
    }

    /**
     *  Returns all of the identifiers.
     *
     *  @return An ArrayList<String> of the identifiers.
     */
    public ArrayList<String> getIdentifiers() {
        return fasta.getIdentifiers();
    }

    /**
     *  Returns all of the sequences.
     *
     *  @return An ArrayList<String> of the sequences.
     */
    public ArrayList<String> getSequences() {
        return fasta.getSequences();
    }

    /**
     *  Returns the NewickTree object containing the tree.
     *
     *  @return A NewickTree object containing the tree.
     */
    public NewickTree getNewickTree() {
        return newickTree;
    }

    /**
     *  Puts an ID and Sequence into the Fasta object.
     *
     *  @param id The ID of the sequence to add.
     *  @param sequence The sequence to add.
     */
    public void put (String id, String sequence) {
        fasta.put(id, sequence);
    }

    /**
     *  Save the Fasta data in this object to an Ecotype Simulation formatted
     *  Fasta file.
     *
     *  @param fileName File name to write the Fasta data to.
     *  @return True if the save was a success, False otherwise.
     */
    public boolean saveFasta (String fileName) {
        return fasta.save(fileName);
    }

    /**
     *  Save the Fasta data in this object to an Ecotype Simulation formatted
     *  Fasta file.
     *
     *  @param file File to write the Fasta data to.
     *  @return True if the save was a success, False otherwise.
     */
    public boolean saveFasta (File file) {
        return fasta.save(file);
    }

    /**
     *  Save the Newick data in this object to an Ecotype Simulation formatted
     *  Newick file.
     *
     *  @param fileName File name to write the Newick data to.
     *  @return True if the save was a success, False otherwise.
     */
    public boolean saveNewick (String fileName) {
        return newickTree.save(fileName);
    }

    /**
     *  Save the Newick data in this object to an Ecotype Simulation formatted
     *  Newick file.
     *
     *  @param file File to write the Newick data to.
     *  @return True if the save was a success, False otherwise.
     */
    public boolean saveNewick (File file) {
        return newickTree.save(file);
    }

    /**
     *  Returns true if Phylogeny has been run, false otherwise.
     *
     *  @return True if Phylogeny has been run, false otherwise.
     */
    public boolean hasRun() {
        return hasRun;
    }

    /**
     *  Changes the value of hasRun.
     *
     *  @param hasRun The new value of hasRun.
     */
    public void setHasRun(boolean hasRun) {
        this.hasRun = hasRun;
    }

    /**
     *  Load a sequence file into this Phylogeny object.
     *
     *  @param sequenceFile The sequence file to load.
     */
    public void loadSequenceFile(File sequenceFile) {
        // Load the sequence file.
        fasta = new Fasta(sequenceFile);
    }

    /**
     *  Load a Newick Tree string into this Phylogeny object.
     *
     *  @param tree The tree string to load.
     */
    public void loadTree(String tree) {
        try {
            newickTree = new NewickTree(tree);
        }
        catch (InvalidNewickException e) {
            e.printStackTrace();
        }
    }

    /**
     *  Load a Newick Tree file into this Phylogeny object.
     *
     *  @param treeFile The tree file to load.
     */
    public void loadTreeFile(File treeFile) {
        try {
            newickTree = new NewickTree(treeFile);
        }
        catch (InvalidNewickException e) {
            e.printStackTrace();
        }
    }

    /**
     *  Verify that the newick tree has the same identifiers as the
     *  sequence file.
     *
     *  @return The validity of this Phylogeny object.
     */
    private boolean verifyPhylogeny() {
        boolean verified = true;
        ArrayList<String> fastaIds = fasta.getIdentifiers();
        ArrayList<NewickTreeNode> newickNodes = newickTree.getDescendants();
        // Verify the sequences in the fasta file.
        int fastaRemoved = 0;
        for (int i = 0; i < fastaIds.size(); i ++) {
            String id = fastaIds.get(i);
            boolean found = false;
            for (int j = 0; j < newickNodes.size(); j ++) {
                if (id.equals(newickNodes.get(j).getName())) {
                    found = true;
                    break;
                }
            }
            if (! found) {
                fasta.remove(id);
                fastaRemoved ++;
            }
        }
        if (fastaRemoved > 0) {
            masterVariables.getLog().append(String.format(
                "  %d sequences in the fasta file were ignored because " +
                "they were not found in the newick file.\n",
                fastaRemoved
            ));
        }
        // Verify the sequences in the newick file.
        int newickRemoved = 0;
        for (int i = 0; i < newickNodes.size(); i ++) {
            String id = newickNodes.get(i).getName();
            boolean found = false;
            for (int j = 0; j < fastaIds.size(); j ++) {
                if (id.equals(fastaIds.get(j))) {
                    found = true;
                    break;
                }
            }
            if (! found) {
                newickTree.removeDescendant(id);
                newickRemoved ++;
            }
        }
        if (newickRemoved > 0) {
            masterVariables.getLog().append(String.format(
                "  %d sequences in the newick file were ignored because " +
                "they were not found in the fasta file.\n",
                newickRemoved
            ));
        }
        if (fastaIds.size() == 0 || fastaIds.size() != newickNodes.size()) {
            verified = false;
        }
        return verified;
    }

    /**
     *  Generate a tree using Phylip's Neibor-Joining or Parsimony methods.
     *
     *  @param method The method to use to generate the tree.
     */
    public File generateTree(String method) {
        Execs execs = masterVariables.getExecs();
        String workingDirectory = masterVariables.getWorkingDirectory();
        File infile = new File (workingDirectory + "infile");
        File outfile = new File (workingDirectory + "outfile");
        File intreeFile = new File (workingDirectory + "intree");
        String outtreeFileName = workingDirectory + "outtree";
        File outtreeFile = new File(outtreeFileName);
        // Create the Phylip infile.
        writeTreeInputFile(infile);
        // Create the tree using the parsimony or neighbor-joining methods.
        switch (method) {
            case "Parsimony":
                execs.runDNAPars();
                break;
            case "Neighbor-Joining":
                execs.runDNADist();
                infile.renameTo(new File (workingDirectory + "infile.dnadist"));
                outfile.renameTo(infile);
                execs.runNJ();
                break;
        }
        // Run retree on the tree to designate the outgroup.
        outtreeFile.renameTo(intreeFile);
        execs.runRetree(fasta.size());
        // Load the Phylip outtree.
        readTreeOutputFile(outtreeFile);
        return outtreeFile;
    }

    /**
     *  Write the numbers file.
     *
     *  @param numbers The numbers.dat file.
     *  @param size The number of environmental sequences.
     *  @param length The length of the environmental sequences.
     */
    private void writeNumbersFile(File numbers, int size, int length) {
        BufferedWriter writer = null;
        try {
            writer = new BufferedWriter(new FileWriter(numbers));
            writer.write(String.format("%d, %d\n", size, length));
        }
        catch (IOException e) {
            System.out.println("Error writing the numbers.dat file for the " +
                               "removegaps program.");
        }
        finally {
            if (writer != null) {
                try {
                    writer.close();
                }
                catch (IOException e) {
                    System.out.println("Error closing the input file for the " +
                                       "removegaps program.");
                }
            }
        }
    }

    /**
     *  Write the pcrerror file.
     *  
     *  @param PCRErrorFile The pcrerror.dat file.
     *  @param PCRError The PCR error.
     */
    private void writePCRErrorFile(File PCRErrorFile, double PCRError) {
        // Create the random number seed; an odd less than 9 digits long
        long randValue = (long)(100000000 * Math.random());
        if (randValue % 2 == 0) {
            randValue ++;
        }
        try {
            BufferedWriter writer = new BufferedWriter(
                new FileWriter(PCRErrorFile)
            );
            writer.write("" + PCRError);
            writer.newLine();
            writer.write("" + randValue);
            writer.newLine();
            writer.close();
        }
        catch (IOException e) {
            e.printStackTrace();
        }
    }

    /**
     *  Private method to read the output file from the correctpcr
     *  program.
     *
     *  @param outputFile The file to read from.
     */
    private void readCorrectPCROutputFile(File outputFile) {
        BufferedReader reader = null;
        try {
            reader = new BufferedReader(new FileReader(outputFile));
            String nextLine = reader.readLine(); // nu, lengthsequence
            int i = 0;
            nextLine = reader.readLine();
            while (nextLine != null) {
                StringTokenizer st = new StringTokenizer(nextLine);
                // Modify the sequence stored in the fasta file.
                fasta.setSequence(i, st.nextToken());
                i ++;
                nextLine = reader.readLine();
            }
        }
        catch (IOException e) {
            System.out.println("Error reading the output file from the " +
                               "correctpcr program.");
        }
        finally {
            if (reader != null) {
                try {
                    reader.close();
                }
                catch (IOException e) {
                    System.out.println("Error closing the output file from " +
                                       "the correctpcr program.");
                }
            }
        }
    }

    /**
     *  Read the tree output file.
     *
     *  @param outputFile The File containing the tree to read.
     */
    private void readTreeOutputFile(File outputFile) {
        newickTree = null;
        try {
            // Load the tree file.
            newickTree = new NewickTree(outputFile);
            // Rename the leaves to match the sequence names.
            ArrayList<NewickTreeNode> leaves = newickTree.getDescendants();
            for (int i = 0; i < leaves.size(); i ++) {
                String name = leaves.get(i).getName();
                int index = new Integer(name).intValue();
                leaves.get(i).setName(fasta.getIdentifier(index));
            }
            newickTree.save(outputFile);
        }
        catch (InvalidNewickException e) {
            e.printStackTrace();
        }
    }

    /**
     *  Writes an infile for the Phylip programs.
     *
     *  @param inputFile The File to write to.
     */
    private void writeTreeInputFile(File inputFile) {
        try {
            BufferedWriter writer = new BufferedWriter(
                new FileWriter(inputFile)
            );
            writer.write(String.format(
                "%d    %d\n",
                fasta.size(), fasta.length()
            ));
            for (int i = 0; i < fasta.size(); i ++) {
                writer.write(
                    String.format("%-10d    %s\n", i, fasta.getSequence(i))
                );
            }
            writer.close();
        }
        catch (IOException e) {
            e.printStackTrace();
        }
    }

    private MasterVariables masterVariables;
    private Fasta fasta;
    private NewickTree newickTree;

    private String sequencesFileName;
    private String numbersFileName;
    private String rgFileName;
    private String populationFileName;
    private String nameofstrainsFileName;
    private String pcrerrorFileName;
    private String correctpcrFileName;

    private boolean hasRun;

}
