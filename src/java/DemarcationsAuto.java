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
import java.io.FileReader;
import java.io.FileWriter;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.IOException;
import java.io.FileNotFoundException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Scanner;
import java.util.StringTokenizer;
import java.text.StringCharacterIterator;
import java.text.CharacterIterator;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import javax.swing.GroupLayout;
import javax.swing.JFrame;
import javax.swing.JMenu;
import javax.swing.JMenuBar;
import javax.swing.JMenuItem;
import javax.swing.JScrollPane;
import javax.swing.JTextArea;
import javax.swing.JOptionPane;

/**
 *  The Demarcations Auto class.
 *
 *  @author Carlo Francisco
 */
public class DemarcationsAuto extends JFrame {

    /**
     *  Creates new form Demarcations
     *
     *  @param hClimbResult FredOutVal containing the HillClimbing result.
     *  @param values The values.
     *  @param inputFasta The input Fasta.
     *  @param treeFile The tree file.
     *  @param masterVariables The MasterVariables.
     */
    public DemarcationsAuto(FredOutVal hClimbResult, BinningAndFred values, Fasta inputFasta, File treeFile, MasterVariables masterVariables) {
        this.hClimbResult = hClimbResult;
        this.masterVariables = masterVariables;
        this.values = new BinningAndFred(masterVariables, values.getSeqVals(), values.getBins());
        this.inputFasta = inputFasta;
        recombs = new ArrayList<String>();
        ecotypes = new ArrayList<ArrayList<String>>();
        narr = masterVariables.getNarrator();
        log = new JTextArea();
        workingDirectory = masterVariables.getWorkingDirectory();
        input = new File(workingDirectory + "demarcationsIn.dat");
        output = new File(workingDirectory + "demarcationsOut.dat");
        fastaOutput = new File(workingDirectory + "fasta.dat");
        numbers = new File(workingDirectory + "numbers.dat");
        narrOut = new File(workingDirectory + "narrDemarcAuto.txt");
        // Find name of outgroup
        outgroup = inputFasta.getIds().get(0);
        BinningFasta binningFasta = new BinningFasta(masterVariables);
        binningFasta.readFile(inputFasta, fastaOutput, numbers);
        fastaRGCopy = new File(workingDirectory + "fastaRGCopy.dat");
        copyFasta();
        boolean generateTree = false;
        if (treeFile == null) {
            generateTree = true;
            treeFile = runTree();
        }
        try {
            NewickTree tree = new NewickTree(treeFile);
            if (generateTree == false) {
                findRecombs(tree);
            }
            getEcotypes(tree.getRoot());
        }
        catch (InvalidNewickException e) {
            System.err.println(e);
        }


        initComponents();
        Iterator iter = ecotypes.iterator();
        for (int i = 1; iter.hasNext(); i++) {
            log.append("Ecotype " + i + ": " + iter.next() + "\n");
        }
        if (recombs.size() > 0) {
            log.append("\nThe following were found to be recombinants and thus ignored: ");
            Iterator recombIter = recombs.iterator();
            log.append("" + recombIter.next());
            while (recombIter.hasNext()) log.append(", " + recombIter.next());
        }
        log.append("\nThe sequence " + outgroup + " was assumed to be the outgroup, and so it was ignored.");
    }

    /**
     *  Get the ecotypes.
     *
     *  @return The ecotypes.
     */
    public ArrayList<ArrayList<String>> getEcotypes() {
        return ecotypes;
    }

    /**
     *  Get the recombinants.
     *
     *  @return The recombinants.
     */
    public ArrayList<String> getRecombs() {
        return recombs;
    }

    /**
     *  This method is called from within the constructor to initialize the form.
     */
    private void initComponents() {
        JScrollPane jScrollPane3 = new JScrollPane();
        JMenuBar jMenuBar1 = new JMenuBar();
        JMenu fileMenu = new JMenu();
        exitItem = new JMenuItem();
        JMenu jMenu1 = new JMenu();
        JMenuItem jMenuItem1 = new JMenuItem();
        addWindowListener(new WindowAdapter() {
            public void windowClosing(WindowEvent we) {
                dispose();
            }
        });
        setTitle("Demarcations");
        log.setColumns(20);
        log.setEditable(false);
        log.setRows(5);
        log.setDoubleBuffered(true);
        jScrollPane3.setViewportView(log);
        fileMenu.setText("File");
        exitItem.setText("Close");
        exitItem.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                exitItemActionPerformed(evt);
            }
        });
        fileMenu.add(exitItem);
        jMenuBar1.add(fileMenu);
        jMenu1.setText("Export");
        jMenuItem1.setText("To CSV");
        jMenuItem1.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                jMenuItem1ActionPerformed(evt);
            }
        });
        jMenu1.add(jMenuItem1);
        jMenuBar1.add(jMenu1);
        setJMenuBar(jMenuBar1);
        GroupLayout layout = new GroupLayout(getContentPane());
        getContentPane().setLayout(layout);
        layout.setHorizontalGroup(
            layout.createParallelGroup(GroupLayout.Alignment.LEADING)
            .addGroup(layout.createSequentialGroup()
                .addContainerGap()
                .addComponent(jScrollPane3, GroupLayout.DEFAULT_SIZE, 670, Short.MAX_VALUE)
                .addContainerGap())
        );
        layout.setVerticalGroup(
            layout.createParallelGroup(GroupLayout.Alignment.LEADING)
            .addGroup(GroupLayout.Alignment.TRAILING, layout.createSequentialGroup()
                .addContainerGap()
                .addComponent(jScrollPane3, GroupLayout.DEFAULT_SIZE, 367, Short.MAX_VALUE)
                .addContainerGap())
        );
        pack();
    }

    private void exitItemActionPerformed(java.awt.event.ActionEvent evt) {
        if (evt.getSource() == exitItem) {
            dispose();
        }
    }

    private void jMenuItem1ActionPerformed(java.awt.event.ActionEvent evt) {
        FileChooser fc = new FileChooser("csv");
        int returnVal = fc.showSaveDialog(this);
        if (returnVal == FileChooser.APPROVE_OPTION) {
            File userFile = fc.getSelectedFile();
            String path = userFile.getPath();
            // If it doesn't already have the csv extension, add it.
            if (! ((path.substring(path.length() - 4, path.length())).equals(".csv"))) {
                userFile = new File(userFile.getPath()+".csv");
            }
            narr.println("Saving to: " + userFile.getName());
            try {
                FileWriter writer = new FileWriter(userFile);
                writer.append("Ecotype Number");
                for (int i = 1; i <= highestSequenceNum(); i++) {
                    writer.append(',' + "Sequence " + i);
                }
                Iterator<ArrayList<String>> ecotypesIterator = ecotypes.iterator();
                for (int j = 1; ecotypesIterator.hasNext(); j++) {
                    writer.append('\n');
                    writer.append("" + j);
                    ArrayList<String> currentEcotype = ecotypesIterator.next();
                    Iterator seqIterator = currentEcotype.iterator();
                    while (seqIterator.hasNext()) {
                        writer.append("," + seqIterator.next());
                    }
                }
                writer.append('\n');
                writer.append("Outgroup");
                writer.append(outgroup);
                writer.append('\n');
                writer.append("Recombinants");
                Iterator<String> iter = recombs.iterator();
                while (iter.hasNext()) {
                    writer.append("," + iter.next());
                }
                writer.flush();
                writer.close();
            }
            catch(IOException e) {
                e.printStackTrace();
            }
        }
    }

    private int highestSequenceNum() {
        int n = 0;
        Iterator<ArrayList<String>> iter = ecotypes.iterator();
        while (iter.hasNext()) {
            ArrayList<String> current = iter.next();
            if (current.size() > n) n = current.size();
        }
        return n;
    }

    private File runTree() {
        File treeFile = new File(workingDirectory + "outtree");
        try {
            // Create empty outtree.
            BufferedWriter output = new BufferedWriter(new FileWriter(treeFile));
            output.close();
        }
        catch (IOException e) {
            e.printStackTrace();
        }
        // Ask for nj or parsimony tree, only if we need to make a tree.
        if (chooseTreeType() != -1) {
            TreeFinder tf = new TreeFinder(inputFasta, treeType, masterVariables);
        }
        return treeFile;
    }

    private int chooseTreeType() {
        Object[] options = {"Dnapars", "Neighbor"};
        int useNJ = JOptionPane.showOptionDialog(DemarcationsAuto.this,
                        "Use Neighbor (neighbor joining tree) or Dnapars (parsimony tree)?",
                        "Tree Type",
                        JOptionPane.YES_NO_OPTION,
                        JOptionPane.QUESTION_MESSAGE,
                        null, options, options[0]);
        if (useNJ == 1) {
            treeType = "nj";
        }
        else {
            treeType = "pars";
        }
        return useNJ;
    }

    private void copyFasta() {
        try {
            BufferedReader input = new BufferedReader(new FileReader(workingDirectory + "fasta.dat"));
            BufferedWriter output = new BufferedWriter(new FileWriter(fastaRGCopy));
            String line = input.readLine();
            while (line != null) {
                output.write(line);
                output.newLine();
                line = input.readLine();
            }
            input.close();
            output.close();
        }
        catch (IOException e) {
            e.printStackTrace();
        }
    }

    private void findRecombs(NewickTree tree) {
        ArrayList<String> seqsInFastaFile = inputFasta.getIds();
        ArrayList<String> seqsInTree = new ArrayList<String>();
        ArrayList<NewickTreeNode> leaves = tree.getDescendants();
        for (int i = 0; i < leaves.size(); i ++) {
            seqsInTree.add(leaves.get(i).getName());
        }
        seqsInTree.removeAll(seqsInFastaFile);
        recombs = seqsInTree;
    }

    private void getEcotypes(NewickTreeNode node) {
        ArrayList<String> sequences = new ArrayList<String>();
        if (node.isLeafNode()) {
            String name = node.getName();
            if (! (name.equals(outgroup) || recombs.contains(name))) {
                sequences.add(name);
                ecotypes.add(sequences);
            }
        }
        else {
            ArrayList<NewickTreeNode> leaves = node.getDescendants();
            for (int i = 0; i < leaves.size(); i ++) {
                String name = leaves.get(i).getName();
                if (! (name.equals(outgroup) || recombs.contains(name))) {
                    sequences.add(name);
                }
            }
            if (sequences.size() == 0) {
                return;
            }
            SelectSeqBins selector = new SelectSeqBins(new Fasta(fastaRGCopy), sequences, masterVariables);
            ArrayList<String> bins = selector.getBins();
            values.setBins(bins);
            values.setSeqVals(selector.getSeqVals());
            DemarcationConfidence demarcConf = new DemarcationConfidence(hClimbResult, values, masterVariables, input, output);
            // If 1 is the lower bound of the confidence interval, add the list of sequences to the list of ecotypes
            if (demarcConf.demarcations()[1] == 1) {
                ecotypes.add(sequences);
            }
            else {
                ArrayList<NewickTreeNode> children = node.getChildren();
                for (int i = 0; i < children.size(); i ++) {
                    getEcotypes(children.get(i));
                }
            }
        }
    }

    private String treeType;
    private MasterVariables masterVariables;
    private String workingDirectory;
    private ArrayList<String> recombs;
    private String outgroup;
    private ArrayList<ArrayList<String>> ecotypes;
    private NarrWriter narr;
    private BinningAndFred values;
    private File input;
    private static File output;
    private File fastaRGCopy;
    private Fasta inputFasta;
    private File fastaOutput;
    private File numbers;
    private File narrOut;
    private FredOutVal hClimbResult;
    private JMenuItem exitItem;
    private JTextArea log;
}
