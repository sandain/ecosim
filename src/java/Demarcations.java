/*
 *    Ecotype Simulation models the sequence diversity within a bacterial clade as
 *    the evolutionary result of net ecotype formation, periodic selection,
 *    and drift, yielding a certain number of ecotypes.
 *
 *    Copyright (C) 2009  Andrew Warner, Wesleyan University
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


package ecosim;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Vector;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import javax.swing.JOptionPane;
import javax.swing.ListModel;

/*
 *  Demarcations.java
 *
 *  @author Andrew
 */
public class Demarcations extends javax.swing.JFrame {

    /**
     *  Creates new form Demarcations.
     *
     *  @param hClimbResult The result from HillClimbing.
     *  @param values The values from BinningAndFred.
     *  @param masterVariables The MasterVariables.
     */
    public Demarcations(FredOutVal hClimbResult, BinningAndFred values, MasterVariables masterVariables) {
        this.hClimbResult = hClimbResult;
        this.currentValue = hClimbResult;
        this.masterVariables = masterVariables;
        workingDirectory = masterVariables.getWorkingDirectory();
        // Create a new binningandfred object so that we can change the bins later, and a user could still
        // get accurate results in the main paramater solutions project.
        this.values = new BinningAndFred(masterVariables, values.getSeqVals(), values.getBins());
        fastaRGCopy = new File(workingDirectory + "fastaRGCopy.dat");
        narr = new NarrWriter(new File(workingDirectory + "narrDemarc.txt"));
        initComponents();
    }

    /** 
     *  This method is called from within the constructor to initialize the form.
     */
    private void initComponents() {
        jScrollPane1 = new javax.swing.JScrollPane();
        allGenes = new javax.swing.JList();
        jScrollPane2 = new javax.swing.JScrollPane();
        selectedGenes = new javax.swing.JList();
        addButton = new javax.swing.JButton();
        removeButton = new javax.swing.JButton();
        removeAllButton = new javax.swing.JButton();
        jScrollPane3 = new javax.swing.JScrollPane();
        log = new javax.swing.JTextArea();
        runAll = new javax.swing.JButton();
        allSeqLabel = new javax.swing.JLabel();
        seqRunLabel = new javax.swing.JLabel();
        jScrollPane4 = new javax.swing.JScrollPane();
        cladeList = new javax.swing.JList();
        cladeLabel = new javax.swing.JLabel();
        defCladeButton = new javax.swing.JButton();
        retCladeButton = new javax.swing.JButton();
        runFullAnalysisButton = new javax.swing.JButton();
        removeCladeButton = new javax.swing.JButton();
        jMenuBar1 = new javax.swing.JMenuBar();
        fileMenu = new javax.swing.JMenu();
        openItem = new javax.swing.JMenuItem();
        loadTreeItem = new javax.swing.JMenuItem();
        saveFastaItem = new javax.swing.JMenuItem();
        exitItem = new javax.swing.JMenuItem();
        /**
         * need this part to get the narrator to close!
         */
        addWindowListener(new WindowAdapter() {
            public void windowClosing(WindowEvent we) {
                if (narr != null)
                    narr.close();
                dispose();
            }
        });
        setTitle("Demarcations");
        allGenes.setModel(new javax.swing.AbstractListModel() {
            String[] strings = { };
            public int getSize() { return strings.length; }
            public Object getElementAt(int i) { return strings[i]; }
        });
        jScrollPane1.setViewportView(allGenes);
        selectedGenes.setModel(new javax.swing.AbstractListModel() {
            String[] strings = {  };
            public int getSize() { return strings.length; }
            public Object getElementAt(int i) { return strings[i]; }
        });
        jScrollPane2.setViewportView(selectedGenes);
        addButton.setText("Add ->");
        addButton.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                addButtonActionPerformed(evt);
            }
        });
        removeButton.setText("Remove");
        removeButton.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                removeButtonActionPerformed(evt);
            }
        });
        removeAllButton.setText("Remove All");
        removeAllButton.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                removeAllButtonActionPerformed(evt);
            }
        });
        log.setColumns(20);
        log.setEditable(false);
        log.setRows(5);
        log.setDoubleBuffered(true);
        jScrollPane3.setViewportView(log);
        runAll.setText("Run Selected Sequences");
        runAll.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                runAllActionPerformed(evt);
            }
        });
        allSeqLabel.setText("All Sequences");
        seqRunLabel.setText("Selected Sequences");
        cladeList.setModel(new javax.swing.AbstractListModel() {
            String[] strings = {  };
            public int getSize() { return strings.length; }
            public Object getElementAt(int i) { return strings[i]; }
        });
        cladeList.setSelectionMode(javax.swing.ListSelectionModel.SINGLE_SELECTION);
        jScrollPane4.setViewportView(cladeList);
        cladeLabel.setText("Clades");
        defCladeButton.setText("Define Clade");
        defCladeButton.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                defCladeButtonActionPerformed(evt);
            }
        });
        retCladeButton.setText("Retrieve Clade");
        retCladeButton.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                retCladeButtonActionPerformed(evt);
            }
        });
        runFullAnalysisButton.setText("Run Full Analysis");
        runFullAnalysisButton.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                runFullAnalysisButtonActionPerformed(evt);
            }
        });
        removeCladeButton.setText("Remove Clade");
        removeCladeButton.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                removeCladeButtonActionPerformed(evt);
            }
        });
        fileMenu.setText("File");
        openItem.setText("Open Fasta File");
        openItem.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                openItemActionPerformed(evt);
            }
        });
        fileMenu.add(openItem);
        loadTreeItem.setText("Load Newick Tree File");
        loadTreeItem.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                loadTreeItemActionPerformed(evt);
            }
        });
        fileMenu.add(loadTreeItem);
        saveFastaItem.setText("Save Sequences as Fasta File");
        saveFastaItem.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                saveFastaItemActionPerformed(evt);
            }
        });
        fileMenu.add(saveFastaItem);
        exitItem.setText("Close");
        exitItem.addActionListener(new java.awt.event.ActionListener() {
            public void actionPerformed(java.awt.event.ActionEvent evt) {
                exitItemActionPerformed(evt);
            }
        });
        fileMenu.add(exitItem);
        jMenuBar1.add(fileMenu);
        setJMenuBar(jMenuBar1);
        javax.swing.GroupLayout layout = new javax.swing.GroupLayout(getContentPane());
        getContentPane().setLayout(layout);
        layout.setHorizontalGroup(
            layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(layout.createSequentialGroup()
                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addGroup(layout.createSequentialGroup()
                        .addGap(35, 35, 35)
                        .addComponent(allSeqLabel)
                        .addGap(231, 231, 231)
                        .addComponent(seqRunLabel)
                        .addGap(97, 97, 97)
                        .addComponent(cladeLabel))
                    .addGroup(layout.createSequentialGroup()
                        .addContainerGap()
                        .addComponent(jScrollPane3, javax.swing.GroupLayout.PREFERRED_SIZE, 449, javax.swing.GroupLayout.PREFERRED_SIZE)
                        .addGap(45, 45, 45)
                        .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.TRAILING, false)
                            .addComponent(runAll, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                            .addComponent(removeCladeButton, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                            .addComponent(retCladeButton, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                            .addComponent(runFullAnalysisButton, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)))
                    .addGroup(layout.createSequentialGroup()
                        .addContainerGap()
                        .addComponent(jScrollPane1, javax.swing.GroupLayout.PREFERRED_SIZE, 133, javax.swing.GroupLayout.PREFERRED_SIZE)
                        .addGap(43, 43, 43)
                        .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING, false)
                            .addComponent(removeAllButton, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                            .addComponent(removeButton, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                            .addComponent(addButton, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE))
                        .addGap(50, 50, 50)
                        .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                            .addGroup(layout.createSequentialGroup()
                                .addGap(10, 10, 10)
                                .addComponent(defCladeButton))
                            .addGroup(layout.createSequentialGroup()
                                .addComponent(jScrollPane2, javax.swing.GroupLayout.PREFERRED_SIZE, 134, javax.swing.GroupLayout.PREFERRED_SIZE)
                                .addGap(49, 49, 49)
                                .addComponent(jScrollPane4, javax.swing.GroupLayout.PREFERRED_SIZE, 137, javax.swing.GroupLayout.PREFERRED_SIZE)))))
                .addGap(35, 35, 35))
        );
        layout.setVerticalGroup(
            layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
            .addGroup(layout.createSequentialGroup()
                .addContainerGap()
                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.BASELINE)
                    .addComponent(seqRunLabel)
                    .addComponent(allSeqLabel)
                    .addComponent(cladeLabel))
                .addGap(16, 16, 16)
                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addGroup(layout.createSequentialGroup()
                        .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                            .addGroup(layout.createSequentialGroup()
                                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.TRAILING, false)
                                    .addComponent(jScrollPane1, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE)
                                    .addGroup(javax.swing.GroupLayout.Alignment.LEADING, layout.createSequentialGroup()
                                        .addComponent(addButton)
                                        .addGap(20, 20, 20)
                                        .addComponent(removeButton)
                                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                                        .addComponent(removeAllButton)
                                        .addGap(10, 10, 10))))
                            .addComponent(jScrollPane2, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                        .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED)
                        .addComponent(defCladeButton))
                    .addComponent(jScrollPane4, javax.swing.GroupLayout.PREFERRED_SIZE, javax.swing.GroupLayout.DEFAULT_SIZE, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addPreferredGap(javax.swing.LayoutStyle.ComponentPlacement.RELATED, javax.swing.GroupLayout.DEFAULT_SIZE, Short.MAX_VALUE)
                .addGroup(layout.createParallelGroup(javax.swing.GroupLayout.Alignment.LEADING)
                    .addGroup(layout.createSequentialGroup()
                        .addComponent(retCladeButton)
                        .addGap(17, 17, 17)
                        .addComponent(removeCladeButton)
                        .addGap(27, 27, 27)
                        .addComponent(runAll)
                        .addGap(22, 22, 22)
                        .addComponent(runFullAnalysisButton))
                    .addComponent(jScrollPane3, javax.swing.GroupLayout.PREFERRED_SIZE, 172, javax.swing.GroupLayout.PREFERRED_SIZE))
                .addContainerGap())
        );
        pack();
    }

    private void loadTreeItemActionPerformed(java.awt.event.ActionEvent evt) {
        if (evt.getSource() == loadTreeItem) {
            FileChooser fc = new FileChooser("newick");
            int returnVal = fc.showOpenDialog(this);
            if (returnVal == FileChooser.APPROVE_OPTION) {
                treeFile = fc.getSelectedFile();
                try {
                    NewickTree tree = new NewickTree(treeFile);
                }
                catch (InvalidNewickException e) {
                    log.append ("That is not a valid tree file, please choose a valid newick tree file.");
                    log.append (System.getProperty ("line.separator"));
                    System.out.println(e);
                    return;
                }
            }
            else {
                return;
            }
            String message = "You must now choose the fasta file corresponding to the tree you just loaded, continue?";
            int option = JOptionPane.showConfirmDialog(null, message);
            if (option != JOptionPane.YES_OPTION) {
                return;
            }
            // Load the fasta file.
            fc = new FileChooser("fasta");
            returnVal = fc.showOpenDialog(this);
            if (returnVal == FileChooser.APPROVE_OPTION) {
                inputFasta = new Fasta(fc.getSelectedFile());
                if (! inputFasta.isValid()) {
                    log.append ("That is not a valid fasta file, please choose a properly formatted fasta file.");
                    log.append (System.getProperty ("line.separator"));
                    return;
                }
                Thread thread = new Thread() {
                    public void run() {
                        runTree();
                    }
                };
                thread.start();
            }
        }
    }

    private void saveFastaItemActionPerformed(java.awt.event.ActionEvent evt) {
        if (evt.getSource() == saveFastaItem) {
            ArrayList<String> data = new ArrayList<String>();
            ListModel model = selectedGenes.getModel();
            for (int i = 0; i < model.getSize(); i ++) {
                data.add((String)model.getElementAt(i));
            }
            if (data.size() == 0) {
                JOptionPane.showMessageDialog(null, "Please add sequences first");
                return;
            }
            FileChooser fc = new FileChooser("fasta");
            int returnVal = fc.showSaveDialog(this);
            if (returnVal == FileChooser.APPROVE_OPTION) {
                File output = fc.getSelectedFile();
                SelectSeqBins fileMaker = new SelectSeqBins(new Fasta(fastaRGCopy), data, output, masterVariables);
                log.append ("Saved to file: " + output.getPath ());
                log.append (System.getProperty ("line.separator"));
            }
        }
    }

    private void openItemActionPerformed(java.awt.event.ActionEvent evt) {
        if (evt.getSource() == openItem) {
            FileChooser fc = new FileChooser("fasta");
            int returnVal = fc.showOpenDialog(this);
            treeFile = null;
            if (returnVal == FileChooser.APPROVE_OPTION) {
                inputFasta = new Fasta(fc.getSelectedFile());
                if (! inputFasta.isValid()) {
                    log.append ("That is not a valid fasta file, please choose a properly formatted fasta file.");
                    log.append (System.getProperty ("line.separator"));
                    return;
                }
                Thread thread = new Thread() {
                    public void run() { runTree(); }
                };
                thread.start();
            }
        }
    }

    private void exitItemActionPerformed(java.awt.event.ActionEvent evt) {
        if (evt.getSource() == exitItem) {
            dispose();
        }
    }

    private void runAllActionPerformed(java.awt.event.ActionEvent evt) {
        if (evt.getSource() == runAll) {
            Thread thread = new Thread() {
                public void run() {
                    runDemarc();
                }
            };
            thread.start();
        }
    }

    private void retCladeButtonActionPerformed(java.awt.event.ActionEvent evt) {
        if (evt.getSource() == retCladeButton) {
            int [] indices = cladeList.getSelectedIndices();
            if (indices.length == 0) {
                JOptionPane.showMessageDialog(null, "Please select a clade first");
                return;
            }
            ListModel model = cladeList.getModel();
            Vector<String> dataTemp = clades.get(model.getElementAt(indices[0]));
            Vector<String> data = new Vector<String>();
            // Copy the data from the cladeList.
            for (int i = 0; i < dataTemp.size(); i ++) {
                data.add(dataTemp.get(i));
            }
            // See if that clade has had homogeneity run on it successfully.
            String message = "Now using the following values for demarcations:";
            log.append (message);
            log.append (System.getProperty ("line.separator"));
            log.append (currentValue.toString ());
            log.append (System.getProperty ("line.separator"));
            narr.println(message);
            narr.println(currentValue.toString());
            selectedGenes.setListData(data);
            cladeList.clearSelection();
        }
    }

    private void removeCladeButtonActionPerformed(java.awt.event.ActionEvent evt) {
        if (evt.getSource() == removeCladeButton) {
            int [] indices = cladeList.getSelectedIndices();
            if (indices.length == 0) {
                JOptionPane.showMessageDialog(null, "Please select a clade first");
                return;
            }
            Vector<String> data = new Vector<String>();
            ListModel model = cladeList.getModel();
            int j;
            for (int i = 0; i < model.getSize(); i ++) {
                for (j = 0; j < indices.length; j ++) {
                    if (indices[j] == i) {
                        break;
                    }
                }
                if (j == indices.length) {
                    data.add((String)model.getElementAt(i));
                }
            }
            cladeList.setListData(data);
        }
    }

    private void runFullAnalysisButtonActionPerformed(java.awt.event.ActionEvent evt) {
        if (evt.getSource()==runFullAnalysisButton) {
            int [] indices = cladeList.getSelectedIndices();
            if (indices.length == 0) {
                JOptionPane.showMessageDialog(null, "Please select a clade first");
                return;
            }
            Thread thread = new Thread() {
                public void run() {
                    runHomoGeneity();
                }
            };
            thread.start();
        }
    }

    private void defCladeButtonActionPerformed(java.awt.event.ActionEvent evt) {
        if (evt.getSource() == defCladeButton) {
            String name = JOptionPane.showInputDialog("Please enter the name of the clade");
            Vector<String> data = new Vector<String>();
            ListModel model = cladeList.getModel();
            for (int i = 0; i < model.getSize(); i ++) {
                // Delete the old clade name if there used to be a clade with this name.
                if (! model.getElementAt(i).equals(name)) {
                    data.add((String)model.getElementAt(i));
                }
            }
            data.add(name);
            cladeList.setListData(data);
            // Add this clade to the cladelist.
            Vector<String> data2 = new Vector<String>();
            ListModel model2 = selectedGenes.getModel();
            for (int j = 0; j < model2.getSize(); j ++) {
                data2.add((String)model2.getElementAt(j));
            }
            if (clades.containsKey(name)) {
                clades.remove(name);
            }
            clades.put(name, data2);
            // Remove all data from the selected genes section.
            selectedGenes.setListData(new Vector<String>());
        }
    }

    private void removeAllButtonActionPerformed(java.awt.event.ActionEvent evt) {
        if (evt.getSource() == removeAllButton) {
            selectedGenes.setListData(new Vector<String>());
        }
    }

    private void removeButtonActionPerformed(java.awt.event.ActionEvent evt) {
        if (evt.getSource() == removeButton) {
            Vector<String> data = new Vector<String>();
            ListModel model = selectedGenes.getModel();
            int [] indices = selectedGenes.getSelectedIndices();
            int j;
            for (int i = 0; i < model.getSize(); i ++) {
                for (j = 0; j < indices.length; j ++) {
                    if (indices[j] == i) {
                        break;
                    }
                }
                if (j == indices.length) {
                    data.add((String)model.getElementAt(i));
                }
            }
            selectedGenes.setListData(data);
            allGenes.clearSelection();
        }
    }

    private void addButtonActionPerformed(java.awt.event.ActionEvent evt) {
        if (evt.getSource() == addButton) {
            Vector<String> data = new Vector<String>();
            ListModel model = selectedGenes.getModel();
            for (int i = 0; i < model.getSize(); i ++) {
                data.add((String)model.getElementAt(i));
            }
            int [] indices = allGenes.getSelectedIndices();
            ListModel model2 = allGenes.getModel();
            for (int i = 0; i < indices.length; i ++) {
                data.add((String)model2.getElementAt(indices[i]));
            }
            selectedGenes.setListData(data);
            allGenes.clearSelection();
        }
    }

    /**
     *  Runs hillclimbing
     *
     *  @param value The value to run hillclimbing on.
     *  @return FredOutVal containing the output value - the optimized set of parameters from hillclimbing.
     */
    private FredOutVal hillClimbing(FredOutVal value) {
        FredOutVal bestFred = value;
        ArrayList<String> bins = values.getBins();
        int[] sequenceVals = values.getSeqVals();
        HillClimbing hillOne = new HillClimbing(bestFred, masterVariables, bins, sequenceVals, masterVariables.NUM_SUCCESSES);
        // Run and time hillclimbing.
        long startTime, stopTime, runTime;
        double hourTime;
        log.append ("Running hillclimbing... ");
        log.append (System.getProperty ("line.separator"));
        narr.println("Running hillclimbing ");
        startTime = System.currentTimeMillis();
        hillOne.run();
        stopTime = System.currentTimeMillis();
        runTime = stopTime - startTime;
        hourTime = (double)runTime / 3600000;
        FredOutVal hClimbResult = hillOne.getValue();
        log.append (System.getProperty ("line.separator"));
        log.append("The values from hill climbing:");
        log.append (System.getProperty ("line.separator"));
        log.append("omega: " + hClimbResult.getOmega ());
        log.append (System.getProperty ("line.separator"));
        log.append("sigma: " + hClimbResult.getSigma ());
        log.append (System.getProperty ("line.separator"));
        log.append("npop: " + hClimbResult.getNpop ());
        log.append (System.getProperty ("line.separator"));
        narr.println();
        narr.println("The values from hill climbing: ");
        narr.println("omega: " + hClimbResult.getOmega());
        narr.println("sigma: " + hClimbResult.getSigma());
        narr.println("npop: " + hClimbResult.getNpop());
        hClimbResult = values.fullLike(hClimbResult);
        log.append ("The full likelihood from hill climbing: ");
        log.append (System.getProperty ("line.separator"));
        log.append (hClimbResult.toString ());
        log.append (System.getProperty ("line.separator"));
        narr.println("The full likelihood from hill climbing: ");
        narr.println(hClimbResult.toString());
        return hClimbResult;
    }

    /**
     *  Reads in the given fasta file, runs binning on it, formats it, then
     *  makes a tree out of using dnapars, and then adds all the sequences
     *  to the gui.
     */
    private void runTree() {
        // Ask for nj or parsimony tree, only if we need to make a tree.
        String type = null;
        if (treeFile == null) {
            String [] choices = new String[2];
            choices[0] = "Neighbor Joining Tree";
            choices[1] = "Parsimony Tree";
            Object returnVal = JOptionPane.showInputDialog(null, "Choose a type of tree:", "Tree Type", JOptionPane.INFORMATION_MESSAGE, null, choices, 0);
            if (returnVal == null) {
                return;
            }
            if (((String)returnVal).equals(choices[0])) {
                type = "nj";
            }
            else {
                type = "pars";
            }
        }
        int sortPer = -1;
        BinningAndFred binner = new BinningAndFred(inputFasta, masterVariables, sortPer);
        binner.runBinningOnly();
        // Check if we already have a tree file or we need to make a tree file.
        Vector<String> values = new Vector<String>();
        if (treeFile == null) {
            TreeFinder tf = new TreeFinder(inputFasta, type, masterVariables);
            values.addAll(tf.getValues());
        }
        else {
            try {
                ArrayList<NewickTreeNode> leaves = new NewickTree(treeFile).getDescendants();
                for (int i = 0; i < leaves.size(); i ++) {
                    values.add(leaves.get(i).getName());
                }
            }
            catch (InvalidNewickException e) {
                System.err.println(e);
            }
        }
        // Fill the allGenes list.
        allGenes.setListData(values);
        // Clear the other lists.
        selectedGenes.setListData(new Vector<String>());
        // Modify currentValue.
        currentValue = hClimbResult;
        // Copy the full fasta.dat file (with gaps removed) to a seperate file.
        copyFasta();
    }

    /**
     *  Copies the contents of fasta.dat to the fasta RG copy file, removing
     *  the first line in the process.
     */
    private void copyFasta() {
        try {
            BufferedReader input = new BufferedReader(new FileReader(workingDirectory + "fasta.dat"));
            BufferedWriter output = new BufferedWriter(new FileWriter(fastaRGCopy));
            String line = input.readLine();
            line = input.readLine();
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

    /**
     *  Runs the homogeneity test.  It bins the sequences in the given clade,
     *  runs the old npop confidence interval on them, finding the best npop
     *  and likelihood, and then runs hillclimbing on that value with the omega
     *  sigma and drift from the first hillclimbing.  If that comes up with a
     *  better likelihood, it makes note of it.
     */
    private void runHomoGeneity() {
        int [] indices = cladeList.getSelectedIndices();
        ListModel model = cladeList.getModel();
        String name = (String)model.getElementAt(indices[0]);
        Vector<String> data = clades.get(model.getElementAt(indices[0]));
        if (data.size() < 2) {
            log.append ("Clade " + name + " is too small, add more more sequences!");
            log.append (System.getProperty ("line.separator"));
            return;
        }
        log.append ("Running full analysis on the following sequences: ");
        log.append (System.getProperty ("line.separator"));
        for (int i = 0; i < data.size(); i ++) {
            log.append ((String)data.get(i));
            log.append (System.getProperty ("line.separator"));
        }
        ArrayList<String> sequences = new ArrayList<String>();
        for (int j = 0; j < data.size(); j ++) {
            sequences.add(data.get(j));
        }
        SelectSeqBins selector = new SelectSeqBins(new Fasta(fastaRGCopy), sequences, masterVariables);
        ArrayList<String> bins = selector.getBins();
        values.setBins(bins);
        // Set the sequencevals so that we can do hillclimbing.
        values.setSeqVals(selector.getSeqVals());
        // Run the demarc confidence interval.
        File input = new File(workingDirectory + "demarcationIn.dat");
        File output = new File(workingDirectory + "demarcationOut.dat");
        DemarcationConfidence demarcConf = new DemarcationConfidence(hClimbResult, values, masterVariables, input, output);
        int [] interval = demarcConf.demarcations();
        double bestLike = demarcConf.getBestLike();
        narr.writeInput(output);
        // Run hillClimbing on the output values.
        int sortPer = masterVariables.getSortPercentage();
        double [] percentages = new double[6];
        percentages[sortPer] = bestLike;
        int npop = interval[0]; // The best npop value.
        FredOutVal homGen = new FredOutVal(hClimbResult.getOmega(),hClimbResult.getSigma(), npop, hClimbResult.getDrift(), percentages, masterVariables);
        log.append (System.getProperty ("line.separator"));
        log.append ("Running hillclimbing with the following value:");
        log.append (System.getProperty ("line.separator"));
        log.append (homGen.toString ());
        log.append (System.getProperty ("line.separator"));
        narr.println();
        narr.println("Running hillclimbing with the following value:");
        narr.println(homGen.toString());
        /**
         * Full likelihood would go here:
         * homoGen = values.fullLike(homoGen);
         */
        homGen = hillClimbing(homGen);
        double [] hClimbPers = hClimbResult.getPercentages();
        double [] homGenPers = homGen.getPercentages();
        String message;
        if ((homGenPers[sortPer] / hClimbPers[sortPer]) >  masterVariables.HOMGEN_COEFF) {
            message = "The selected clade had a better likelihood than the whole set together.";
        }
        else {
            message = "The selected clade had a likelihood that was not significantly better than the likelihood for the whole set";
        }
        log.append (message);
        log.append (System.getProperty ("line.separator"));
        narr.println(message);
    }

    /**
     *  Runs demarctions.  It bins the selected sequences, runs the demarcations
     *  confidence interval, finds the best value from that, and then finds
     *  the confidence interval, making sure to point out if 1 is in that
     *  confidence interval or not.
     */
    private void runDemarc() {
        String message;
        ArrayList<String> selected = new ArrayList<String>();
        ListModel model = selectedGenes.getModel();
        for (int i = 0; i < model.getSize(); i ++) {
            selected.add((String)model.getElementAt(i));
        }
        if (selected.size() == 0) {
            JOptionPane.showMessageDialog(null, "Please add sequences first");
            return;
        }
        log.append ("Running demarcation on the following sequences:");
        log.append (System.getProperty ("line.separator"));
        for (int i = 0; i < selected.size(); i ++) {
            log.append (selected.get (i));
            log.append (System.getProperty ("line.separator"));
        }
        SelectSeqBins selector = new SelectSeqBins(new Fasta(fastaRGCopy), selected, masterVariables);
        ArrayList<String> bins = selector.getBins();
        values.setBins(bins);
        values.setSeqVals(selector.getSeqVals());
        // Run the demarc confidence interval.
        File input = new File(workingDirectory + "demarcationIn.dat");
        File output = new File(workingDirectory + "demarcationOut.dat");
        DemarcationConfidence demarcConf = new DemarcationConfidence(currentValue, values, masterVariables, input, output);
        int [] interval = demarcConf.demarcations();
        narr.writeInput(output);
        // Write the output to the log and narrator.
        message = "The optimal value was at " + interval[0] + " and the confidence interval stretched from " + interval[1] + " to " + interval[2];
        log.append (message);
        log.append (System.getProperty ("line.separator"));
        narr.println(message);
        if (interval[1] == 1) {
            message = "The confidence interval included 1";
            log.append (message);
            log.append (System.getProperty ("line.separator"));
            narr.println(message);
        }
    }

    private String workingDirectory;
    private MasterVariables masterVariables;
    private File treeFile;
    private File fastaRGCopy;
    private BinningAndFred values;
    private FredOutVal currentValue;
    private FredOutVal hClimbResult;
    private NarrWriter narr;
    private Fasta inputFasta;
    private HashMap<String, Vector<String>> clades = new HashMap<String, Vector<String>>();
    private javax.swing.JButton addButton;
    private javax.swing.JList allGenes;
    private javax.swing.JLabel allSeqLabel;
    private javax.swing.JLabel cladeLabel;
    private javax.swing.JList cladeList;
    private javax.swing.JButton defCladeButton;
    private javax.swing.JMenuItem exitItem;
    private javax.swing.JMenu fileMenu;
    private javax.swing.JMenuBar jMenuBar1;
    private javax.swing.JScrollPane jScrollPane1;
    private javax.swing.JScrollPane jScrollPane2;
    private javax.swing.JScrollPane jScrollPane3;
    private javax.swing.JScrollPane jScrollPane4;
    private javax.swing.JMenuItem loadTreeItem;
    private javax.swing.JTextArea log;
    private javax.swing.JMenuItem openItem;
    private javax.swing.JButton removeAllButton;
    private javax.swing.JButton removeButton;
    private javax.swing.JButton removeCladeButton;
    private javax.swing.JButton retCladeButton;
    private javax.swing.JButton runAll;
    private javax.swing.JButton runFullAnalysisButton;
    private javax.swing.JMenuItem saveFastaItem;
    private javax.swing.JList selectedGenes;
    private javax.swing.JLabel seqRunLabel;
}
