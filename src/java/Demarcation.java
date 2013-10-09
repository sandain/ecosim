/*
 *    Ecotype Simulation models the sequence diversity within a bacterial clade
 *    as the evolutionary result of net ecotype formation and periodic
 *    selection, yielding a certain number of ecotypes.
 *
 *    Copyright (C) 2009       Carlo Francisco, Wesleyan University
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

import java.util.ArrayList;

/**
 *  Demarcates ecotypes based on the hillclimbing values and the phylogeny of
 *  the sequences.
 *
 *  @author Carlo Francisco
 *  @author Jason M. Wood
 *  @copyright GNU General Public License
 */
public class Demarcation implements Runnable {

    /**
     *  Creates new form Demarcations
     *
     *  @param masterVariables The MasterVariables.
     *  @param phylogeny The phylogeny data.
     *  @param binning The binning results.
     *  @param hillclimb The hillclimbing results.
     */
    public Demarcation(MasterVariables masterVariables,
        Phylogeny phylogeny, Binning binning, Hillclimb hillclimb) {
        this.masterVariables = masterVariables;
        this.phylogeny = phylogeny;
        this.binning = binning;
        this.hillclimb = hillclimb;
        hasRun = false;
        ecotypes = new ArrayList<ArrayList<String>>();
        workingDirectory = masterVariables.getWorkingDirectory();
        log = masterVariables.getLog();
    }

    /**
     *  Run the demarcation program.
     */
    public void run () {
        // Find the ecotypes.
        NewickTree tree = phylogeny.getNewickTree();
        findEcotypes(tree.getRoot(), hillclimb.getResult(), 0);
        // Set the flag stating that the demarcation program has run.
        hasRun = true;
    }

    /**
     *  Returns the demarcation result as a String.
     *
     *  @return the demarcation result.
     */
    public String toString() {
        String string = "";
        // Append the ecotypes to the string.
        for (int i = 0; i < ecotypes.size(); i++) {
            string += String.format(
                "  Ecotype %4d: %s\n", (i + 1), ecotypes.get(i)
            );
        }
        return string;
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
     *  Set the ecotypes.
     *
     *  @param ecotypes The ecotypes.
     */
    public void setEcotypes(ArrayList<ArrayList<String>> ecotypes) {
        this.ecotypes = ecotypes;
    }

    /**
     *  Returns true if demarcation has been run, false otherwise.
     *
     *  @return True if demarcation has been run, false otherwise.
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
     *  Find the ecotypes using a recursive algorithm.  If the subclade of the
     *  tree represented by node has a demarcation confidence interval with a
     *  lower bound of 1, demarcate the subclade as an ecotype.  Otherwise,
     *  recurse on the children of node.
     *
     *  @param node The current node representing the subclade.
     *  @param parentHillclimbResult The Hillclimb result of the parent node.
     *  @param iteration The number of nodes already visited.
     */
    private void findEcotypes(NewickTreeNode node,
        ParameterSet parentHillclimbResult, int iteration) {
        ArrayList<String> sequences = new ArrayList<String>();
        if (node.isLeafNode()) {
            String name = node.getName();
            if (! name.equals(phylogeny.getOutgroupIdentifier())) {
                sequences.add(name);
                ecotypes.add(sequences);
            }
        }
        else {
            ArrayList<NewickTreeNode> leaves = node.getDescendants();
            for (int i = 0; i < leaves.size(); i ++) {
                String name = leaves.get(i).getName();
                if (! name.equals(phylogeny.getOutgroupIdentifier())) {
                    sequences.add(name);
                }
            }
            if (sequences.size() == 0) {
                return;
            }
            // Append a suffix to all file names used by Demarcation.
            String suffix = "-demarcation-" + iteration;
            // Create a new Phylogeny object containing just the sequences to
            // be tested.
            Phylogeny samplePhylogeny = new Phylogeny(masterVariables, suffix);
            // Add the sequences from the subtree to the new Phylogeny object.
            for (int i = 0; i < sequences.size(); i ++) {
                String id = sequences.get(i);
                samplePhylogeny.put(
                    id,
                    phylogeny.getSequence(id)
                );
            }
            // Add the subtree to the new Phylogeny object.
            samplePhylogeny.loadTree(node.toString());
            // Save a copy of the sequence data.
            samplePhylogeny.saveFasta(
                workingDirectory + "sequence" + suffix + ".dat"
            );
            samplePhylogeny.setHasRun(true);
            // Run the binning program.
            Binning sampleBinning = new Binning(
                masterVariables, samplePhylogeny, suffix
            );
            sampleBinning.run();
            // Run the hillclimb program.
            Hillclimb sampleHillclimb = new Hillclimb(
                masterVariables, samplePhylogeny, sampleBinning,
                parentHillclimbResult, suffix
            );
            sampleHillclimb.run();
            // Run the demarcation confidence interval program.
            DemarcationConfidenceInterval demarcConf =
                new DemarcationConfidenceInterval(
                masterVariables, phylogeny, samplePhylogeny,
                sampleBinning, sampleHillclimb, suffix
            );
            demarcConf.run();
            // If 1 is the lower bound of the confidence interval, add the
            // list of sequences to the list of ecotypes
            if (demarcConf.getResult() == 1) {
                ecotypes.add(sequences);
            }
            else {
                ArrayList<NewickTreeNode> children = node.getChildren();
                for (int i = 0; i < children.size(); i ++) {
                    iteration ++;
                    findEcotypes(
                        children.get(i), sampleHillclimb.getResult(), iteration
                    );
                }
            }
        }
    }

    private boolean hasRun;
    private String workingDirectory;
    private ArrayList<ArrayList<String>> ecotypes;
    private MasterVariables masterVariables;
    private Logger log;
    private Phylogeny phylogeny;
    private Binning binning;
    private Hillclimb hillclimb;

}
