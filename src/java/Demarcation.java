/*
 *    Ecotype Simulation models the sequence diversity within a bacterial
 *    clade as the evolutionary result of net ecotype formation and periodic
 *    selection, yielding a certain number of ecotypes.
 *
 *    Copyright (C) 2009       Carlo Francisco, Wesleyan University
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


package ecosim;

import ecosim.tree.Node;
import ecosim.tree.InvalidTreeException;
import ecosim.tree.Tree;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.text.NumberFormat;
import java.text.ParseException;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.Locale;
import java.util.StringTokenizer;

/**
 *  Demarcates ecotypes based on the hillclimbing values and the phylogeny of
 *  the sequences.
 *
 *  @author Carlo Francisco
 *  @author Jason M. Wood
 *  @copyright GNU General Public License
 */
public class Demarcation extends Tree {

    public static final int DEMARCATION_METHOD_MONOPHYLY = 1501;
    public static final int DEMARCATION_METHOD_PARAPHYLY = 1502;

    /**
     *  Creates new form Demarcations
     *
     *  @param mainVariables The MainVariables.
     *  @param execs The Execs object.
     *  @param nu The number of environmental sequences.
     *  @param length The length of the environmental sequences.
     *  @param outgroup The name of the outgroup.
     *  @param tree The phylogeny data.
     *  @param hclimbResult The result from hillclimbing.
     *  @param method The method to use for demarcation.
     */
    public Demarcation (MainVariables mainVariables, Execs execs,
        Integer nu, Integer length, String outgroup, Tree tree,
        ParameterSet hclimbResult, int method)
        throws InvalidTreeException {
        super (tree);
        this.mainVariables = mainVariables;
        this.execs = execs;
        this.nu = nu;
        this.length = length;
        this.outgroup = outgroup;
        this.hclimbResult = hclimbResult;
        this.method = method;
        hasRun = false;
        ecotypes = new ArrayList<ArrayList<String>> ();
        workingDirectory = mainVariables.getWorkingDirectory ();
    }

    /**
     *  Run the demarcation program.
     */
    public void run () throws InvalidTreeException {
        iteration = 0;
        // Find the ecotypes.
        findEcotypes (this);
        // Set the flag stating that the demarcation program has run.
        hasRun = true;
    }

    /**
     *  Returns the demarcation result as a String.
     *
     *  @return the demarcation result.
     */
    public String toString () {
        String string = "";
        // Append the ecotypes to the string.
        for (int i = 0; i < ecotypes.size (); i++) {
            string += String.format (
                "  Ecotype %4d: %s\n", (i + 1), ecotypes.get (i)
            );
        }
        return string;
    }

    /**
     *  Get the ecotypes.
     *
     *  @return The ecotypes.
     */
    public ArrayList<ArrayList<String>> getEcotypes () {
        return ecotypes;
    }

    /**
     *  Get the demarcation method.
     *
     *  @return The demarcation method.
     */
    public int getMethod () {
        return method;
    }

    /**
     *  Set the ecotypes.
     *
     *  @param ecotypes The ecotypes.
     */
    public void setEcotypes (ArrayList<ArrayList<String>> ecotypes) {
        for (ArrayList<String> ecotype: ecotypes) {
            addEcotype (ecotype);
        }
        hasRun = true;
    }

    /**
     *  Add an ecotype.
     *
     *  @param ecotype The ecotype to add.
     */
    public void addEcotype (ArrayList<String> ecotype) {
        Node node = lastCommonAncestor (ecotype);
        ecotypes.add (ecotype);
        String name = String.format (
            "Ecotype%04d-%.4f",
            ecotypes.size (),
            node.maximumDistanceBetweenLeafNodes ()
        );
        if (node.isLeafNode ()) {
            node.addChild (new Node (node.getName (), 0.0d));
        }
        node.setName (name);
        node.collapse ();
    }

    /**
     *  Returns true if demarcation has been run, false otherwise.
     *
     *  @return True if demarcation has been run, false otherwise.
     */
    public boolean hasRun () {
        return hasRun;
    }

    /**
     *  Changes the value of hasRun.
     *
     *  @param hasRun The new value of hasRun.
     */
    public void setHasRun (boolean hasRun) {
        this.hasRun = hasRun;
    }

    /**
     *  Find ecotypes using the provided phylogeny data.
     *
     *  @param tree The phylogeny data.
     */
    private void findEcotypes (Tree tree) throws InvalidTreeException {
        switch (method) {
            case DEMARCATION_METHOD_MONOPHYLY:
                findMonophylyEcotypes (tree.getRoot ());
                break;
            case DEMARCATION_METHOD_PARAPHYLY:
                findParaphylyEcotypes (tree);
                break;
            default:
                System.err.println ("Demarcation: Invalid method.");
                System.exit (1);
        }
    }

    /**
     *  Find paraphyletic ecotypes using the provided phylogeny data. Starting
     *  with the most divergent leaf node, work backwards in time through the
     *  tree until the subclade no longer has an optimal npop value of 1. Once
     *  the last common ancestor node of the subclade is known (where npop
     *  still equals 1), demarcate the subclade and remove the entire clade
     *  from the tree. Repeat until there are no leaf node descendants left.
     *
     *  @param tree The phylogeny data.
     */
    private void findParaphylyEcotypes (Tree tree) throws InvalidTreeException {
        // Make a copy of the tree to avoid destroying it.
        Tree sampleTree;
        try {
            sampleTree = new Tree (tree);
        }
        catch (InvalidTreeException e) {
            System.err.println ("Error creating subtree.");
            return;
        }
        // Make a list of leaf node descendants.
        ArrayList<Node> leaves = sampleTree.getDescendants ();
        // Sort the leaves by their distance using a custom Comparator.
        Comparator<Node> comparator = new Comparator<Node> () {
            public int compare (Node nodeA, Node nodeB) {
                Double a = nodeA.getDistance ();
                Double b = nodeB.getDistance ();
                return a.compareTo (b);
            }
        };
        Heapsorter<Node> sorter = new Heapsorter<Node> (comparator);
        sorter.sort (leaves);
        // Find the ecotypes.
        while (leaves.size () > 0) {
            // Get the first leaf on the list.
            Node leaf = leaves.get (0);
            // Skip the outgroup.
            if (outgroup.equals (leaf.getName ())) {
                leaves.remove (0);
                continue;
            }
            // Find the ancestor node of the leaf node whose descendants make
            // up a single ecotype.
            Node node = leaf;
            while (true) {
                Node parent = node.getParent ();
                // Exit the loop if the parent node is the root node.
                if (parent.isRootNode ()) break;
                // Predict the number of ecotypes using the parent node and
                // exit the loop if the result is greater than one.
                NpopValue result = runSample (parent);
                if (result.npop > 1L) break;
                // Move the node pointer to the parent node.
                node = parent;
            }
            // Demarcate the ecotype.
            ArrayList<String> ecotype = new ArrayList<String> ();
            if (node.isLeafNode ()) {
                // Demarcate a singleton ecotype.
                ecotype.add (node.getName ());
                // Remove the node from the tree.
                leaves.remove (node);
                sampleTree.removeDescendant (node);
            }
            else {
                // Demarcate an ecotype with multiple representatives.
                for (Node descendant: node.getDescendants ()) {
                    // Add the descendant to the ecotype.
                    ecotype.add (descendant.getName ());
                    // Remove the descendant from the tree.
                    leaves.remove (descendant);
                }
                // Remove the node from the tree.
                Node parent = node.getParent ();
                parent.removeChild (node);
            }
            ecotypes.add (ecotype);
        }
    }

    /**
     *  Find monophyletic ecotypes using the provided phylogeny data.  If the
     *  subclade of the tree represented by node has an optimal npop value of
     *  1, demarcate the subclade as an ecotype. Otherwise, recurse on the
     *  node's children.
     *
     *  @param node The current node representing the subclade.
     */
    private void findMonophylyEcotypes (Node node) throws InvalidTreeException {
        ArrayList<String> sample = new ArrayList<String> ();
        String ecotype = String.format (
            "Ecotype%04d-%.4f",
            ecotypes.size () + 1,
            node.maximumDistanceBetweenLeafNodes ()
        );
        if (node.isLeafNode ()) {
            String name = node.getName ();
            if (! name.equals (outgroup)) {
                sample.add (name);
                ecotypes.add (sample);
                node.setName (ecotype);
                node.addChild (new Node (name, 0.0d));
                node.collapse ();
            }
        }
        else {
            ArrayList<Node> leaves = node.getDescendants ();
            for (int i = 0; i < leaves.size (); i ++) {
                String name = leaves.get (i).getName ();
                if (! name.equals (outgroup)) {
                    sample.add (name);
                }
            }
            if (sample.size () == 0) {
                return;
            }
            // Predict the npop value for the sample.
            NpopValue result = runSample (node);
            // If npop = 1, demarcate the list of sequences as a new ecotype.
            if (result.npop == 1L) {
                ecotypes.add (sample);
                node.setName (ecotype);
                node.collapse ();
            }
            else {
                // Npop > 1, recurse on children nodes.
                ArrayList<Node> children = node.getChildren ();
                for (int i = 0; i < children.size (); i ++) {
                    findMonophylyEcotypes (children.get (i));
                }
            }
        }
    }

   /**
    *  A private helper method to run a sample through the demarcation
    *  program.
    *
    *  @param node The Node describing the sample to run.
    *  @return The npop value tested and its likelihood
    */
   private NpopValue runSample (Node node) throws InvalidTreeException {
        // Increment the iteration variable used in the file names.
        iteration ++;
        File inputFile = new File (
            workingDirectory + "demarcationIn-" + iteration + ".dat"
        );
        File outputFile = new File (
            workingDirectory + "demarcationOut-" + iteration + ".dat"
        );
        File newickFile = new File (
            workingDirectory + "demarcationTree-" + iteration + ".dat"
        );
        // Create a new Tree containing just the sequences to
        // be tested.
        Tree sampleTree = new Tree (node.toString ());
        sampleTree.toNewick (newickFile);
        Integer sampleNu = numberOfDescendants (node);
        // Run the binning program on the sample tree.
        Binning sampleBinning = new Binning (sampleTree);
        sampleBinning.run ();
        // Use the omega and sigma values from hillclimbing.
        Double omega = hclimbResult.getOmega ();
        Double sigma = hclimbResult.getSigma ();
        // Estimate the value of npop by multiplying the npop value found
        // in hillclimbing by the ratio of the sample size to the total
        // number of environmental sequences.
        Long npop = hclimbResult.getNpop () * sampleNu / nu;
        if (npop < 1L) {
            npop = 1L;
        }
        // Write the input values for the demarcation program.
        writeInputFile (
            inputFile, sampleBinning, sampleNu, omega, sigma, npop,
            hclimbResult.getLikelihood ()
        );
        // Run the demarcation program.
        execs.runDemarcation (inputFile, outputFile);
        // Get the output provided by the demarcation program.
        // [0] npop=1
        // [1] most likely npop
        NpopValue[] results = readOutputFile (outputFile);
        return results[1];
    }

    /**
     *  Private method to write the input file for the demarcation program.
     *
     *  @param inputFile The file to write to.
     *  @param binning The binning results.
     *  @param sampleNu The number of environmental sequences.
     *  @param omega The omega estimate.
     *  @param sigma The sigma estimate.
     *  @param npop The npop estimate.
     *  @param likelihood The likelihood of the omega, sigma, npop estimates.
     */
    private void writeInputFile (File inputFile, Binning binning,
        Integer sampleNu, Double omega, Double sigma, Long npop,
        Double likelihood) {
        ArrayList<BinLevel> bins = binning.getBins ();
        BufferedWriter writer = null;
        try {
            writer = new BufferedWriter (new FileWriter (inputFile));
            writer.write (String.format (
                Locale.US,
                "%-20d numcrit\n",
                bins.size ()
            ));
            // Output the crit levels and the number of bins.
            for (int j = 0; j < bins.size (); j ++) {
                writer.write (String.format (
                    Locale.US,
                    "%-20.6f %-20d\n",
                    bins.get (j).getCrit (),
                    bins.get (j).getLevel ()
                ));
            }
            // Write the omega value.
            writer.write (String.format (
                Locale.US,
                "%-20.5f omega\n",
                omega
            ));
            // Write the sigma value.
            writer.write (String.format (
                Locale.US,
                "%-20.5f sigma\n",
                sigma
            ));
            // Write the npop value.
            writer.write (String.format (
                Locale.US,
                "%-20d npop\n",
                npop
            ));
            // Write the step value.
            writer.write (String.format (
                Locale.US,
                "%-20d step\n",
                step
            ));
            // Write the nu value.
            writer.write (String.format (
                Locale.US,
                "%-20d nu\n",
                sampleNu
            ));
            // Write the nrep value.
            writer.write (String.format (
                Locale.US,
                "%-20d nrep\n",
                nrep
            ));
            // Create the random number seed; an odd integer less than nine
            // digits long.
            long iii = (long)(100000000 * Math.random ());
            if (iii % 2 == 0) {
                iii ++;
            }
            // Write the random number seed.
            writer.write (String.format (
                Locale.US,
                "%-20d iii (random number seed)\n",
                iii
            ));
            // Write the length of the sequences.
            writer.write (String.format (
                Locale.US,
                "%-20d lengthseq (after deleting gaps, etc.)\n",
                length
            ));
            // Write the whichavg value.
            int whichavg = mainVariables.getCriterion ();
            writer.write (String.format (
                Locale.US,
                "%-20d whichavg\n",
                whichavg
            ));
        }
        catch (IOException e) {
            System.out.println ("Error writing the input file.");
        }
        finally {
            if (writer != null) {
                try {
                    writer.close ();
                }
                catch (IOException e) {
                    System.out.println ("Error closing the input file.");
                }
            }
        }
    }

    /**
     *  Private method to read the output file from the demarcation program.
     *
     *  @param outputFile The file to read from.
     *  @return The npop values tested and their likelihood.
     */
    private NpopValue[] readOutputFile (File outputFile) {
        BufferedReader reader = null;
        NumberFormat format = NumberFormat.getInstance (Locale.US);
        NpopValue result[] = {
            new NpopValue (0L, 0.0d),
            new NpopValue (0L, 0.0d)
        };
        try {
            reader = new BufferedReader (new FileReader (outputFile));
            String nextLine = reader.readLine ();
            Integer i = 0;
            // Each line of the output contains the tested npop value and
            // its likelihood.
            while (nextLine != null) {
                StringTokenizer st = new StringTokenizer (nextLine);
                st.nextToken (); // "npop".
                result[i].npop = format.parse (st.nextToken ()).longValue ();
                st.nextToken (); // "likelihood".
                result[i].likelihood = format.parse (st.nextToken ()).doubleValue ();
                nextLine = reader.readLine ();
                i ++;
            }
        }
        catch (IOException e) {
            System.out.println ("Error reading the output file.");
        }
        catch (ParseException e) {
            System.out.println ("Error parsing a number.");
        }
        finally {
            if (reader != null) {
                try {
                    reader.close ();
                }
                catch (IOException e) {
                    System.out.println ("Error closing the output file.");
                }
            }
        }
        return result;
    }

    /**
     *  A private class to store a npop value and its likelihood.
     */
    private class NpopValue {
        public NpopValue (Long npop, Double likelihood) {
            this.npop = npop;
            this.likelihood = likelihood;
        }
        public Long npop;
        public Double likelihood;
    }

    private boolean hasRun;
    private String workingDirectory;
    private ArrayList<ArrayList<String>> ecotypes;
    private MainVariables mainVariables;
    private Execs execs;
    private String outgroup;
    private Integer length;
    private Integer nu;
    private ParameterSet hclimbResult;

    private Integer method;

    private Integer nrep = 1000;
    private Integer step = 1;

    private int iteration;

}
