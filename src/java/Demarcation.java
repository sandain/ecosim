/*
 *    Ecotype Simulation models the sequence diversity within a bacterial
 *    clade as the evolutionary result of net ecotype formation and periodic
 *    selection, yielding a certain number of ecotypes.
 *
 *    Copyright (C) 2009       Carlo Francisco, Wesleyan University
 *    Copyright (C) 2009-2016  Jason M. Wood, Montana State University
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

import ecosim.tree.Node;
import ecosim.tree.InvalidTreeException;
import ecosim.tree.Tree;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Comparator;
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

    public static final int DEMARCATION_PRECISION_COARSE_SCALE = 2501;
    public static final int DEMARCATION_PRECISION_FINE_SCALE = 2502;

    /**
     *  Creates new form Demarcations
     *
     *  @param masterVariables The MasterVariables.
     *  @param execs The Execs object.
     *  @param nu The number of environmental sequences.
     *  @param length The length of the environmental sequences.
     *  @param outgroup The name of the outgroup.
     *  @param tree The phylogeny data.
     *  @param hclimbResult The result from hillclimbing.
     *  @param precision The precision to use for demarcation.
     */
    public Demarcation (MasterVariables masterVariables, Execs execs,
        Integer nu, Integer length, String outgroup, Tree tree,
        ParameterSet hclimbResult, int precision)
        throws InvalidTreeException {
        super (tree);
        this.masterVariables = masterVariables;
        this.execs = execs;
        this.nu = nu;
        this.length = length;
        this.outgroup = outgroup;
        this.hclimbResult = hclimbResult;
        this.precision = precision;
        hasRun = false;
        ecotypes = new ArrayList<ArrayList<String>> ();
        workingDirectory = masterVariables.getWorkingDirectory ();
    }

    /**
     *  Run the demarcation program.
     */
    public void run () throws InvalidTreeException {
        iteration = 0;
        // Find the ecotypes.
        findEcotypes (root);
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
     *  Set the ecotypes.
     *
     *  @param ecotypes The ecotypes.
     */
    public void setEcotypes (ArrayList<ArrayList<String>> ecotypes) {
        this.ecotypes = ecotypes;
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
     *  Find the ecotypes using a recursive algorithm.  If the subclade of
     *  the tree represented by node has an optimal npop value of 1, demarcate
     *  the subclade as an ecotype. Otherwise, recurse on the node's children.
     *
     *  @param node The current node representing the subclade.
     */
    private void findEcotypes (Node node) throws InvalidTreeException {
        ArrayList<String> sample = new ArrayList<String> ();
        String ecotype = String.format ("Ecotype%03d", ecotypes.size () + 1);
        if (node.isLeafNode ()) {
            String name = node.getName ();
            if (! name.equals (outgroup)) {
                sample.add (name);
                ecotypes.add (sample);
                node.setName (ecotype);
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
                node.collapse ();
                node.setName (ecotype);
            }
            else {
                // Npop > 1, recurse on children nodes.
                ArrayList<Node> children = node.getChildren ();
                for (int i = 0; i < children.size (); i ++) {
                    findEcotypes (children.get (i));
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
        Integer sampleNu = node.numberOfDescendants ();
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
        NpopValue result = new NpopValue (0L, 0.0D);
        // Determine which result to return.
        switch (precision) {
            case DEMARCATION_PRECISION_FINE_SCALE:
                // Fine scale uses the most likely result.
                result = results[1];
                break;
            case DEMARCATION_PRECISION_COARSE_SCALE:
                // Coarse scale uses the result of npop=1 if likelihood > 0.
                Double likelihood = results[0].likelihood;
                if (likelihood.compareTo (MasterVariables.EPSILON) > 0) {
                    result = results[0];
                }
                else {
                    result = results[1];
                }
                break;
            default:
                System.err.println ("Error, unknown precision level.");
        }
        return result;
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
            writer.write (String.format ("%-20d numcrit\n", bins.size ()));
            // Output the crit levels and the number of bins.
            for (int j = 0; j < bins.size (); j ++) {
                writer.write (String.format (
                    "%-20.6f %-20d\n",
                    bins.get (j).getCrit (),
                    bins.get (j).getLevel ()
                ));
            }
            // Write the omega value.
            writer.write (String.format ("%-20.5f omega\n", omega));
            // Write the sigma value.
            writer.write (String.format ("%-20.5f sigma\n", sigma));
            // Write the npop value.
            writer.write (String.format ("%-20d npop\n", npop));
            // Write the step value.
            writer.write (String.format ("%-20d step\n", step));
            // Write the nu value.
            writer.write (String.format ("%-20d nu\n", sampleNu));
            // Write the nrep value.
            writer.write (String.format ("%-20d nrep\n", nrep));
            // Create the random number seed; an odd integer less than nine
            // digits long.
            long iii = (long)(100000000 * Math.random ());
            if (iii % 2 == 0) {
                iii ++;
            }
            // Write the random number seed.
            writer.write (
                String.format ("%-20d iii (random number seed)\n", iii)
            );
            // Write the length of the sequences.
            writer.write (
                String.format (
                    "%-20d lengthseq (after deleting gaps, etc.)\n",
                    length
                )
            );
            // Write the whichavg value.
            int whichavg = masterVariables.getCriterion ();
            writer.write (String.format ("%-20d whichavg\n", whichavg));
            // Write the likelihoodsolution value.
            writer.write (String.format (
                "%-20.5f likelihoodsolution\n", likelihood
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
                result[i].npop = new Long (st.nextToken ());
                st.nextToken (); // "likelihood".
                result[i].likelihood = new Double (st.nextToken ());
                nextLine = reader.readLine ();
                i ++;
            }
        }
        catch (IOException e) {
            System.out.println ("Error reading the output file.");
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
    private MasterVariables masterVariables;
    private Execs execs;
    private String outgroup;
    private Integer length;
    private Integer nu;
    private ParameterSet hclimbResult;

    private Integer precision;

    private Integer nrep = 1000;
    private Integer step = 1;

    private int iteration;

}
