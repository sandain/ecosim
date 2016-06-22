/*
 *    Ecotype Simulation models the sequence diversity within a bacterial
 *    clade as the evolutionary result of net ecotype formation and periodic
 *    selection, yielding a certain number of ecotypes.
 *
 *    Copyright (C) 2013-2016  Jason M. Wood, Montana State University
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
import ecosim.tree.Tree;

import java.util.ArrayList;

/**
 *  Object to to estimate the number of bins in a provided tree using a
 *  complete-linkage function.
 *
 *  @author Jason M. Wood
 *  @copyright GNU General Public License
 */
public class Binning implements Runnable {

    /**
     *  A default constructor for Binning objects.  All crit levels
     *  will be initialized with zeros.
     */
    public Binning () {
        bins = new ArrayList<BinLevel> ();
    }

    /**
     *  Object to estimate the number of bins in a provided tree.
     *
     *  @param tree The Tree object.
     */
    public Binning (Tree tree) {
        bins = new ArrayList<BinLevel> ();
        this.tree = tree;
    }

    /**
     *  Run the binning program.
     */
    public void run () {
        // Only run binning if a tree has been loaded.
        if (tree == null) return;
        // Run complete-linkage binning on the provided tree for each
        // sequence criterion level.
        for (Double crit: binLevels) {
            Integer level = getNumberBins (crit, tree.getRoot ());
            bins.add (new BinLevel (crit, level));
        }
    }

    /**
     *  Returns an ArrayList<BinLevel> containing the bin levels.
     *
     *  @return The bin levels.
     */
    public ArrayList<BinLevel> getBins () {
        return bins;
    }

    /**
     *  Change the bin levels stored in this object.
     *
     *  @param bins The new bin levels.
     */
    public void setBins (ArrayList<BinLevel> bins) {
        this.bins = bins;
    }

    /**
     *  Add a bin level to this object.
     *
     *  @param binLevel The new BinLevel to add.
     */
    public void addBinLevel (BinLevel binLevel) {
        bins.add (binLevel);
    }

    /**
     *  Returns the binning result as a String.
     *
     *  @return the binning result.
     */
    public String toString () {
        String str = "";
        for (int i = 0; i < bins.size (); i ++) {
            str += bins.get (i).toString ();
            if (i < bins.size () - 1) {
                str += ", ";
            }
        }
        return str;
    }

    /**
     *  A private recursive method to estimate the number of bins in a
     *  provided tree using a complete-linkage method.
     *
     *  @param crit The sequence criterion level.
     *  @param node The current node in the tree to examine.
     *  @return The number of bins.
     */
    private Integer getNumberBins (Double crit, Node node) {
        Integer num = 0;
        // Return zero bins if the node is the outgroup.
        if (node.isOutgroup ()) return 0;
        // Return one bin if the node is a leaf node.
        if (node.isLeafNode ()) return 1;
        // Calculate the maximum distance between the leaf node ancestors of
        // the current node.
        Double distance = node.maximumDistanceBetweenLeafNodes ();
        // Recurse on the child nodes if the maximum distance exceeds the
        // crit value, otherwise colapse this branch into one bin.
        if (distance > 1.000d - crit - MasterVariables.EPSILON) {
            for (Node child: node.getChildren ()) {
                num += getNumberBins (crit, child);
            }
        }
        else {
            num = 1;
        }
        return num;
    }

    private ArrayList<BinLevel> bins;
    private Tree tree;

    /**
     *  The default bin levels.
     */
    public static Double[] binLevels = {
        0.800d, 0.810d, 0.820d, 0.830d, 0.840d, 0.850d, 0.860d,
        0.870d, 0.880d, 0.890d, 0.900d, 0.910d, 0.920d, 0.930d,
        0.940d, 0.950d, 0.955d, 0.960d, 0.965d, 0.970d, 0.975d,
        0.980d, 0.981d, 0.982d, 0.983d, 0.984d, 0.985d, 0.986d,
        0.987d, 0.988d, 0.989d, 0.990d, 0.991d, 0.992d, 0.993d,
        0.994d, 0.995d, 0.996d, 0.997d, 0.998d, 0.999d, 1.000d
    };

}
