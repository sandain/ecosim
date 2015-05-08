/*
 *    Ecotype Simulation models the sequence diversity within a bacterial
 *    clade as the evolutionary result of net ecotype formation and periodic
 *    selection, yielding a certain number of ecotypes.
 *
 *    Copyright (C) 2013-2014  Jason M. Wood, Montana State University
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
 *  Object to to estimate the number of bins in a provided tree using a
 *  complete-linkage function.
 *
 *  @author Jason M. Wood
 *  @copyright GNU General Public License
 */
public class Binning {

    /**
     *  Object to estimate the number of bins in a provided tree.
     *
     *  @param masterVariables The MasterVariables object.
     *  @param tree The NewickTree object.
     */
    public Binning (MasterVariables masterVariables, NewickTree tree) {
        this.masterVariables = masterVariables;
        bins = new ArrayList<BinLevel> ();
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
    private Integer getNumberBins (Double crit, NewickTreeNode node) {
        Integer num = 0;
        // Return zero bins if the node is the outgroup.
        if (node.isOutgroup ()) return 0;
        // Return one bin if the node is a leaf node.
        if (node.isLeafNode ()) return 1;
        // There can only be two children.
        ArrayList<NewickTreeNode> children = node.getChildren ();
        NewickTreeNode a = children.get (0);
        NewickTreeNode b = children.get (1);
        // Calculate the maximum distance between the leaf node ancestors
        // of the two child nodes.
        Double distance = 0.0d;
        distance += a.maximumDistanceFromLeafNode () + a.getDistance ();
        distance += b.maximumDistanceFromLeafNode () + b.getDistance ();
        // Recurse on the child nodes if the maximum distance exceeds the
        // crit value, otherwise colapse this branch into one bin.
        if (distance > 1.000d - crit - masterVariables.EPSILON) {
            for (NewickTreeNode child: children) {
                num += getNumberBins (crit, child);
            }
        }
        else {
            num = 1;
        }
        return num;
    }

    private MasterVariables masterVariables;

    private ArrayList<BinLevel> bins;

    /**
     *  The default bin levels.
     */
    private Double[] binLevels = {
        0.800d, 0.850d, 0.900d, 0.950d, 0.960d, 0.970d, 0.980d, 0.990d, 0.995d, 1.000d
    };

}
