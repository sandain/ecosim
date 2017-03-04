/*
 *    Ecotype Simulation models the sequence diversity within a bacterial clade as
 *    the evolutionary result of net ecotype formation, periodic selection,
 *    and drift, yielding a certain number of ecotypes.
 * 
 *    Copyright (C) 2009  Andrew Warner, Wesleyan University
 *                        Jason M. Wood, Montana State University
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

/**
 *  Reads in a Newick tree from a file and provides options to traverse it.
 *
 *  @author Andrew Warner & Jason M. Wood
 */
public class NewickTree {

    /**
     *  Constructor for objects of class NewickTree.
     *
     *  @param tree Newick formated tree to read.
     */
    public NewickTree(File tree) throws InvalidNewickException {
        root = loadTree(tree);
    }

    /**
     *  Returns the root node of this tree.
     *
     *  @return NewickTreeNode containing the root node.
     */
    public NewickTreeNode getRoot() {
        return root;
    }

    /**
     *  Returns an ArrayList of NewickTreeNodes that makes of the leaves of this tree.
     *
     *  @return ArrayList containing the leaves.
     */
    public ArrayList<NewickTreeNode> getDescendants() {
        return root.getDescendants();
    }

    /**
     *  Returns this Newick tree formatted as a String.
     *
     *  @return String containing the Newick tree.
     */
    public String toString() {
        return root.toString();
    }

    /**
     *  Save the tree data in this object to a Newick formatted file.
     *
     *  @param fileName File name to write the Newick tree to.
     *  @return True if the save was a success, False otherwise.
     */
    public boolean save (String fileName) {
        return save(new File(fileName));
    }

    /**
     *  Save the tree data in this object to a Newick formatted file.
     *
     *  @param file File to write the Newick tree to.
     *  @return True if the save was a success, False otherwise.
     */
    public boolean save (File file) {
        boolean success = false;
        try {
            BufferedWriter out = new BufferedWriter(new FileWriter(file));
            out.write (
                root.toString() + System.getProperty ("line.separator")
            );
            out.close();
            success = true;
        }
        catch (IOException e) {
            System.out.println("Error writing to output file.");
        }
        return success;
    }

    /**
     *  Reads the tree in from an input file.
     *
     *  @param treeFile The Newick formatted tree file to load.
     *  @return NewickTreeNode containing the root of the Newick tree.
     */
    private NewickTreeNode loadTree(File treeFile) throws InvalidNewickException {
        NewickTreeNode node = new NewickTreeNode();
        try {
            // Read in the tree file.
            BufferedReader input = new BufferedReader(new FileReader(treeFile));
            StringBuffer buffer = new StringBuffer();
            String line = input.readLine();
            while (line != null){
                buffer.append(line);
                line = input.readLine();
            }
            // Check for the semicolon at the end of the Newick tree.
            String bufferString = buffer.toString();
            if (bufferString.contains(";")) {
                // Phylip can return tree files that contain multiple trees, handle only the first one.
                String[] trees = bufferString.split(";");
                // Remove any spaces from the Newick tree before continuing.
                trees[0] = trees[0].replace(" ", "");
                // Generate the tree from the Newick tree.
                node = getTree(trees[0]);
            }
            else {
                throw new InvalidNewickException("Malformed Newick tree, no semicolon found.");
            }
        }
        catch (java.io.FileNotFoundException e) {
            throw new InvalidNewickException("File not found.");
        }
        catch (java.io.IOException e) {
            throw new InvalidNewickException("IO error.");
        }
        // Make sure that we actually have a tree.
        if (node.getChildren().size() <= 1) {
            throw new InvalidNewickException("Malformed Newick tree, not enough leaves found.");
        }
        return node;
    }

    /**
     *  The recursive method used to parse through a Newick tree.
     *
     *  @param tree The subtree to parse.
     *  @return A NewickTreeNode containing the subtree 
     */
    private NewickTreeNode getTree(String tree) throws InvalidNewickException {
        NewickTreeNode node = new NewickTreeNode();
        String metaString;
        // Check if the current node has a subTree.
        if (tree.length() > 0 && tree.charAt(0) == '(') {
            // Split the sub tree from the meta data.
            int length = getSpaceInBetweenParens(0, tree);
            String subTree = tree.substring(1, length);
            metaString = tree.substring(length + 1, tree.length());
            // Parse this node's sub tree.
            StringBuffer buffer = new StringBuffer();
            int i = 0;
            while (i < subTree.length()) {
                // Check to see if this child node has a subtree, append it to the buffer.
                if (subTree.charAt(i) == '(') {
                    int j = getSpaceInBetweenParens(i, subTree);
                    buffer.append(subTree.substring(i, j));
                    i = j;
                }
                // Keep appending the child node's meta data to the buffer until we find
                // a new child node.  If we do find a new child node, add it to the current
                // node by recursion with the buffer as the subTree, then clear the buffer.
                if (subTree.charAt(i) != ',') {
                    buffer.append(subTree.charAt(i));
                }
                else {
                    node.addChild(getTree(buffer.toString()));
                    buffer = new StringBuffer();
                }
                i ++;
            }
            // Add a new child to the current node by recursion with the buffer as the subTree.
            node.addChild(getTree(buffer.toString()));
        }
        else {
            metaString = tree;
        }
        // Parse this node's meta data.
        if (metaString.length() > 0) {
            String[] meta = metaString.split(":", 2);
            if (meta.length > 0 && meta[0].length() > 0) {
                node.setName(meta[0]);
            }
            if (meta.length > 1 && meta[1].length() > 0) {
                try {
                    Double distance = new Double(meta[1]);
                    node.setDistance(distance.doubleValue());
                }
                catch (NumberFormatException e) {
                    throw new InvalidNewickException("Malformed Newick tree, expected to find a number.");
                }
            }
        }
        // Each child of this node needs to know who its parent is.
        for (int i = 0; i < node.getChildren().size(); i ++) {
            node.getChildren().get(i).setParent(node);
        }
        return node;
    }

    /**
     *  Get the space in between a set of parentheses, ie (   ), which could contain
     *  other open and close parentheses within it.
     *
     *  @param index The index of the first open parenthesis.
     *  @param tree The tree, or the tree segment.
     *  @return The index of the close parenthesis corresponding to the given index of the open parenthesis.
     */
    private int getSpaceInBetweenParens(int index, String tree) throws InvalidNewickException {
        int numOpenParensSeen = 0;
        for (int i = index; i < tree.length(); i ++) {
            if (tree.charAt(i) == '(') {
                numOpenParensSeen ++;
            }
            else if (tree.charAt(i) == ')') {
                numOpenParensSeen --;
                if (numOpenParensSeen == 0) {
                    return i;
                }
            }
        }
        // Throw an exception if there is there are unmatched parentheses.
        if (numOpenParensSeen != 0) {
            throw new InvalidNewickException("Malformed Newick tree, unmatched parentheses.");
        }
        return 0;
    }

    private NewickTreeNode root;

}
