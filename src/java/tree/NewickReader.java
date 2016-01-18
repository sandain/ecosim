/*
 *    Ecotype Simulation models the sequence diversity within a bacterial
 *    clade as the evolutionary result of net ecotype formation and periodic
 *    selection, yielding a certain number of ecotypes.
 *
 *    Copyright (C) 2009       Andrew Warner, Wesleyan University
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


package ecosim.tree;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.Reader;
import java.util.ArrayList;

/**
 *  Read a Newick formated tree.
 *
 *  @author Andrew Warner
 *  @author Jason M. Wood
 *  @copyright GNU General Public License
 */
public class NewickReader extends BufferedReader {

    /**
     *  Read a Newick formated tree. 
     */
    public NewickReader (Reader reader) {
        super (reader);
    }

    /**
     *  Read a Newick formated tree.
     *
     *  @return Node containing the root of the Newick tree.
     */
    public Node readTree () throws InvalidTreeException {
        Node node = null;
        try {
            // Create a buffer to hold the tree.
            StringBuffer buffer = new StringBuffer ();
            // Start reading the tree.
            String line = readLine ();
            while (line != null) {
                // Replace spaces in the line.
                line = line.replace (" ", "");
                // Look for the end of the Newick tree.
                if (line.contains (";")) {
                    buffer.append (line.substring (0, line.indexOf (";")));
                    break;
                }
                // Append the line to the buffer.
                buffer.append (line);
                // Read the next line of the file.
                line = readLine ();
            }
            // If the buffer contains data, attempt to parse the tree.
            if (buffer.length () > 0) {
                node = parseTree (buffer.toString ());
            }
        }
        catch (IOException e) {
            throw new InvalidTreeException ("Unable to read from file.");
        }
        return node;
    }

    /**
     *  The recursive method used to parse through a Newick tree.
     *
     *  @param tree The subtree to parse.
     *  @return A Node containing the subtree.
     */
    private Node parseTree (String tree) throws InvalidTreeException {
        Node node = new Node ();
        String metaString;
        // Check if the current node has a subTree.
        if (tree.length () > 0 && tree.charAt (0) == '(') {
            // Split the sub tree from the meta data.
            int length = getSpaceInBetweenParens (0, tree);
            String subTree = tree.substring (1, length);
            metaString = tree.substring (length + 1, tree.length ());
            // Parse this node's sub tree.
            StringBuffer buffer = new StringBuffer ();
            int i = 0;
            while (i < subTree.length ()) {
                // Check to see if this child node has a subtree, append it to
                // the buffer.
                if (subTree.charAt (i) == '(') {
                    int j = getSpaceInBetweenParens (i, subTree);
                    buffer.append (subTree.substring (i, j));
                    i = j;
                }
                // Keep appending the child node's meta data to the buffer
                // until we find a new child node.  If we do find a new child
                // node, add it to the current node by recursion with the
                // buffer as the subTree, then clear the buffer.
                if (subTree.charAt (i) != ',') {
                    buffer.append (subTree.charAt (i));
                }
                else {
                    node.addChild (parseTree (buffer.toString ()));
                    buffer = new StringBuffer ();
                }
                i ++;
            }
            // Add a new child to the current node by recursion with the
            // buffer as the subTree.
            node.addChild (parseTree (buffer.toString ()));
        }
        else {
            metaString = tree;
        }
        // Parse this node's meta data.
        if (metaString.length () > 0) {
            String[] meta = metaString.split (":", 2);
            if (meta.length > 0 && meta[0].length () > 0) {
                node.setName (meta[0]);
            }
            if (meta.length > 1 && meta[1].length () > 0) {
                try {
                    Double distance = new Double (meta[1]);
                    node.setDistance (distance.doubleValue ());
                }
                catch (NumberFormatException e) {
                    throw new InvalidTreeException (
                        "Malformed Newick tree, expected a number." + e
                    );
                }
            }
        }
        return node;
    }

    /**
     *  Get the space in between a set of parentheses, ie (   ), which could
     *  contain other open and close parentheses within it.
     *
     *  @param index The index of the first open parenthesis.
     *  @param tree The tree, or the tree segment.
     *  @return The index of the close parenthesis corresponding to the given
     *  index of the open parenthesis.
     */
    private int getSpaceInBetweenParens (int index, String tree)
        throws InvalidTreeException {
        int numOpenParensSeen = 0;
        for (int i = index; i < tree.length (); i ++) {
            if (tree.charAt (i) == '(') {
                numOpenParensSeen ++;
            }
            else if (tree.charAt (i) == ')') {
                numOpenParensSeen --;
                if (numOpenParensSeen == 0) {
                    return i;
                }
            }
        }
        // Throw an exception if there is there are unmatched parentheses.
        if (numOpenParensSeen != 0) {
            throw new InvalidTreeException (
                "Malformed Newick tree, unmatched parentheses."
            );
        }
        return 0;
    }

}
