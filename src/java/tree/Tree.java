/*
 *    Ecotype Simulation models the sequence diversity within a bacterial
 *    clade as the evolutionary result of net ecotype formation and periodic
 *    selection, yielding a certain number of ecotypes.
 *
 *    Copyright (C) 2009       Andrew Warner, Wesleyan University
 *    Copyright (C) 2009-2015  Jason M. Wood, Montana State University
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

import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.Reader;
import java.io.StringReader;

import java.util.ArrayList;

/**
 *  Reads in a Newick tree from a file and provides options to traverse it.
 *
 *  @author Andrew Warner
 *  @author Jason M. Wood
 *  @copyright GNU General Public License
 */
public class Tree {

    /**
     *  The default constructor for objects of class Tree.
     */
    public Tree () {
        root = new Node ();
    }

    /**
     *  Constructor for objects of class Tree.
     *
     *  @param tree String containing a Newick formated tree to read.
     */
    public Tree (String tree) throws InvalidTreeException {
        readTree (new StringReader (tree));
    }

    /**
     *  Constructor for objects of class Tree.
     *
     *  @param file File containing a Newick formated tree to read.
     */
    public Tree (File file) throws InvalidTreeException {
        FileReader reader = null;
        try {
            reader = new FileReader (file);
            readTree (reader);
        }
        catch (FileNotFoundException e) {
            throw new InvalidTreeException ("File not found.");
        }
        finally {
            try {
                reader.close ();
            }
            catch (IOException e) {
                throw new InvalidTreeException ("IO error.");
            }
        }
    }

    /**
     *  Constructor for objects of class Tree.
     *
     *  @param tree Tree to use.
     */
    public Tree (Tree tree) throws InvalidTreeException {
        readTree (new StringReader (tree.toString ()));
    }

    /**
     *  Returns the root node of this tree.
     *
     *  @return Node containing the root node.
     */
    public Node getRoot () {
        return root;
    }

    /**
     *  Returns the descendant with the given name.
     *
     *  @param name The name of the descendant.
     *  @return Node of the descendant.
     */
    public Node getDescendant (String name) {
        Node namedDescendant = new Node ();
        for (Node descendant: root.getDescendants ()) {
            if (descendant.getName ().equals (name)) {
                namedDescendant = descendant;
                break;
            }
        }
        return namedDescendant;
    }

    /**
     *  Returns an array of nodes that are descendants of the root node.
     *
     *  @return The descendants of the root node.
     */
    public ArrayList<Node> getDescendants () {
        return root.getDescendants ();
    }

    /**
     *  Remove a leaf node descendant from the tree.
     *
     *  @param name The name of the descendant to remove.
     */
    public void removeDescendant (String name) {
        Node descendant = getDescendant (name);
        removeDescendant (descendant);
    }

    /**
     *  Remove a leaf node descendant from the tree.
     *
     * Since we assume that each internal node has exactly two children,
     * the parent node of the descendant being destroyed will also have to
     * be destroyed.  The other child will have its grandparent become its
     * new parent, and its distance will become the sum of its distance and
     * its parents distance.
     *
     *  @param descendant The descendant to remove.
     */
    public void removeDescendant (Node descendant) {
        if (! descendant.isLeafNode ()) return;
        // Grab the parent node of the descendant.
        Node parent = descendant.getParent ();
        // Remove the descendant as a child of the parent node.
        parent.removeChild (descendant);
        // Grab the other child node of the parent.
        Node otherChild = parent.getChildren ().get (0);
        // Calculate a new distance for the other child node.
        Double distance = otherChild.getDistance () + parent.getDistance ();
        otherChild.setDistance (distance);
        // Setup the relationships of the other child node.
        if (parent.isRootNode ()) {
            // The parent is the root node, the other child becomes the new 
            // root node.
            otherChild.setParent (null);
            root = otherChild;
        }
        else {
            // The grandparent of the child node is the new parent.
            Node grandparent = parent.getParent ();
            grandparent.removeChild (parent);
            grandparent.addChild (otherChild);
        }
    }

    /**
     *  Reroot the tree so that the outgroup descendant is next to the root.
     *
     *  @param outgroup The outgroup descendant.
     */
    public void reroot (Node outgroup) {
        Node newRoot = new Node ();
        Node oldParent = outgroup.getParent ();
        // Add the outgroup to the new root, and remove it from its old
        // parents children.
        oldParent.removeChild (outgroup);
        newRoot.addChild (outgroup);
        double distance = outgroup.getDistance () * 0.5;
        double oldDistance = oldParent.getDistance ();
        outgroup.setDistance (distance);
        outgroup.setOutgroup (true);
        oldParent.setDistance (distance);
        Node newParent = newRoot;
        while (oldParent != root) {
            Node node = oldParent;
            oldParent = node.getParent ();
            oldParent.removeChild (node);
            newParent.addChild (node);
            newParent = node;
            node.setDistance (distance);
            distance = oldDistance;
            oldDistance = oldParent.getDistance ();
        }
        for (Node child: oldParent.getChildren ()) {
          newParent.addChild (child);
          distance = child.getDistance () + distance;
          child.setDistance (distance);
        }
        // Save the the new root.
        root = newRoot;
    }


    /**
     *  Return the number of living descendants of the root node.
     *
     *  @return The number of descendants.
     */
    public int size () {
        return root.numberOfDescendants ();
    }

    /**
     *  Returns this tree as a Newick formatted String.
     *
     *  @return Newick formatted String containing the tree.
     */
    public String toString () {
        return root.toString ();
    }

    /**
     *  Save the tree data in this object to a Newick formatted file.
     *
     *  @param file File to write the Newick tree to.
     *  @return True if the save was a success, False otherwise.
     */
    public boolean toNewick (File file) {
        boolean success = false;
        NewickWriter out = null;
        try {
            out = new NewickWriter (new FileWriter (file));
            out.write (this);
        }
        catch (IOException e) {
            System.out.println ("Error writing to output file.");
        }
        finally {
            try {
                out.close ();
                success = true;
            }
            catch (IOException e) {
                System.out.println ("Error closing output file.");
            }
        }
        return success;
    }

    /**
     *  Check if this Tree object is valid.
     *
     *  @return True if this is a valid Tree object, False if not.
     */
    public boolean isValid () {
        boolean valid = false;
        if (root.getDescendants ().size () > 0) {
            valid = true;
        }
        return valid;
    }

    /**
     *  Read a Newick formatted tree.
     *
     *  @param tree Reader to use.
     */
    private void readTree (Reader reader) throws InvalidTreeException {
        // Read in the tree file.
        NewickReader newick = new NewickReader (reader);
        root = newick.readTree ();
        // Make sure that we actually have a tree.
        if (root.getChildren ().size () <= 1) {
            throw new InvalidTreeException (
                "Malformed Newick tree, not enough leaves found."
            );
        }
    }

    private Node root;

}
