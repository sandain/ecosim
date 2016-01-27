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

import ecosim.api.Painter;

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
     *  Compare this tree with another.
     *
     *  @param other The tree compare with.
     *  @return 0 if the trees are equal, -1 or 1 otherwise.
     */
    public int compareTo (Tree other) {
        return root.compareTo (other.getRoot ());
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
     *  Returns the maximum width of the tree from root node to the most
     *  divergent leaf node.
     *
     *  @return Node maximum width of the tree.
     */
    public Double maximumWidth () {
        return root.maximumDistanceFromLeafNode ();
    }

    /**
     *  Returns the maximum X value of the tree.
     *
     *  The maximum X value is the distance from the root node to the most
     *  divergent leaf node.
     *
     *  @return Node maximum X value of the tree.
     */
    public Double maximumX () {
        return root.maximumDistanceFromLeafNode ();
    }

    /**
     *  Returns the maximum Y value of the tree.
     *
     *  The maximum Y value is the distance from the root node to the most
     *  divergent leaf node.
     *
     *  @return Node maximum Y value of the tree.
     */
    public Double maximumY () {
        return 1.0d * root.numberOfDescendants ();
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
     *  Reroot the tree so that the outgroup descendant is next to the root.
     *
     *  @param name The name of the outgroup descendant.
     */
    public void reroot (String name) {
        Node descendant = getDescendant (name);
        reroot (descendant);
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
        // Recalculate the XY location of all nodes.
        calculateNodeXY (root, 0.0d);
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
     *  Paint the tree using the given Painter.
     *
     *  @param painter The Painter to use.
     */
    public void paintTree (Painter painter) {
        int fontHeight = painter.fontHeight ();
        int fontWidth = painter.fontWidth ();
        int xModifier = 1000;
        int xSpacer = fontWidth / 2;
        int height = fontHeight * (size () + 1);
        int max = 0;
        for (Node node: getDescendants ()) {
            String name = node.getName ();
            int labelWidth = (name.length () + 1) * fontWidth;
            int x = fontWidth + labelWidth + Math.round (
                node.getX ().floatValue () * xModifier
            );
            if (x > max) max = x;
        }
        int width = max;
        painter.start (width, height);
        paintNode (painter, root);
        painter.end ();
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
     *  Recursive method to paint a node and all of its descendants.
     *
     *  @param painter The Painter to use.
     *  @param node The Node to paint.
     */
    protected void paintNode (Painter painter, Node node) {
        int fontHeight = painter.fontHeight ();
        int fontWidth = painter.fontWidth ();
        int xModifier = 1000;
        int yModifier = fontHeight;
        int xSpacer = (int)Math.floor (0.5d * fontWidth);
        int ySpacer = (int)Math.floor (0.5d * fontHeight);
        int nodeX = fontWidth + Math.round (
            node.getX ().floatValue () * xModifier
        );
        int nodeY = fontHeight + Math.round (
            node.getY ().floatValue () * yModifier
        );
        if (node.isLeafNode () || node.isCollapsed ()) {
            // Paint the name of the node.
            painter.drawString (
                node.getName (), nodeX + xSpacer, nodeY + ySpacer - 2
            );
        }
        else {
            for (Node child: node.getChildren ()) {
                int childX = fontWidth + Math.round (
                    child.getX ().floatValue () * xModifier
                );
                int childY = fontHeight + Math.round (
                    child.getY ().floatValue () * yModifier
                );
                // Paint a vertical line connecting the node to it's parent.
                painter.drawLine (nodeX, nodeY, nodeX, childY);
                // Paint a triangle if the child node is collapsed, otherwise
                // draw a horizontal line.
                if (child.isCollapsed () && child.numberOfDescendants () > 1) {
                    int a = childY - ySpacer + 1;
                    int b = childY + ySpacer - 1;
                    // Paint a triangle.
                    painter.drawLine (nodeX, childY, childX, a);
                    painter.drawLine (nodeX, childY, childX, b);
                    painter.drawLine (childX, a, childX, b);
                }
                else {
                    // Paint a horizontal line.
                    painter.drawLine (nodeX, childY, childX, childY);
                }
                // Paint the child node.
                paintNode (painter, child);
            }
        }
    }

    /**
     *  Read a Newick formatted tree.
     *
     *  @param reader The Reader to use.
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
        // Calculate the XY location of all nodes.
        calculateNodeXY (root, 0.0d);
    }

    /**
     *  Private recursive method to calculate the XY location of the node
     *  and all of its descendants.
     *
     *  @param node The current Node.
     *  @param height The current height.
     */
    protected void calculateNodeXY (Node node, double height) {
        Node parent = node.getParent ();
        // The X coordinate is based on the node's distance from its parent.
        double x = node.getDistance ();
        // Add the parent's X coordinate.
        if (parent != null) x += parent.getX ();
        // If the node is collapsed, add the descendants distance as well.
        if (node.isCollapsed () && node.numberOfDescendants () > 1) {
            double max = node.maximumDistanceFromLeafNode ();
            if (max < 0.01d) max = 0.01d;
            x += max;
        }
        node.setX (x);
        // The Y coordinate is calculated differently for leaf and internal
        // nodes.
        if (node.isLeafNode () || node.isCollapsed ()) {
            node.setY (height);
        }
        else {
            double minY = Double.MAX_VALUE;
            double maxY = 0.0d;
            double h = height;
            for (Node child: node.getChildren ()) {
                calculateNodeXY (child, h);
                int numLeaf = 1;
                if (! (child.isLeafNode () || child.isCollapsed ())) {
                    numLeaf = child.numberOfDescendants ();
                }
                h += numLeaf;
                double childY = child.getY ();
                if (childY < minY) minY = childY;
                if (childY > maxY) maxY = childY;
            }
            node.setY ((minY + maxY) / 2);
        }
    }

    protected Node root;

}
