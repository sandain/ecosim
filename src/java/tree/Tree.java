/*
 *    Ecotype Simulation models the sequence diversity within a bacterial
 *    clade as the evolutionary result of net ecotype formation and periodic
 *    selection, yielding a certain number of ecotypes.
 *
 *    Copyright (C) 2009       Andrew Warner, Wesleyan University
 *    Copyright (C) 2009-2018  Jason M. Wood, Montana State University
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

    public static final int PAINT_METHOD_NORMAL     = 1010;
    public static final int PAINT_METHOD_COLLAPSED  = 1020;
    public static final int PAINT_METHOD_DEMARCATED = 1030;

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
        root = new Node (tree.getRoot ());
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
     *  Returns the maximum diversity of this tree.
     *
     *  @return The maximum diversity of this tree.
     */
    public Double getDiversity () {
        Double diversity = 0.0d;
        for (Node child: root.getChildren ()) {
            if (child.isOutgroup ()) continue;
            diversity = 1.0d - child.maximumDistanceBetweenLeafNodes ();
        }
        return diversity;
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
     *  Get the paint method.
     *
     *  @return The paint method.
     */
    public int getPaintMethod () {
        return paintMethod;
    }

    /**
     *  Set the paint method.
     *
     *  @param paintMethod The paint method.
     */
    public void setPaintMethod (int paintMethod) {
        this.paintMethod = paintMethod;
    }

    /**
     *  Set the scale used for tree display.
     *
     *  @param scale The scale to use.
     */
    public void setScale (int scale) {
        xModifier = scale;
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
        return 1.0d * numberOfDescendants (root);
    }

    /**
     *  Returns the number of living descendants of the Node.
     *
     *  @return The number of living descendants of the Node.
     */
    public int numberOfDescendants (Node node) {
        int descendants = 0;
        boolean paintCollapsed = (paintMethod == PAINT_METHOD_COLLAPSED);
        for (Node child: node.getChildren ()) {
            if (child.isLeafNode () ||
                (child.isCollapsed () && paintCollapsed)) {
                descendants ++;
            }
            else {
                descendants += numberOfDescendants (child);
            }
        }
        return descendants;
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
     *  Returns an array of nodes that are collapsed descendants of the
     *  root node.
     *
     *  @return The collapsed descendants of the root node.
     */
    public ArrayList<Node> getCollapsed () {
        return root.getCollapsed ();
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
     *  Make this tree binary.
     */
    public void makeBinary () {
        makeBinary (root);
    }

    /**
     *  Private recursive method to make this tree binary.
     *
     *  @param node The current node being analyzed.
     */
    private void makeBinary (Node node) {
        ArrayList<Node> children = new ArrayList<Node> (node.getChildren ());
        // If this node has more than 2 children, split the node in two.
        if (children.size () > 2) {
            // Skip the first child.
            children.remove (0);
            // Move all remaining children to a new parent.
            Node parent = new Node ();
            for (Node child: children) {
                node.removeChild (child);
                parent.addChild (child);
            }
            // Add the new parent as a child to this node.
            node.addChild (parent);
        }
        // Make sure each child is binary as well.
        for (Node child: node.getChildren ()) {
            makeBinary (child);
        }
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
    }

    /**
     *  Return the number of living descendants of the root node.
     *
     *  @return The number of descendants.
     */
    public int size () {
        return numberOfDescendants (root);
    }

    /**
     *  Returns this tree as a Newick formatted String.
     *
     *  @return Newick formatted String containing the tree.
     */
    public String toString () {
        validateTree (root);
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
        int xSpacer = fontWidth / 2;
        // Calculate the XY location of all nodes.
        calculateNodeXY (root, 0.0d);
        // Calculate the max X value.
        int max = 0;
        for (Node node: getDescendants ()) {
            String name = node.getName ();
            int labelWidth = painter.stringWidth (name);
            int x = fontWidth + labelWidth + Math.round (
                node.getX ().floatValue () * xModifier
            );
            if (x > max) max = x;
        }
        // Calculate the height and width needed for the tree.
        int height = fontHeight * (size () + 3);
        int width = max;
        // Make room for the demarcation line if needed.
        if (paintMethod == PAINT_METHOD_DEMARCATED) {
            // Add room for the demarcation line.
            width += 100;
            // Add room for the label.
            int maxLabel = 0;
            for (Node node: getCollapsed ()) {
                String name = node.getName ();
                int labelWidth = painter.stringWidth (name);
                if (labelWidth > maxLabel) maxLabel = labelWidth;
            }
            width += xSpacer + maxLabel;
        }
        // Paint the tree.
        painter.start (width, height);
        paintNode (painter, root);
        paintScaleBar (painter, 25, height - fontHeight);
        if (paintMethod == PAINT_METHOD_DEMARCATED) {
            paintDemarcation (painter, root, max + 10);
        }
        painter.end ();
    }

    /**
     *  Check if this Tree object is valid.
     *
     *  @return True if this is a valid Tree object, False if not.
     */
    public boolean isValid () {
        boolean valid = false;
        if (root != null && root.getDescendants ().size () > 0) {
            valid = true;
        }
        return valid;
    }

    /**
     *  Find the last common ancestor node for the named descendants.
     *
     *  @param names The names of the descendants.
     *  @return The last common ancestor node of the descendants.
     */
    public Node lastCommonAncestor (ArrayList<String> names) {
        ArrayList<Node> descendants = new ArrayList<Node> ();
        for (String name: names) {
            descendants.add (getDescendant (name));
        }
        if (descendants.isEmpty ()) return null;
        return lastCommonAncestor (descendants, descendants.get (0));
    }

    /**
     *  A private recursive method to find the last common ancestor node for
     *  the given descendants.
     *
     *  @param descendants The list of descendants.
     *  @param node The current node being checked.
     *  @return The last common ancestor node of the descendants.
     */
    private Node lastCommonAncestor (ArrayList<Node> descendants, Node node) {
        boolean has = true;
        for (Node descendant: descendants) {
            has &= node == descendant || node.hasDescendant (descendant);
        }
        if (has) {
            return node;
        }
        else {
           return lastCommonAncestor (descendants, node.getParent ());
        }
    }

    /**
     *  Recursive method to paint a node and all of its descendants.
     *
     *  @param painter The Painter to use.
     *  @param node The Node to paint.
     */
    private void paintNode (Painter painter, Node node) {
        int fontHeight = painter.fontHeight ();
        int fontWidth = painter.fontWidth ();
        int stroke = 1;
        int demarcationStroke = 10;
        int yModifier = fontHeight;
        int xSpacer = (int)Math.floor (0.5d * fontWidth);
        int ySpacer = (int)Math.floor (0.5d * fontHeight);
        int nodeX = fontWidth + Math.round (
            node.getX ().floatValue () * xModifier
        );
        int nodeY = fontHeight + Math.round (
            node.getY ().floatValue () * yModifier
        );
        boolean paintCollapsed = (paintMethod == PAINT_METHOD_COLLAPSED);
        boolean paintDemarcated = (paintMethod == PAINT_METHOD_DEMARCATED);
        // Is this a leaf node or a collapsed node?
        if (node.isLeafNode () || (node.isCollapsed () && paintCollapsed)) {
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
                // Paint a vertical line connecting the node to its
                // parent.
                painter.drawLine (nodeX, nodeY, nodeX, childY, stroke);
                // Paint a triangle if the child node is collapsed,
                // otherwise draw a horizontal line.
                int num = numberOfDescendants (child);
                if (child.isCollapsed () && paintCollapsed && num > 1) {
                    int a = childY - ySpacer + 1;
                    int b = childY + ySpacer - 1;
                    // Paint a triangle.
                    painter.drawLine (nodeX, childY, childX, a, stroke);
                    painter.drawLine (nodeX, childY, childX, b, stroke);
                    painter.drawLine (childX, a, childX, b, stroke);
                }
                else {
                    // Paint a horizontal line.
                    painter.drawLine (
                        nodeX, childY, childX, childY, stroke
                    );
                }
                // Paint the child node.
                paintNode (painter, child);
            }
        }
    }

    /**
     *  Recursive method to paint demarcation bars for the tree.
     *
     *  @param painter The Painter to use.
     *  @param node The current Node to paint.
     *  @param x The X location to start the demarcation bars.
     */
    private void paintDemarcation (Painter painter, Node node, int x) {
        if (node.isCollapsed ()) {
            int fontHeight = painter.fontHeight ();
            int fontWidth = painter.fontWidth ();
            int demarcationStroke = 10;
            int ySpacer = (int)Math.floor (0.5d * fontHeight);
            int minY = Integer.MAX_VALUE;
            int maxY = 0;
            for (Node descendant: node.getDescendants ()) {
                int y = fontHeight + Math.round (
                    descendant.getY ().floatValue () * fontHeight
                );
                if (y < minY) minY = y;
                if (y > maxY) maxY = y;
            }
            int a = minY - ySpacer + 2;
            int b = maxY + ySpacer - 2;
            int c = Math.round (0.5f * (minY + maxY)) + ySpacer - 2;
            painter.drawLine (x, a, x, b, demarcationStroke);
            painter.drawString (node.getName (), x + fontWidth, c);
        }
        else {
            for (Node child: node.getChildren ()) {
                paintDemarcation (painter, child, x);
            }
        }
    }

    /**
     *  Paint a scale bar for the tree.
     *
     *  @param painter The Painter to use.
     *  @param x The X location to start the scale bar.
     *  @param y The Y location to start the scale bar.
     */
    private void paintScaleBar (Painter painter, int x, int y) {
        Float size = 0.01f;
        if (xModifier < 3000) size = 0.05f;
        int stroke = 1;
        int fontHeight = painter.fontHeight ();
        int height = fontHeight / 2;
        int x2 = x + Math.round (size * xModifier);
        int y1 = y - height;
        int y2 = y + height;
        int x3 = (x + x2 - painter.stringWidth (size.toString ())) / 2;
        painter.drawLine (x, y, x2, y, stroke);
        painter.drawLine (x, y1, x, y2, stroke);
        painter.drawLine (x2, y1, x2, y2, stroke);
        painter.drawString (size.toString (), x3, y + fontHeight);
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
        if (root == null || root.getChildren ().size () <= 1) {
            throw new InvalidTreeException (
                "Malformed Newick tree, not enough leaves found."
            );
        }
    }

    /**
     *  Private recursive method to calculate the XY location of the node
     *  and all of its descendants.
     *
     *  @param node The current Node.
     *  @param height The current height.
     */
    private void calculateNodeXY (Node node, double height) {
        boolean paintCollapsed = (paintMethod == PAINT_METHOD_COLLAPSED);
        boolean paintDemarcated = (paintMethod == PAINT_METHOD_DEMARCATED);
        Node parent = node.getParent ();
        double x = 0.0d;
        // The X coordinate is based on the node's distance from its parent.
            x += node.getDistance ();
        // Add the parent's X coordinate.
        if (parent != null) x += parent.getX ();
        // If the node is collapsed, add the descendants distance as well.
        int num = numberOfDescendants (node);
        if (node.isCollapsed () && paintCollapsed && num > 1) {
            double max = node.maximumDistanceFromLeafNode ();
            if (max < 0.01d) max = 0.01d;
            x += max;
        }
        node.setX (x);
        // The Y coordinate is calculated differently for leaf and internal
        // nodes.
        if (node.isLeafNode () || (node.isCollapsed () && paintCollapsed)) {
            node.setY (height);
        }
        else {
            double minY = Double.MAX_VALUE;
            double maxY = 0.0d;
            double h = height;
            for (Node child: node.getChildren ()) {
                calculateNodeXY (child, h);
                int numLeaf = 1;
                if (! (child.isLeafNode () ||
                    (child.isCollapsed () && paintCollapsed))) {
                    numLeaf = numberOfDescendants (child);
                }
                h += numLeaf;
                double childY = child.getY ();
                if (childY < minY) minY = childY;
                if (childY > maxY) maxY = childY;
            }
            node.setY ((minY + maxY) / 2);
        }
    }

    /**
     *  A private recursive method to validate all nodes descended from the
     *  node requested. Valid internal nodes will have two child nodes, valid
     *  leaf nodes will have no children.
     *
     *  @param node The current node being validated.
     */
    private void validateTree (Node node) {
        ArrayList<Node> children = new ArrayList<Node> (node.getChildren ());
        Node parent = node.getParent ();
        if (children.size () == 1) {
            Node child = children.get (0);
            child.setDistance (child.getDistance () + node.getDistance ());
            parent.addChild (child);
            parent.removeChild (node);
            validateTree (child);
        }
        else {
            for (Node child: children) {
                validateTree (child);
            }
        }
    }

    private Node root;
    private int paintMethod = PAINT_METHOD_NORMAL;
    private int xModifier = 5000;

}
