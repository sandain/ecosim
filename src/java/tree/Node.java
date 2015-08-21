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

import ecosim.Heapsorter;
import ecosim.MasterVariables;

import java.util.ArrayList;

/**
 *  A node potentially contains a name, a parent node, a distance from the
 *  parent node, and a list of all child nodes.
 *
 *  @author Andrew Warner
 *  @author Jason M. Wood
 *  @copyright GNU General Public License
 */
public class Node implements Comparable<Node> {

    /**
     *  Constructor for objects of class Node.
     *
     *  @param name The name of this Node.
     *  @param distance The distance from the parent node.
     *  @param parent The parent Node of this node.
     *  @param children A list of TreeNodes which are the children of this
     *  node.
     */
    public Node (String name, Double distance,
        Node parent, ArrayList<Node> children) {
        this.name = name;
        this.distance = distance;
        this.parent = parent;
        this.children = children;
        outgroup = false;
        x = 0;
        y = 0;
    }

    /**
     *  Default constructor for object of class Node.  This node
     *  will have no name, distance, parent, or children unless defined
     *  later.
     */
    public Node () {
        this ("", 0.0d, null, new ArrayList<Node> ());
    }

    /**
     *  Compare nodes.
     */
    public int compareTo (Node other) {
        // Compare the names of the two nodes.
        int c = other.getName ().compareTo (name);
        if (c != 0) return c;
        // Check if the two nodes are root nodes.
        c = other.isRootNode ().compareTo (isRootNode ());
        if (c != 0) return c;
        if (! isRootNode ()) {
            // Check if the parents of the two nodes are root nodes.
            Node otherParent = other.getParent ();
            c = otherParent.isRootNode ().compareTo (parent.isRootNode ());
            if (c != 0) return c;
            // Compare the distances of the two nodes from their parents.
            if (parent.isRootNode ()) {
                // When the parent is the root node, compare the sums of the
                // distance of all the root node's children. 
                Double a = 0.0d;
                for (Node child: parent.getChildren ()) {
                    a += child.getDistance ();
                }
                Double b = 0.0d;
                for (Node child: otherParent.getChildren ()) {
                    b += child.getDistance ();
                }
                c = 0;
                if (b > a + MasterVariables.EPSILON) c = -1;
                if (a > b + MasterVariables.EPSILON) c = 1;
                if (c != 0) return c;
            }
            else {
                // When the parent is not the root node, compare just the
                // distances.
                Double otherDist = other.getDistance ();
                c = 0;
                if (distance > otherDist + MasterVariables.EPSILON) c = -1;
                if (otherDist > distance + MasterVariables.EPSILON) c = 1;
                if (c != 0) return c;
            }
        }
        // Compare the children of the two nodes.
        for (Node a: children) {
            // Compare each child with the children of the other node.
            for (Node b: other.getChildren ()) {
                c = b.compareTo (a);
                if (c == 0) break;
            }
            if (c != 0) return c;
        }
        // Nodes are equal.
        return 0;
    }

    /**
     *  Set the name of this node.
     *
     *  @param name The name of this node.
     */
     public void setName (String name) {
         this.name = name;
     }

    /**
     *  Set the distance from the parent node to this node.
     *
     *  @param distance The distance to the parent.
     */
    public void setDistance (Double distance) {
        this.distance = distance;
    }

    /**
     *  Let this node know whether or not it is the outgroup.
     *
     *  @param outgroup Whether or not this node is the outgroup.
     */
    public void setOutgroup (Boolean outgroup) {
        this.outgroup = outgroup;
    }

    /**
     *  Set the parent of this node.
     *
     *  @param parent The parent of this node.
     */
    public void setParent (Node parent) {
        this.parent = parent;
    }

    /**
     *  Add a child to this node.
     *
     *  @param child The child to add.
     */
    public void addChild (Node child) {
        child.setParent (this);
        children.add (child);
    }

    /**
     *  Sort the children of this node using the default sorter.
     */
    public void sortChildren () {
        sortChildren (new Heapsorter<Node> ());
    }

    /**
     *  Sort the children of this node using the provided sorter.
     *
     *  @param sorter The sorter to use.
     */
    public void sortChildren (Heapsorter<Node> sorter) {
        for (Node child: children) {
            child.sortChildren (sorter);
        }
        sorter.sort (children);
    }

    /**
     *  Returns the name of this node.
     *
     *  @return String containing the name of this node.
     */
    public String getName () {
        return name;
    }

    /**
     *  Get the distance from the parent node to this node.
     *
     *  @return Double containing the distance to the parent.
     */
    public Double getDistance () {
        return distance;
    }

    /**
     *  Returns the parent of this node.
     *
     *  @return Node containing the parent of this node.
     */
    public Node getParent () {
        return parent;
    }

    /**
     *  Returns the children of this node.
     *
     *  @return ArrayList<Node> with the children of this node.
     */
    public ArrayList<Node> getChildren () {
        return children;
    }

    /**
     *  Removes a child of this node.
     *
     *  @param child The child Node to remove.
     */
    public void removeChild (Node child) {
        children.remove (child);
        child.setParent (null);
    }

    /**
     *  Returns a Boolean stating whether this node is a leaf or not.
     *
     *  @return A Boolean stating whether this node is a leaf or not.
     */
    public Boolean isLeafNode () {
        return children.isEmpty ();
    }

    /**
     *  Returns a Boolean stating whether this node is the outgroup or not.
     *
     *  @return A Boolean stating whether this node is the outgroup or not.
     */
    public Boolean isOutgroup () {
        return outgroup;
    }

    /**
     *  Returns whether this node is the root node or not.
     *
     *  @return True if this node is the root node.
     */
    public Boolean isRootNode () {
        return parent == null;
    }

    /**
     *  Returns the distance of this node from the root node.
     *
     *  @return The distance of this node from the root node.
     */
    public Double distanceFromRootNode () {
        Double distanceFromRoot = 0.0d;
        if (parent != null) {
            distanceFromRoot = distance + parent.distanceFromRootNode ();
        }
        return distanceFromRoot;
    }

    /**
     *  Returns the maximum distance of this node from a leaf node.
     *
     *  @return The maximum distance of this node from a leaf node.
     */
    public Double maximumDistanceFromLeafNode () {
        Double maximumDistance = 0.0d;
        for (Node child: children) {
            Double childDistance =
                child.maximumDistanceFromLeafNode () + child.getDistance ();
            if (childDistance > maximumDistance) {
                maximumDistance = childDistance;
            }
        }
        return maximumDistance;
    }

    /**
     *  Returns the minimum distance of this node from a leaf node.
     *
     *  @return The minimum distance of this node from a leaf node.
     */
    public Double minimumDistanceFromLeafNode () {
        Double minimumDistance = Double.MAX_VALUE;
        for (Node child: children) {
            Double childDistance =
                child.minimumDistanceFromLeafNode () + child.getDistance ();
            if (childDistance < minimumDistance) {
                minimumDistance = childDistance;
            }
        }
        return minimumDistance;
    }

    /**
     *  Returns an array of nodes that are descendants of this node.
     *
     *  @return The descendants of this node.
     */
    public ArrayList<Node> getDescendants () {
        ArrayList<Node> descendants = new ArrayList<Node> ();
        for (Node child: children) {
            if (child.isLeafNode ()) {
                descendants.add (child);
            }
            else {
                descendants.addAll (child.getDescendants ());
            }
        }
        return descendants;
    }

    /**
     *  Returns the number of living descendants of this node.
     *
     *  @return The number of living descendants of this node.
     */
    public int numberOfDescendants () {
        int descendants = 0;
        for (Node child: children) {
            if (child.isLeafNode ()) {
                descendants ++;
            }
            else {
                descendants += child.numberOfDescendants ();
            }
        }
        return descendants;
    }

    /**
     *  Returns this Node as a Newick formatted String.
     *
     *  @return A Newick formatted String representing this Node.
     */
    public String toString () {
        String newick = "";
        for (int i = 0; i < children.size (); i ++) {
            if (i != 0) {
                newick += ",";
            }
            newick += children.get (i);
        }
        if (newick.length () > 0) {
            newick = "(" + newick + ")";
        }
        if (isLeafNode ()) {
            newick += String.format ("%s:%.5f", name, distance);
        }
        else {
            newick += String.format (":%.5f", distance);
        }
        if (parent == null) {
            newick += ";";
        }
        return newick;
    }

    public int getX () {
        return x;
    }

    public int getY () {
        return y;
    }

    public void setX (int x) {
        this.x = x;
    }

    public void setY (int y) {
        this.y = y;
    }

    /**
     *  The name of this node.
     */
    private String name;

    /**
     *  The distance of this node from it's parent.
     */
    private Double distance;

    /**
     *  The parent of this node.
     */
    private Node parent;

    /**
     *  A list of children of this node.
     */
    private ArrayList<Node> children;

    /**
     *  Whether or not this node is the outgroup.
     */
    private boolean outgroup;

    /**
     *  The X coordinate of this node.  Used only when creating a graphical
     *  representation of this node.
     */
    private int x;

    /**
     *  The Y coordinate of this node.  Used only when creating a graphical
     *  representation of this node.
     */
    private int y;

}
