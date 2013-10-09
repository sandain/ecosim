/*
 *    Ecotype Simulation models the sequence diversity within a bacterial clade
 *    as the evolutionary result of net ecotype formation and periodic
 *    selection, yielding a certain number of ecotypes.
 *
 *    Copyright (C) 2009       Andrew Warner, Wesleyan University
 *    Copyright (C) 2009-2013  Jason M. Wood, Montana State University
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
 *  A Newick tree node potentially contains a name, a parent node,
 *  a distance from the parent node, and a list of all child nodes.
 *
 *  @author Andrew Warner
 *  @author Jason M. Wood
 *  @copyright GNU General Public License
 */
public class NewickTreeNode {

    /**
     *  Constructor for objects of class NewickTreeNode.
     *
     *  @param name The name of this treenode.
     *  @param distance The distance from the parent node.
     *  @param parent The parent NewickTreeNode of this node.
     *  @param children A list of NewickTreeNodes which are the children of this
     *  node.
     */
    public NewickTreeNode(String name, double distance, NewickTreeNode parent,
        ArrayList<NewickTreeNode> children) {
        this.name = name;
        this.distance = distance;
        this.parent = parent;
        this.children = children;
    }

    /**
     *  Default constructor for object of class NewickTreeNode.  This node
     *  will have no name, distance, parent, or children unless defined later.
     */
    public NewickTreeNode() {
        this ("", 0.0, null, new ArrayList<NewickTreeNode>());
    }

    /**
     *  Set the name of this node.
     *
     *  @param name The name of this node.
     */
     public void setName(String name) {
         this.name = name;
     }

    /**
     *  Set the distance from the parent node to this node.
     *
     *  @param distance The distance to the parent.
     */
    public void setDistance(double distance) {
        this.distance = distance;
    }

    /**
     *  Set the parent of this node.
     *
     *  @param parent The parent of this node.
     */
    public void setParent(NewickTreeNode parent) {
        this.parent = parent;
    }

    /**
     *  Add a child to this node.
     *
     *  @param child The child to add.
     */
    public void addChild(NewickTreeNode child) {
        child.setParent(this);
        children.add(child);
    }

    /**
     *  Returns the name of this node.
     *
     *  @return String containing the name of this node.
     */
    public String getName() {
        return name;
    }

    /**
     *  Get the distance from the parent node to this node.
     *
     *  @return double containing the distance to the parent.
     */
    public double getDistance() {
        return distance;
    }

    /**
     *  Returns the parent of this node.
     *
     *  @return NewickTreeNode containing the parent of this node.
     */
    public NewickTreeNode getParent() {
        return parent;
    }

    /**
     *  Returns the children of this node.
     *
     *  @return ArrayList<NewickTreeNode> containing the children of this node.
     */
    public ArrayList<NewickTreeNode> getChildren() {
        return children;
    }

    /**
     *  Removes a child of this node.
     *
     *  @param child The child to remove.
     */
    public void removeChild(NewickTreeNode child) {
        children.remove(child);
    }

    /**
     *  Returns a boolean stating whether this node is a leaf or not.
     *
     *  @return A boolean stating whether this node is a leaf or not.
     */
    public boolean isLeafNode() {
        return children.isEmpty();
    }

    /**
     *  Returns the distance of this node from the root node.
     *
     *  @return The distance of this node from the root node.
     */
    public double distanceFromRootNode() {
        double distanceFromRoot = 0.0;
        if (parent != null) {
            distanceFromRoot = distance + parent.distanceFromRootNode();
        }
        return distanceFromRoot;
    }

    /**
     *  Returns whether this node is the root node or not.
     *
     *  @return True if this node is the root node.
     */
    public boolean isRootNode() {
          return parent == null;
    }


    /**
     *  Returns an array of nodes that are descendants of this node.
     *
     *  @return The descendants of this node.
     */
    public ArrayList<NewickTreeNode> getDescendants() {
        ArrayList<NewickTreeNode> descendants = new ArrayList<NewickTreeNode>();
        if (children.size() > 0) {
            for (int i = 0; i < children.size(); i ++) {
                NewickTreeNode child = children.get(i);
                if (child.isLeafNode()) {
                    descendants.add(child);
                }
                else {
                    descendants.addAll(child.getDescendants());
                }
            }
        }
        return descendants;
    }

    /**
     *  Returns the number of living descendants of this node.
     *
     *  @return The number of living descendants of this node.
     */
    public int numberOfDescendants() {
        int descendants = 0;
        // If this node has children then the children are the descendants,
        // otherwise this node is a desendant.
        if (children.size() > 0) {
            for (int i = 0; i < children.size(); i ++) {
                if (children.get(i).isLeafNode()) {
                    descendants ++;
                }
                else {
                    descendants += children.get(i).numberOfDescendants();
                }
            }
        }
        else {
            descendants = 1;
        }
        return descendants;
    }

    /**
     *  Returns this node formatted as a String.
     *
     *  @return String containing the name and distance of this node.
     */
    public String toString() {
        String newick = "";
        for (int i = 0; i < children.size(); i ++) {
            if (i != 0) {
                newick += ",";
            }
            newick += children.get(i);
        }
        if (newick.length() > 0) {
            newick = "(" + newick + ")";
        }
        newick += name + ":" + distance;
        if (parent == null) {
            newick += ";";
        }
        return newick;
    }

    /**
     *  The name of this node.
     */
    private String name;

    /**
     *  The distance of this node from it's parent.
     */
    private double distance;

    /**
     *  The parent of this node.
     */
    private NewickTreeNode parent;

    /**
     *  A list of children of this node.
     */
    private ArrayList<NewickTreeNode> children;

}
