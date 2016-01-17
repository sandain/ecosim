/*
 *    Ecotype Simulation models the sequence diversity within a bacterial
 *    clade as the evolutionary result of net ecotype formation and periodic
 *    selection, yielding a certain number of ecotypes.
 *
 *    Copyright (C) 2016  Jason M. Wood, Montana State University
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


package ecosim.gui;

import ecosim.MasterVariables;
import ecosim.tree.Node;
import ecosim.tree.Tree;

import java.awt.Color;
import java.awt.Dimension;
import java.awt.Font;
import java.awt.FontMetrics;
import java.awt.Graphics;
import java.awt.Graphics2D;
import java.awt.Rectangle;
import java.awt.geom.Rectangle2D;
import java.awt.image.BufferedImage;
import javax.swing.JPanel;
import javax.swing.Scrollable;

/**
 *  Create a JPanel to display a phylogenetic tree.
 *
 *  @author Jason M. Wood
 *  @copyright GNU General Public License
 */
public class TreePane extends JPanel implements Scrollable {

    /**
     *  A custom JPanel used to draw a Tree.
     *
     *  This object creates and renders a tree to a tiled surface made up
     *  of a 2D BufferedImage grid stored in an array named tiles.
     *
     *  @param tree The tree to draw.
     */
    public TreePane (Tree tree) {
        this.tree = tree;
        // Calculate the X,Y location of all nodes in the tree.
        tree.calculateXY ();
        // Calculate the number of tiles needed to accommodate the width
        // of the tree.
        double width = 0.0d;
        BufferedImage tile = new BufferedImage (
            TILE_SIZE, TILE_SIZE, IMAGE_TYPE
        );
        for (Node node: tree.getDescendants ()) {
            Rectangle2D bounds = metrics.getStringBounds (
                node.getName (), tile.getGraphics ()
            );
            double w = bounds.getWidth () + Math.round (
                node.getX ().floatValue () * xModifier
            );
            if (w > width) width = w;
        }
        // Add room for the whitespace.
        width += maxAdvance * 3;
        numTilesWidth = (int)Math.ceil (width / TILE_SIZE);
        // Calculate the number of tiles needed to accommodate the height
        // of the tree.
        double height = maxAscent * tree.size ();
        // Add room for the whitespace.
        height += maxAscent * 2;
        numTilesHeight = (int)Math.ceil (height / TILE_SIZE);
        // Create the tiles array.
        tiles = new BufferedImage[numTilesWidth][numTilesHeight];
        // Set the preferred size.
        setPreferredSize (new Dimension ((int)width, (int)height));
        // Draw the tree, starting with the root node.
        drawNode (tree.getRoot ());
    }

    /**
     *  Paint the tiles.
     *
     *  @param g The Graphics object used to paint the tiles.
     */
    @Override
    public void paintComponent (Graphics g) {
        super.paintComponent (g);
        Rectangle clip = g.getClipBounds ();
        int startX = clip.x - (clip.x % TILE_SIZE);
        int startY = clip.y - (clip.y % TILE_SIZE);
        int stopX = clip.x + clip.width;
        int stopY = clip.y + clip.height;
        for (int x = startX; x < stopX; x += TILE_SIZE) {
            for (int y = startY; y < stopY; y += TILE_SIZE) {
                BufferedImage tile = getTile (x / TILE_SIZE, y / TILE_SIZE);
                g.drawImage (tile, x, y, null);
            }
        }
    }

    /**
     *  Calculate the preferred scrollable viewport size.
     *
     *  @return The preferred scrollable viewport size.
     */
    @Override
    public Dimension getPreferredScrollableViewportSize () {
        return new Dimension (
            numTilesWidth * TILE_SIZE, numTilesHeight * TILE_SIZE
        );
    }

    /**
     *  Calculate the scrollable block increment.
     *
     *  @param r The visible rectangle.
     *  @param o The orientation.
     *  @param d The direction.
     *  @return The block increment.
     */
    @Override
    public int getScrollableBlockIncrement (Rectangle r, int o, int d) {
        return TILE_SIZE * Math.max (1, numTilesHeight - 1);
    }

    /**
     *  The max ascent is used as the scrollable unit increment.
     *
     *  @param r The visible rectangle.
     *  @param o The orientation.
     *  @param d The direction.
     *  @return The unit increment.
     */
    @Override
    public int getScrollableUnitIncrement (Rectangle r, int o, int d) {
        return maxAscent;
    }

    /**
     *  Don't force the height of this JPanel to track the viewport height.
     *
     *  @return False
     */
    @Override
    public boolean getScrollableTracksViewportHeight () {
        return false;
    }

    /**
     *  Don't force the width of this JPanel to track the viewport width.
     *
     *  @return False
     */
    @Override
    public boolean getScrollableTracksViewportWidth () {
        return false;
    }

    /**
     *  Get a tile and initialize it if it is empty.
     *
     *  @param x The X position of the tile.
     *  @param y The Y position of the tile.
     *  @return The tile.
     */
    private BufferedImage getTile (int x, int y) {
        // Initialize the tile if it is empty.
        if (tiles[x][y] == null) {
            tiles[x][y] = new BufferedImage (
                TILE_SIZE, TILE_SIZE, IMAGE_TYPE
            );
            Graphics2D g2 = (Graphics2D)tiles[x][y].getGraphics ();
            g2.setBackground (Color.WHITE);
            g2.clearRect (0, 0, TILE_SIZE, TILE_SIZE);
        }
        return tiles[x][y];
    }

    /**
     *  Check whether or not the third point three is between the other two.
     *
     *  @param x1 X value for the first point.
     *  @param y1 Y value for the first point.
     *  @param x2 X value for the second point.
     *  @param y2 Y value for the second point.
     *  @param x3 X value for the third point.
     *  @param y3 Y value for the third point.
     *  @return True if point three is between the other two.
     */
    private boolean isBetween (
        int x1, int y1, int x2, int y2, int x3, int y3
    ) {
        boolean between = false;
        if (x3 != Integer.MAX_VALUE && y3 != Integer.MAX_VALUE) {
            boolean x = (x1 < x3 && x2 > x3) || (x2 < x3 && x1 > x3);
            boolean y = (y1 < y3 && y2 > y3) || (y2 < y3 && y1 > y3);
            if (x1 == x2 && y1 != y2) between = y;
            else if (y1 == y2 && x1 != x2) between = x;
            else between = x && y;
        }
        return between;
    }

    /**
     *  Calculate the euclidean distance between two points using the
     *  Pythagorean theorem.
     *
     *  @param x1 X value for the first point.
     *  @param y1 Y value for the first point.
     *  @param x2 X value for the second point.
     *  @param y2 Y value for the second point.
     *  @return The distance between the two points.
     */
    private double distanceBetween (int x1, int y1, int x2, int y2) {
        double a2 = (x1 - x2) * (x1 - x2);
        double b2 = (y1 - y2) * (y1 - y2);
        return Math.sqrt (a2 + b2);
    }

    /**
     *  Calculate X for the given value of Y on the line defined by its
     *  slope and y-intercept.
     *
     *  @param y The Y value.
     *  @param b The Y-intercept of the line.
     *  @param m The slope of the line.
     *  @return The X value.
     */
    private int calculateX (int y, double b, double m) {
        Long x = Math.round ((y - b) / m);
        return x.intValue ();
    }

    /**
     *  Calculate Y for the given value of X on the line defined by its
     *  slope and Y-intercept.
     *
     *  @param x The X value.
     *  @param b The Y-intercept of the line.
     *  @param m The slope of the line.
     *  @return The Y value.
     */
    private int calculateY (int x, double b, double m) {
        Long y = Math.round (x * m + b);
        return y.intValue ();
    }

    /**
     *  A private recursive method to draw a line between two points.
     *
     *  @param x1 X value for the first point.
     *  @param y1 Y value for the first point.
     *  @param x2 X value for the second point.
     *  @param y2 Y value for the second point.
     */
    private void drawLine (int x1, int y1, int x2, int y2) {
        // Get the tile containing the first point.
        int tileX = x1 / TILE_SIZE;
        int tileY = y1 / TILE_SIZE;
        BufferedImage tile = getTile (tileX, tileY);
        // Calculate the bounding box.
        int startX = tileX * TILE_SIZE - 1;
        int startY = tileY * TILE_SIZE - 1;
        int endX = startX + TILE_SIZE + 1;
        int endY = startY + TILE_SIZE + 1;
        // Calculate the slope and y-intercept of the line.
        Double m = (double)(y1 - y2) / (x1 - x2);
        Double b = y1 - m * x1;
        // Calculate possible places the line could cross the tile bounds.
        int northX = Integer.MAX_VALUE;
        int northY = Integer.MAX_VALUE;
        int southX = Integer.MAX_VALUE;
        int southY = Integer.MAX_VALUE;
        int westX = Integer.MAX_VALUE;
        int westY = Integer.MAX_VALUE;
        int eastX = Integer.MAX_VALUE;
        int eastY = Integer.MAX_VALUE;
        if (m.isInfinite ()) {
            // Line is vertical.
            northX = x1;
            northY = startY;
            southX = x1;
            southY = endY;
        }
        else if (Math.abs (m) < MasterVariables.EPSILON) {
            // Line is horizontal.
            westX = startX;
            westY = y1;
            eastX = endX;
            eastY = y1;
        }
        else {
            // Line has a slope.
            northX = calculateX (startY, b, m);
            northY = startY;
            southX = calculateX (endY, b, m);
            southY = endY;
            westX = startX;
            westY = calculateY (startX, b, m);
            eastX = endX;
            eastY = calculateY (endX, b, m);
        }
        // Calculate the euclidean distance between points one and two.
        double distance = distanceBetween (x1, y1, x2, y2);
        // Find the point (north, south, west, east) that is closest to the
        // first point.
        int pointX = x2;
        int pointY = y2;
        if (isBetween (x1, y1, x2, y2, northX, northY)) {
            double d = distanceBetween (x1, y1, northX, northY);
            if (d < distance) {
                pointX = northX;
                pointY = northY;
                distance = d;
            }
        }
        if (isBetween (x1, y1, x2, y2, southX, southY)) {
            double d = distanceBetween (x1, y1, southX, southY);
            if (d < distance) {
                pointX = southX;
                pointY = southY;
                distance = d;
            }
        }
        if (isBetween (x1, y1, x2, y2, westX, westY)) {
            double d = distanceBetween (x1, y1, westX, westY);
            if (d < distance) {
                pointX = westX;
                pointY = westY;
                distance = d;
            }
        }
        if (isBetween (x1, y1, x2, y2, eastX, eastY)) {
            double d = distanceBetween (x1, y1, eastX, eastY);
            if (d < distance) {
                pointX = eastX;
                pointY = eastY;
                // no reason to save the distance.
            }
        }
        // If a point was found between p1 and p2, recurse on the line
        // segment (point -> p2).
        if (pointX != x2 || pointY != y2) {
            drawLine (pointX, pointY, x2, y2);
        }
        // Translate the x,y coordinates into tile space and draw the line.
        int x1pr = x1 - tileX * TILE_SIZE;
        int y1pr = y1 - tileY * TILE_SIZE;
        int x2pr = pointX - tileX * TILE_SIZE;
        int y2pr = pointY - tileY * TILE_SIZE;
        Graphics2D g2 = (Graphics2D)tile.getGraphics ();
        g2.setColor (Color.BLACK);
        // Draw each line twice, from x1pr,1ypr -> x2pr,y2pr and
        // x2pr,y2pr -> x1pr,1ypr to work around a bug in Java.
        g2.drawLine (x1pr, y1pr, x2pr, y2pr);
        g2.drawLine (x2pr, y2pr, x1pr, y1pr);
    }

    /**
     *  Private recursive method to draw a String.
     *
     *  @param str The String to draw.
     *  @param tileX The X location in tile space.
     *  @param tileY The Y location in tile space.
     *  @param x The X location to draw the String.
     *  @param y The Y location to draw the String.
     */
    private void drawString (String str, int tileX, int tileY, int x, int y) {
        BufferedImage tile = getTile (tileX, tileY);
        Graphics2D g2 = (Graphics2D)tile.getGraphics ();
        g2.setColor (Color.BLACK);
        g2.drawString (str, x, y);
        // Get the bounds of the string.
        Rectangle2D bounds = metrics.getStringBounds (str, g2);
        int strHeight = (int)bounds.getHeight ();
        int strWidth = (int)bounds.getWidth ();
        // Catch strings being drawn on the north border of a tile.
        if (y >= 0 && y <= strHeight) {
            drawString (str, tileX, tileY - 1, x, TILE_SIZE + y);
        }
        // Catch strings being drawn on the east border of a tile.
        if (x + maxAdvance + strWidth > TILE_SIZE) {
            drawString (str, tileX + 1, tileY, x - TILE_SIZE, y);
        }
    }

    /**
     *  Private recursive method to draw a node and all of its descendants.
     *
     *  @param node The node to draw.
     */
    private void drawNode (Node node) {
        int nodeX = maxAdvance + Math.round (
            node.getX ().floatValue () * xModifier
        );
        int nodeY = maxAscent + Math.round (
            node.getY ().floatValue () * yModifier
        );
        if (node.isLeafNode ()) {
            // Add some space before the node's name.
            nodeX += maxAdvance;
            // Center the node's name.
            nodeY += Math.floor (maxAscent / 2);
            // Translate the X,Y coordinates of the node into tile space and
            // draw it's name.
            int tileX = nodeX / TILE_SIZE;
            int tileY = nodeY / TILE_SIZE;
            drawString (
                node.getName (), tileX, tileY,
                nodeX - tileX * TILE_SIZE,
                nodeY - tileY * TILE_SIZE
            );
        }
        else {
            for (Node child: node.getChildren ()) {
                int childX = maxAdvance + Math.round (
                    child.getX ().floatValue () * xModifier
                );
                int childY = maxAscent + Math.round (
                    child.getY ().floatValue () * yModifier
                );
                // Add the horizontal line.
                drawLine (nodeX, childY, childX, childY);
                // Add the vertical line.
                drawLine (nodeX, nodeY, nodeX, childY);
                // Draw the child node.
                drawNode (child);
            }
        }
    }

    private int numTilesWidth = 1;
    private int numTilesHeight = 1;

    private final BufferedImage[][] tiles;

    private final Font font = new Font ("monospaced", Font.PLAIN, 12);
    private final FontMetrics metrics = getFontMetrics (font);

    private final int maxAdvance = metrics.getMaxAdvance ();
    private final int maxAscent = metrics.getMaxAscent ();

    private int xModifier = 1000;
    private int yModifier = maxAscent;

    private final int TILE_SIZE = 1000;

    private final int IMAGE_TYPE = BufferedImage.TYPE_BYTE_BINARY;

    private Tree tree;
    
}
