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
import ecosim.api.Painter;

import java.awt.BasicStroke;
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
 *  Create a tiled JPanel with a custom Painter.
 *
 *  @author Jason M. Wood
 *  @copyright GNU General Public License
 */
public class TiledPainter extends JPanel implements Scrollable, Painter {

    /**
     *  A custom tile-based JPanel using a custom Painter.
     *
     *  This object creates a tiled surface made up of a 2D BufferedImage grid
     *  stored in an array named tiles.
     */
    public TiledPainter () {
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
     *  The font height is used as the scrollable unit increment.
     *
     *  @param r The visible rectangle.
     *  @param o The orientation.
     *  @param d The direction.
     *  @return The unit increment.
     */
    @Override
    public int getScrollableUnitIncrement (Rectangle r, int o, int d) {
        return fontHeight;
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
     *  Start a Tiled Painter of the requested size.
     *
     *  @param width The width of the Painter.
     *  @param height The height of the Painter.
     */
    public void start (int width, int height) {
        numTilesWidth = (int)Math.ceil ((float)width / TILE_SIZE);
        numTilesHeight = (int)Math.ceil ((float)height / TILE_SIZE);
        // Create the tiles array.
        tiles = new BufferedImage[numTilesWidth][numTilesHeight];
        for (int i = 0; i < numTilesWidth; i ++) {
            for (int j = 0; j < numTilesHeight; j ++) {
                // Initialize the tile.
                tiles[i][j] = new BufferedImage (
                    TILE_SIZE, TILE_SIZE, IMAGE_TYPE
                );
                Graphics2D g2 = (Graphics2D)tiles[i][j].getGraphics ();
                g2.setBackground (Color.WHITE);
                g2.clearRect (0, 0, TILE_SIZE, TILE_SIZE);
            }
        }
        // Set the preferred size.
        setPreferredSize (new Dimension (width, height));
    }

    /**
     *  End the Tiled Painter.
     */
    public void end () {
    }

    /**
     *  A recursive method to draw a line between two points.
     *
     *  @param x1 X value for the first point.
     *  @param y1 Y value for the first point.
     *  @param x2 X value for the second point.
     *  @param y2 Y value for the second point.
     *  @param stroke The stroke width of the line.
     */
    public void drawLine (int x1, int y1, int x2, int y2, int stroke) {
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
            drawLine (pointX, pointY, x2, y2, stroke);
        }
        // Translate the x,y coordinates into tile space and draw the line.
        int x1pr = x1 - tileX * TILE_SIZE;
        int y1pr = y1 - tileY * TILE_SIZE;
        int x2pr = pointX - tileX * TILE_SIZE;
        int y2pr = pointY - tileY * TILE_SIZE;
        Graphics2D g2 = (Graphics2D)tile.getGraphics ();
        g2.setColor (Color.BLACK);
        g2.setStroke (new BasicStroke (stroke, BasicStroke.CAP_BUTT, BasicStroke.JOIN_MITER));
        // Draw each line twice, from x1pr,1ypr -> x2pr,y2pr and
        // x2pr,y2pr -> x1pr,1ypr to work around a bug in Java.
        g2.drawLine (x1pr, y1pr, x2pr, y2pr);
        g2.drawLine (x2pr, y2pr, x1pr, y1pr);
    }

    /**
     *  Draw a String at the provided X,Y location.
     *
     *  @param str The String to paint.
     *  @param x The X location to paint the String.
     *  @param y The Y location to paint the String.
     */
    public void drawString (String str, int x, int y) {
        // Translate the X,Y coordinates into tile space and draw the String.
        int tileX = (int)Math.floor ((float)x / TILE_SIZE);
        int tileY = (int)Math.floor ((float)y / TILE_SIZE);
        int xpr = x - tileX * TILE_SIZE;
        int ypr = y - tileY * TILE_SIZE;
        drawString (str, tileX, tileY, xpr, ypr);
    }

    /**
     *  Return the width of the font.
     *
     *  @return The font width.
     */
    public int fontWidth () {
        return fontWidth;
    }

    /**
     *  Return the height of the font.
     *
     *  @return The font height.
     */
    public int fontHeight () {
        return fontHeight;
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
        if (x + fontWidth + strWidth > TILE_SIZE) {
            drawString (str, tileX + 1, tileY, x - TILE_SIZE, y);
        }
    }

    /**
     *  Get a tile and initialize it if it is empty.
     *
     *  @param x The X position of the tile.
     *  @param y The Y position of the tile.
     *  @return The tile.
     */
    private BufferedImage getTile (int x, int y) {
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

    private int numTilesWidth = 1;
    private int numTilesHeight = 1;

    private BufferedImage[][] tiles;

    private final Font font = new Font ("monospaced", Font.PLAIN, 12);
    private final FontMetrics metrics = getFontMetrics (font);

    private final int fontHeight = metrics.getMaxAscent ();
    private final int fontWidth = metrics.getMaxAdvance ();

    private final int TILE_SIZE = 1000;

    private final int IMAGE_TYPE = BufferedImage.TYPE_BYTE_BINARY;
    
}
