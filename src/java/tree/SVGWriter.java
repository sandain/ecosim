/*
 *    Ecotype Simulation models the sequence diversity within a bacterial
 *    clade as the evolutionary result of net ecotype formation and periodic
 *    selection, yielding a certain number of ecotypes.
 *
 *    Copyright (C) 2015  Jason M. Wood, Montana State University
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

import java.io.BufferedWriter;
import java.io.IOException;
import java.io.Writer;
import java.util.ArrayList;

/**
 *  Converts a Node based tree into a SVG formatted image.
 *
 *  @author Jason M. Wood
 *  @copyright GNU General Public License
 */

public class SVGWriter extends BufferedWriter {

    public SVGWriter (Writer writer) {
        super (writer);
    }

    /**
     *  Save the Node based tree data in this object to a SVG formatted file.
     *
     *  @param tree The tree to write.
     *  @return True if the save was a success, False otherwise.
     */
    public void write (Tree tree) throws IOException {
        // Calculate the X,Y location of all nodes in the tree.
        tree.calculateXY ();
        // Calculate the height and width.
        int height = fontHeight * tree.size () + yOffset;
        int max = 0;
        for (Node node: tree.getDescendants ()) {
            String name = node.getName ();
            int labelWidth = (name.length () + 1) * fontWidth;
            int x = labelXOffset + labelWidth + Math.round (
                node.getX ().floatValue () * xModifier
            );
            if (x > max) max = x;
        }        
        int width = max;
        // Add the XML information.
        writeln (
            "<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"no\"?>"
        );
        // Add the SVG doctype.
        writeln (
            "<!DOCTYPE svg PUBLIC \"-//W3C//DTD SVG 1.1//EN\" " +
            "\"http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd\">"
        );
        // Start the SVG document.
        writeln (
            "<svg width=\"" + width + "\" height=\"" + height + "\" " +
            "version=\"1.1\" xmlns=\"http://www.w3.org/2000/svg\">"
        );
        // Add some CSS style to the SVG document.
        writeln ("<defs>");
        writeln ("<style type=\"text/css\"><![CDATA[");
        writeln ("  line {");
        writeln ("    stroke: black;");
        writeln ("    stroke-width: 1;");
        writeln ("  }");
        writeln ("  text {");
        writeln ("    font-family: monospace;");
        writeln ("    font-size: " + fontHeight + "px;");
        writeln ("    stroke-width: 0;");
        writeln ("    fill: blue;");
        writeln ("  }");
        writeln ("]]></style>");
        writeln ("</defs>");
        // Add the tree to the SVG document.
        nodeToSVG (tree.getRoot ());
        // End the SVG document.
        writeln ("</svg>");
    }

    /**
     *  Private recursive method to convert a Node and all of its descendants
     *  into a SVG representation.
     */
    private void nodeToSVG (Node node) throws IOException {
        int nodeX = xOffset + Math.round (
            node.getX ().floatValue () * xModifier
        );
        int nodeY = yOffset + Math.round (
            node.getY ().floatValue () * yModifier
        );
        if (node.isLeafNode () || node.isCollapsed ()) {
            // Add the name of the node to the SVG as a text element.
            String name = node.getName ();
            writeln (String.format (
                "<text x=\"%d\" y=\"%d\" textLength=\"%d\">%s</text>",
                nodeX + labelXOffset,
                nodeY + labelYOffset,
                name.length () * fontWidth,
                name
            ));
        }
        else {
            for (Node child: node.getChildren ()) {
                int childX = xOffset + Math.round (
                    child.getX ().floatValue () * xModifier
                );
                int childY = yOffset + Math.round (
                    child.getY ().floatValue () * yModifier
                );
                // Draw a triangle if the child node is collapsed, otherwise
                // draw a horizontal line.
                if (child.isCollapsed ()) {
                    int a = childY - fontHeight / 2 + 1;
                    int b = childY + fontHeight / 2 - 1;
                    // Draw a triangle.
                    writeln (String.format (
                        "<line x1=\"%d\" y1=\"%d\" x2=\"%d\" y2=\"%d\"/>",
                        nodeX, childY, childX, a
                    ));
                    writeln (String.format (
                        "<line x1=\"%d\" y1=\"%d\" x2=\"%d\" y2=\"%d\"/>",
                        nodeX, childY, childX, b
                    ));
                    writeln (String.format (
                        "<line x1=\"%d\" y1=\"%d\" x2=\"%d\" y2=\"%d\"/>",
                        childX, a, childX, b
                    ));
                }
                else {
                    // Draw a horizontal line.
                    writeln (String.format (
                        "<line x1=\"%d\" y1=\"%d\" x2=\"%d\" y2=\"%d\"/>",
                        nodeX, childY, childX, childY
                    ));
                }
                // Draw a vertical line connecting the node to it's parent.
                writeln (String.format (
                    "<line x1=\"%d\" y1=\"%d\" x2=\"%d\" y2=\"%d\"/>",
                    nodeX, nodeY, nodeX, childY
                ));
                // Add the SVG data from the child node.
                nodeToSVG (child);
            }
        }
    }

    /**
     *  Private method to write a string of text with a line separator.
     *
     *  @param line The line of text.
     */
    private void writeln (String line) throws IOException {
        write (line + System.getProperty ("line.separator"));
    }

    private int fontHeight = 12;
    private int fontWidth = 7;
    private int xOffset = 5;
    private int yOffset = 5;
    private int labelXOffset = 3;
    private int labelYOffset = 5;
    private int xModifier = 1000;
    private int yModifier = fontHeight;

}
