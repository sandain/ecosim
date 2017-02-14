/*
 *    Ecotype Simulation models the sequence diversity within a bacterial
 *    clade as the evolutionary result of net ecotype formation and periodic
 *    selection, yielding a certain number of ecotypes.
 *
 *    Copyright (C) 2016  Jason M. Wood, Montana State University
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

import java.io.IOException;
import java.io.File;
import java.io.FileWriter;

/**
 *  Paints to a SVG formatted image.
 *
 *  @author Jason M. Wood
 *  @copyright GNU General Public License
 */

public class SVGPainter implements Painter {

    /**
     *  Create a Painter for painting to SVG formatted images.
     *
     *  @param file The file to paint the SVG image to.
     */
    public SVGPainter (File file) {
        try {
            writer = new FileWriter (file);
        }
        catch (IOException e) {
            System.out.println ("Error opening SVG file: " + e);
        }
    }

    /**
     *  Initialize the SVG Painter.
     *
     *  @param width The width of the SVG document.
     *  @param height The height of the SVG document.
     */
    public void start (int width, int height) {
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
        writeln ("    stroke: " + strokeColor + ";");
        writeln ("    stroke-width: " + strokeWidth + ";");
        writeln ("  }");
        writeln ("  text {");
        writeln ("    font-family: monospace;");
        writeln ("    font-size: " + fontHeight + "px;");
        writeln ("    stroke-width: 0;");
        writeln ("    fill: " + fontColor + ";");
        writeln ("  }");
        writeln ("]]></style>");
        writeln ("</defs>");
    }

    /**
     *  End the SVG Painter.
     */
    public void end () {
        try {
            try {
                // End the SVG document.
                writer.write ("</svg>");
                writer.write (System.getProperty ("line.separator"));
            }
            finally {
                // Close the file IO.
                writer.close ();
            }
        }
        catch (IOException e) {
            System.out.println ("Error closing SVG file: " + e);
            e.printStackTrace ();
        }
    }

    /**
     *  Draw a line between two provided points.
     *
     *  @param x1 The start X location of the line to draw.
     *  @param y1 The start Y location of the line to draw.
     *  @param x2 The end X location of the line to draw.
     *  @param y2 The end Y location of the line to draw.
     *  @param stroke The stroke width of the line.
     */
    public void drawLine (int x1, int y1, int x2, int y2, int stroke) {
        if (stroke == strokeWidth) {
            writeln (String.format (
                "<line x1=\"%d\" y1=\"%d\" x2=\"%d\" y2=\"%d\"/>",
                x1, y1, x2, y2
            ));
        }
        else {
            writeln (String.format (
                "<line x1=\"%d\" y1=\"%d\" x2=\"%d\" y2=\"%d\" " +
                "stroke-width=\"%d\"/>",
                x1, y1, x2, y2, stroke
            ));
        }
    }

    /**
     *  Draw a String at the provided X,Y location.
     *
     *  @param str The String to paint.
     *  @param x The X location to paint the String.
     *  @param y The Y location to paint the String.
     */
    public void drawString (String str, int x, int y) {
        writeln (String.format (
            "<text x=\"%d\" y=\"%d\">%s</text>", x, y, str
        ));
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
     *  Return the width of the font.
     *
     *  @return The font width.
     */
    public int fontWidth () {
        return fontWidth;
    }

    /**
     *  Return the width of the string as drawn.
     *
     *  @return The width of the string as drawn.
     */
    public int stringWidth (String str) {
        return (str.length () + 1 ) * fontWidth;
    }

    /**
     *  Private method to write a string of text with a line separator.
     *
     *  @param line The line of text.
     */
    private void writeln (String line) {
        try {
            writer.write (line + System.getProperty ("line.separator"));
        }
        catch (IOException e) {
            System.out.println ("Error writing to SVG file: " + line);
            e.printStackTrace ();
        }
    }

    private FileWriter writer;

    private int fontHeight = 12;
    private int fontWidth = 7;
    private String fontColor = "blue";

    private int strokeWidth = 1;
    private String strokeColor = "black";
}
