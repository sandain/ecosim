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


package ecosim.api;

public interface Painter {

    /**
     *  Start a Painter of the requested size.
     *
     *  @param width The width of the Painter.
     *  @param height The height of the Painter.
     */
    public void start (int width, int height);

    /**
     *  End a Painter.
     */
    public void end ();

    /**
     *  Draw a line between two provided points.
     *
     *  @param x1 The start X location of the line to draw.
     *  @param y1 The start Y location of the line to draw.
     *  @param x2 The end X location of the line to draw.
     *  @param y2 The end Y location of the line to draw.
     *  @param stroke The stroke width of the line.
     */
    public void drawLine (int x1, int y1, int x2, int y2, int stroke);

    /**
     *  Draw a String at the provided X,Y location.
     *
     *  @param str The String to paint.
     *  @param x The X location to paint the String.
     *  @param y The Y locaiton to paint the String.
     */
    public void drawString (String str, int x, int y);

    /**
     *  Return the width of the font.
     *
     *  @return The font width.
     */
    public int fontWidth ();

    /**
     *  Return the height of the font.
     *
     *  @return The font height.
     */
    public int fontHeight ();

    /**
     *  Return the width of the string as drawn by the Painter.
     *
     *  @return The width of the string as drawn by the Painter.
     */
    public int stringWidth (String str);

}
