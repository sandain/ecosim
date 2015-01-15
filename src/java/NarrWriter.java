/*
 *    Ecotype Simulation models the sequence diversity within a bacterial clade as
 *    the evolutionary result of net ecotype formation, periodic selection,
 *    and drift, yielding a certain number of ecotypes.
 * 
 *    Copyright (C) 2009  Fred Cohan, Wesleyan University
 *                        Carlo Francisco, Wesleyan University
 *                        Danny Krizanc, Wesleyan University
 *                        Andrew Warner, Wesleyan University
 *                        Jason Wood, Montana State University
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

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;

/**
 *  Writes a bare-bones narrative file for review.
 *
 *  @author Andrew Warner
 */
public class NarrWriter {

    /**
     * Creates a NarrWriter object.
     *
     * @param output The file to write the narrative to.
     */
    public NarrWriter(File output) {
        try {
            this.output = new BufferedWriter(new FileWriter(output));
        }
        catch (IOException e) {
            e.printStackTrace();
        }
    }

    /**
     *  Writes a string to the narrative file.
     *
     *  @param s The String to be written.
     */
    public void print(String s) {
        try {
            output.write(s);
            allText.append(s);
        }
        catch (IOException e) {
            e.printStackTrace();
        }
    }

    /**
     *  Writes a line to the narrative file.
     *
     *  @param s The String to be written.
     */
    public void println(String s) {
        try {
            output.write(s);
            output.newLine();
            allText.append(s + "\n");
        }
        catch(IOException e) {
            e.printStackTrace();
        }
    }

    /**
     *  Overloaded the println method to print a blank line.
     */
    public void println() {
        try {
            output.newLine();
            allText.append("\n");
        }
        catch (IOException e) {
            e.printStackTrace();
        }
    }

    /**
     *  Closes the narrator file.
     */
    public void close() {
        try {
            output.close();
        }
        catch(IOException e) {
            e.printStackTrace();
        }
    }

    /**
     *  Takes a fred method / hillclimbing input or output file and writes it to
     *  the narrative file via the narrator.
     *
     *  @param inputFile The file to write to the narrator.
     */
    public void writeInput(File inputFile) {
        BufferedReader input;
        try {
            input = new BufferedReader(new FileReader(inputFile));
            String nextLine = input.readLine();
            while (nextLine != null) {
                println(nextLine);
                nextLine = input.readLine();
            }
            input.close();
        }
        catch (IOException e) {
            e.printStackTrace();
        }
    }

    /**
     *  Return all text stored in this Narrator.
     *
     *  @return String containing the text.
     */
    public String getText() {
       return allText.toString();
    }

    /**
     *  The output writer.
     */
    private BufferedWriter output;

    /**
     *  The text from the narrator.
     */
    private StringBuffer allText = new StringBuffer();
}
