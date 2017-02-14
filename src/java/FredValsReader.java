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
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.StringTokenizer;

/**
 *  Reads the values from the output of Fred's method, usually contained in the file acinas.out
 *
 *  @author Andrew Warner
 */
public class FredValsReader {

    /**
     *  Constructor for FredValsReader.
     *
     *  @param masterVariables The MasterVariables.
     */
    public FredValsReader(MasterVariables masterVariables) {
        this.masterVariables = masterVariables;
    }

    /**
     *  Read the values of the Fred's method output file (usually output.dat) into an ArrayList.
     *
     *  pre: filename is the name of a valid input file.
     *  post: the returned arraylist is an array of FredOutVals from the input file.
     *
     *  @param filename The file to read the FredOutVals from.
     *  @return The values of Fred's method output.
     */
    public ArrayList<FredOutVal> readFredOutVals(File filename) {
        // The list of returned FredOutVals.
        ArrayList<FredOutVal> returnVals = new ArrayList<FredOutVal>();
        // A reader for the Fred's method output file
        try {
            BufferedReader input = new BufferedReader(new FileReader(filename));
            StringTokenizer reader;
            /* Rather than creating seperate variables for each value of the output line
             * it seems simpler just to create an array of double since npop will be
             * cast into an integer in the FredOutVal constructor.
             */
            double[] vals = new double[4];
            double[] percentages = new double[6];
            // A necessary object needed to cast strings to doubles.
            Double convertedVal;
            String nextLine = input.readLine();
            String token;
            // Create a FredOutVal corresponding to each line of output.dat and add it to the arraylist.
            while (nextLine != null) {
                reader = new StringTokenizer(nextLine);
                // Read the appropriate tokens and place them in the array.
                for (int i = 0; i < vals.length; i++) {
                    if (! reader.hasMoreTokens()) {
                        throw new IndexOutOfBoundsException("Invalid input file for Fred sorting");
                    }
                    token = reader.nextToken();
                    // Delete the comma at the end.
                    token = token.substring(0, token.length() - 1);
                    convertedVal = new Double(token);
                    vals[i] = convertedVal.doubleValue();
                }
                // Read the percentage values and place them in a seperate array.
                for (int i = 0; i < percentages.length; i++) {
                    if (! reader.hasMoreTokens()) {
                        throw new IndexOutOfBoundsException("Invalid input file for Fred sorting");
                    }
                    token = reader.nextToken();
                    // Delete the comma at the end, the last value does not have a comma so we must check that this token is not the last value.
                    if (i != percentages.length - 1) {
                        token = token.substring(0, token.length() - 1);
                    }
                    convertedVal = new Double(token);
                    percentages[i] = convertedVal.doubleValue();
                }
                // Place the values from the "vals" array into a FredOutVal object.
                FredOutVal nextVal = new FredOutVal(vals[0], vals[1], (int)vals[2], vals[3], percentages.clone(), masterVariables);
                returnVals.add(nextVal);
                nextLine = input.readLine();
            }
            input.close();
        }
        catch (IOException e) {
            System.out.println("Error reading file.");
            e.printStackTrace();
            return null;
        }
        return returnVals;
    }

    /**
     *  Reads the values of omega, sigma, and npop from the hillclimbing output file,
     *  exponentiates them, and outputs them as a FredOutVal with likelihoods as 0.
     *
     *  pre: hClimbOut is a valid hill climbing output file.
     *  post: FredOutVal contains the final data points from hill climbing.
     *
     *  @param hClimbOut A valid hill climbing output file.
     *  @return FredOutVal containing the final omega, sigma, and npop for hillclimbing.
     */
    public FredOutVal readHillClimbOutput(File hClimbOut) {
        FredOutVal returnVal = null;
        try {
            BufferedReader input = new BufferedReader(new FileReader(hClimbOut));
            String nextLine = input.readLine();
            StringTokenizer tok;
            while (nextLine!=null) {
                tok = new StringTokenizer(nextLine);
                if (tok.hasMoreTokens()) {
                    if (tok.nextToken().equals("Minimum")) {
                        break;
                    }
                }
                nextLine = input.readLine();
            }
            if (nextLine == null) {
                System.out.println("Bad hillclimb output file");
                return null;
            }
            nextLine = input.readLine();
            nextLine = input.readLine();
            tok = new StringTokenizer(nextLine);
            double omega = (new Double(tok.nextToken())).doubleValue();
            double sigma = (new Double(tok.nextToken())).doubleValue();
            // Round npop to an int.
            int npop = (int)((new Double(tok.nextToken())).doubleValue() + 0.5);
            omega = Math.exp(omega);
            sigma = Math.exp(sigma);
            // Placeholder for percentages; all values 0
            double[] percentages = new double[6];
            returnVal = new FredOutVal(omega, sigma, npop, 1.0e25, percentages, masterVariables);
            input.close();
        }
        catch (IOException e) {
            e.printStackTrace();
        }
        return returnVal;
    }

    private MasterVariables masterVariables;
}
