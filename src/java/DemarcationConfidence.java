/*
 *    Ecotype Simulation models the sequence diversity within a bacterial clade as
 *    the evolutionary result of net ecotype formation, periodic selection,
 *    and drift, yielding a certain number of ecotypes.
 * 
 *    Copyright (C) 2009  Andrew Warner, Wesleyan University
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


package ecosim;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.StringTokenizer;
import javax.swing.JTextArea;

/**
 *  Runs the confidence interval needed in demarcations.
 *
 *  @author Andrew Warner
 */
public class DemarcationConfidence extends ConfidenceInterval {

    /**
     *  Constructor for objects of class DemarcationConfidence.
     *
     *  @param value The inputted value from hillclimbing that includes theoptimized omega, sigma, and npop.
     *  @param results Results of this confidence interval.
     *  @param masterVariables The MasterVariables.
     *  @param input File containing the values fred method to examine.
     *  @param output File to write to.
     */
    public DemarcationConfidence(FredOutVal value, BinningAndFred results, MasterVariables masterVariables, File input, File output) {
        super("demarcation", value, results, masterVariables, input, output);
    }

    /**
     *  Runs demarcations using the old npop confidence interval program.
     *
     *  @return Array of int such that [0] is the optimal npop, [1] is
     *  the lower bound of the confidence interval, and [2] is the upper bound.
     */
    public int[] demarcations() {
        writeInputDemarcations();
        runCI(); // Run the confidence interval, we don't care about the
                 // last value as we do in normal confidence interval
                 // runs, we will read the output seperately.
        int [] interval = getDemarcOutput();
        return interval;
    }

    /**
     *  Reads the output from demarcations from the output file.
     *
     *  @return Array of int such that [0] is the optimal npop, [1] is
     *  the lower bound of the confidence interval, and [2] is the upper bound.
     */
    protected int [] getDemarcOutput() {
        try {
            BufferedReader input = new BufferedReader(new FileReader(output));
            bestLike = 0.0;
            double nextLike;
            String line = input.readLine();
            StringTokenizer tk;
            int [] demarcValues = new int[3];
            int nextNpop;
            while (line != null) {
                tk = new StringTokenizer(line);
                nextNpop = (int)(new Double(tk.nextToken())).doubleValue();
                tk.nextToken();
                nextLike = -1 * (new Double(tk.nextToken())).doubleValue();
                if (nextLike > bestLike) {
                    bestLike = nextLike;
                    demarcValues[0] = nextNpop;
                }
                line = input.readLine();
            }
            input.close();
            // Read the same file again, only this time to find the confidence interval.
            double probThresh = bestLike / masterVariables.CI_NUMBER;
            input = new BufferedReader(new FileReader(output));
            line = input.readLine();
            while (line != null) {
                tk = new StringTokenizer(line);
                nextNpop = (int)(new Double(tk.nextToken())).doubleValue();
                tk.nextToken();
                nextLike = -1 * (new Double(tk.nextToken())).doubleValue();
                // If we have reached the border of the confidence interval, set the value appropriately.
                if (demarcValues[1] == 0) {
                    if (nextLike > probThresh) {
                        demarcValues[1] = nextNpop;
                    }
                }
                // Otherwise check if we are out of the confidence interval yet.
                else if (nextLike < probThresh) {
                    demarcValues[2] = nextNpop - 1;
                    break;
                }
                line = input.readLine();
                if (line == null && ! (nextLike < probThresh))
                    demarcValues[2] = nextNpop;
            }
            input.close();
            return demarcValues;
        }
        catch (IOException e) {
            e.printStackTrace();
        }
        return null;
    }

    /**
     *  @return The best likelihood associated with the best npop from the demarcation confidence interval
     */
    public double getBestLike() {
        return bestLike;
    }

    /**
     *  Writes the input file for a demarcations run.
     */
    private void writeInputDemarcations() {
        // Declare default values.
        if (value.getNpop() < sequenceVals[0]) {
            upperBound = value.getNpop();
        }
        else {
            upperBound = sequenceVals[0];
        }
        double[] npopRange = {1, upperBound};
        // Npop range declared as a double - see ConfidenceInterval.writeInput documentation.
        int[] xnumics = {0, 0, 1, 0};
        int sortPer = masterVariables.getSortPercentage();
        double [] percentages = value.getPercentages();
        double probThresh = percentages[sortPer] / masterVariables.CI_NUMBER;
        super.writeInput(input, npopRange, xnumics, probThresh);
        narr.println();
        narr.println("The input for the demarcations run: ");
        narr.writeInput(input);
    }

    /**
     * The best likelihood from the demarcation confidence interval.
     */
    private double bestLike;
}
