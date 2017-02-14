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

import java.io.File;
import javax.swing.JTextArea;

/**
 *  Finds the lower and upper bounds for the npop confidence interval.
 *
 *  @author Andrew Warner
 */
public class NpopConfidence extends ConfidenceInterval {

    /**
     *  Constructor for NpopConfidence.
     *
     *  @param value The inputted value from hillclimbing that includes theoptimized omega, sigma, and npop.
     *  @param results Results of this confidence interval.
     *  @param masterVariables The master variables.
     *  @param input File containing the values fred method to examine.
     *  @param output Narrative to write to.
     */
    public NpopConfidence(FredOutVal value, BinningAndFred results, MasterVariables masterVariables, File input, File output) {
        super("npop", value, results, masterVariables, input, output);
    }

    /**
     *  Finds the upper bound of the confidence interval for npop.
     *
     *  @return FredOutVal containing the upper bound value of the confidence interval of npop.
     */
    public FredOutVal upperBound() {
        writeInputUpperBound();
        FredOutVal result = runCI();
        return result;
    }

    /**
     *  Finds the lower bound of the confidence interval for npop.
     *
     *  @return FredOutVal containing the lower bound of the confidenc interval of npop.
     */
    public FredOutVal lowerBound() {
        writeInputLowerBound();
        FredOutVal result = runCI();
        return result;
    }

    /**
     *  Writes in the input file for the upper bound of the Npop confidence interval.
     */
    private void writeInputUpperBound() {
        // Declare default values.
        double[] npopRange = {value.getNpop(), sequenceVals[0]};
        // Npop range declared as a double - see ConfidenceInterval.writeInput documentation.
        int[] xnumics = {0, 0, 1, 0};
        int sortPer = masterVariables.getSortPercentage();
        double [] percentages = value.getPercentages();
        double probThresh = percentages[sortPer] / masterVariables.CI_NUMBER;
        super.writeInput(input, npopRange, xnumics, probThresh);
        narr.println();
        narr.println("The input for upper bound npop confidence interval: ");
        narr.writeInput(input);
    }

    /**
     * Writes the input file, npopIn.dat for the lower bound confidence interval.
     */
    private void writeInputLowerBound() {
        // Declare default values.
        double[] npopRange = {value.getNpop(), 1};
        // Npop range declared as a double - see ConfidenceInterval.writeInput documentation.
        int[] xnumics = {0, 0, -1, 0};
        int sortPer = masterVariables.getSortPercentage();
        double [] percentages = value.getPercentages();
        double probThresh = percentages[sortPer] / masterVariables.CI_NUMBER;
        super.writeInput(input, npopRange, xnumics, probThresh);
        narr.println();
        narr.println("The input for lower bound npop confidence interval: ");
        narr.writeInput(input);
    }
}
