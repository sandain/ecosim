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
 *  Finds the lower and upper bounds for the omega confidence interval
 *
 *  @author Andrew Warner
 */
public class OmegaConfidence extends ConfidenceInterval {

    /**
     *  Constructor for OmegaConfidence.
     *
     *  @param value The inputted value from hillclimbing that includes theoptimized omega, sigma, and npop.
     *  @param results Results of this confidence interval.
     *  @param masterVariables The MasterVariables.
     *  @param input File containing the values fred method to examine.
     *  @param output Narrative to write to.
     */
    public OmegaConfidence(FredOutVal value, BinningAndFred results,
        MasterVariables masterVariables, File input, File output) {
        super("omega", value, results, masterVariables, input, output);
        this.incPerOM = 6;
        this.upperBound = 100.0;
        this.incUpper = (int)(incPerOM * Math.abs(Math.log10(value.getOmega() / upperBound)));
        this.incLower = (int)(incPerOM * Math.abs(Math.log10(value.getOmega() / lowerBound)));
    }

    /**
     *  Finds the upper bound of the confidence interval for omega.
     *
     *  @return FredOutVal containing the upper bound values of the confidence interval of omega.
     */
    public FredOutVal upperBound() {
        writeInputUpperBound();
        FredOutVal result = runCI();
        return result;
    }

    /**
     *  Finds the lower bound of the confidence interval for omega.
     *
     *  @return FredOutVal containing the lower bound values of the confidence interval of omega.
     */
    public FredOutVal lowerBound() {
        writeInputLowerBound();
        FredOutVal result = runCI();
        return result;
    }

    /**
     *  Writes in the input file for the upper bound of the omega confidence interval.
     */
    private void writeInputUpperBound() {
        // Declare default values.
        double[] omegaRange = {value.getOmega(), upperBound};
        // Omega range declared as a double - see ConfidenceInterval.writeInput documentation.
        int[] xnumics = {incUpper,0,0,0};
        int sortPer = masterVariables.getSortPercentage();
        double [] percentages = value.getPercentages();
        double probThresh = percentages[sortPer] / masterVariables.CI_NUMBER;
        super.writeInput(input, omegaRange, xnumics, probThresh);
        narr.println();
        narr.println("The input for upper bound omega confidence interval: ");
        narr.writeInput(input);
    }

    /**
     * Writes the input file, omegaIn.dat for the lower bound confidence interval.
     */
    private void writeInputLowerBound() {
        // Declare default values.
        double[] omegaRange = {value.getOmega(), lowerBound};
        // Omega range declared as a double - see ConfidenceInterval.writeInput documentation.
        int[] xnumics = {incLower, 0, 0, 0};
        int sortPer = masterVariables.getSortPercentage();
        double [] percentages = value.getPercentages();
        double probThresh = percentages[sortPer] / masterVariables.CI_NUMBER;
        super.writeInput(input, omegaRange, xnumics, probThresh);
        narr.println();
        narr.println("The input for lower bound omega confidence interval: ");
        narr.writeInput(input);
    }
}
