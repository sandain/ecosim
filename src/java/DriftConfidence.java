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
 *  The drift confidence interval.
 */
public class DriftConfidence extends ConfidenceInterval {

    /**
     *  Constructor for objects of class DriftConfidence.
     *
     *  @param value The inputted value from hillclimbing that includes theoptimized omega, sigma, and npop.
     *  @param results Results of this confidence interval.
     *  @param masterVariables The MasterVariables.
     *  @param input File containing the values fred method to examine.
     *  @param output Narrative to write to.
     */
    public DriftConfidence(FredOutVal value, BinningAndFred results, MasterVariables masterVariables, File input, File output) {
        super("drift", value, results, masterVariables, input, output);
        this.incPerOM = 3;
        this.incLower = (int)(incPerOM * Math.abs(Math.log10(value.getSigma() / lowerBound)));
    }

    /**
     *  Finds the lower bound of the confidence interval for drift.
     *
     *  @return FredOutVal containing the lower bound of the confidence interval for drift.
     */
    public FredOutVal lowerBound() {
        writeInputLowerBound();
        FredOutVal result = runCI();
        return result;
    }

    /**
     *  Writes the input file, driftIn.dat for the lower bound confidence interval.
     */
    private void writeInputLowerBound() {
        // Declare default values.
        double[] driftRange = {masterVariables.DRIFT_CI_START, lowerBound};
        // Drift range declared as a double - see ConfidenceInterval.writeInput documentation.
        int[] xnumics = {0, 0, 0, incLower};
        int sortPer = masterVariables.getSortPercentage();
        double [] percentages = value.getPercentages();
        double probThresh = percentages[sortPer] / masterVariables.ONETAIL_CI_NUMBER;
        super.writeInput(input, driftRange, xnumics, probThresh);
        narr.println();
        narr.println("The input for lower bound drift confidence interval: ");
        narr.writeInput(input);
    }
}
