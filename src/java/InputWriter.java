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

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.StringTokenizer;

/**
 *  Writes input files for Fred's method.
 * 
 *  @author Andrew Warner
 */
public class InputWriter {

    /**
     *  No constructor; used only from static context.
     */
    public InputWriter() {
    }
    
    /**
     *  Write the input file for a run of Fred's method, Hill climbing, or confidence intervals.
     *
     *  pre - All the input values are formatted properly.
     *  post - The input file (usually acinas.dat) for Fred method is written.
     *
     *  @param output The output file to write to.
     *  @param bins The top part of the input file, ie, the percentage differences
     *  such that the sequences fit into more than one bin.
     *  @param omegaRange Must be an array of size 2, containing the range of omega; IE if the
     *  range is .001 to 100.0, then omegaRange = {.001, 100.0}.
     *  @param sigmaRange Same as omegaRange.
     *  @param npopRange The range for npop.
     *  @param driftRange The drift range.
     *  @param xnumincs An array of size 4 where the values correspond to those of the xnumics line of input.
     *  @param nu The value of nu.
     *  @param nrep The number of reps.
     *  @param lenSequence The length of the gene sequences.
     *  @param jwhichxavg Used for everything except fred method, the program will not write it if jwhichxavg is -1.
     *  @return True if writing the output file was successful, false otherwise.
     */
    public static boolean writeFile(File output, ArrayList<String> bins, double[] omegaRange, double[] sigmaRange, int[] npopRange, 
        double[] driftRange, int[] xnumincs, int nu, int nrep, int lenSequence, int jwhichxavg) {
        return writeFile(output, bins, omegaRange, sigmaRange, npopRange, driftRange, xnumincs, nu, nrep, lenSequence, jwhichxavg, 1.0);
    }
    
    /**
     *  Write the input file for a run of Fred's method, Hill climbing, or confidence intervals.
     *
     *  pre - All the input values are formatted properly
     *  post - The input file (usually acinas.dat) for Fred method is written
     *
     *  @param output The output file to write to.
     *  @param bins The top part of the input file, ie, the percentage differences such that the sequences fit into more than one bin.
     *  @param omegaRange An array of size 2, containing the range of omega; IE if the range is .001 to 100.0, then omegaRange = {.001, 100.0}.
     *  @param sigmaRange The same as omegaRange.
     *  @param npopRange The range for npop.
     *  @param driftRange The drift range.
     *  @param xnumincs An array of size 4 where the values correspond to those of the xnumics line of input.
     *  @param nu The value of nu.
     *  @param nrep The number of reps.
     *  @param lenSequence The length of the gene sequences.
     *  @param jwhichxavg Used for everything except fred method, the program will not write it if jwhichxavg is -1.
     *  @param probThres The probability threshold at which to continue. If the value gets a better probability
     *  or to stop if the best hill climbing value at that value of the given parameter is outside of the confidence interval.
     *  @return True if writing the output file was successful, false otherwise
     */
    public static boolean writeFile(File output, ArrayList<String> bins, double[] omegaRange, double[] sigmaRange, int[] npopRange, 
        double[] driftRange, int[] xnumincs, int nu, int nrep, int lenSequence, int jwhichxavg, double probThres) {
        try {
            BufferedWriter out;
            out = new BufferedWriter(new FileWriter(output));
            out.write("" + bins.size() + "        numcrit"); // Write the first line of the output file.
            out.newLine();
            StringTokenizer a;
            for (int i = 0; i < bins.size(); i ++) {
                a = new StringTokenizer(bins.get(i));
                a.nextToken();
                out.write(a.nextToken());
                out.newLine();
            } // Write the bin data; has to be written twice, first without percentages.
            for (int j = 0; j < bins.size(); j ++) {
                out.write(bins.get(j));
                out.newLine();
            }
            // Now write the range of omega.
            out.write("" + omegaRange[0] + "," + omegaRange[1] + "      omega");
            out.newLine();
            // Range of sigma.
            out.write("" + sigmaRange[0] + "," + sigmaRange[1] + "       sigma");
            out.newLine();
            // Npop.
            out.write("" + npopRange[0] + "," + npopRange[1] + "         npop");
            out.newLine();
            // Drift.
            out.write("" + driftRange[0] + "," + driftRange[1] + "       drift");
            out.newLine();
            // Xnumics values.
            out.write("" + xnumincs[0] + "," + xnumincs[1] + "," + xnumincs[2] + "," + xnumincs[3] + "      xnumincs");
            out.newLine();
            // Write nu.
            out.write("" + nu + "       nu");
            out.newLine();
            // Write nrep.
            out.write("" + nrep + "        nrep (number of replicate runs per quartet  of parameter values)");
            out.newLine();
            // Create the random number seed; an odd less than 9 digits long.
            long randNumSeed = (long)(100000000 * Math.random());
            if (randNumSeed % 2 == 0) {
                randNumSeed ++;
            }
            // Write the random number seed.
            out.write("" + randNumSeed + "           iii (random number seed)");
            out.newLine();
            out.write("" + lenSequence + "      lengthseq (after deleting gaps, etc.)");
            if (jwhichxavg > 0) {
                out.newLine();
                out.write("" + jwhichxavg + "  jwhichxavg (which criterion)  (#4 = 1.25x criterion)");
            }
            out.newLine();
            out.write(probThres + "     probthresh");
            out.newLine();
            out.close();
        }
        catch (IOException e) {
            System.out.println("Invalid output file.");
            e.printStackTrace();
            return false;
        }
        return true;
    }
}
