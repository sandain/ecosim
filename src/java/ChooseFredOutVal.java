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

/**
 *  Sorts the output values from Fred's method and chooses the best one.
 *
 *  @author Andrew Warner
 */
public class ChooseFredOutVal {

    /**
     *  This class is only used from a static context so it has no constructor
     */
    public ChooseFredOutVal() {
    }

    /**
     *  Takes in a list of FredOutVal objects, chooses the first column from
     *  the right with a non-zero likelihood (ie first it tries 1.05x, 1.1x,
     *  then 1.25x, etc, and then sorts the objects by that value, using the
     *  greater likelihoods as tiebreakers, and outputs the best Fred Value
     *  from that sort
     *
     *  pre - values is a list of appropriate fredoutvals and debug is a valid output filename
     *  post - the best FredOutVal is returned and the sorted list is written to debug
     *
     *  @param values The output from a run of Fred's method.
     *  @param debug An output file to write the sorted list of values to, in case the chosen value was
     *      not optimal and the user wants to go back and see what other values might be used.
     *  @param masterVariables The MasterVariable.
     *  @return FredOutVal containing the optimal choice for hill climbing.
     */
    public static FredOutVal chooseBest(ArrayList<FredOutVal> values, File debug, MasterVariables masterVariables) {
        try {
            BufferedWriter out = new BufferedWriter(new FileWriter(debug));
            int startPer; // The starting percentage to sort by (5 is 1.05x, etc).
            double [] storage; // Storage for the array of percentages.
            int i; // A counter.
            for (startPer = 5; startPer > 0; startPer --) {
                for (i = 0; i < values.size(); i ++){
                    storage = values.get(i).getPercentages();
                    if (storage[startPer] > 0) {
                        break;
                    }
                }
                if (i != values.size()) {
                    // If we found a non-zero percentage before traversing the list, then we have found the percentage to start sorting by.
                    break;
                }
            }
            masterVariables.setSortPercentage(startPer);
            Heapsorter<FredOutVal> h = new Heapsorter<FredOutVal>();
            h.heapSort(values);
            for (int j = 0; j < values.size(); j ++) {
                out.write(values.get(j).toString());
                out.newLine();
            }
            out.close();
        }
        catch (IOException e) {
            System.out.println("Error creating output file.");
            e.printStackTrace();
            return null;
        }
        return values.get(0);
    }
}
