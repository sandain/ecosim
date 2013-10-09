/*
 *    Ecotype Simulation models the sequence diversity within a bacterial clade
 *    as the evolutionary result of net ecotype formation and periodic
 *    selection, yielding a certain number of ecotypes.
 *
 *    Copyright (C) 2013  Jason M. Wood, Montana State University
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

/**
 *  Stores the likelihood values for the Bruteforce object.
 *
 *  @author Jason M. Wood
 *  @copyright GNU General Public License
 */
public class Likelihood implements Comparable<Likelihood> {

    /**
     *  Constructor for an array of likelihood values.
     *
     *  @param masterVariables The MasterVariables object.
     *  @param likelihood An array of likelihood values.
     */
    public Likelihood(MasterVariables masterVariables, Double [] likelihood) {
        this.masterVariables = masterVariables;
        this.likelihood = likelihood;
    }

    /**
     *  Constructor for likelihood values stored as a coma separated list.
     *
     *  @param masterVariables The MasterVariables object.
     *  @param str A comma separated String of likelihood values.
     */
    public Likelihood(MasterVariables masterVariables, String str) {
        this.masterVariables = masterVariables;
        likelihood = new Double[6];
        String[] strings = str.split(",\\s*");
        for (int i = 0; i < 6; i ++) {
            likelihood[i] = new Double(strings[i]).doubleValue();
        }
    }

    /**
     *  Returns the value stored at the given index.
     *
     *  @param index The index of the value wanted.
     *  @return The value at the given index.
     */
    public Double get(int index) {
        return likelihood[index];
    }

    /**
     *  Compares the likelihood values stored in this object with those
     *  stored in another Likelihood object.
     *
     *  @param other The other set of likelihood values to compare to.
     *  @return -1 if less than, 0 if equal, and 1 if more than.
     */
    public int compareTo(Likelihood other) {
        for (int i = masterVariables.getCriterion() - 1; i > -1; i --) {
            // If the percentage is less for this output value, return -1.
            if (likelihood[i] < other.get(i)) {
                return -1;
            }
            // If it is greater, return 1.
            if (likelihood[i] > other.get(i)) {
                return 1;
            }
        }
        // If they are equal for the starting percentage, compare the next
        // highest likelihood, (ie if we started at 1.5x, compare 2x now).
        return 0;
    }

    /**
     *  Returns the likelihood values as a comma separated list.
     *
     *  @return The likelihood values as a comma separated list in String from.
     */
    public String toString() {
        String str = "";
        for (int i = 0; i < 6; i ++) {
            str += String.format("%.5g", likelihood[i]);
            if (i < 5) {
                str += ", ";
            }
        }
        return str;
    }

    private MasterVariables masterVariables;
    private Double [] likelihood;
}

