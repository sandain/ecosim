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

/**
 *  This Class holds output values from the results of a run of Fred's method,
 *  including values for omega, sigma, npop, drift, 5x, 2x, 1.5x, 1.25x, 1.1x,
 *  and 1.05x.
 *
 *  @author Andrew Warner
 */
public class FredOutVal implements Comparable<FredOutVal> {

    /**
     *  Default constructor for a Fred's method output value, takes in values
     *  for omega, sigma, etc and sets the appropriate value for this object
     *  percentages is an array of values where percentages[0]=5x, [1]=2x,
     *  [2]=1.5x, [3]=1.25x, [4]=1.1x, [5]=1.05x ; this structure makes it
     *  easier for a user to set which values to sort by.
     * 
     *  pre: all the given values are of the appropriate form and percentages is an array of size 6.
     *
     *  @param omega The value of omega.
     *  @param sigma The value of sigma.
     *  @param npop The value of npop.
     *  @param drift The value of drift.
     *  @param percentages The percentages array.
     *  @param masterVariables The MasterVariables.
     */
    public FredOutVal(double omega, double sigma, int npop, double drift, double [] percentages, MasterVariables masterVariables) {
        this.omega = omega;
        this.sigma = sigma;
        this.npop = npop;
        this.drift = drift;
        this.percentages = percentages;
        this.masterVariables = masterVariables;
    }

    /**
     *  Compares the given FredOutVal to this one.
     *
     *  pre: other is another object of type FredOutVal.
     *  post: -1 is returned if the set percentages are less in this object than the other object, otherwise.
     *
     *  @param other The other FredOutVal to be compared.
     */
    public int compareTo(FredOutVal other) {
        double [] otherVals = other.getPercentages();
        for (int i = masterVariables.getSortPercentage(); i > -1; i --) {
            // If the percentage is less for this output value, return -1.
            if (percentages[i] < otherVals[i]) {
                return -1;
            }
            // If it is greater, return 1.
            if (percentages[i] > otherVals[i]) {
                return 1;
            }
        }
        // If they are equal for the starting percentage, compare the next
        // highest percentages, (ie if we started at 1.5x, compare 2x now).
        return 0;
    }

    /**
     *  Returns the array "percentages" from this object.
     *
     *  @return double array containing the percentages.
     */
    public double [] getPercentages() {
        return percentages.clone();
    }

    /**
     *  Sets the percentages field.
     *
     *  @param percentages The percentages.
     */
    public void setPercentages(double[] percentages) {
        this.percentages = percentages;
    }

    /**
     *  Returns the FredOutVal as a string, in the form omega, sigma, npop,
     *  drift, 5x, 2x, 1.5x, 1.25x, 1.1x, 1.05x
     *
     *  @return String containing the values of omega, sigma, npop, drift, and the 6 percentages.
     */
    public String toString() {
        String s = String.format ("%.3g", omega);
        s += String.format (",   %.3g", sigma);
        s += String.format (",   %1d", npop);
        s += String.format (",   %.3g", drift);
        s += String.format (",   %.2g", percentages[0]);
        s += String.format (",   %.2g", percentages[1]);
        s += String.format (",   %.2g", percentages[2]);
        s += String.format (",   %.2g", percentages[3]);
        s += String.format (",   %.2g", percentages[4]);
        s += String.format (",   %.2g", percentages[5]);
        return s;
    }

    /**
     *  Set drift.
     *
     *  @param drift The value of drift.
     */  
 
    public void setDrift(double drift) {
        this.drift = drift;
    }

    /**
     *  Get omega.
     *
     *  @return double containing the value of omega.
     */  
    public double getOmega() {
        return omega;
    }

    /**
     *  Get sigma.
     *
     *  @return double containing the value of sigma.
     */  
    public double getSigma() {
        return sigma;
    }

    /**
     *  Get npop.
     *
     *  @return int containing the value of npop.
     */  
    public int getNpop() {
        return npop;
    }

    /**
     *  Get drift.
     *
     *  @return double containing the value of drift.
     */  
    public double getDrift() {
        return drift;
    }

    private double omega;
    private double sigma;
    private int npop;
    private double drift;
    private double[] percentages;
    private MasterVariables masterVariables;
}
