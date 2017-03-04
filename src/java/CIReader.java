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
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.StringTokenizer;

/**
 *  A confidence Interval reader, reads through the output from the confidence
 *  interval program and returns the values of omega, sigma, npop, drift, and
 *  percentage likelihood
 *
 *  @author Andrew Warner
 */
public class CIReader extends Thread {

    /**
     *  Constructor for objects of class CIReader.
     *
     *  @param is The inputstream from the confidence interval program running.
     *  @param p The process currently running (whatever confidence interval program is running.
     *  @param lastVal The output value from hillclimbing.
     *  @param type The type of confidence interval we are running, either "omega", "sigma", "npop", or "drift".
     *  @param isNew Whether to use the new CI or the old CI program.
     *  @param masterVariables The MasterVariables.
     */
    public CIReader(InputStream is, Process p, FredOutVal lastVal, String type, boolean isNew, MasterVariables masterVariables) {
        this.is = is;
        this.p = p;
        this.lastVal = lastVal;
        this.type = type;
        this.isNew = isNew;
        this.masterVariables = masterVariables;
        double [] percentages = lastVal.getPercentages();
        likelihood = percentages[masterVariables.getSortPercentage()];
        if (type.equals("npop") || type.equals("sigma") || type.equals("omega")) {
            ciNumber = masterVariables.CI_NUMBER;
        }
        else {
            ciNumber = masterVariables.ONETAIL_CI_NUMBER;
        }
        hadBetterValue = false;
    }

    /**
     *  Run this CIReader.
     */
    public void run() {
        double omega, sigma, drift, yvalue;
        int npop;
        double[] percentages;
        int sortPer = masterVariables.getSortPercentage();
        try {
            BufferedReader input = new BufferedReader(new InputStreamReader(is));
            String line = input.readLine();
            StringTokenizer tk;
            while (line != null) {
                tk = new StringTokenizer(line);
                while (! tk.nextToken().equals("starting")) {
                    System.out.println(line);
                    line = input.readLine();
                    tk = new StringTokenizer(line);
                }
                System.out.println(line);
                // Read omega.
                line = input.readLine();
                System.out.println(line);
                tk = new StringTokenizer(line);
                tk.nextToken();
                if (! isNew) {
                    tk.nextToken();
                    tk.nextToken();
                }
                // If we are reading the output from the old confidence interval
                // program, omega is the 4th token of the line, rather than the second.
                omega = (new Double(tk.nextToken())).doubleValue();
                // Read sigma.
                line = input.readLine();
                System.out.println(line);
                tk = new StringTokenizer(line);
                tk.nextToken();
                if (! isNew) {
                    tk.nextToken();
                    tk.nextToken();
                }
                sigma = (new Double(tk.nextToken())).doubleValue();
                // Read npop.
                line = input.readLine();
                System.out.println(line);
                tk = new StringTokenizer(line);
                tk.nextToken();
                if (! isNew) {
                    tk.nextToken();
                    tk.nextToken();
                }
                npop = (new Integer(tk.nextToken())).intValue();
                // Read drift.
                line = input.readLine();
                System.out.println(line);
                tk = new StringTokenizer(line);
                tk.nextToken();
                if (! isNew) {
                    tk.nextToken();
                    tk.nextToken();
                }
                drift = (new Double(tk.nextToken())).doubleValue();
                // Read the seed line.
                line = input.readLine();
                System.out.println(line);
                // Read the likelihood.
                line = input.readLine();
                System.out.println(line);
                tk = new StringTokenizer(line);
                tk.nextToken();
                yvalue = (new Double(tk.nextToken())).doubleValue() * (-1.0);
                percentages = new double[6];
                percentages[sortPer] = yvalue;
                // Find out if a hillclimbing has been completed, if so
                // compare the likelihood from the last value of the confidence
                // interval program to the value we got in hillclimbing.
                // If it is less than the likelihood we got from hillclimbing
                // divided by 6.83 we are done.
                // If the likelihood is better than the value we got from
                // hillclimbing by 20% we need to later restart hillclimbing
                // using that value as a seed.
                if (nextVal != null) {
                    double thisRunVal, compValue;
                    if (type.equals("omega")) {
                        thisRunVal = omega;
                        compValue = nextVal.getOmega();
                    }
                    else if (type.equals("sigma")) {
                        thisRunVal = sigma;
                        compValue = nextVal.getSigma();
                    }
                    else if (type.equals("npop")) {
                        thisRunVal = npop;
                        compValue = nextVal.getNpop();
                    }
                    else {
                        thisRunVal = drift;
                        compValue = nextVal.getDrift();
                    }
                    // If we have finished the confidence interval and the
                    // end value is out of confidence or more than 20% better
                    // than our likelihood from hillclimbing, we are done.
                    if (thisRunVal != compValue) {
                        double [] tempPercentages = nextVal.getPercentages();
                        if (tempPercentages[sortPer] < likelihood / ciNumber) {
                             p.destroy();
                             break;
                        }
                        // If we are running one of the new confidence interval
                        // programs we would like to check if the value from this
                        // hillclimbing is better than the one we came into
                        // the confidence interval with.
                        if (isNew) {
                            if (tempPercentages[sortPer] > (likelihood * 1.2)) {
                                hadBetterValue = true;
                                p.destroy();
                                break;
                            }
                        }
                        lastVal = nextVal;
                    }
                }
                // Set the last value that we saw and attempt to read a new line of output.
                nextVal = new FredOutVal(omega, sigma, npop, drift, percentages, masterVariables);
                line = input.readLine();
                System.out.println(line);
                line = input.readLine();
                if (line != null) {
                    tk = new StringTokenizer(line);
                }
                while (line != null && ! tk.nextToken().equals("starting")) {
                    System.out.println(line);
                    line = input.readLine();
                    if (line != null) {
                        tk = new StringTokenizer(line);
                    }
                }
                System.out.println("---------- END ONE REP ---------");
            }
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }

    /**
     *  Returns the last value still in confidence.
     *
     *  @return The last FredOutVal still in confidence.
     */
    public FredOutVal lastValue() {
        return lastVal;
    }

    /**
     *  returns the most recent value (usually out of confidence but there
     *  will be some cases when it is needed because the likelihood is better
     *  than the likelihood we got from hillclimbing).
     *
     *  @return The most recent FredOutVal.
     */
    public FredOutVal nextValue() {
        return nextVal;
    }

    /**
     *  Returns whether or not the given confidence interval got a value
     *  with a better likelihood for hillclimbing.
     *
     *  @return True if the given confidence interval has a better likelihood, false otherwise.
     */
    public boolean hadBetterValue() {
        return hadBetterValue;
    }

    private InputStream is;
    private double likelihood;
    private Process p;
    private double ciNumber;
    private FredOutVal lastVal;
    private FredOutVal nextVal;
    private String type;
    private boolean isNew;
    private boolean hadBetterValue;
    private MasterVariables masterVariables;
}
