/*
 *    Ecotype Simulation models the sequence diversity within a bacterial
 *    clade as the evolutionary result of net ecotype formation and periodic
 *    selection, yielding a certain number of ecotypes.
 *
 *    Copyright (C) 2015  Jason M. Wood, Montana State University
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

import java.util.List;
import java.util.ArrayList;
import java.util.Collections;

/**
 *  Object to estimate the parameter values from the binning results.
 *
 *  @author Jason M. Wood
 *  @copyright GNU General Public License
 */
public class ParameterEstimate implements Runnable {

    /**
     *  Run the parameter estimate program.
     *
     *  @param masterVariables The MasterVariables object.
     *  @param nu The number of environmental sequences.
     *  @param length The length of the sequences being analyzed.
     *  @param binning The Binning object.
     */
    public ParameterEstimate (MasterVariables masterVariables,
        Execs execs, Integer nu, Integer length, Binning binning) {
        this.masterVariables = masterVariables;
        this.execs = execs;
        this.nu = nu;
        this.length = length;
        this.binning = binning;
        estimate = new ParameterSet ();
        sigma = new Line (new ArrayList<Point> ());
        omega = new Line (new ArrayList<Point> ());
        // The threshold needs to be modified by the number of environmental
        // sequences.
        threshold *= Math.log (nu) / logTwo;
        hasRun = false;
    }

    /**
     *  Run the parameter estimate program.
     */
    public void run () {
        // Calculate the list of points that the sigma and omega lines will
        // be fitted to.
        List<Point> points = getPoints (length, binning);
        // Fit the sigma line to the points.
        Integer [] sigmaBounds = fitLinePoints (points, 1);
        sigma = new Line (points.subList (sigmaBounds[0], sigmaBounds[1]));
        // Fit the omega line to the points.
        Integer [] omegaBounds = fitLinePoints (points, sigmaBounds[1] + 1);
        omega = new Line (points.subList (omegaBounds[0], omegaBounds[1]));
        // Omega is estimated from the slope of the omega line.
        Double omegaEstimate = -1.0d * omega.m;
        // Sigma is estimated from the slope of the sigma line.
        Double sigmaEstimate = -1.0d * sigma.m;
        // Npop is estimated by calculating the intersection of the omega
        // and sigma lines.
        Long npopEstimate = Math.round (Math.pow (
            2, omega.m * (sigma.b - omega.b) / (omega.m - sigma.m) + omega.b
        ));
        // Store the estimated parameter values.
        estimate = new ParameterSet (
            npopEstimate, omegaEstimate, sigmaEstimate, 0.0d
        );
        // Setup the fredmethod.
        FredMethod fredmethod = new FredMethod (
            masterVariables, execs, nu, length, binning, estimate
        );
        fredmethod.run ();
        estimate.setLikelihood (fredmethod.getResult ().getLikelihood ());
        // Set the flag stating that the parameter estimate program has
        // been run.
        hasRun = true;
    }

    /**
     *  Returns true if parameter estimate has been run, false otherwise.
     *
     *  @return True if parameter estimate has been run, false otherwise.
     */
    public boolean hasRun () {
        return hasRun;
    }

    /**
     *  Return the parameter estimate result.
     *
     *  @return The parameter estimate result.
     */
    public ParameterSet getResult () {
        return estimate;
    }

    /**
     *  Changes the value of hasRun.
     *
     *  @param hasRun The new value of hasRun.
     */
    public void setHasRun (boolean hasRun) {
        this.hasRun = hasRun;
    }

    /**
     *  Set the parameter estimate.
     */
    public void setResult (ParameterSet result) {
        estimate = result;
    }

    /**
     *  Return the parameter estimate as a String.
     *
     *  @return The parameter estimate as a String.
     */
    public String toString () {
        return estimate.toString ();
    }

    /**
     *  Return the slope and y-intercept calculated for omega.
     *
     *  @return An array of form [slope, intercept].
     */
    public double[] getOmega () {
        return new double[] { omega.m, omega.b };
    }

    /**
     *  Return the slope and y-intercept calculated for sigma.
     *
     *  @return An array of form [slope, intercept].
     */
    public double[] getSigma () {
        return new double[] { sigma.m, sigma.b };
    }

    /**
     *  Find all of the points that can fit a line.
     *
     *  @param points All of the points.
     *  @param start The index of points to start the line calculation.
     *  @return The bounds of the points array that fit a line.
     */
    private Integer [] fitLinePoints (List<Point> points, Integer start) {
        Integer [] bounds = { start, start + 2 };
        // Catch errors before they happen.
        if (bounds[0] > points.size () || bounds[1] > points.size ()) {
            throw new ArrayIndexOutOfBoundsException (
                "Bounds exceed the size of the points list while fitting the line."
            );
        }
        // Calculate the line using the current set of points.
        Line line = new Line (
            points.subList (bounds[0], bounds[1])
        );
        // Slurp up any points that are close to the line.
        for (int i = bounds[1]; i < points.size (); i ++) {
            Double error = squaredError (points.get (i), line);
            if (error > threshold) break;
            bounds[1] = i;
        }
        return bounds;
    }

    /**
     *  Calculate the squared error of the point compared to the line.
     *
     *  @param point The point.
     *  @param line The line.
     *  @return The squared error of the point compared to the line.
     */
    private Double squaredError (Point point, Line line) {
        // Calculate the distance of the point from the line as the error.
        Double error = Math.abs (-1 * line.m * point.x + point.y - line.b) /
            Math.sqrt (line.m * line.m + 1);
        // Square the error.
        return error * error;
    }

    /**
     *  Calculate the list of points from the binning results.
     *
     *  @param length The length of the sequences.
     *  @param binning The binning results to transform.
     *  @return The transformed binning results as a list of XY points.
     */
    private List<Point> getPoints (Integer length, Binning binning) {
        List<Point> points = new ArrayList<Point> ();
        Integer previous = -1;
        for (BinLevel bin: binning.getBins ()) {
            Double crit = bin.getCrit ();
            Integer level = bin.getLevel ();
            // Don't include binning results == 1.
            if (level == 1) continue;
            // Don't include duplicate levels.
            if (level == previous) continue;
            // Transform the sequence criterion value into the number of SNPs.
            Double x = (1.0d - crit) * length;
            // Transform the number of sequence clusters (bins) into log base 2 scale.
            Double y = Math.log (level) / Math.log (2);
            // Add the XY point to the list.
            points.add (new Point (x, y));
            // Save the level in previous.
            previous = level;
        }
        // Return the points in reverse order.
        Collections.reverse (points);
        return points;
    }

    private MasterVariables masterVariables;
    private Execs execs;
    private Integer nu;
    private Integer length;
    private Binning binning;

    private Line sigma;
    private Line omega;

    private Double threshold = 0.1d;

    private final Double logTwo = Math.log (2);

    /**
     * The estimate for the parameter values.
     */
    private ParameterSet estimate;

    private boolean hasRun;

    /**
     *  A private class to calculate the best fit line from from a selection
     *  of points.
     */
    private class Line {
        public Line (List<Point> points) {
            Double sumX = 0.0d;
            Double sumY = 0.0d;
            Double sumXY = 0.0d;
            Double sumX2 = 0.0d;
            for (Point point: points) {
                sumX += point.x;
                sumY += point.y;
                sumXY += point.x * point.y;
                sumX2 += point.x * point.x;
            }
            // Calculate the slope and Y intercept of the line.
            int n = points.size ();
            b = (sumY * sumX2 - sumX * sumXY) / (n * sumX2 - sumX * sumX);
            m = (n * sumXY - sumX * sumY) / (n * sumX2 - sumX * sumX);
        }
        public Double b;
        public Double m;
    }

    /**
     *  A private class to store a simple XY point value.
     */
    private class Point {
        public Point (Double x, Double y) {
            this.x = x;
            this.y = y;
        }
        public Double x;
        public Double y;
    }

}
