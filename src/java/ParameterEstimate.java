/*
 *    Ecotype Simulation models the sequence diversity within a bacterial
 *    clade as the evolutionary result of net ecotype formation and periodic
 *    selection, yielding a certain number of ecotypes.
 *
 *    Copyright (C) 2016  Jason M. Wood, Montana State University
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
     *  The parameter estimate program.
     */
    public ParameterEstimate () {
        this (0, 0, new Binning ());
    }

    /**
     *  The parameter estimate program.
     *
     *  @param nu The number of environmental sequences.
     *  @param length The length of the sequences being analyzed.
     *  @param binning The Binning object.
     */
    public ParameterEstimate (Integer nu, Integer length, Binning binning) {
        this.nu = nu;
        this.length = length;
        this.binning = binning;
        result = new ParameterSet ();
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
        List<Point> points = getPoints (length);
        // Estimate the bounds of the sigma and omega lines.
        Integer[] bounds = getBounds (points);
        Integer[] sigmaBounds = {
            bounds[0],
            (int)Math.floor ((bounds[1] - bounds[0]) / 2 + bounds[0])
        };
        Integer[] omegaBounds = {
            (int)Math.ceil ((bounds[1] - bounds[0]) / 2 + bounds[0]),
            bounds[1]
        };
        // Optimize the sigma and omega lines using a custom variant of
        // Lloyd's Algorithm.
        Double error = 0.0D;
        Double previousError = 0.0D;
        Double deltaError = 1.0D;
        // Keep track of the bounds that produce the smallest error.
        Integer[] minSigmaBounds = sigmaBounds;
        Integer[] minOmegaBounds = omegaBounds;
        int num = 0;
        while (deltaError > MasterVariables.EPSILON) {
            // Calculate the sigma and omega lines for the current guess.
            sigma = new Line (points.subList (sigmaBounds[0], sigmaBounds[1]));
            omega = new Line (points.subList (omegaBounds[0], omegaBounds[1]));
            // Save the previous error.
            previousError = error;
            // Calculate the current error.
            error = 0.0D;
            for (int i = 1; i < points.size (); i ++) {
                Double sigmaError = squaredError (points.get (i), sigma);
                Double omegaError = squaredError (points.get (i), omega);
                if (sigmaError < threshold && sigmaError < omegaError) {
                    // Save the bounds of the points for the sigma line.
                    if (i < sigmaBounds[0]) {
                        sigmaBounds[0] = i;
                    }
                    if (i > sigmaBounds[1]) {
                        sigmaBounds[1] = i;
                        omegaBounds[0] = i + 1;
                    }
                    // Add to the total error.
                    error += sigmaError;
                }
                else if (omegaError < threshold) {
                    // Save the bounds of the points for the omega line.
                    if (i < omegaBounds[0]) {
                        omegaBounds[0] = i;
                        sigmaBounds[1] = i - 1;
                    }
                    if (i > omegaBounds[1]) {
                        omegaBounds[1] = i;
                    }
                    // Add to the total error.
                    error += omegaError;
                }
            }
            // Calculate the delta error.
            deltaError = Math.abs (error - previousError);
            // Make sure we don't get stuck in an infinite loop.
            num ++;
            if (error < previousError) {
                minSigmaBounds = sigmaBounds;
                minOmegaBounds = omegaBounds;
            }
            if (num > 10) {
                sigmaBounds = minSigmaBounds;
                omegaBounds = minOmegaBounds;
                break;
            }
        }
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
        result = new ParameterSet (
            npopEstimate, omegaEstimate, sigmaEstimate, null
        );
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
        return result;
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
     *  Set the parameter estimate result.
     */
    public void setResult (ParameterSet result) {
        this.result = result;
    }

    /**
     *  Set the slope and intercept of the omega line.
     *
     *  @param slope The slope of the omega line.
     *  @param intercept The Y intercept of the omega line.
     */
    public void setOmega (double slope, double intercept) {
        omega.m = slope;
        omega.b = intercept;
    }

    /**
     *  Set the slope and intercept of the sigma line.
     *
     *  @param slope The slope of the sigma line.
     *  @param intercept The Y intercept of the sigma line.
     */
    public void setSigma (double slope, double intercept) {
        sigma.m = slope;
        sigma.b = intercept;
    }

    /**
     *  Return the parameter estimate as a String.
     *
     *  @return The parameter estimate as a String.
     */
    public String toString () {
        return result.toString ();
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
     *  Find the first non-repeating point value on either end of the list
     *  and return the index values for the two bounds.
     *
     *  @param points The list of points.
     *  @return The bounds of the point list.
     */
    private Integer [] getBounds (List<Point> points) {
        Integer [] bounds = new Integer[2];
        // Look for the lower bound.
        Double previous = points.get (0).y;
        for (int i = 1; i < points.size (); i ++) {
            if (previous - points.get (i).y > MasterVariables.EPSILON) {
                bounds[0] = i - 1;
                break;
            }
            previous = points.get (i).y;
        }
        // Look for the upper bound.
        previous = points.get (points.size () - 1).y;
        for (int i = points.size () - 1; i > 0; i --) {
            if (points.get (i).y - previous > MasterVariables.EPSILON) {
                bounds[1] = i + 1;
                break;
            }
            previous = points.get (i).y;
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
        Double error = Double.MAX_VALUE;
        // If the line has a non-zero slope, calculate the squared distance of
        // the point from the line as the error.
        if (Math.abs (line.m) < MasterVariables.EPSILON) {
            Double a = Math.abs (-1 * line.m * point.x + point.y - line.b);
            Double b = Math.sqrt (line.m * line.m + 1);
            Double dist = a / b;
            // The error is the square of the distance.
            error = dist * dist;
        }
        return error;
    }

    /**
     *  Calculate the list of points from the binning results.
     *
     *  @param length The length of the sequences.
     *  @return The transformed binning results as a list of XY points.
     */
    private List<Point> getPoints (Integer length) {
        List<Point> points = new ArrayList<Point> ();
        for (BinLevel bin: binning.getBins ()) {
            Double crit = bin.getCrit ();
            Integer level = bin.getLevel ();
            // Don't include binning results == 1.
            if (level == 1) continue;
            // Transform the sequence criterion value into the number of SNPs.
            Double x = (1.0d - crit) * length;
            // Transform the number of sequence clusters (bins) into
            // log base 2 scale.
            Double y = Math.log (level) / logTwo;
            // Add the XY point to the list.
            points.add (new Point (x, y));
        }
        // Return the points in reverse order.
        Collections.reverse (points);
        return points;
    }

    private Integer nu;
    private Integer length;
    private Binning binning;

    private Line sigma;
    private Line omega;

    private Double threshold = 0.1d;

    private final Double logTwo = Math.log (2);

    /**
     * The resulting estimate for the parameter values.
     */
    private ParameterSet result;

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
