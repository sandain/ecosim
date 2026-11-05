/*
 *    Ecotype Simulation models the sequence diversity within a bacterial
 *    clade as the evolutionary result of net ecotype formation and periodic
 *    selection, yielding a certain number of ecotypes.
 *
 *    Copyright (C) 2016-2026  Jason M. Wood
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

import java.util.List;
import java.util.ArrayList;
import java.util.Collections;

/**
 *  Object to estimate the parameter values from the binning results.
 *
 *  The binning results are plotted as the number of sequence differences
 *  against the log base 2 of the number of bins. The flat ends of the curve
 *  are excluded, and the remaining points are split in two with a
 *  breakpoint search: every split that leaves at least two points on each
 *  side is tried, a straight line is fitted to each side, and the split with
 *  the smallest total squared error is kept. Only splits where both lines
 *  fall and the low divergence line is the steeper are accepted. The slope
 *  of the steeper line gives sigma, the slope of the shallower line gives
 *  omega, and the intersection of the two lines gives npop.
 *
 *  If no split is accepted, the estimate is not valid: omega and sigma are
 *  set to NaN and npop to 1.
 *
 *  @author Jason M. Wood
 *  @copyright GNU General Public License
 */
public class ParameterEstimate implements Runnable {

    /**
     *  The parameter estimate program.
     */
    public ParameterEstimate () {
        this (0, new Binning ());
    }

    /**
     *  The parameter estimate program.
     *
     *  @param length The length of the sequences being analyzed.
     *  @param binning The Binning object.
     */
    public ParameterEstimate (Integer length, Binning binning) {
        this.length = length;
        this.binning = binning;
        result = new ParameterSet ();
        sigma = new Line (new ArrayList<Point> ());
        omega = new Line (new ArrayList<Point> ());
        hasRun = false;
    }

    /**
     *  Run the parameter estimate program.
     */
    public void run () {
        // Calculate the list of points that the sigma and omega lines will
        // be fitted to.
        List<Point> points = getPoints (length);
        // Find the bounds of the points between the flat ends of the curve.
        Integer[] bounds = getBounds (points);
        // Search for the split of the points into a sigma line and an omega
        // line with the smallest total squared error.
        Line bestSigma = null;
        Line bestOmega = null;
        Double bestError = Double.MAX_VALUE;
        int first = bounds[0] + MIN_POINTS;
        int last = bounds[1] - MIN_POINTS;
        for (int split = first; split <= last; split ++) {
            List<Point> sigmaPoints = points.subList (bounds[0], split);
            List<Point> omegaPoints = points.subList (split, bounds[1]);
            Line sigmaLine = new Line (sigmaPoints);
            Line omegaLine = new Line (omegaPoints);
            // Only accept splits where both lines fall and the sigma line is
            // the steeper of the two.
            if (! (sigmaLine.m < omegaLine.m && omegaLine.m < 0.0d)) continue;
            Double error = sigmaLine.squaredError (sigmaPoints) +
                omegaLine.squaredError (omegaPoints);
            if (error < bestError) {
                bestError = error;
                bestSigma = sigmaLine;
                bestOmega = omegaLine;
            }
        }
        // Without an accepted split, the estimate is not valid.
        if (bestSigma == null) {
            sigma.m = Double.NaN;
            sigma.b = Double.NaN;
            omega.m = Double.NaN;
            omega.b = Double.NaN;
            result = new ParameterSet (1L, Double.NaN, Double.NaN, null);
            hasRun = true;
            return;
        }
        sigma = bestSigma;
        omega = bestOmega;
        // Omega is estimated from the slope of the omega line.
        Double omegaEstimate = -1.0d * omega.m;
        // Sigma is estimated from the slope of the sigma line.
        Double sigmaEstimate = -1.0d * sigma.m;
        // Npop is estimated by calculating the intersection of the omega
        // and sigma lines.
        Long npopEstimate = Math.round (Math.pow (
            2, omega.m * (sigma.b - omega.b) / (omega.m - sigma.m) + omega.b
        ));
        // Don't let the npop estimate fall below 1.
        if (npopEstimate < 1L) npopEstimate = 1L;
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
        // Initialize the bounds with that of the points array.
        Integer [] bounds = new Integer[] {
          Integer.valueOf (0), Integer.valueOf (points.size () - 1)
        };
        // Look for the lower bound.
        Double previous = points.get (0).y;
        for (int i = 1; i < points.size (); i ++) {
            if (previous - points.get (i).y > MainVariables.EPSILON) {
                bounds[0] = i - 1;
                break;
            }
            previous = points.get (i).y;
        }
        // Look for the upper bound.
        previous = points.get (points.size () - 1).y;
        for (int i = points.size () - 1; i > bounds[0]; i --) {
            if (points.get (i).y - previous > MainVariables.EPSILON) {
                bounds[1] = i + 1;
                break;
            }
            previous = points.get (i).y;
        }
        return bounds;
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

    private Integer length;
    private Binning binning;

    private Line sigma;
    private Line omega;

    private final Double logTwo = Math.log (2);

    /**
     *  The minimum number of points each line is fitted to.
     */
    private static final int MIN_POINTS = 2;

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

        /**
         *  Calculate the sum of the squared vertical distances of the points
         *  from this line.
         *
         *  @param points The points.
         *  @return The sum of the squared errors.
         */
        public Double squaredError (List<Point> points) {
            Double error = 0.0d;
            for (Point point: points) {
                Double residual = point.y - (m * point.x + b);
                error += residual * residual;
            }
            return error;
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
