/*
 *    Ecotype Simulation models the sequence diversity within a bacterial
 *    clade as the evolutionary result of net ecotype formation and periodic
 *    selection, yielding a certain number of ecotypes.
 *
 *    Copyright (C) 2013-2015  Jason M. Wood, Montana State University
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
 *  An object to store the parameter values.
 *
 *  @author Jason M. Wood
 *  @copyright GNU General Public License
 */
public class ParameterSet implements Comparable<ParameterSet> {

    /**
     *  A basic ParameterSet with all values set to zero.
     */
    public ParameterSet () {
        this (0L, 0.0D, 0.0D, 0.0D);
    }

    /**
     *  ParameterSet stores npop, omega, and sigma values, plus their
     *  likelihood.
     *
     *  @param npop The value of npop.
     *  @param omega The value of omega.
     *  @param sigma The value of sigma.
     *  @param likelihood Likelihood of the (npop, omega, sigma) parameter set.
     */
    public ParameterSet (Long npop, Double omega, Double sigma, Double likelihood) {
        this.npop = npop;
        this.omega = omega;
        this.sigma = sigma;
        this.likelihood = likelihood;
    }

    /**
     *  Compare the value of this parameter set with that of another.
     *
     *  @param other The other parameter set to comapre.
     *  @return 0 if the value of this parameter set is within EPSILON of
     *  being equal to the other parameter set's value. -1 if the
     *  the value is less.  1 if the value is more.
     */
    public int compareTo (ParameterSet other) {
        return likelihood.compareTo (other.getLikelihood ());
    }

    /**
     *  Get the npop value.
     *
     *  @return The npop value.
     */
    public Long getNpop () {
        return npop;
    }

    /**
     *  Get the omega value.
     *
     *  @return The omega value.
     */
    public Double getOmega () {
        return omega;
    }

    /**
     *  Get the sigma value.
     *
     *  @return The sigma value.
     */
    public Double getSigma () {
        return sigma;
    }

    /**
     *  Get the likelihood value.
     *
     *  @return The likelihood value.
     */
    public Double getLikelihood () {
        return likelihood;
    }

    /**
     *  Set the npop value.
     *
     *  @param npop The npop value.
     */
    public void setNpop (Long npop) {
        this.npop = npop;
    }

    /**
     *  Set the omega value.
     *
     *  @param omega The omega value.
     */
    public void setOmega (Double omega) {
        this.omega = omega;
    }

    /**
     *  Set the sigma value.
     *
     *  @param sigma The sigma value.
     */
    public void setSigma (Double sigma) {
        this.sigma = sigma;
    }

    /**
     *  Set the likelihood value.
     *
     *  @param likelihood The likelihood value.
     */
    public void setLikelihood (Double likelihood) {
        this.likelihood = likelihood;
    }

    /**
     *  Returns the parameter set as a string.
     *
     *  @return String containing the values of omega, sigma, and npop, and
     *  the likelihood of those values.
     */
    public String toString () {
        return String.format (
            "  npop:        %d\n" +
            "  omega:       %.4f\n" +
            "  sigma:       %.4f\n" +
            "  likelihood:  %.4f", npop, omega, sigma, likelihood
        );
    }

    private Long npop;
    private Double omega;
    private Double sigma;
    private Double likelihood;

}
