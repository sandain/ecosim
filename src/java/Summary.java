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

import ecosim.tree.Tree;

import java.util.ArrayList;
import java.util.Observable;

/**
 *  The summary class is used to display summary information to the user.
 *  This class utilizes the observer pattern to update listeners when the
 *  summary information has changed.
 *
 *  @author Jason M. Wood
 *  @copyright GNU General Public License
 */
public class Summary extends Observable {

    /**
     *  An object to store summary information about the simulation.
     */
    public Summary () {
        bins = new ArrayList<BinLevel> ();
        nu = 0;
        length = 0;
        outgroup = "";
        estimate = new ParameterEstimate ();
        hillclimbing = new ParameterSet ();
        confidenceInterval = new ParameterSet[] { 
            new ParameterSet (), new ParameterSet ()
        };
        demarcationList = new ArrayList<Demarcation> ();
    }

    /**
     *  Get the phylogenetic tree.
     *
     *  @return The tree.
     */
    public Tree getTree () {
        return tree;
    }

    /**
     *  Get the bins.
     *
     *  @return The bins.
     */
    public ArrayList<BinLevel> getBins () {
        return bins;
    }

    /**
     *  Get the number of environmental sequences.
     *
     *  @return The number of environmental sequences.
     */
    public Integer getNu () {
        return nu;
    }

    /**
     *  Get the length of the environmental sequences.
     *
     *  @return The length of the environmental sequences.
     */
    public Integer getLength () {
        return length;
    }

    /**
     *  Get the name of the outgroup.
     *
     *  @return The name of the outgroup.
     */
    public String getOutgroup () {
        return outgroup;
    }

    /**
     *  Get the parameter estimate results.
     *
     *  @return The parameter estimate results.
     */
    public ParameterEstimate getEstimate () {
        return estimate;
    }

    /**
     *  Get the hillclimbing results.
     *
     *  @return The hillclimbing results.
     */
    public ParameterSet getHillclimbing () {
        return hillclimbing;
    }

    /**
     *  Get the confidence interval.
     */
    public ParameterSet[] getConfidenceInterval () {
        return confidenceInterval;
    }

    /**
     *  Get the number of demarcation results.
     *
     *  @return The number of demarcation results.
     */
    public int getDemarcationCount () {
        return demarcationList.size ();
    }

    /**
     *  Get a specific demarcation result.
     *
     *  @param index The index of the demarcation result to return.
     *  @return The demarcation result requested.
     */
    public Demarcation getDemarcation (int index) {
        Demarcation last = null;
        if (index < demarcationList.size ()) {
            last = demarcationList.get (index);
        }
        return last;
    }

    /**
     *  Add the phylogenetic tree to the summary.
     *
     *  @param tree The phylogenetic tree.
     */
    public void setTree (Tree tree) {
        this.tree = tree;
        // Notify observers of the change.
        setChanged ();
        notifyObservers (this);
    }

    /**
     *  Add the binning results to the summary.
     *
     *  @param bins The binning results.
     */
    public void setBins (ArrayList<BinLevel> bins) {
        this.bins = bins;
        // Notify observers of the change.
        setChanged ();
        notifyObservers (this);
    }

    /**
     *  Add the number of environmental sequences to the summary.
     *
     *  @param nu The number of environmental sequences.
     */
    public void setNu (Integer nu) {
        this.nu = nu;
        // Notify observers of the change.
        setChanged ();
        notifyObservers (this);
    }

    /**
     *  Add the length of the environmental sequences to the summary.
     *
     *  @param length The length of the environmental sequences.
     */
    public void setLength (Integer length) {
        this.length = length;
        // Notify observers of the change.
        setChanged ();
        notifyObservers (this);
    }

    /**
     *  Add the name of the outgroup to the summary.
     *
     *  @param outgroup The name of the outgroup.
     */
    public void setOutgroup (String outgroup) {
        this.outgroup = outgroup;
        // Notify observers of the change.
        setChanged ();
        notifyObservers (this);
    }

    /**
     *  Add the parameter estimate result to the summary.
     *
     *  @param estimate The parameter estimate result.
     */
    public void setEstimate (ParameterEstimate estimate) {
        this.estimate = estimate;
        // Notify observers of the change.
        setChanged ();
        notifyObservers (this);
    }

    /**
     *  Add the hillclimbing result to the summary.
     *
     *  @param hillclimbing The hillclimbing result.
     */
    public void setHillclimbing (ParameterSet hillclimbing) {
        this.hillclimbing = hillclimbing;
        // Notify observers of the change.
        setChanged ();
        notifyObservers (this);
    }

    /**
     *  Add the confidence interval result to the summary.
     *
     *  @param confidenceInterval The confidence interval result.
     */
    public void setConfidenceInterval (ParameterSet[] confidenceInterval) {
        this.confidenceInterval = confidenceInterval;
        // Notify observers of the change.
        setChanged ();
        notifyObservers (this);
    }

    /**
     *  Add a demarcation result to the summary.
     *
     *  @param demarcation The demarcation result.
     */
    public void addDemarcation (Demarcation demarcation) {
        demarcationList.add (demarcation);
        // Notify observers of the change.
        setChanged ();
        notifyObservers (this);
    }

    private Tree tree;
    private Integer nu;
    private Integer length;
    private String outgroup;
    private ArrayList<BinLevel> bins;
    private ParameterEstimate estimate;
    private ParameterSet hillclimbing;
    private ParameterSet[] confidenceInterval;
    private ArrayList<Demarcation> demarcationList;

}
