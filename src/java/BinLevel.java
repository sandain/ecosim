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
 *  Stores the bin levels for the Binning object.
 *
 *  @author Jason M. Wood
 *  @copyright GNU General Public License
 */
public class BinLevel {

    /**
     *  Constructor for the BinLevel object.
     *
     *  @param crit The crit level.
     *  @param level The number of clusters.
     */
    public BinLevel(Float crit, Integer level) {
        this.crit = crit;
        this.level = level;
    }

    /**
     *  Get the crit value.
     *
     *  @return The crit value.
     */
    public Float getCrit() {
        return crit;
    }

    /**
     *  Get the number of clusters for this crit level.
     *
     *  @return The number of clusters.
     */
    public Integer getLevel() {
        return level;
    }

    /**
     *  Returns this BinLevel as a String.
     *
     *  @return This BinLevel as a String.
     */
    public String toString() {
        return String.format("%5.3f: %d", crit, level);
    }

    private Float crit;
    private Integer level;
}
