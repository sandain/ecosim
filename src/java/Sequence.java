/*
 *    Ecotype Simulation models the sequence diversity within a bacterial
 *    clade as the evolutionary result of net ecotype formation and periodic
 *    selection, yielding a certain number of ecotypes.
 *
 *    Copyright (C) 2014-2019  Jason M. Wood, Montana State University
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


public class Sequence {

    /**
     *  Constructor for a Sequence object.
     */
    public Sequence () {
        this ("", "", "");
    }

    /**
     *  Constructor for a Sequence object.
     *
     *  @param identifier The identifier of this Sequence.
     *  @param sequence The sequence of this Sequence.
     */
    public Sequence (String identifier, String sequence) {
        this (identifier, "", sequence);
    }

    /**
     *  Constructor for a Sequence object.
     *
     *  @param identifier The identifier of this Sequence.
     *  @param description The description of this Sequence.
     *  @param sequence The sequence of this Sequence.
     */
    public Sequence (String identifier, String description,
        String sequence) {
        this.identifier = identifier;
        this.description = description;
        this.sequence = sequence;
    }

    /**
     *  Set the identifier of this Sequence object.
     *
     *  @param identifier The identifier of this Sequence.
     */
    public void setIdentifier (String identifier) {
        this.identifier = identifier;
    }

    /**
     *  Set the description of this Sequence object.
     *
     *  @param description The description of this Sequence.
     */
    public void setDescription (String description) {
        this.description = description;
    }

    /**
     *  Set the sequence of this Sequence object.
     *
     *  @param sequence The sequence of this Sequence.
     */
    public void setSequence (String sequence) {
        this.sequence = sequence;
    }

    /**
     *  Return the identifier of this Sequence object.
     *
     *  @return The identifier of this Sequence.
     */
    public String getIdentifier () {
        return identifier;
    }

    /**
     *  Return the description of this Sequence object.
     *
     *  @return The description of this Sequence.
     */
    public String getDescription () {
        return description;
    }

    /**
     *  Return the sequence of this Sequence object.
     *
     *  @return The sequence of this Sequence.
     */
    public String getSequence () {
        return sequence;
    }

    /**
     *  Return the length of the sequence of this Sequence object.
     *
     *  @return The length of the sequence of this Sequence.
     */
    public int length () {
        return sequence.length ();
    }

    /**
     *  Converts this Sequence object into a fasta formatted String.
     *
     *  @return This Sequence object as a fasta formatted String.
     */
    public String toString () {
        String fasta = ">" + identifier;
        if (description.length () > 0) {
            fasta += " " + description;
        }
        fasta += "\n" + sequence;
        return fasta;
    }

    private String identifier;
    private String description;
    private String sequence;

}
