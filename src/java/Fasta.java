/*
 *    Ecotype Simulation models the sequence diversity within a bacterial
 *    clade as the evolutionary result of net ecotype formation and periodic
 *    selection, yielding a certain number of ecotypes.
 *
 *    Copyright (C) 2009-2026 Jason M. Wood
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
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.HashSet;

/**
 *  Reads the information needed from an aligned Fasta formatted file: the
 *  identifier of the outgroup, which is the first sequence, the length of
 *  the alignment, and the number of sequences.
 *
 *  The file is read once, line by line, and is closed before the constructor
 *  returns; the sequences themselves are not kept. While reading, the file
 *  is verified to start with a header line, to give every sequence a unique
 *  identifier and some sequence data, and to give every sequence the same
 *  length as the outgroup.
 *
 *  @author Jason M. Wood
 *  @copyright GNU General Public License
 */
public class Fasta {

    /**
     *  Read and verify a Fasta formatted file.
     *
     *  @param fastaFile A File containing Fasta formatted data.
     *  @throws InvalidFastaException If the file can not be read, or is not
     *    an aligned Fasta formatted file with unique identifiers.
     */
    public Fasta (File fastaFile) throws InvalidFastaException {
        try (BufferedReader reader = new BufferedReader (
            new FileReader (fastaFile)
        )) {
            read (reader);
        }
        catch (FileNotFoundException e) {
            throw new InvalidFastaException ("Fasta file not found: " + e);
        }
        catch (IOException e) {
            throw new InvalidFastaException ("Fasta IO error: " + e);
        }
    }

    /**
     *  Get the identifier of the outgroup, the first sequence in the file.
     *
     *  @return The identifier of the outgroup.
     */
    public String getOutgroup () {
        return outgroup;
    }

    /**
     *  Get the length of the alignment.
     *
     *  @return The length of the sequences.
     */
    public int getLength () {
        return length;
    }

    /**
     *  Get the number of sequences in the file, including the outgroup.
     *
     *  @return The number of sequences.
     */
    public int size () {
        return size;
    }

    /**
     *  Private method to read the Fasta formatted data, recording the
     *  outgroup, the length of the alignment, and the number of sequences.
     *
     *  @param reader The reader to read the Fasta formatted data from.
     */
    private void read (BufferedReader reader)
        throws IOException, InvalidFastaException {
        HashSet<String> identifiers = new HashSet<String> ();
        // The identifier and sequence length of the current sequence.
        String identifier = null;
        int sequenceLength = 0;
        String line;
        while ((line = reader.readLine ()) != null) {
            line = line.trim ();
            // Skip empty lines.
            if (line.length () == 0) continue;
            if (line.charAt (0) == '>') {
                // A header line ends the previous sequence.
                if (identifier != null) {
                    endSequence (identifier, sequenceLength);
                }
                // The identifier is the header up to the first white space.
                identifier = line.substring (1).split ("\\s+", 2)[0];
                if (identifier.length () == 0) {
                    throw new InvalidFastaException (
                        "Fasta file contains a sequence without an " +
                        "identifier."
                    );
                }
                if (! identifiers.add (identifier)) {
                    throw new InvalidFastaException (
                        "Fasta file contains duplicate sequence " +
                        "identifier: " + identifier
                    );
                }
                sequenceLength = 0;
            }
            else {
                // Sequence data must follow a header line.
                if (identifier == null) {
                    throw new InvalidFastaException (
                        "Not a Fasta formatted file: sequence data before " +
                        "the first header line."
                    );
                }
                sequenceLength += line.length ();
            }
        }
        if (identifier == null) {
            throw new InvalidFastaException (
                "Fasta file contains no sequences."
            );
        }
        endSequence (identifier, sequenceLength);
    }

    /**
     *  Private method to verify a sequence once all of its data has been
     *  read, and to record it. The first sequence is the outgroup and sets
     *  the length of the alignment.
     *
     *  @param identifier The identifier of the sequence.
     *  @param sequenceLength The length of the sequence.
     */
    private void endSequence (String identifier, int sequenceLength)
        throws InvalidFastaException {
        if (sequenceLength == 0) {
            throw new InvalidFastaException (
                "Fasta file contains a sequence without data: " + identifier
            );
        }
        if (outgroup == null) {
            outgroup = identifier;
            length = sequenceLength;
        }
        else if (sequenceLength != length) {
            throw new InvalidFastaException (String.format (
                "Sequence %s has length %d, but the outgroup %s has length " +
                "%d. The sequences must be aligned.",
                identifier, sequenceLength, outgroup, length
            ));
        }
        size ++;
    }

    private String outgroup;
    private int length;
    private int size;

}
