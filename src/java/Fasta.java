/*
 *    Ecotype Simulation models the sequence diversity within a bacterial
 *    clade as the evolutionary result of net ecotype formation and periodic
 *    selection, yielding a certain number of ecotypes.
 *
 *    Copyright (C) 2009-2014  Jason M. Wood, Montana State University
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

import java.io.File;
import java.io.RandomAccessFile;
import java.io.IOException;
import java.io.FileNotFoundException;
import java.lang.IndexOutOfBoundsException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Random;

/**
 *  Handles the input and output of fasta formatted text files.
 *
 *  @author Jason M. Wood
 *  @copyright GNU General Public License
 */
public class Fasta {

    /**
     *  Constructor for a Fasta object.
     */
    public Fasta () {
        file = null;
        current = 0;
        length = 0;
        size = 0;
        ids = new ArrayList<String> ();
        sequenceLengths = new HashMap<String, Integer> ();
        sequenceOffsets = new HashMap<String, Long> ();
    }

    /**
     *  Constructor for a Fasta object using a user supplied File.
     *
     *  @param fastaFile A File containing Fasta formatted data.
     */
    public Fasta (File fastaFile) throws InvalidFastaException {
        try {
            file = new RandomAccessFile (fastaFile, "rw");
            current = 0;
            length = 0;
            size = 0;
            ids = new ArrayList<String> ();
            sequenceLengths = new HashMap<String, Integer> ();
            sequenceOffsets = new HashMap<String, Long> ();
            long offset = 0;
            Sequence seq;
            while ((seq = parseSequence ()) != null) {
                String id = seq.getIdentifier ();
                int length = seq.length ();
                ids.add (id);
                sequenceOffsets.put (id, offset);
                sequenceLengths.put (id, length);
                if (length > this.length) this.length = length;
                size ++;
                // Update the file offset.
                offset = file.getFilePointer ();
            }
            file.seek (0);
        }
        catch (FileNotFoundException e) {
            throw new InvalidFastaException ("Fasta file not found: " + e);
        }
        catch (IOException e) {
            throw new InvalidFastaException ("IO error: " + e);
        }
    }

    /**
     *  Close this Fasta.
     */
    public void close () throws InvalidFastaException {
        if (file == null) return;
        try {
            file.close ();
        }
        catch (IOException e) {
            throw new InvalidFastaException ("IO error: " + e);
        }
    }

    /**
     *  Check if this Fasta object is valid.
     *
     *  @return True if this is a valid Fasta object, False if not.
     */
    public boolean isValid () {
        boolean valid = false;
        if (size > 0) {
            valid = true;
        }
        return valid;
    }

    /**
     *  Retrieves the next sequence from the Fasta formated file.
     *
     *  @return The next sequence in the buffer.
     */
    public Sequence nextSequence () throws InvalidFastaException {
        Sequence seq = new Sequence ();
        if (file == null) return seq;
        // Make sure there are sequences available.
        if (size > current) {
            // Seek to the location of the current sequence in the file.
            long offset = sequenceOffsets.get (ids.get (current));
            try {
                file.seek (offset);
            }
            catch (IOException e) {
                throw new InvalidFastaException ("IO Error: " + e);
            }
            // Parse the current sequence.
            seq = parseSequence ();
            // Increment the current sequence.
            current ++;
        }
        return seq;
    }

    /**
     *  Gets the Sequence with the provided index.
     *
     *  @param index The index of the Sequence to return.
     *  @return The Sequence linked to the provided index.
     */
    public Sequence get (int index) throws InvalidFastaException {
        Sequence seq = new Sequence ();
        if (index < size) {
            seq = get (ids.get (index));
        }
        return seq;
    }

    /**
     *  Gets the Sequence with the provided identifier.
     *
     *  @param id The identifier of the Sequence to return.
     *  @return The Sequence linked to the provided identifier.
     */
    public Sequence get (String id) throws InvalidFastaException {
        Sequence seq = new Sequence ();
        if (file == null) return seq;
        if (sequenceOffsets.containsKey (id)) {
            long offset = sequenceOffsets.get (id);
            try {
                file.seek (offset);
            }
            catch (IOException e) {
                throw new InvalidFastaException ("IO Error: " + e);
            }
            seq = parseSequence ();
        }
        return seq;
    }

    /**
     *  Gets the identifier with the provided index.
     *
     *  @param index The index of the identifier to return.
     *  @return The identifier linked to the provided index.
     */
    public String getIdentifier (int index) {
        String id = "";
        if (index < size) {
            id = ids.get (index);
        }
        return id;
    }

    /**
     *  Gets the sequence with the provided identifier.
     *
     *  @param id The ID of the sequence to return.
     *  @return The sequence linked to the provided identifier.
     */
    public String getSequence (String id) throws InvalidFastaException {
        Sequence seq = get (id);
        return seq.getSequence ();
    }

    /**
     *  Gets the sequence with the provided index.
     *
     *  @param index The index of the sequence to return.
     *  @return The sequence linked to the provided index.
     */
    public String getSequence (int index) throws InvalidFastaException {
        Sequence seq = get (index);
        return seq.getSequence ();
    }

    /**
     *  Gets the description with the provided identifier.
     *
     *  @param id The ID of the description to return.
     *  @return The description linked to the provided identifier.
     */
    public String getDescription (String id) throws InvalidFastaException {
        Sequence seq = get (id);
        return seq.getDescription ();
    }

    /**
     *  Gets the description with the provided index.
     *
     *  @param index The index of the description to return.
     *  @return The description linked to the provided index.
     */
    public String getDescription (int index) throws InvalidFastaException {
        Sequence seq = get (index);
        return seq.getDescription ();
    }

    /**
     *  Get the Identifiers stored in this Fasta.
     *
     *  @return ArrayList of Strings containing the identifiers.
     */
    public ArrayList<String> getIdentifiers () {
        ArrayList<String> idsCopy = new ArrayList<String> ();
        for (int i = 0; i < size; i ++) {
            idsCopy.add (ids.get (i));
        }
        return idsCopy;
    }

    /**
     *  Get the sequences stored in this Fasta.
     *
     *  @return ArrayList of Strings containing the sequences.
     */
    public ArrayList<String> getSequences () throws InvalidFastaException {
        ArrayList<String> sequences = new ArrayList<String> ();
        int oldCurrent = current;
        Sequence seq = null;
        while ((seq = nextSequence ()) != null) {
            sequences.add (seq.getSequence ());
        }
        current = oldCurrent;
        return sequences;
    }


    /**
     *  Write the next sequence from the Fasta formated file.
     *
     *  @return The next sequence in the buffer.
     */
    public void writeSequence (Sequence seq) throws InvalidFastaException {
        String id = seq.getIdentifier ();
        int length = seq.length ();
        if (file == null) return;
        try {
            // Seek to the end of the file.
            long offset = file.length ();
            file.seek (offset);
            // Write the Sequence to the file.
            file.writeChars (seq.toString ());
            // Add the sequence to this Fasta object.
            ids.add (id);
            sequenceLengths.put (id, length);
            sequenceOffsets.put (id, offset);
            if (length > this.length) this.length = length;
            size ++;
        }
        catch (IOException e) {
            throw new InvalidFastaException ("IO Error: " + e);
        }
    }

    /**
     *  Returns the length of the longest sequence.
     */
    public int length () {
        return length;
    }

    /**
     *  Returns the number of sequences stored in this Fasta object.
     */
    public int size () {
        return size;
    }

    /**
     *  Private method to parse a Fasta formatted sequence.
     *
     *  @return A Sequence object.
     */
    private Sequence parseSequence () throws InvalidFastaException {
        Sequence seq = null;
        try {
            // Don't try to read past the end of the file.
            if (file.getFilePointer () == file.length ()) {
                return seq;
            }
            // Read the line at the current location into a buffer.
            String buffer = file.readLine ();
            // Make sure there is a sequence at the current location.
            if (buffer.charAt (0) != '>') {
                throw new InvalidFastaException (
                    "Not a fasta formated sequence."
                );
            }
            // Grab the sequence identifier and description from the buffer.
            String[] header = buffer.split ("\\s+", 2);
            String id = header[0].substring (1);
            String description = "";
            if (header.length == 2) {
                description = header[1];
            }
            // Save the current file offset.
            long offset = file.getFilePointer ();
            // Grab the sequence data from the buffer.
            String sequence = "";
            while ((buffer = file.readLine ()) != null) {
                // Stop when the next sequence is found.
                if (buffer.charAt (0) == '>') {
                    file.seek (offset);
                    break;
                }
                // Store the sequence data.
                sequence += buffer;
                // Update the file offset.
                offset = file.getFilePointer ();
            }
            // Store the sequence.
            if (! id.equals ("null")) {
                seq = new Sequence (id, description, sequence);
            }
        }
        catch (IOException e) {
            throw new InvalidFastaException ("IO Error: " + e);
        }
        return seq;
    }

    private RandomAccessFile file;
    private ArrayList<String> ids;
    private HashMap<String, Integer> sequenceLengths;
    private HashMap<String, Long> sequenceOffsets;
    private int current;
    private int length;
    private int size;

}
