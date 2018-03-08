/*
 *    Ecotype Simulation models the sequence diversity within a bacterial
 *    clade as the evolutionary result of net ecotype formation and periodic
 *    selection, yielding a certain number of ecotypes.
 *
 *    Copyright (C) 2009-2018 Jason M. Wood, Montana State University
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

import java.io.File;
import java.io.RandomAccessFile;
import java.io.IOException;
import java.io.FileNotFoundException;
import java.util.HashMap;

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
        outgroup = null;
    }

    /**
     *  Constructor for a Fasta object using a user supplied File.
     *
     *  @param fastaFile A File containing Fasta formatted data.
     */
    public Fasta (File fastaFile) throws InvalidFastaException {
        try {
            file = new RandomAccessFile (fastaFile, "rw");
            // Load the first sequence as the outgroup.
            outgroup = parseSequence ();
            if (outgroup == null) {
                throw new InvalidFastaException (
                    "Fasta file contains no sequences"
                );
            }
        }
        catch (FileNotFoundException e) {
            throw new InvalidFastaException ("Fasta file not found: " + e);
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
        return outgroup != null;
    }

    /**
     *  Get the outgroup Sequence.
     *
     *  @return The outgroup Sequence.
     */
    public Sequence getOutgroup () {
        return outgroup;
    }

    /**
     *  Retrieves the next sequence from the Fasta formated file.
     *
     *  @return The next sequence in the buffer.
     */
    public Sequence nextSequence () throws InvalidFastaException {
        if (file == null) return null;
        // Parse the current sequence.
        return parseSequence ();
    }

    /**
     *  Private method to parse a Fasta formatted sequence.
     *
     *  @return A Sequence object.
     */
    private Sequence parseSequence () throws InvalidFastaException {
        Sequence seq = null;
        try {
            // Read the line at the current location into a buffer.
            String buffer = "";
            // Ignore empty lines at the start of the file.
            while (buffer.length () == 0) {
                // Don't try to read past the end of the file.
                if (file.getFilePointer () == file.length ()) {
                    return seq;
                }
                // Read the line at the current location into a buffer.
                buffer = file.readLine ();
            }
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
                if (buffer.length () > 0 && buffer.charAt (0) == '>') {
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
    private Sequence outgroup;

}
