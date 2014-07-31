/*
 *    Ecotype Simulation models the sequence diversity within a bacterial
 *    clade as the evolutionary result of net ecotype formation and periodic
 *    selection, yielding a certain number of ecotypes.
 *
 *    Copyright (C) 2009-2013  Jason M. Wood, Montana State University
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

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
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
        length = 0;
        size = 0;
        ids = new ArrayList<String> ();
        descriptionHash = new HashMap<String, String> ();
        sequenceHash = new HashMap<String, String> ();
    }

    /**
     *  Constructor for a Fasta object using a user supplied File.
     *
     *  @param fastaFile A File containing Fasta formatted data.
     */
    public Fasta (File fastaFile) throws InvalidFastaException {
        file = fastaFile;
        length = 0;
        size = 0;
        parseFasta (fastaFile);
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
     *  Puts an ID and Sequence into this Fasta object.
     *
     *  @param id The ID of the sequence to add.
     *  @param sequence The sequence to add.
     */
    public void put (String id, String sequence) {
        put (id, "", sequence);
    }

    /**
     *  Puts an ID and Sequence into this Fasta object.
     *
     *  @param id The ID of the sequence to add.
     *  @param description The descripiton of the sequence to add.
     *  @param sequence The sequence to add.
     */
    public void put (String id, String description, String sequence) {
        ids.add (id);
        descriptionHash.put (id, description);
        sequenceHash.put (id, sequence.toLowerCase ());
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
        String value = null;
        if (id != null) {
            value = sequenceHash.get (id);
        }
        return value;
    }

    /**
     *  Gets the sequence with the provided index.
     *
     *  @param index The index of the sequence to return.
     *  @return The sequence linked to the provided index.
     */
    public String getSequence (int index) throws InvalidFastaException {
        String value = null;
        try {
            String id = ids.get (index);
            value = getSequence (id);
        }
        catch (IndexOutOfBoundsException e) {
            e.printStackTrace ();
        }
        return value;
    }

    /**
     *  Changes the sequence for the provided identifier.
     *
     *  @param id The identifier of the sequence to change.
     *  @param sequence The new sequence.
     */
    public void setSequence (String id, String sequence) {
        if (id != null) {
            sequenceHash.put (id, sequence.toLowerCase ());
        }
    }

    /**
     *  Changes the sequence at the provided index.
     *
     *  @param index The index of the sequence to change.
     *  @param sequence The new sequence.
     */
    public void setSequence (int index, String sequence) {
        try {
            String id = ids.get (index);
            setSequence (id, sequence);
        }
        catch (IndexOutOfBoundsException e) {
            e.printStackTrace ();
        }
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
        for (int i = 0; i < size; i ++) {
            sequences.add (sequenceHash.get (ids.get (i)));
        }
        return sequences;
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
     *  Save the sequence data in this object to a Fasta formatted file.
     *
     *  @param fileName File name to write the sequences to.
     *  @return True if the save was a success, False otherwise.
     */
    public boolean save (String fileName) throws InvalidFastaException {
        return save (new File (fileName));
    }

    /**
     *  Save the sequence data in this object to a Fasta formatted file.
     *
     *  @param fastaFile File to write the sequences to.
     *  @return True if the save was a success, False otherwise.
     */
    public boolean save (File fastaFile) throws InvalidFastaException {
        boolean success = false;
        BufferedWriter out = null;
        try {
            out = new BufferedWriter (new FileWriter (fastaFile));
            for (int i = 0; i < size; i ++) {
                String key = ids.get (i);
                out.write (">" + key + " ");
                out.write (descriptionHash.get (key) + "\n");
                out.write (sequenceHash.get (key) + "\n");
            }
        }
        catch (IOException e) {
            System.out.println ("Error writing to output file.");
        }
        finally {
            if (out != null) {
                try {
                    out.close ();
                    this.file = fastaFile;
                    success = true;
                }
                catch (IOException e) {
                    System.out.println ("Error closing output file.");
                }
            }
        }
        return success;
    }

    /**
     *  Create a new Fasta object containing a random subset of the
     *  sequences contained in this Fasta object.
     *
     *  @param number The number of sequences to select.
     *  @return The random subset of the sequences.
     */
    public Fasta randomSubset (int number) {
        Fasta randomFasta = new Fasta ();
        ArrayList<String> oldIds = getIdentifiers ();
        Random random = new Random ();
        for (int i = 0; i < number; i ++) {
            // Choose a random sequence to add to the random fasta.
            int rand = random.nextInt (oldIds.size ());
            String id = oldIds.get (rand);
            randomFasta.put (
                id,
                descriptionHash.get (id),
                sequenceHash.get (id)
            );
            // Remove that sequence from being picked again.
            oldIds.remove (id);
        }
        return randomFasta;
    }

    /**
     *  Private method to parse through a Fasta formatted text file.
     */
    private void parseFasta (File fastaFile) {
        BufferedReader input = null;
        sequenceHash = new HashMap<String, String> ();
        descriptionHash = new HashMap<String, String> ();
        ids = new ArrayList<String> ();
        try {
            input = new BufferedReader (new FileReader (fastaFile));
            String id = "null";
            String description = "";
            String sequence = "";
            String line = null;
            while ((line = input.readLine ()) != null) {
                if (! line.isEmpty () && line.charAt (0) == '>') {
                    // This line contains a header.
                    if (! id.equals ("null")) {
                        // Save the previous id/sequence before parsing the
                        // line.
                        sequenceHash.put (id, sequence.toLowerCase ());
                        descriptionHash.put (id, description);
                        int length = sequence.length ();
                        ids.add (id);
                        if (length > this.length) this.length = length;
                        size ++;
                        sequence = "";
                    }
                    // Grab the id out of the header.
                    String[] header = line.split ("\\s+", 2);
                    id = header[0].substring (1);
                    if (header.length == 2) {
                        description = header[1];
                    }
                }
                else {
                    // This line contains sequence (or it is empty and we can
                    // concatenate it anyway).
                    sequence = sequence + line;
                }
            }
            // Save the last sequence.
            if (! id.equals ("null")) {
                sequenceHash.put (id, sequence.toLowerCase ());
                descriptionHash.put (id, description);
                ids.add (id);
                size ++;
            }
        }
        catch (IOException e) {
            System.out.println ("Error reading from input file.");
        }
        finally {
            if (input != null) {
                try {
                    input.close ();
                }
                catch (IOException e) {
                    System.out.println ("Error closing input file.");
                }
            }
        }
    }

    private File file;
    private ArrayList<String> ids;
    private HashMap<String, String> descriptionHash;
    private HashMap<String, String> sequenceHash;
    private int length;
    private int size;

}
