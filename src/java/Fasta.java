/*
 *    Ecotype Simulation models the sequence diversity within a bacterial clade as
 *    the evolutionary result of net ecotype formation, periodic selection,
 *    and drift, yielding a certain number of ecotypes.
 * 
 *    Copyright (C) 2009  Fred Cohan, Wesleyan University
 *                        Carlo Francisco, Wesleyan University
 *                        Danny Krizanc, Wesleyan University
 *                        Andrew Warner, Wesleyan University
 *                        Jason Wood, Montana State University
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

/**
 *  Handles the input and output of Fasta formatted text files.
 *
 *  @author Jason Wood
 */
public class Fasta {

    /**
     *  Constructor for a Fasta object.
     */
    public Fasta () {
        length = 0;
        fastaHash = new HashMap<String, String>();
        ids = new ArrayList<String>();
    }

    /**
     *  Constructor for a Fasta object using a user supplied file name.
     *
     *  @param fileName The path and file name of the file containing the Fasta data.
     */
    public Fasta (String fileName) {
        parseFasta(new File(fileName));
    }

    /**
     *  Constructor for a Fasta object using a user supplied File.
     *
     *  @param file A File containing the Fasta data.
     */
    public Fasta (File file) {
        parseFasta(file);
    }

    /**
     *  Constructor for a Fasta object, used for making a copy of a Fasta object.
     *
     *  @param fasta The Fasta object that contains the data for this new Fasta object.
     */
    public Fasta (Fasta fasta) {
        length = fasta.length();
        fastaHash = new HashMap<String, String>(fasta.getHash());
        ids = new ArrayList<String>(fasta.getIds());
    }

    /**
     *  Constructor for a Fasta object, used for making a copy of a Fasta object, using
     *  only the selected IDs.
     *
     *  @param fasta The Fasta object that contains the data for this new Fasta object.
     *  @param selectedIds The IDs to keep in this new Fasta.
     */
    public Fasta (Fasta fasta, ArrayList<String> selectedIds) {
        ids = selectedIds;
        length = fasta.length();
        HashMap<String, String> hash = new HashMap<String, String>(fasta.getHash());
        fastaHash = new HashMap<String, String>();
        for (int i = 0; i < ids.size(); i ++) {
            fastaHash.put(ids.get(i), hash.get(ids.get(i)));
        }
    }

    /**
     *  Check if this Fasta object is valid.
     *
     *  @return True if this is a valid Fasta object, False if not.
     */
    public boolean isValid () {
        boolean valid = false;
        if (ids.size() > 0) {
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
        fastaHash.put(id, sequence);
        ids.add(id);
    }

    /**
     *  Gets a sequence from this Fasta object that has the provided ID.
     *
     *  @param id The ID of the sequence to return.
     *  @return The sequence linked to the provided ID.
     */
    public String get (String id) {
        String value = null;
        if (id != null) {
            value = fastaHash.get(id);
        }
        return value;
    }

    /**
     *  Gets a sequence from this Fasta object that has the provided index.
     *
     *  @param index The index of the sequence to return.
     *  @return The sequence linked to the provided index.
     */
    public String get (int index) {
        String value = null;
        try {
            String id = ids.get(index);
            value = get(id);
        }
        catch (IndexOutOfBoundsException e) {
            e.printStackTrace();
        }
        return value;
    }

    /**
     *  Get the HashMap containing the IDs and Sequences stored in this Fasta.
     *
     *  @return The HashMap containing the IDs and Sequences stored in this Fasta.
     */
    public HashMap<String, String> getHash() {
        return fastaHash;
    }

    /**
     *  Get the IDs stored in this Fasta.
     *
     *  @return ArrayList of Strings containing the IDs.
     */
    public ArrayList<String> getIds() {
        return ids;
    }

    /**
     *  Get the sequences stored in this Fasta.
     *
     *  @return ArrayList of Strings containing the sequences.
     */
    public ArrayList<String> getSequences() {
        ArrayList<String> sequences = new ArrayList<String>();
        for (int i = 0; i < ids.size(); i ++) {
            sequences.add(fastaHash.get(ids.get(i)));
        }
        return sequences;
    }

    /**
     *  Removes a sequences from this Fasta object.
     *
     *  @param id The ID of the sequence to remove.
     *  @return True if the sequence was successfully removed, False if not.
     */
    public boolean remove (String id) {
        boolean success = false;
        if (id != null) {
            String value = fastaHash.remove(id);
            ids.remove(id); 
            if (value != null) {
                success = true;
            }
        }
        return success;
    }

    /**
     *  Removes a sequence from this Fasta object.
     *
     *  @param index The index of the sequence to remove.
     *  @return True if the sequences was successfully removed, False if not.
     */
    public boolean remove (int index) {
        boolean success = false;
        try {
            String id = ids.get(index);
            success = remove(id);
        }
        catch (IndexOutOfBoundsException e) {
            e.printStackTrace();
        }
        return success;
    }

    /**
     *  Returns the length of the longest sequence stored in this Fasta object.
     */
    public int length () {
        return length;
    }

    /**
     *  Returns the number of sequences stored in this Fasta object.
     */
    public int size () {
        return ids.size();
    }

    /**
     *  Loads the given Fasta formatted file into this object.
     *
     *  @param fileName The path and file name to load.
     */
    public void load (String fileName) {
        parseFasta(new File(fileName));
    }

    /**
     *  Loads the given Fasta formatted file into this object.
     *
     *  @param file The file to load.
     */
    public void load (File file) {
        parseFasta(file);
    }

    /**
     *  Save the Fasta data in this object to an Ecotype Simulation formatted Fasta file.
     *
     *  @param fileName File name to write the Fasta data to.
     *  @return True if the save was a success, False otherwise.
     */
    public boolean save (String fileName) {
        return save(new File(fileName));
    }

    /**
     *  Save the Fasta data in this object to an Ecotype Simulation formatted Fasta file.
     *
     *  @param file File to write the Fasta data to.
     *  @return True if the save was a success, False otherwise.
     */
    public boolean save (File file) {
        boolean success = false;
        try {
            BufferedWriter out = new BufferedWriter(new FileWriter(file));
            for (int i = 0; i < ids.size(); i ++) {
                String key = ids.get(i);
                out.write (
                    ">" + key + System.getProperty ("line.separator")
                );
                out.write (
                    fastaHash.get (key) +
                    System.getProperty ("line.separator")
                );
            }
            out.close();
            success = true;
        }
        catch (IOException e) {
            System.out.println("Error writing to output file.");
        }
        return success;
    }

    /**
     *  Private method to parse through a Fasta formatted text file.
     */
    private void parseFasta (File fastaFile) {
        fastaHash = new HashMap<String, String>();
        ids = new ArrayList<String>();
        length = 0;
        try {
            BufferedReader input = new BufferedReader(new FileReader(fastaFile));
            String id = "null";
            String sequence = "";
            String line = null;
            while ((line = input.readLine()) != null) {
                if (! line.isEmpty() && line.charAt(0) == '>') {
                    // This line contains a header.
                    if (! id.equals("null")) {
                        // Save the previous id/sequence before parsing the line.
                        fastaHash.put(id, sequence);
                        ids.add(id);
                        if (sequence.length() > length) {
                            length = sequence.length();
                        }
                        sequence = "";
                    }
                    // Grab the id out of the header.
                    String[] header = line.split("\\s+", 2);
                    id = header[0].substring(1);
                }
                else {
                    // This line contains sequence (or it is empty and we can concatenate it anyway).
                    sequence = sequence + line;
                }
            }
            input.close();
            // Save the last sequence.
            if (! id.equals("null")) {
                fastaHash.put(id, sequence);
                ids.add(id);
                if (sequence.length() > length) {
                    length = sequence.length();
                }
            }
        }
        catch (IOException e) {
            System.out.println("Error reading from input file.");
        }
    }

    private int length;
    private HashMap<String, String> fastaHash;
    private ArrayList<String> ids;

}
