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

import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import org.xml.sax.Attributes;
import org.xml.sax.InputSource;
import org.xml.sax.XMLReader;
import org.xml.sax.helpers.DefaultHandler;
import org.xml.sax.helpers.XMLReaderFactory;

/**
 *  Perform Input and Output of the options in XML format.
 *
 *  @author Jason Wood
 */
public class OptionsIO {

    /**
     *  No constructor needed, only used in a static context.
     */
    public OptionsIO() {
    }

    /**
     *  Saves the options in XML format.
     *
     *  @param options The options to save.
     */
    public static void save(Options options) {
        try {
            BufferedWriter out = new BufferedWriter(new FileWriter(options.OPTIONS_FILE));
            // Output the XML header.
            out.write("<?xml version=\"1.0\" encoding=\"UTF-8\"?>");
            out.newLine();
            // Output the options element.
            out.write("<options>");
            out.newLine();
            // Output each tab element.
            String[] tabs = options.getTabs();
            for (int i = 0; i < tabs.length; i ++) {
                out.write("  <tab name=\"" + tabs[i] + "\">");
                out.newLine();
                // Output all of the options for this tab.
                String[] keys = options.getKeys(tabs[i]);
                for (int j = 0; j < keys.length; j ++) {
                    String value = options.getValue(tabs[i], keys[j]);
                    out.write("    <option name=\"" + keys[j] + "\" value=\"" + value + "\" />");
                    out.newLine();
                }
                out.write("  </tab>");
                out.newLine();
            }
            out.write("</options>");
            out.newLine();
            out.close();
        }
        catch (IOException e) {
            e.printStackTrace();
        }
    }
  
    /**
     *  Loads the options file.
     *
     *  @return Options.
     */
    public static Options load() {
        Options options;
        XMLHandler handler = new XMLHandler();
        try {
            XMLReader xr = XMLReaderFactory.createXMLReader();
            xr.setContentHandler(handler);
            xr.setErrorHandler(handler);
            xr.parse(new InputSource(new FileReader(Options.OPTIONS_FILE)));
            options = handler.getOptions();
        }
        catch (Exception e) {
            // File not found, use the default values.
            options = new Options();
        }
        return options;
    }

    /**
     *  Create a DefaultHandler for parsing the options file.
     */
    private static class XMLHandler extends DefaultHandler  {

        public XMLHandler() {
            options = new Options();
            activeTab = "none";
        }

        /**
         *  Override the startElement method in DefaultHandler to grab the key/value pairs
         *  stored in the XML document.
         */
        public void startElement (String uri, String localName, String qName, Attributes attrs) {
            // Figure out which tab we are working with.
            if (localName.equals("tab")) {
                activeTab = attrs.getValue(uri, "name");
            }
            else if (localName.equals("option")) {
                // Grab the directory options.
                if (activeTab.equals(Options.DIRECTORY_OPTIONS)) {
                    options.setDirectoryOption(attrs.getValue(uri, "name"), attrs.getValue(uri, "value"));
                }
                // Grab the helper programs.
                else if (activeTab.equals(Options.HELPER_PROGRAMS)) {
                    options.setHelperProgram(attrs.getValue(uri, "name"), attrs.getValue(uri, "value"));
                }
            }
        }

        /**
         *  Override the endElement method in DefaultHandler to help keep track of our location in
         *  the XML document.
         */
        public void endElement (String uri, String localName, String qName) {
            // Look for the tab end element.
            if (localName.equals("tab")) {
                activeTab = "none";
            }
        }

        /**
         *  Returns the options stored in the XML file.
         *
         *  @return Options - The options.
         */
        public Options getOptions() {
            return options;
        }

        private String activeTab;

        /**
         *  The options.
         */
        private Options options;
    }
}
