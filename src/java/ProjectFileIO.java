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

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import org.xml.sax.Attributes;
import org.xml.sax.InputSource;
import org.xml.sax.XMLReader;
import org.xml.sax.helpers.DefaultHandler;
import org.xml.sax.helpers.XMLReaderFactory;

/**
 *  Perform input and output operations for the XML formatted project file.
 *
 *  @author Jason M. Wood
 *  @copyright GNU General Public License
 */
public class ProjectFileIO {

    /**
     *  Load and save the current project.
     *
     *  @param masterVariables The master variables.
     */
    public ProjectFileIO(MasterVariables masterVariables) {
        this.masterVariables = masterVariables;
        // Create new objects for each item in the project file.
        phylogeny = new Phylogeny(masterVariables);
        binning = new Binning(
            masterVariables, phylogeny
        );
        bruteforce = new Bruteforce(
            masterVariables, phylogeny, binning
        );
        hillclimb = new Hillclimb(
            masterVariables, phylogeny, binning, bruteforce.getBestResult()
        );
        omegaCI = new OmegaConfidenceInterval(
            masterVariables, phylogeny, binning, hillclimb
        );
        sigmaCI = new SigmaConfidenceInterval(
            masterVariables, phylogeny, binning, hillclimb
        );
        npopCI = new NpopConfidenceInterval(
            masterVariables, phylogeny, binning, hillclimb
        );
        demarcation = new Demarcation(
            masterVariables, phylogeny, binning, hillclimb
        );
    }

    /**
     *  Load and save the current project.
     *
     *  @param masterVariables The master variables.
     *  @param phylogeny The Phylogeny object.
     *  @param binning The Binning object.
     *  @param bruteforce The Bruteforce object.
     *  @param hillclimb The Hillclimb object.
     *  @param omegaCI The OmegaConfidenceInterval object.
     *  @param sigmaCI The SigmaConfidenceInterval object.
     *  @param npopCI The NpopConfidenceInterval object.
     *  @param demarcation The Demarcation object.
     */
    public ProjectFileIO(MasterVariables masterVariables, Phylogeny phylogeny,
        Binning binning, Bruteforce bruteforce, Hillclimb hillclimb,
        OmegaConfidenceInterval omegaCI, SigmaConfidenceInterval sigmaCI,
        NpopConfidenceInterval npopCI, Demarcation demarcation) {
        this.masterVariables = masterVariables;
        this.phylogeny = phylogeny;
        this.binning = binning;
        this.bruteforce = bruteforce;
        this.hillclimb = hillclimb;
        this.omegaCI = omegaCI;
        this.sigmaCI = sigmaCI;
        this.npopCI = npopCI;
        this.demarcation = demarcation;
    }

    /**
     *  Save the project file.
     *
     *  @param projectFile The project file to save.
     */
    public void save(File projectFile) {
        try {
            BufferedWriter out = new BufferedWriter(
                new FileWriter(projectFile)
            );
            // Output the XML header.
            out.write("<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n");
            out.write(String.format(
                "<ecosim type=\"savefile\" version=\"%s\">\n",
                masterVariables.getVersion()
            ));
            // Output the current criterion.
            out.write(
                "  <criterion value=\"" + 
                masterVariables.getCriterion() + "\"/>\n"
            );
            // Output the current PCR error.
            out.write(
                "  <pcrerror value=\"" + 
                masterVariables.getPCRError() + "\"/>\n"
            );
            // Output the phylogeny data.
            if (phylogeny != null && phylogeny.hasRun()) {
                out.write(String.format(
                    "  <phylogeny size=\"%d\" length=\"%d\">\n",
                    phylogeny.getNu(), phylogeny.length()
                ));
                out.write(String.format(
                    "    <tree value=\"%s\"/>\n",
                    phylogeny.getNewickTree().toString()
                ));
                ArrayList<String> identifiers = phylogeny.getIdentifiers();
                for (int i = 0; i < identifiers.size(); i ++) {
                    out.write(String.format(
                        "    <sequence identifier=\"%s\" sequence=\"%s\"/>\n",
                        phylogeny.getIdentifier(i),
                        phylogeny.getSequence(i)
                    ));
                }
                out.write(
                    "  </phylogeny>\n"
                );
            }
            // Output the binning data.
            if (binning != null && binning.hasRun()) {
                ArrayList<BinLevel> bins = binning.getBinLevels();
                out.write("  <binning>\n");
                out.write("    <bins size=\"" + bins.size() + "\">\n");
                // Output the crit levels and the number of bins.
                for (int i = 0; i < bins.size(); i ++) {
                    out.write(
                        "      <bin crit=\"" +
                        String.format("%.6f", bins.get(i).getCrit()) +
                        "\" value=\"" + bins.get(i).getLevel() + "\"/>\n"
                    );
                }
                out.write("    </bins>\n");
                out.write("  </binning>\n");
            }
            // Output the bruteforce data.
            if (bruteforce != null && bruteforce.hasRun()) {
                ArrayList<ParameterSet<Likelihood>> results =
                    bruteforce.getResults();
                double [] omegaRange = bruteforce.getOmegaRange();
                double [] sigmaRange = bruteforce.getSigmaRange();
                int [] npopRange = bruteforce.getNpopRange();
                out.write("  <bruteforce>\n");
                out.write(String.format(
                    "    <omegaRange low=\"%.5f\" high=\"%.5f\"/>\n",
                    omegaRange[0],
                    omegaRange[1]
                ));
                out.write(String.format(
                    "    <sigmaRange low=\"%.5f\" high=\"%.5f\"/>\n",
                    sigmaRange[0],
                    sigmaRange[1]
                ));
                out.write(String.format(
                    "    <npopRange low=\"%d\" high=\"%d\"/>\n",
                    npopRange[0],
                    npopRange[1]
                ));
                out.write(String.format(
                    "    <results size=\"%d\">\n", results.size()
                ));
                for (int i = 0; i < results.size(); i ++) {
                    out.write(String.format(
                        "      <result omega=\"%.5f\" sigma=\"%.5f\" " +
                        "npop=\"%d\" likelihood=\"%s\"/>\n",
                        results.get(i).getOmega(),
                        results.get(i).getSigma(),
                        results.get(i).getNpop(),
                        results.get(i).getValue().toString()
                    ));
                }
                out.write("    </results>\n");
                out.write("  </bruteforce>\n");
            }
            // Output the hillclimb data.
            if (hillclimb != null && hillclimb.hasRun()) {
                ParameterSet result = hillclimb.getResult();
                out.write("  <hillclimb>\n");
                out.write(String.format(
                    "    <result omega=\"%.5f\" sigma=\"%.5f\" " +
                    "npop=\"%d\" likelihood=\"%.5g\"/>\n",
                    result.getOmega(),
                    result.getSigma(),
                    result.getNpop(),
                    result.getValue()
                ));
                out.write("  </hillclimb>\n");
            }
            // Output the OmegaCI data.
            if (omegaCI != null && omegaCI.hasRun()) {
                double[] result = omegaCI.getResult();
                double[] likelihood = omegaCI.getLikelihood();
                out.write("  <omegaCI>\n");
                out.write(String.format(
                    "    <lower value=\"%.5f\" likelihood=\"%.5g\"/>\n",
                    result[0],
                    likelihood[0]
                ));
                out.write(String.format(
                    "    <upper value=\"%.5f\" likelihood=\"%.5g\"/>\n",
                    result[1],
                    likelihood[1]
                ));
                out.write("  </omegaCI>\n");
            }
            // Output the SigmaCI data.
            if (sigmaCI != null && sigmaCI.hasRun()) {
                double[] result = sigmaCI.getResult();
                double[] likelihood = sigmaCI.getLikelihood();
                out.write("  <sigmaCI>\n");
                out.write(String.format(
                    "    <lower value=\"%.5f\" likelihood=\"%.5g\"/>\n",
                    result[0],
                    likelihood[0]
                ));
                out.write(String.format(
                    "    <upper value=\"%.5f\" likelihood=\"%.5g\"/>\n",
                    result[1],
                    likelihood[1]
                ));
                out.write("  </sigmaCI>\n");
            }
            // Output the NpopCI data.
            if (npopCI != null && npopCI.hasRun()) {
                int[] result = npopCI.getResult();
                double[] likelihood = npopCI.getLikelihood();
                out.write("  <npopCI>\n");
                out.write(String.format(
                    "    <lower value=\"%d\" likelihood=\"%.5g\"/>\n",
                    result[0],
                    likelihood[0]
                ));
                out.write(String.format(
                    "    <upper value=\"%d\" likelihood=\"%.5g\"/>\n",
                    result[1],
                    likelihood[1]
                ));
                out.write("  </npopCI>\n");
            }
            // Output the Demarcation data.
            if (demarcation != null && demarcation.hasRun()) {
                ArrayList<ArrayList<String>> ecotypes =
                    demarcation.getEcotypes();
                out.write("  <demarcation>\n");
                out.write(String.format(
                    "    <ecotypes size=\"%d\">\n",
                    ecotypes.size()
                ));
                for (int i = 0; i < ecotypes.size(); i ++) {
                    ArrayList<String> ecotype = ecotypes.get(i);
                    out.write(String.format(
                        "      <ecotype number=\"%d\" size=\"%d\">\n",
                        (i + 1), ecotype.size()
                    ));
                    for (int j = 0; j < ecotype.size(); j ++) {
                        out.write(String.format(
                            "        <member name=\"%s\"/>\n",
                            ecotype.get(j)
                        ));
                    }
                    out.write("      </ecotype>\n");
                }
                out.write("    </ecotypes>\n");
                out.write("  </demarcation>\n");
            }
            out.write("</ecosim>\n");
            out.close();
        }
        catch (IOException e) {
            e.printStackTrace();
        }
    }

    /**
     *  Load the project file.
     *
     *  @param projectFile The project file to load.
     */
    public void load(File projectFile) {
        XMLHandler handler = new XMLHandler();
        try {
            XMLReader xr = XMLReaderFactory.createXMLReader();
            xr.setContentHandler(handler);
            xr.setErrorHandler(handler);
            xr.parse(new InputSource(new FileReader(projectFile)));
        }
        catch (Exception e) {
            e.printStackTrace();
        }
    }

    /**
     *  Get the Phylogeny object.
     *
     *  @return The Phylogeny object.
     */
    public Phylogeny getPhylogeny() {
        return phylogeny;
    }

    /**
     *  Get the Binning object.
     *
     *  @return The Binning object.
     */
    public Binning getBinning() {
        return binning;
    }

    /**
     *  Get the Bruteforce object.
     *
     *  @return The Bruteforce object.
     */
    public Bruteforce getBruteforce() {
        return bruteforce;
    }

    /**
     *  Get the Hillclimb object.
     *
     *  @return The Hillclimb object.
     */
    public Hillclimb getHillclimb() {
        return hillclimb;
    }

    /**
     *  Get the OmegaConfidenceInterval object.
     *
     *  @return The OmegaConfidenceInterval object.
     */
    public OmegaConfidenceInterval getOmegaCI() {
        return omegaCI;
    }

    /**
     *  Get the SigmaConfidenceInterval object.
     *
     *  @return The SigmaConfidenceInterval object.
     */
    public SigmaConfidenceInterval getSigmaCI() {
        return sigmaCI;
    }

    /**
     *  Get the NpopConfidenceInterval object.
     *
     *  @return The NpopConfidenceInterval object.
     */
    public NpopConfidenceInterval getNpopCI() {
        return npopCI;
    }

    /**
     *  Get the Demarcation object.
     *
     *  @return The Demarcation object.
     */
    public Demarcation getDemarcation() {
        return demarcation;
    }

    private MasterVariables masterVariables;
    private Phylogeny phylogeny;
    private Binning binning;
    private Bruteforce bruteforce;
    private Hillclimb hillclimb;
    private OmegaConfidenceInterval omegaCI;
    private SigmaConfidenceInterval sigmaCI;
    private NpopConfidenceInterval npopCI;
    private Demarcation demarcation;

    /**
     *  Create a DefaultHandler for parsing the options file.
     */
    private class XMLHandler extends DefaultHandler  {

        /**
         *  The XMLHandler constructor.
         */
        public XMLHandler() {
            elements = new ArrayList<String>();
            elements.add("phylogeny");
            elements.add("binning");
            elements.add("bruteforce");
            elements.add("hillclimb");
            elements.add("omegaCI");
            elements.add("sigmaCI");
            elements.add("npopCI");
            elements.add("demarcation");
            activeElement = "none";
            isProjectFile = false;
        }

        /**
         *  Override the startElement method in DefaultHandler to grab the
         *  key/value pairs stored in the XML document.
         */
        public void startElement (String uri, String localName, String qName,
            Attributes attrs) {
            // Make sure that this is a save file.
            if (localName.equals("ecosim") &&
                attrs.getValue(uri, "type").equals("savefile")) {
                isProjectFile = true;
            }
            if (isProjectFile) {
                // Update the active element.
                if (elements.contains(localName)) {
                    activeElement = localName;
                }
                // Look for the criterion element.
                if (localName.equals("criterion")) {
                    masterVariables.setCriterion(
                        new Integer(attrs.getValue(uri, "value")).intValue()
                    );
                }
                // Look for the pcrerror element.
                if (localName.equals("pcrerror")) {
                    masterVariables.setPCRError(
                        new Double(attrs.getValue(uri, "value")).doubleValue()
                    );
                }
                // Look for the elements within phylogeny.
                if (activeElement.equals("phylogeny")) {
                    // Look for the tree element.
                    if (localName.equals("tree")) {
                         phylogeny.loadTree(
                            attrs.getValue(uri, "value")
                        );
                    }
                    // Look for the sequence elements.
                    if (localName.equals("sequence")) {
                        phylogeny.put(
                            attrs.getValue(uri, "identifier"),
                            attrs.getValue(uri, "sequence")
                        );
                    }
                }
                // Look for elements within binning.
                if (activeElement.equals("binning")) {
                    if (localName.equals("bin")) {
                        binning.addBinLevel(new BinLevel(
                            new Float(attrs.getValue(uri, "crit")),
                            new Integer(attrs.getValue(uri, "value"))
                        ));
                    }
                }
                // Look for elements within bruteforce.
                if (activeElement.equals("bruteforce")) {
                    // Grab the omega range.
                    if (localName.equals("omegaRange")) {
                        double [] omegaRange = new double [2];
                        omegaRange[0] = new Double(
                            attrs.getValue(uri, "low")
                        ).doubleValue();
                        omegaRange[1] = new Double(
                            attrs.getValue(uri, "high")
                        ).doubleValue();
                        bruteforce.setOmegaRange(omegaRange);
                    }
                    // Grab the sigma range.
                    if (localName.equals("sigmaRange")) {
                        double [] sigmaRange = new double [2];
                        sigmaRange[0] = new Double(
                            attrs.getValue(uri, "low")
                        ).doubleValue();
                        sigmaRange[1] = new Double(
                            attrs.getValue(uri, "high")
                        ).doubleValue();
                        bruteforce.setSigmaRange(sigmaRange);
                    }
                    // Grab the sigma range.
                    if (localName.equals("npopRange")) {
                        int [] npopRange = new int [2];
                        npopRange[0] = new Integer(
                            attrs.getValue(uri, "low")
                        ).intValue();
                        npopRange[1] = new Integer(
                            attrs.getValue(uri, "high")
                        ).intValue();
                        bruteforce.setNpopRange(npopRange);
                    }
                    if (localName.equals("result")) {
                        bruteforce.addResult(new ParameterSet<Likelihood>(
                            new Double(attrs.getValue(uri, "omega")),
                            new Double(attrs.getValue(uri, "sigma")),
                            new Integer(attrs.getValue(uri, "npop")),
                            new Likelihood(
                                masterVariables,
                                attrs.getValue(uri, "likelihood")
                            )
                        ));
                    }
                }
                // Look for elements within hillclimb.
                if (activeElement.equals("hillclimb")) {
                    if (localName.equals("result")) {
                        hillclimb.setResult(
                            new ParameterSet<Double>(
                                new Double(attrs.getValue(uri, "omega")),
                                new Double(attrs.getValue(uri, "sigma")),
                                new Integer(attrs.getValue(uri, "npop")),
                                new Double(attrs.getValue(uri, "likelihood"))
                            )
                        );
                    }
                }
                // Look for elements within omegaCI.
                if (activeElement.equals("omegaCI")) {
                    if (localName.equals("lower")) {
                        double value = new Double(
                            attrs.getValue(uri, "value")
                        ).doubleValue();
                        double likelihood = new Double(
                            attrs.getValue(uri, "likelihood")
                        ).doubleValue();
                        omegaCI.setLowerResult(value, likelihood);
                    }
                    if (localName.equals("upper")) {
                        double value = new Double(
                            attrs.getValue(uri, "value")
                        ).doubleValue();
                        double likelihood = new Double(
                            attrs.getValue(uri, "likelihood")
                        ).doubleValue();
                        omegaCI.setUpperResult(value, likelihood);
                    }
                }
                // Look for elements within sigmaCI.
                if (activeElement.equals("sigmaCI")) {
                    if (localName.equals("lower")) {
                        double value = new Double(
                            attrs.getValue(uri, "value")
                        ).doubleValue();
                        double likelihood = new Double(
                            attrs.getValue(uri, "likelihood")
                        ).doubleValue();
                        sigmaCI.setLowerResult(value, likelihood);
                    }
                    if (localName.equals("upper")) {
                        double value = new Double(
                            attrs.getValue(uri, "value")
                        ).doubleValue();
                        double likelihood = new Double(
                            attrs.getValue(uri, "likelihood")
                        ).doubleValue();
                        sigmaCI.setUpperResult(value, likelihood);
                    }
                }
                // Look for elements within npopCI.
                if (activeElement.equals("npopCI")) {
                    if (localName.equals("lower")) {
                        int value = new Integer(
                            attrs.getValue(uri, "value")
                        ).intValue();
                        double likelihood = new Double(
                            attrs.getValue(uri, "likelihood")
                        ).doubleValue();
                        npopCI.setLowerResult(value, likelihood);
                    }
                    if (localName.equals("upper")) {
                        int value = new Integer(
                            attrs.getValue(uri, "value")
                        ).intValue();
                        double likelihood = new Double(
                            attrs.getValue(uri, "likelihood")
                        ).doubleValue();
                        npopCI.setUpperResult(value, likelihood);
                    }
                }
                // Look for elements within demarcation.
                if (activeElement.equals("demarcation")) {
                    if (localName.equals("ecotypes")) {
                        int size = new Integer(
                            attrs.getValue(uri, "size")
                        ).intValue();
                        ecotypes = new ArrayList<ArrayList<String>>(size);
                    }
                    if (localName.equals("ecotype")) {
                        int size = new Integer(
                            attrs.getValue(uri, "size")
                        ).intValue();
                        ecotypes.add(new ArrayList<String>(size));
                        ecotypeNumber = new Integer(
                            attrs.getValue(uri, "number")
                        ).intValue();
                    }
                    if (localName.equals("member") && ecotypes != null) {
                        ecotypes.get(ecotypeNumber - 1).add(
                            attrs.getValue(uri, "name")
                        );
                    }
                }
            }
        }

        /**
         *  Override the endElement method in DefaultHandler to help keep track
         *  of our location in the XML document.
         */
        public void endElement (String uri, String localName, String qName) {
            if (isProjectFile) {
                // Look for the end of the ecosim save file.
                if (localName.equals("ecosim")) {
                    isProjectFile = false;
                }
                // Look for the end of the phylogeny element.
                if (localName.equals("phylogeny")) {
                    phylogeny.setHasRun(true);
                }
                // Look for the end of the binning element.
                if (localName.equals("binning")) {
                    binning.setHasRun(true);
                }
                // Look for the end of the bruteforce element.
                if (localName.equals("bruteforce")) {
                    bruteforce.setHasRun(true);
                }
                // Look for the end of the hillclimb element.
                if (localName.equals("hillclimb")) {
                    hillclimb.setHasRun(true);
                }
                // Look for the end of the omegaCI element.
                if (localName.equals("omegaCI")) {
                    omegaCI.setHasRun(true);
                }
                // Look for the end of the sigmaCI element.
                if (localName.equals("sigmaCI")) {
                    sigmaCI.setHasRun(true);
                }
                // Look for the end of the npopCI element.
                if (localName.equals("npopCI")) {
                    npopCI.setHasRun(true);
                }
                // Look for the end of the demarcation element.
                if (localName.equals("demarcation")) {
                    demarcation.setEcotypes(ecotypes);
                    demarcation.setHasRun(true);
                }
                // Update the active element.
                if (elements.contains(localName)) {
                    activeElement = "none";
                }
            }
        }

        private String activeElement;

        private int ecotypeNumber;
        private ArrayList<ArrayList<String>> ecotypes;

        private boolean isProjectFile;

        private ArrayList<String> elements;
    }
}
