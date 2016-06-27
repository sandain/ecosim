/*
 *    Ecotype Simulation models the sequence diversity within a bacterial
 *    clade as the evolutionary result of net ecotype formation and periodic
 *    selection, yielding a certain number of ecotypes.
 *
 *    Copyright (C) 2013-2016  Jason M. Wood, Montana State University
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

import ecosim.tree.InvalidTreeException;
import ecosim.tree.Tree;

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
     *  @param execs The Execs object.
     */
    public ProjectFileIO (MasterVariables masterVariables,
        Execs execs) {
        this.masterVariables = masterVariables;
        this.execs = execs;
        Fasta fasta = new Fasta ();
        // Create new objects for each item in the project file.
        nu = 0;
        length = 0;
        outgroup = null;
        tree = null;
        binning = null;
        estimate = null;
        hillclimb = null;
        npopCI = null;
        omegaCI = null;
        sigmaCI = null;
        demarcation = null;
    }

    /**
     *  Load and save the current project.
     *
     *  @param masterVariables The master variables.
     *  @param execs The Execs object.
     *  @param nu The number of environmental sequences.
     *  @param length the length of the environmental sequences.
     *  @param outgroup The identifier of the outgroup sequence.
     *  @param tree The Tree object.
     *  @param binning The Binning object.
     *  @param estimate The ParameterEstimate object.
     *  @param hillclimb The Hillclimb object.
     *  @param npopCI The NpopConfidenceInterval object.
     *  @param omegaCI The OmegaConfidenceInterval object.
     *  @param sigmaCI The SigmaConfidenceInterval object.
     *  @param demarcation The Demarcation object.
     */
    public ProjectFileIO (MasterVariables masterVariables, Execs execs,
        Integer nu, Integer length, String outgroup, Tree tree,
        Binning binning, ParameterEstimate estimate, Hillclimb hillclimb,
        NpopConfidenceInterval npopCI, OmegaConfidenceInterval omegaCI,
        SigmaConfidenceInterval sigmaCI, Demarcation demarcation) {
        this.masterVariables = masterVariables;
        this.execs = execs;
        this.nu = nu;
        this.length = length;
        this.outgroup = outgroup;
        this.tree = tree;
        this.binning = binning;
        this.estimate = estimate;
        this.hillclimb = hillclimb;
        this.npopCI = npopCI;
        this.omegaCI = omegaCI;
        this.sigmaCI = sigmaCI;
        this.demarcation = demarcation;
    }

    /**
     *  Save the project file.
     *
     *  @param projectFile The project file to save.
     */
    public void save (File projectFile) {
        try {
            BufferedWriter out = new BufferedWriter (
                new FileWriter (projectFile)
            );
            // Output the XML header.
            out.write ("<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n");
            out.write (String.format (
                "<ecosim type=\"savefile\" version=\"%s\">\n",
                masterVariables.getVersion ()
            ));
            // Output the current criterion.
            out.write (
                "  <criterion value=\"" +
                masterVariables.getCriterion () + "\"/>\n"
            );
            // Output the phylogeny data.
            if (tree != null && tree.isValid ()) {
                out.write (String.format (
                    "  <phylogeny size=\"%d\" length=\"%d\">\n", nu, length
                ));
                out.write (String.format (
                    "    <outgroup value=\"%s\"/>\n", outgroup
                ));
                out.write (String.format (
                    "    <tree value=\"%s\"/>\n", tree.toString ()
                ));
                out.write (
                    "  </phylogeny>\n"
                );
            }
            // Output the binning data.
            if (binning != null) {
                ArrayList<BinLevel> bins = binning.getBins ();
                out.write ("  <binning>\n");
                out.write ("    <bins size=\"" + bins.size () + "\">\n");
                // Output the crit levels and the number of bins.
                for (int i = 0; i < bins.size (); i ++) {
                    out.write (
                        "      <bin crit=\"" +
                        String.format ("%.6f", bins.get (i).getCrit ()) +
                        "\" value=\"" + bins.get (i).getLevel () + "\"/>\n"
                    );
                }
                out.write ("    </bins>\n");
                out.write ("  </binning>\n");
            }
            // Output the parameter estimate data.
            if (estimate != null && estimate.hasRun ()) {
                ParameterSet result = estimate.getResult ();
                out.write ("  <estimate>\n");
                out.write (String.format ("    " +
                    "<result npop=\"%d\" omega=\"%.5f\" sigma=\"%.5f\"/>\n",
                    result.getNpop (), result.getOmega (), result.getSigma ()
                ));
                double omega[] = estimate.getOmega ();
                out.write (String.format ("    " +
                    "<omega slope=\"%.5f\" intercept=\"%.5f\"/>\n",
                    omega[0], omega[1]
                ));
                double sigma[] = estimate.getSigma ();
                out.write (String.format ("    " +
                    "<sigma slope=\"%.5f\" intercept=\"%.5f\"/>\n",
                    sigma[0], sigma[1]
                ));
                out.write ("  </estimate>\n");
            }
            // Output the hillclimb data.
            if (hillclimb != null && hillclimb.hasRun ()) {
                ParameterSet result = hillclimb.getResult ();
                out.write ("  <hillclimb>\n");
                out.write (String.format ("    " +
                    "<result npop=\"%d\" omega=\"%.5f\" sigma=\"%.5f\" " +
                    "likelihood=\"%.5g\"/>\n",
                    result.getNpop (), result.getOmega (), result.getSigma (),
                    result.getLikelihood ()
                ));
                out.write ("  </hillclimb>\n");
            }
            // Output the NpopCI data.
            if (npopCI != null && npopCI.hasRun ()) {
                Long [] result = npopCI.getResult ();
                Double [] likelihood = npopCI.getLikelihood ();
                out.write ("  <npopCI>\n");
                out.write (String.format (
                    "    <lower value=\"%d\" likelihood=\"%.5g\"/>\n",
                    result[0],
                    likelihood[0]
                ));
                out.write (String.format (
                    "    <upper value=\"%d\" likelihood=\"%.5g\"/>\n",
                    result[1],
                    likelihood[1]
                ));
                out.write ("  </npopCI>\n");
            }
            // Output the OmegaCI data.
            if (omegaCI != null && omegaCI.hasRun ()) {
                Double[] result = omegaCI.getResult ();
                Double[] likelihood = omegaCI.getLikelihood ();
                out.write ("  <omegaCI>\n");
                out.write (String.format (
                    "    <lower value=\"%.5f\" likelihood=\"%.5g\"/>\n",
                    result[0],
                    likelihood[0]
                ));
                out.write (String.format (
                    "    <upper value=\"%.5f\" likelihood=\"%.5g\"/>\n",
                    result[1],
                    likelihood[1]
                ));
                out.write ("  </omegaCI>\n");
            }
            // Output the SigmaCI data.
            if (sigmaCI != null && sigmaCI.hasRun ()) {
                Double[] result = sigmaCI.getResult ();
                Double[] likelihood = sigmaCI.getLikelihood ();
                out.write ("  <sigmaCI>\n");
                out.write (String.format (
                    "    <lower value=\"%.5f\" likelihood=\"%.5g\"/>\n",
                    result[0],
                    likelihood[0]
                ));
                out.write (String.format (
                    "    <upper value=\"%.5f\" likelihood=\"%.5g\"/>\n",
                    result[1],
                    likelihood[1]
                ));
                out.write ("  </sigmaCI>\n");
            }
            // Output the Demarcation data.
            if (demarcation != null && demarcation.hasRun ()) {
                ArrayList<ArrayList<String>> ecotypes =
                    demarcation.getEcotypes ();
                String precision = "";
                switch (demarcation.getPrecision ()) {
                    case Demarcation.DEMARCATION_PRECISION_FINE_SCALE:
                        precision = "fine-scale";
                        break;
                    case Demarcation.DEMARCATION_PRECISION_COARSE_SCALE:
                        precision = "coarse-scale";
                        break;
                }
                out.write ("  <demarcation>\n");
                out.write ("    <precision value=\"" + precision + "\"/>\n");
                out.write (String.format (
                    "    <ecotypes size=\"%d\">\n",
                    ecotypes.size ()
                ));
                for (int i = 0; i < ecotypes.size (); i ++) {
                    ArrayList<String> ecotype = ecotypes.get (i);
                    out.write (String.format (
                        "      <ecotype number=\"%d\" size=\"%d\">\n",
                        (i + 1), ecotype.size ()
                    ));
                    for (int j = 0; j < ecotype.size (); j ++) {
                        out.write (String.format (
                            "        <member name=\"%s\"/>\n",
                            ecotype.get (j)
                        ));
                    }
                    out.write ("      </ecotype>\n");
                }
                out.write ("    </ecotypes>\n");
                out.write ("  </demarcation>\n");
            }
            out.write ("</ecosim>\n");
            out.close ();
        }
        catch (IOException e) {
            e.printStackTrace ();
        }
    }

    /**
     *  Load the project file.
     *
     *  @param projectFile The project file to load.
     */
    public void load (File projectFile) {
        XMLHandler handler = new XMLHandler ();
        try {
            XMLReader xr = XMLReaderFactory.createXMLReader ();
            xr.setContentHandler (handler);
            xr.setErrorHandler (handler);
            xr.parse (new InputSource (new FileReader (projectFile)));
        }
        catch (Exception e) {
            e.printStackTrace ();
        }
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
     *  Get the length of the sequences.
     *
     *  @return The length of the sequences.
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
     *  Get the Tree object.
     *
     *  @return The Tree object.
     */
    public Tree getTree () {
        return tree;
    }

    /**
     *  Get the Binning object.
     *
     *  @return The Binning object.
     */
    public Binning getBinning () {
        return binning;
    }

    /**
     *  Get the ParamterEstimate object.
     *
     *  @return The ParamterEstimate object.
     */
    public ParameterEstimate getEstimate () {
        return estimate;
    }

    /**
     *  Get the Hillclimb object.
     *
     *  @return The Hillclimb object.
     */
    public Hillclimb getHillclimb () {
        return hillclimb;
    }

    /**
     *  Get the NpopConfidenceInterval object.
     *
     *  @return The NpopConfidenceInterval object.
     */
    public NpopConfidenceInterval getNpopCI () {
        return npopCI;
    }

    /**
     *  Get the OmegaConfidenceInterval object.
     *
     *  @return The OmegaConfidenceInterval object.
     */
    public OmegaConfidenceInterval getOmegaCI () {
        return omegaCI;
    }

    /**
     *  Get the SigmaConfidenceInterval object.
     *
     *  @return The SigmaConfidenceInterval object.
     */
    public SigmaConfidenceInterval getSigmaCI () {
        return sigmaCI;
    }

    /**
     *  Get the Demarcation object.
     *
     *  @return The Demarcation object.
     */
    public Demarcation getDemarcation () {
        return demarcation;
    }

    private MasterVariables masterVariables;
    private Execs execs;
    private Integer nu;
    private Integer length;
    private String outgroup;
    private Tree tree;
    private Binning binning;
    private ParameterEstimate estimate;
    private Hillclimb hillclimb;
    private NpopConfidenceInterval npopCI;
    private OmegaConfidenceInterval omegaCI;
    private SigmaConfidenceInterval sigmaCI;
    private Demarcation demarcation;

    /**
     *  Create a DefaultHandler for parsing the options file.
     */
    private class XMLHandler extends DefaultHandler  {

        /**
         *  The XMLHandler constructor.
         */
        public XMLHandler () {
            elements = new ArrayList<String> ();
            elements.add ("phylogeny");
            elements.add ("binning");
            elements.add ("estimate");
            elements.add ("hillclimb");
            elements.add ("npopCI");
            elements.add ("omegaCI");
            elements.add ("sigmaCI");
            elements.add ("demarcation");
            precision = null;
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
            if (localName.equals ("ecosim") &&
                attrs.getValue (uri, "type").equals ("savefile")) {
                isProjectFile = true;
            }
            if (isProjectFile) {
                // Update the active element.
                if (elements.contains (localName)) {
                    activeElement = localName;
                }
                // Look for the criterion element.
                if (localName.equals ("criterion")) {
                    masterVariables.setCriterion (new Integer (
                        attrs.getValue (uri, "value")
                    ).intValue ());
                }
                // Look for the phylogeny element.
                if (localName.equals ("phylogeny")) {
                    nu = new Integer (attrs.getValue (uri, "size"));
                    length = new Integer (attrs.getValue (uri, "length"));
                }
                // Look for the elements within phylogeny.
                if (activeElement.equals ("phylogeny")) {
                    // Look for the outgroup element.
                    if (localName.equals ("outgroup")) {
                        outgroup = attrs.getValue (uri, "value");
                    }
                    // Look for the tree element.
                    if (localName.equals ("tree")) {
                        try {
                            tree = new Tree (
                                attrs.getValue (uri, "value")
                            );
                        }
                        catch (InvalidTreeException e) {
                            System.err.println (
                                "Invalid Newick formatted tree found."
                            );
                        }
                    }
                }
                // Look for elements within binning.
                if (activeElement.equals ("binning")) {
                    // Initialize the Binning object if need be.
                    if (binning == null) binning = new Binning ();
                    // Look for a bin to add to the Binning.
                    if (localName.equals ("bin")) {
                        binning.addBinLevel (new BinLevel (
                            new Double (attrs.getValue (uri, "crit")),
                            new Integer (attrs.getValue (uri, "value"))
                        ));
                    }
                }
                // Look for elements within estimate.
                if (activeElement.equals ("estimate")) {
                    // Initialize the ParameterEstimate object if need be.
                    if (estimate == null) {
                        if (binning != null) {
                            estimate = new ParameterEstimate (
                                nu, length, binning
                            );
                        }
                        else {
                            System.err.println (
                                "Error in project file: " +
                                "estimated value without Binning."
                            );
                            System.exit (1);
                        }
                    }
                    // Look for a result to add to the the ParameterEstimate.
                    if (localName.equals ("result")) {
                        estimate.setResult (new ParameterSet (
                            new Long (attrs.getValue (uri, "npop")),
                            new Double (attrs.getValue (uri, "omega")),
                            new Double (attrs.getValue (uri, "sigma")),
                            0.0d
                        ));
                    }
                    // Look for the slope and intercept of the omega line.
                    if (localName.equals ("omega")) {
                        estimate.setOmega (
                            new Double (attrs.getValue (uri, "slope")),
                            new Double (attrs.getValue (uri, "intercept"))
                        );
                    }
                    // Look for the slope and intercept of the sigma line.
                    if (localName.equals ("sigma")) {
                        estimate.setSigma (
                            new Double (attrs.getValue (uri, "slope")),
                            new Double (attrs.getValue (uri, "intercept"))
                        );
                    }
                }
                // Look for elements within hillclimb.
                if (activeElement.equals ("hillclimb")) {
                    // Initialize the Hillclimb object if need be.
                    if (hillclimb == null) {
                        if (binning != null && estimate != null) {
                            hillclimb = new Hillclimb (
                                masterVariables,
                                execs,
                                nu,
                                length,
                                binning,
                                estimate.getResult ()
                            );
                        }
                        else if (binning == null) {
                            System.err.println (
                                "Error in project file: " +
                                "hillclimb value without Binning."
                            );
                            System.exit (1);
                        }
                        else {
                            System.err.println (
                                "Error in project file: " +
                                "hillclimb value without ParameterEstimate."
                            );
                            System.exit (1);
                        }
                    }
                    // Look for a result to add to the Hillclimb object.
                    if (localName.equals ("result")) {
                        hillclimb.setResult (new ParameterSet (
                            new Long (attrs.getValue (uri, "npop")),
                            new Double (attrs.getValue (uri, "omega")),
                            new Double (attrs.getValue (uri, "sigma")),
                            new Double (attrs.getValue (uri, "likelihood"))
                        ));
                    }
                }
                // Look for elements within npopCI.
                if (activeElement.equals ("npopCI")) {
                    // Initialize NpopConfidenceInterval if need be.
                    if (npopCI == null) {
                        if (binning != null && hillclimb != null) {
                            npopCI = new NpopConfidenceInterval (
                                masterVariables,
                                execs,
                                nu,
                                length,
                                binning,
                                hillclimb.getResult ()
                            );
                        }
                        else if (binning == null) {
                            System.err.println (
                                "Error in project file: " +
                                "npopCI value without Binning."
                            );
                            System.exit (1);
                        }
                        else {
                            System.err.println (
                                "Error in project file: " +
                                "npopCI value without Hillclimb."
                            );
                            System.exit (1);
                        }
                    }
                    // Look for the lower bound of the confidence interval.
                    if (localName.equals ("lower")) {
                        Long value = new Long (
                            attrs.getValue (uri, "value")
                        );
                        Double likelihood = new Double (
                            attrs.getValue (uri, "likelihood")
                        );
                        npopCI.setLowerResult (value, likelihood);
                    }
                    // Look for the upper bound of the confidence interval.
                    if (localName.equals ("upper")) {
                        Long value = new Long (
                            attrs.getValue (uri, "value")
                        );
                        Double likelihood = new Double (
                            attrs.getValue (uri, "likelihood")
                        );
                        npopCI.setUpperResult (value, likelihood);
                    }
                }
                // Look for elements within omegaCI.
                if (activeElement.equals ("omegaCI") && omegaCI != null) {
                    // Initialize OmegaConfidenceInterval if need be.
                    if (omegaCI == null) {
                        if (binning != null && hillclimb != null) {
                            omegaCI = new OmegaConfidenceInterval (
                                masterVariables,
                                execs,
                                nu,
                                length,
                                binning,
                                hillclimb.getResult ()
                            );
                        }
                        else if (binning == null) {
                            System.err.println (
                                "Error in project file: " +
                                "omegaCI value without Binning."
                            );
                            System.exit (1);
                        }
                        else {
                            System.err.println (
                                "Error in project file: " +
                                "omegaCI value without Hillclimb."
                            );
                            System.exit (1);
                        }
                    }
                    // Look for the lower bound of the confidence interval.
                    if (localName.equals ("lower")) {
                        Double value = new Double (
                            attrs.getValue (uri, "value")
                        );
                        Double likelihood = new Double (
                            attrs.getValue (uri, "likelihood")
                        );
                        omegaCI.setLowerResult (value, likelihood);
                    }
                    // Look for the upper bound of the confidence interval.
                    if (localName.equals ("upper")) {
                        Double value = new Double (
                            attrs.getValue (uri, "value")
                        );
                        Double likelihood = new Double (
                            attrs.getValue (uri, "likelihood")
                        );
                        omegaCI.setUpperResult (value, likelihood);
                    }
                }
                // Look for elements within sigmaCI.
                if (activeElement.equals ("sigmaCI") && sigmaCI != null) {
                    // Initialize SigmaConfidenceInterval if need be.
                    if (sigmaCI == null) {
                        if (binning != null && hillclimb != null) {
                            sigmaCI = new SigmaConfidenceInterval (
                                masterVariables,
                                execs,
                                nu,
                                length,
                                binning,
                                hillclimb.getResult ()
                            );
                        }
                        else if (binning == null) {
                            System.err.println (
                                "Error in project file: " +
                                "sigmaCI value without Binning."
                            );
                            System.exit (1);
                        }
                        else {
                            System.err.println (
                                "Error in project file: " +
                                "sigmaCI value without Hillclimb."
                            );
                            System.exit (1);
                        }
                    }
                    // Look for the lower bound of the confidence interval.
                    if (localName.equals ("lower")) {
                        Double value = new Double (
                            attrs.getValue (uri, "value")
                        );
                        Double likelihood = new Double (
                            attrs.getValue (uri, "likelihood")
                        );
                        sigmaCI.setLowerResult (value, likelihood);
                    }
                    // Look for the upper bound of the confidence interval.
                    if (localName.equals ("upper")) {
                        Double value = new Double (
                            attrs.getValue (uri, "value")
                        );
                        Double likelihood = new Double (
                            attrs.getValue (uri, "likelihood")
                        );
                        sigmaCI.setUpperResult (value, likelihood);
                    }
                }
                // Look for elements within demarcation.
                if (activeElement.equals ("demarcation")) {
                    if (localName.equals ("method")) {
                        switch (attrs.getValue (uri, "value")) {
                            case "monophyly":
                                method = Demarcation.DEMARCATION_METHOD_MONOPHYLY;
                                break;
                            case "paraphyly":
                                method = Demarcation.DEMARCATION_METHOD_PARAPHYLY;
                                break;
                        }
                    }
                    if (localName.equals ("precision")) {
                        switch (attrs.getValue (uri, "value")) {
                            case "fine-scale":
                                precision = Demarcation.DEMARCATION_PRECISION_FINE_SCALE;
                                break;
                            case "coarse-scale":
                                precision = Demarcation.DEMARCATION_PRECISION_COARSE_SCALE;
                                break;
                        }
                    }
                    if (localName.equals ("ecotypes")) {
                        Integer size = new Integer (
                            attrs.getValue (uri, "size")
                        );
                        ecotypes = new ArrayList<ArrayList<String>> (size);
                    }
                    if (localName.equals ("ecotype")) {
                        Integer size = new Integer (
                            attrs.getValue (uri, "size")
                        );
                        ecotypes.add (new ArrayList<String> (size));
                        ecotypeNumber = new Integer (
                            attrs.getValue (uri, "number")
                        );
                    }
                    if (localName.equals ("member") && ecotypes != null) {
                        ecotypes.get (ecotypeNumber - 1).add (
                            attrs.getValue (uri, "name")
                        );
                    }
                }
            }
        }

        /**
         *  Override the endElement method in DefaultHandler to help keep
         *  track of our location in the XML document.
         */
        public void endElement (String uri, String localName, String qName) {
            if (isProjectFile) {
                // Look for the end of the ecosim save file.
                if (localName.equals ("ecosim")) {
                    isProjectFile = false;
                }
                // Look for the end of the estimate element.
                if (localName.equals ("estimate") && estimate != null) {
                    estimate.setHasRun (true);
                }
                // Look for the end of the hillclimb element.
                if (localName.equals ("hillclimb") && hillclimb != null) {
                    hillclimb.setHasRun (true);
                }
                // Look for the end of the npopCI element.
                if (localName.equals ("npopCI") && npopCI != null) {
                    npopCI.setHasRun (true);
                }
                // Look for the end of the omegaCI element.
                if (localName.equals ("omegaCI") && omegaCI != null) {
                    omegaCI.setHasRun (true);
                }
                // Look for the end of the sigmaCI element.
                if (localName.equals ("sigmaCI") && sigmaCI != null) {
                    sigmaCI.setHasRun (true);
                }
                // Look for the end of the demarcation element.
                if (localName.equals ("demarcation") && demarcation != null) {
                    if (demarcation == null) {
                        if (outgroup != null && tree != null && hillclimb != null && precision != null) {
                            try {
                                demarcation = new Demarcation (
                                    masterVariables,
                                    execs,
                                    nu,
                                    length,
                                    outgroup,
                                    tree,
                                    hillclimb.getResult (),
                                    method,
                                    precision
                                );
                            }
                            catch (InvalidTreeException e) {
                                System.err.println (
                                    "Error in project file: invalid tree."
                                );
                            }
                        }
                        else if (outgroup == null) {
                            System.err.println (
                                "Error in project file: " +
                                "demarcation value without an outgroup."
                            );
                            System.exit (1);

                        }
                        else if (tree == null) {
                            System.err.println (
                                "Error in project file: " +
                                "demarcation value without a phylogeny."
                            );
                            System.exit (1);
                        }
                        else {
                            System.err.println (
                                "Error in project file: " +
                                "demarcation value without Hillclimb."
                            );
                            System.exit (1);
                        }
                    }
                    demarcation.setEcotypes (ecotypes);
                    demarcation.setHasRun (true);
                }
                // Update the active element.
                if (elements.contains (localName)) {
                    activeElement = "none";
                }
            }
        }

        private Integer method;
        private Integer precision;

        private String activeElement;

        private Integer ecotypeNumber;
        private ArrayList<ArrayList<String>> ecotypes;

        private boolean isProjectFile;

        private ArrayList<String> elements;
    }
}
