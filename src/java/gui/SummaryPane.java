/*
 *    Ecotype Simulation models the sequence diversity within a bacterial
 *    clade as the evolutionary result of net ecotype formation and periodic
 *    selection, yielding a certain number of ecotypes.
 *
 *    Copyright (C) 2015-2019  Jason M. Wood, Montana State University
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


package ecosim.gui;

import ecosim.Binning;
import ecosim.BinLevel;
import ecosim.MainVariables;
import ecosim.ParameterEstimate;
import ecosim.ParameterSet;
import ecosim.Summary;

import java.awt.BasicStroke;
import java.awt.BorderLayout;
import java.awt.Color;
import java.awt.Component;
import java.awt.FlowLayout;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.Insets;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.Observable;
import java.util.Observer;
import javax.swing.BorderFactory;
import javax.swing.JLayeredPane;
import javax.swing.JPanel;
import javax.swing.JTable;
import javax.swing.JTextArea;
import javax.swing.table.DefaultTableCellRenderer;
import javax.swing.table.JTableHeader;
import javax.swing.table.TableColumn;
import javax.swing.table.TableColumnModel;
import org.jfree.chart.ChartPanel;
import org.jfree.chart.JFreeChart;
import org.jfree.chart.annotations.XYTitleAnnotation;
import org.jfree.chart.axis.LogAxis;
import org.jfree.chart.axis.NumberAxis;
import org.jfree.chart.axis.NumberTickUnit;
import org.jfree.chart.block.LineBorder;
import org.jfree.chart.plot.XYPlot;
import org.jfree.chart.renderer.xy.XYLineAndShapeRenderer;
import org.jfree.chart.title.LegendTitle;
import org.jfree.data.xy.DefaultXYDataset;
import org.jfree.ui.RectangleAnchor;
import org.jfree.ui.RectangleEdge;
import org.jfree.ui.RectangleInsets;

/**
 *  This class is used to display summary information to the user.
 *
 *  @author Jason M. Wood
 *  @copyright GNU General Public License
 */
public class SummaryPane extends JPanel {

    public SummaryPane (Summary summary) {
        this.summary = summary;
        setBorder (BorderFactory.createEmptyBorder (15, 15, 15, 15));
        setLayout (new GridBagLayout ());
        // Setup the constraints for the GridBagLayout.
        GridBagConstraints northWest = new GridBagConstraints (
            0, 0,                           // gridx, gridy
            1, 1,                           // gridwidth, gridheight
            0.9D, 0.9D,                     // weightx, weighty
            GridBagConstraints.CENTER,      // anchor
            GridBagConstraints.BOTH,        // fill
            new Insets (0, 0, 0, 0),        // insets
            0, 0                            // ipadx, ipady
        );
        GridBagConstraints northEast = new GridBagConstraints (
            1, 0,                           // gridx, gridy
            1, 1,                           // gridwidth, gridheight
            0.1D, 0.5D,                     // weightx, weighty
            GridBagConstraints.NORTH,       // anchor
            GridBagConstraints.HORIZONTAL,  // fill
            new Insets (0, 0, 0, 0),        // insets
            0, 0                            // ipadx, ipady
        );
        GridBagConstraints south = new GridBagConstraints (
            0, 1,                           // gridx, gridy
            2, 1,                           // gridwidth, gridheight
            1.0D, 0.1D,                     // weightx, weighty
            GridBagConstraints.SOUTH,       // anchor
            GridBagConstraints.HORIZONTAL,  // fill
            new Insets (0, 0, 0, 0),        // insets
            0, 0                            // ipadx, ipady
        );
        // Add everything to the summary pane.
        add (makeBinningChart (), northWest);
        add (makeTextPane (), northEast);
        add (makeTablePane (), south);
    }

    /**
     *  Private method to build the binning chart.
     *
     *  @return A ChartPanel containing the binning chart.
     */
    private ChartPanel makeBinningChart () {
        final DefaultXYDataset binData = new DefaultXYDataset ();
        final NumberFormat nf = NumberFormat.getInstance ();
        final NumberAxis xAxis = new NumberAxis (
            "Sequence identity criterion"
        );
        nf.setMinimumFractionDigits (2);
        xAxis.setLowerBound (Binning.binLevels[0]);
        xAxis.setUpperBound (1.0D);
        xAxis.setTickUnit (new NumberTickUnit (0.05D, nf));
        LogAxis yAxis = new LogAxis ("Number of bins");
        yAxis.setBase (2.0D);
        yAxis.setNumberFormatOverride (NumberFormat.getInstance ());
        yAxis.setTickUnit (new NumberTickUnit (2.0D));
        yAxis.setMinorTickMarksVisible (true);
        yAxis.setAutoRangeMinimumSize (4.0D);
        yAxis.setSmallestValue (1.0D);
        XYLineAndShapeRenderer renderer = new XYLineAndShapeRenderer (
            true, false
        );
        for (int i = 0; i < seriesColors.length; i++) {
            renderer.setSeriesPaint (i, seriesColors[i]);
            renderer.setSeriesStroke (i, new BasicStroke (seriesStroke[i]));
        }
        XYPlot plot = new XYPlot (binData, xAxis, yAxis, renderer);
        JFreeChart binChart = new JFreeChart (
            null, JFreeChart.DEFAULT_TITLE_FONT, plot, false
        );
        binChart.setPadding (new RectangleInsets (0.0D, 0.0D, 0.0D, 10.0D));
        LegendTitle legend = new LegendTitle (plot);
        legend.setMargin (new RectangleInsets (1.0D, 1.0D, 1.0D, 1.0D));
        legend.setFrame (new LineBorder ());
        legend.setBackgroundPaint (Color.white);
        legend.setPosition (RectangleEdge.BOTTOM);
        plot.addAnnotation (new XYTitleAnnotation (
            0.001D, 0.999D, legend, RectangleAnchor.TOP_LEFT
        ));
        final ChartPanel pane = new ChartPanel (
            binChart, false, true, true, false, false
        );
        // Watch for changes to the Summary object.
        summary.addObserver (new Observer () {
            public void update (Observable o, Object obj) {
                Summary s = (Summary)obj;
                ParameterEstimate estimate = s.getEstimate ();
                ArrayList<BinLevel> bins = s.getBins ();
                if (bins.size () > 0) {
                    double[][] values = new double[2][bins.size ()];
                    Double low = 1.0d;
                    for (int i = 0; i < bins.size (); i++) {
                        BinLevel bin = bins.get (i);
                        values[0][i] = bin.getCrit ();
                        values[1][i] = bin.getLevel ();
                        if (values[0][i] < low) low = values[0][i];
                    }
                    binData.addSeries ("sequences", values);
                    xAxis.setLowerBound (low);
                    if (low > 0.95d - MainVariables.EPSILON) {
                        xAxis.setTickUnit (new NumberTickUnit (0.005D, nf));
                    }
                    else if (low > 0.90d - MainVariables.EPSILON) {
                        xAxis.setTickUnit (new NumberTickUnit (0.010D, nf));
                    }
                    else if (low > 0.80d - MainVariables.EPSILON) {
                        xAxis.setTickUnit (new NumberTickUnit (0.025D, nf));
                    }
                    if (estimate != null) {
                        double[][] omega = new double[2][bins.size ()];
                        double[][] sigma = new double[2][bins.size ()];
                        double[] omegaLine = estimate.getOmega ();
                        double[] sigmaLine = estimate.getSigma ();
                        for (int i = 0; i < bins.size (); i++) {
                            double crit = 1.0D - values[0][i];
                            double snp = s.getLength () * crit;
                            omega[0][i] = values[0][i];
                            sigma[0][i] = values[0][i];
                            omega[1][i] = Math.pow (
                                2.0D, snp * omegaLine[0] + omegaLine[1]
                            );
                            sigma[1][i] = Math.pow (
                                2.0D, snp * sigmaLine[0] + sigmaLine[1]
                            );
                        }
                        if (-1.0D * omegaLine[0] > MainVariables.EPSILON) {
                            binData.addSeries ("omega", omega);
                        }
                        if (-1.0D * sigmaLine[0] > MainVariables.EPSILON) {
                            binData.addSeries ("sigma", sigma);
                        }
                    }
                    // Repaint the summary pane.
                    pane.repaint ();
                }
            }
        });
        return pane;
    }

    /**
     *  Private method to build the text pane.
     *
     *  @return A JLayeredPane containing the text pane.
     */
    private JLayeredPane makeTextPane () {
        final String ls = System.getProperty ("line.separator");
        final String fmt =
            "Outgroup: %s" + ls +
            "Number: %,d" + ls +
            "Length: %,d" + ls +
            "Diversity: %.2f" + ls;

        final JLayeredPane pane = new JLayeredPane ();
        final JTextArea summaryTextArea = new JTextArea (String.format (
            fmt,
            summary.getOutgroup (),
            summary.getNu (),
            summary.getLength (),
            summary.getDiversity ()
        ));
        summaryTextArea.setBackground (getBackground ());
        pane.setBorder (BorderFactory.createTitledBorder ("Sequences"));
        pane.setLayout (new FlowLayout (0));
        pane.add (summaryTextArea);
        // Watch for changes to the Summary object.
        summary.addObserver (new Observer () {
            public void update (Observable o, Object obj) {
                Summary s = (Summary)obj;
                ParameterEstimate estimate = s.getEstimate ();
                summaryTextArea.setText (String.format (
                    fmt,
                    s.getOutgroup (),
                    s.getNu (),
                    s.getLength (),
                    s.getDiversity ()
                ));
                pane.repaint ();
            }
        });
        return pane;
    }

    /**
     *  Private method to build the table pane.
     *
     *  @return A JLayeredPane containing the table pane.
     */
    private JLayeredPane makeTablePane () {
        String[] columnNames = {
            "", "Estimate", "Hillclimbing", "Low", "High"
        };
        Integer[] columnWidths = {
            250, 60, 60, 60, 60
        };
        Object[][] rowData = {
            { "Number of putative ecotypes (npop)", null, null, null, null },
            { "Rate of ecotype formation (omega)", null, null, null, null },
            { "Rate of periodic selection (sigma)", null, null, null, null }
        };
        final JLayeredPane pane = new JLayeredPane ();
        final JTable table = new JTable (rowData, columnNames) {
            public boolean isCellEditable (int row, int column) {
                return false;
            }
        };
        final JTableHeader header = table.getTableHeader ();
        pane.setLayout (new BorderLayout ());
        header.setReorderingAllowed (false);
        TableColumnModel cm = table.getColumnModel ();
        for (int i = 0; i < columnWidths.length; i++) {
            TableColumn column = cm.getColumn (i);
            column.setMinWidth (columnWidths[i]);
            if (i == 0) {
                column.setCellRenderer (new DefaultTableCellRenderer () {
                    public Component getTableCellRendererComponent (
                        JTable table, Object value, boolean isSelected,
                        boolean hasFocus, int row, int column
                    ) {
                        Component cell = super.getTableCellRendererComponent (
                            table, value, isSelected, hasFocus, row, column
                        );
                        cell.setBackground (header.getBackground ());
                        cell.setForeground (header.getForeground ());
                        return cell;
                    }
                });
            }
        }
        pane.add (table.getTableHeader (), "North");
        pane.add (table, "Center");
        // Watch for changes to the Summary object.
        summary.addObserver (new Observer () {
            public void update (Observable o, Object obj) {
                Summary s = (Summary)obj;
                ParameterEstimate estimate = s.getEstimate ();
                ParameterSet hillclimbing = s.getHillclimbing ();
                ParameterSet[] ci = s.getConfidenceInterval ();
                if (estimate != null) {
                    ParameterSet e = estimate.getResult ();
                    table.setValueAt (e.getNpop (), 0, 1);
                    table.setValueAt (
                        String.format ("%.4f", e.getOmega ()), 1, 1
                    );
                    table.setValueAt (
                        String.format ("%.4f", e.getSigma ()), 2, 1
                    );
                }
                if (hillclimbing != null) {
                    table.setValueAt (hillclimbing.getNpop (), 0, 2);
                    table.setValueAt (
                        String.format ("%.4f", hillclimbing.getOmega ()), 1, 2
                    );
                    table.setValueAt (
                        String.format ("%.4f", hillclimbing.getSigma ()), 2, 2
                    );
                }
                if (ci[0].getNpop () != null) {
                    table.setValueAt (ci[0].getNpop (), 0, 3);
                    table.setValueAt (ci[1].getNpop (), 0, 4);
                }
                if (ci[1].getOmega () != null) {
                    table.setValueAt (
                        String.format ("%.4f", ci[0].getOmega ()), 1, 3
                    );

                    String fmt = "%.4f";
                    if (ci[1].getOmega () > 10.0D) {
                        fmt = "%.1f";
                    }
                    else if (ci[1].getOmega () > 1.0D) {
                        fmt = "%.2f";
                    }
                    table.setValueAt (
                        String.format (fmt, ci[1].getOmega ()), 1, 4
                    );
                }
                if (ci[1].getSigma () != null) {
                    table.setValueAt (
                        String.format ("%.4f", ci[0].getSigma ()), 2, 3
                    );
                    if (ci[1].getSigma () > 100.0D - MainVariables.EPSILON) {
                        table.setValueAt ("≥100", 2, 4);
                    }
                    else {
                        String fmt = "%.4f";
                        if (ci[1].getSigma () > 10.0D) {
                            fmt = "%.1f";
                        }
                        else if (ci[1].getSigma () > 1.0D) {
                            fmt = "%.2f";
                        }
                        table.setValueAt (
                            String.format (fmt, ci[1].getSigma ()), 2, 4
                        );
                    }
                }
                pane.repaint ();
            }
        });
        return pane;
    }

    private Summary summary;

    private static final Color[] seriesColors = {
        Color.RED, Color.BLUE, Color.GREEN
    };

    private static final float[] seriesStroke = {
        2.0F, 1.0F, 1.0F
    };

}
