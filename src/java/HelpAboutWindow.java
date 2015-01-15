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

import java.awt.BorderLayout;
import java.awt.Desktop;
import java.awt.Dimension;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import javax.swing.JEditorPane;
import javax.swing.JFrame;
import javax.swing.JScrollPane;
import javax.swing.event.HyperlinkEvent;
import javax.swing.event.HyperlinkListener;
import javax.swing.text.html.HTMLFrameHyperlinkEvent;
import javax.swing.text.html.HTMLDocument;

public class HelpAboutWindow extends JFrame implements Runnable {

    public HelpAboutWindow(MasterVariables masterVariables) {
        this.masterVariables = masterVariables;
        initialize();
    }

    public void run() {
    }

    /**
     *  Set the HelpAboutWindow visible, with the given pane selected.
     *
     *  @param type The type of pane to display.
     */
    public void setVisible(String type) {
        if (type.equals(userGuide)) {
            scrollPane.setViewportView(userGuidePane);
        }
        else if (type.equals(about)) {
            scrollPane.setViewportView(aboutPane);
        }
        setVisible(true);
    }

    /**
     *  Initialize this HelpAboutWindow.
     */
    private void initialize() {
        userGuidePane = makeIndexPane(userGuide);
        aboutPane = makeIndexPane(about);
        scrollPane = new JScrollPane();
        setTitle("Help / About");
        setMinimumSize(new Dimension(640, 480));
        setPreferredSize(new Dimension(1420, 800));
        setLayout(new BorderLayout());
        add(scrollPane, BorderLayout.CENTER);
        // Listen for the window to close.
        addWindowListener(new WindowAdapter() {
            public void windowClosing(WindowEvent evt) {
                exitActionPerformed();
            }
        });
        pack();
    }

    /**
     *  Make the Index page.
     *
     *  @param type The type of index to display.
     */
    private JEditorPane makeIndexPane(String type) {
        JEditorPane editorPane = new JEditorPane();
        String index = new String();
        index += "<!DOCTYPE html>\n";
        index += "<html lang=\"en\" xmlns=\"http://www.w3.org/1999/xhtml\">\n";
        index += "  <head>\n";
        index += "    <title>\n";
        index += "      Help / About\n";
        index += "    </title>\n";
        index += "  </head>\n";
        index += "  <frameset cols=\"20%,80%\">\n";
        index += "    <frame src=\"file:///" + masterVariables.getHelpDirectory() + "menu.html" + "\" name=\"menu\" />\n";
        if (type.equals(userGuide)) {
            index += "    <frame src=\"file:///" + masterVariables.getHelpDirectory() + "userguide.html" + "\" name=\"body\" />\n";
        }
        else if (type.equals(about)) {
            index += "    <frame src=\"file:///" + masterVariables.getHelpDirectory() + "about.html" + "\" name=\"body\" />\n";
        }
        index += "  </frameset>\n";
        index += "</html>\n";
        try {
            File temp = File.createTempFile("index", ".html", new File(masterVariables.getWorkingDirectory()));
            BufferedWriter out = new BufferedWriter(new FileWriter(temp));
            out.write(index);
            out.close();
            editorPane.setPage(temp.toURI().toURL());
            temp.deleteOnExit();
        }
        catch (Exception e) {
            e.printStackTrace();
        }
        editorPane.setContentType("text/html;charset=UTF-8");
        editorPane.setEditable(false);
        editorPane.addHyperlinkListener(new Hyperactive());
        return editorPane;
    }

    /**
     *  Perform the exit action.
     */
    private void exitActionPerformed() {
        dispose();
    }

    /**
     *  Availiable pages.
     */
    public static final String userGuide = "User Guide";
    public static final String about = "About";
    public static final String version = "Version";
    public static final String license = "License";

    /**
     *  Private variables.
     */
    private MasterVariables masterVariables;
    private JEditorPane userGuidePane;
    private JEditorPane aboutPane;
    private JScrollPane scrollPane;

    /**
     *  Private class Hyperactive listens for clicked hyperlinks.
     */
    private class Hyperactive implements HyperlinkListener {
        public void hyperlinkUpdate(HyperlinkEvent evt) {
            if (evt.getEventType() == HyperlinkEvent.EventType.ACTIVATED) {
                JEditorPane pane = (JEditorPane) evt.getSource();
                try {
                    // Open external links in the default browser.
                    if (evt.getURL().toString().startsWith("http")) {
                        Desktop d = Desktop.getDesktop();
                        d.browse(evt.getURL().toURI());
                    }
                    // Open mailto links in the default mailer.
                    else if (evt.getURL().toString().startsWith("mailto")) {
                        Desktop d = Desktop.getDesktop();
                        d.mail(evt.getURL().toURI());
                    }
                    // Open internal links within the proper frame.
                    else if (evt instanceof HTMLFrameHyperlinkEvent) {
                        HTMLFrameHyperlinkEvent  frameEvt = (HTMLFrameHyperlinkEvent) evt;
                        HTMLDocument doc = (HTMLDocument) pane.getDocument();
                        doc.processHTMLFrameHyperlinkEvent(frameEvt);
                    }
                }
                catch (Exception e) {
                    e.printStackTrace();
                }
            }
        }
    }
}
