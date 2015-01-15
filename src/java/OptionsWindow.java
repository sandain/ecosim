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

import java.awt.Dimension;
import java.awt.BorderLayout;
import java.awt.FlowLayout;
import java.awt.GridLayout;
import java.awt.event.WindowAdapter;
import java.awt.event.WindowEvent;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.awt.event.ItemEvent;
import java.awt.event.ItemListener;
import java.io.File;
import java.util.HashMap;
import javax.swing.Box;
import javax.swing.JButton;
import javax.swing.JComponent;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JOptionPane;
import javax.swing.JPanel;
import javax.swing.JTabbedPane;
import javax.swing.JTextField;
import javax.swing.event.DocumentEvent;
import javax.swing.event.DocumentListener;
import javax.swing.text.BadLocationException;

/**
 *  Creates the Option Window.
 *
 *  @author Jason Wood
 */
public class OptionsWindow extends JFrame {

    /**
     *  Constructor for OptionsWindow.
     *
     *  @param masterVariables The MasterVariables.
     */

    public OptionsWindow(MasterVariables masterVariables) {
        this.masterVariables = masterVariables;
        options = masterVariables.getOptions();
        directoryOptions = new HashMap<String, Item>();
        helperPrograms = new HashMap<String, Item>();
        makeOptionsWindow();
        loadActionPerformed();
    }

    /**
     *  Make the Options Window.
     */
    private void makeOptionsWindow() {
        setTitle("Options");
        setMinimumSize(new Dimension(600, 400));
        setMaximumSize(new Dimension(600, 400));
        setLayout(new BorderLayout());
        // Create Tabbed Pane
        add(makeTabbedPane(), BorderLayout.CENTER);
        // Create Button Pane
        add(makeButtonPane(), BorderLayout.SOUTH);
        // Listen for the window to close.
        addWindowListener(new WindowAdapter() {
            public void windowClosing(WindowEvent evt) {
                exitActionPerformed();
            }
        });
        pack();
        setVisible(true);
    }

    /**
     *  Make the Tabbed pane.
     *
     *  @return JComponent - The tabbed pane.
     */
    private JComponent makeTabbedPane() {
        JTabbedPane tabbedPane = new JTabbedPane();
        tabbedPane.addTab(Options.DIRECTORY_OPTIONS, makeDirectoryOptionsTab());
        tabbedPane.addTab(Options.HELPER_PROGRAMS, makeHelperProgramsTab());
        return tabbedPane;
    }

    /**
     *  Make the Button pane.
     *
     *  @return JComponent - The button pane.
     */
    private JComponent makeButtonPane() {
        JPanel buttonPane = new JPanel();
        buttonPane.setLayout(new GridLayout(1, 3));
        buttonPane.add(new Item(Item.BUTTON, exitButton));
        buttonPane.add(new Item(Item.BUTTON, cancelButton));
        buttonPane.add(new Item(Item.BUTTON, saveButton));
        return buttonPane;
    }

    /**
     *  Make the Directory Options tab.
     *
     *  @return JComponent - The directory options tab.
     */
    private JComponent makeDirectoryOptionsTab() {
        JPanel panel = new JPanel(false);
        JPanel directoryOptionsPanel = new JPanel(false);
        JPanel textFieldPanel = new JPanel(false);
        // Add a text field for each directory option.
        String[] directoryOptionList = options.getDirectoryOptionList();
        textFieldPanel.setLayout(new GridLayout(directoryOptionList.length, 1));
        for (int i = 0; i < directoryOptionList.length; i ++) {
            Item item = new Item(Item.DIRECTORY_CHOOSER, directoryOptionList[i]);
            directoryOptions.put(directoryOptionList[i], item);
            textFieldPanel.add(item);
        }
        // Create the directory options panel.
        directoryOptionsPanel.setLayout(new BorderLayout());
        directoryOptionsPanel.add(textFieldPanel, BorderLayout.NORTH);
        // Create a border around the directory options.
        panel.setLayout(new BorderLayout());
        panel.add(Box.createRigidArea(new Dimension(0, 10)), BorderLayout.NORTH);
        panel.add(Box.createRigidArea(new Dimension(0, 10)), BorderLayout.SOUTH);
        panel.add(Box.createRigidArea(new Dimension(10, 0)), BorderLayout.EAST);
        panel.add(Box.createRigidArea(new Dimension(10, 0)), BorderLayout.WEST);
        panel.add(directoryOptionsPanel, BorderLayout.CENTER);
        return panel;
    }

    /**
     *  Make the Helper Programs tab.
     *
     *  @return JComponent - The helper programs tab.
     */
    private JComponent makeHelperProgramsTab() {
        JPanel panel = new JPanel(false);
        JPanel helperProgramsPanel = new JPanel(false);
        JPanel textFieldPanel = new JPanel(false);
        // Add a text field for each helper program.
        String[] helperProgramList = options.getHelperProgramList();
        textFieldPanel.setLayout(new GridLayout(helperProgramList.length, 1));
        for (int i = 0; i < helperProgramList.length; i ++) {
            Item item = new Item(Item.FILE_CHOOSER, helperProgramList[i]);
            helperPrograms.put(helperProgramList[i], item);
            textFieldPanel.add(item);
        }
        // Create the helper programs panel.
        helperProgramsPanel.setLayout(new BorderLayout());
        helperProgramsPanel.add(textFieldPanel, BorderLayout.NORTH);
        // Create a border around the helper programs.
        panel.setLayout(new BorderLayout());
        panel.add(Box.createRigidArea(new Dimension(0, 10)), BorderLayout.NORTH);
        panel.add(Box.createRigidArea(new Dimension(0, 10)), BorderLayout.SOUTH);
        panel.add(Box.createRigidArea(new Dimension(10, 0)), BorderLayout.EAST);
        panel.add(Box.createRigidArea(new Dimension(10, 0)), BorderLayout.WEST);
        panel.add(helperProgramsPanel, BorderLayout.CENTER);
        return panel;
    }

    /**
     *  Perform the action for the Button pushed.
     */
    private void buttonActionPerformed(String label) {
        if (label.equals(saveButton)) {
            saveActionPerformed();
        }
        else if (label.equals(cancelButton)) {
            loadActionPerformed();
        }
        else if (label.equals(exitButton)) {
            exitActionPerformed();
        }
    }

    /**
     *  Save the options.
     */
    private void saveActionPerformed() {
        String[] keys;
        // Save the directory options.
        keys = options.getDirectoryOptionList();
        for (int i = 0; i < keys.length; i ++) {
            options.setDirectoryOption(keys[i], directoryOptions.get(keys[i]).getValue());
        }
        // Save the helper programs.
        keys = options.getHelperProgramList();
        for (int i = 0; i < keys.length; i ++) {
            options.setHelperProgram(keys[i], helperPrograms.get(keys[i]).getValue());
        }
        OptionsIO.save(options);
        saved = true;
        changed = false;
    }

    /**
     *  Load the options file.
     */
    private void loadActionPerformed() {
        options = OptionsIO.load();
        String[] keys;
        // Load the directory options.
        keys = options.getDirectoryOptionList();
        for (int i = 0; i < keys.length; i ++) {
            directoryOptions.get(keys[i]).setValue(options.getDirectoryOption(keys[i]));
        }
        // Load the helper programs.
        keys = options.getHelperProgramList();
        for (int i = 0; i < keys.length; i ++) {
            helperPrograms.get(keys[i]).setValue(options.getHelperProgram(keys[i]));
        }
        // Update the options stored in the MasterVariables.
        masterVariables.setOptions(options);
        // Set the saved and changed flags.
        saved = false;
        changed = false;
    }

    /**
     *  Perform the exit action.  If the user has not saved, verify that he/she
     *  really wants to exit.
     */
    private void exitActionPerformed() {
        if (changed && ! saved) {
            if (exitWithSaveDialog()) {
                saveActionPerformed();
            }
        }
        dispose();
    }

    /**
     *  Display a Confirmation Dialog, asking if the user wants to save their
     *  changes before exiting.
     *
     *  @return boolean - True if Yes, False if No.
     */
    private boolean exitWithSaveDialog() {
        boolean save = false;
        String title = "Save changes?";
        String message = "You have unsaved changes, would you like to save them now?";
        int n = JOptionPane.showConfirmDialog(this, message, title, JOptionPane.YES_NO_OPTION);
        if (n == JOptionPane.YES_OPTION) {
            save = true;
        }
        return save;
    }

    private class Item extends javax.swing.JPanel {

        private Item(int type, String label) {
            this.type = type;
            this.label = label;
            this.value = "";
            initialize();
        }

        private Item(int type, String label, String value) {
            this.type = type;
            this.label = label;
            this.value = value;
            initialize();
        }

        public void setValue(String value) {
            this.value = value;
            switch(type) {
                case BUTTON:
                    break;
                case TEXT:
                case FILE_CHOOSER:
                case DIRECTORY_CHOOSER:
                    valueTextField.setText(value);
                    break;
                default:
                    break;
            }
        }

        public String getValue() {
            return value;
        }

        public static final int BUTTON = 1;
        public static final int TEXT = 2;
        public static final int FILE_CHOOSER = 3;
        public static final int DIRECTORY_CHOOSER = 4;

        private void initialize() {
            JLabel jLabel;
            switch(type) {
                case BUTTON:
                    button = new JButton();
                    button.setText(label);
                    button.addActionListener (
                        new ActionListener() {
                            public void actionPerformed(ActionEvent evt) {
                                buttonActionPerformed(label);
                            }
                        }
                    );
                    add(button);
                    break;
                case TEXT:
                    valueTextField = new JTextField(value, 35);
                    valueTextField.getDocument().addDocumentListener (
                        new DocumentListener() {
                            public void changedUpdate(DocumentEvent e) {
                                changed = true;
                            }
                            public void insertUpdate(DocumentEvent e) {
                                changed = true;
                            }
                            public void removeUpdate(DocumentEvent e) {
                                changed = true;
                            }
                        }
                    );
                    jLabel = new JLabel(label + ": ", JLabel.TRAILING);
                    jLabel.setLabelFor(valueTextField);
                    add(jLabel);
                    add(valueTextField);
                    break;
                case FILE_CHOOSER:
                    // Create the text field.
                    valueTextField = new JTextField(value, 35);
                    valueTextField.getDocument().addDocumentListener (
                        new DocumentListener() {
                            public void changedUpdate(DocumentEvent e) {
                                changed = true;
                            }
                            public void insertUpdate(DocumentEvent e) {
                                changed = true;
                            }
                            public void removeUpdate(DocumentEvent e) {
                                changed = true;
                            }
                        }
                    );
                    // Create the label for the text field.
                    jLabel = new JLabel(label + ": ", JLabel.TRAILING);
                    jLabel.setLabelFor(valueTextField);
                    // Create the button.
                    button = new JButton();
                    button.setText(fileChooserButton);
                    button.addActionListener (
                        new ActionListener() {
                            public void actionPerformed(ActionEvent evt) {
                                FileChooser fc = new FileChooser("file", new File(System.getProperty("user.dir")));
                                int returnVal = fc.showOpenDialog(button);
                                if (returnVal == FileChooser.APPROVE_OPTION) {
                                    value = fc.getSelectedFile().getPath();
                                    valueTextField.setText(value);
                                }
                            }
                        }
                    );
                    // Set the layout.
                    setLayout(new FlowLayout(FlowLayout.RIGHT));
                    // Add all of the components.
                    add(jLabel);
                    add(valueTextField);
                    add(button);
                    break;
                case DIRECTORY_CHOOSER:
                    // Create the text field.
                    valueTextField = new JTextField(value, 35);
                    valueTextField.getDocument().addDocumentListener (
                        new DocumentListener() {
                            public void changedUpdate(DocumentEvent e) {
                                changed = true;
                            }
                            public void insertUpdate(DocumentEvent e) {
                                changed = true;
                            }
                            public void removeUpdate(DocumentEvent e) {
                                changed = true;
                            }
                        }
                    );
                    // Create the label for the text field.
                    jLabel = new JLabel(label + ": ", JLabel.TRAILING);
                    jLabel.setLabelFor(valueTextField);
                    // Create the button.
                    button = new JButton();
                    button.setText(fileChooserButton);
                    button.addActionListener (
                        new ActionListener() {
                            public void actionPerformed(ActionEvent evt) {
                                FileChooser fc = new FileChooser("directory", new File(System.getProperty("user.dir")));
                                int returnVal = fc.showOpenDialog(button);
                                if (returnVal == FileChooser.APPROVE_OPTION) {
                                    value = fc.getSelectedFile().getPath();
                                    valueTextField.setText(value);
                                }
                            }
                        }
                    );
                    // Set the layout.
                    setLayout(new FlowLayout(FlowLayout.RIGHT));
                    // Add all of the components.
                    add(jLabel);
                    add(valueTextField);
                    add(button);
                    break;
                default:
                    break;
            }
        }

        private JButton button;
        private JTextField valueTextField;
        private JTextField labelTextField;

        private int type;
        private String label;
        private String value;
    }

    private MasterVariables masterVariables;
    private Options options;

    private HashMap<String, Item> directoryOptions;
    private HashMap<String, Item> helperPrograms;

    private boolean saved;
    private boolean changed;
    private final String saveButton = "Save";
    private final String cancelButton = "Cancel";
    private final String exitButton = "Exit";
    private final String fileChooserButton = " ";

}
