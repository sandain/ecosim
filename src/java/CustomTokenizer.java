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

import java.util.ArrayList;
import java.util.List;
import java.util.Stack;
import java.text.StringCharacterIterator;
import java.util.Iterator;

public class CustomTokenizer {

    /**
     *  CustomTokenizer
     *
     *  @param s String to tokenize.
     */
    public CustomTokenizer(String s) {
        StringCharacterIterator iter = new StringCharacterIterator(s);
        tokens = new ArrayList<String>();
        for(char c = iter.first(); c != StringCharacterIterator.DONE; c = iter.next()) {
            String st = "";
            if (c == '(') {
                st = "(";
                Stack<String> stack = new Stack<String>();
                stack.push("(");
                while (!stack.empty()) {
                    c = iter.next();
                    if (c == '(')
                        stack.push("(");
                    else if (c == ')')
                        stack.pop();
                    st += c;
                }
            }
            else if (c != ',') {
                while ((c != ',') && (c != StringCharacterIterator.DONE)) {
                    st += c;
                    c = iter.next();
                }
            }
            tokens.add(st);
        }
    }
    
    /**
     *  Get the tokens.
     *
     *  @return List<String> of tokens.
     */
    public List<String> getTokens() {
        return tokens;
    }

    private List<String> tokens;
}
