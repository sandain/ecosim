/*
 *    Ecotype Simulation models the sequence diversity within a bacterial
 *    clade as the evolutionary result of net ecotype formation and periodic
 *    selection, yielding a certain number of ecotypes.
 *
 *    Copyright (C) 2018  Jason M. Wood, Montana State University
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


package ecosim.os;

import ecosim.api.OperatingSystem;

import java.nio.file.Path;

/**
 *  Performs operations specific to the Windows Operating System.
 *
 *  @author Jason M. Wood
 *  @copyright GNU General Public License
 */

public class Windows implements OperatingSystem {

    @Override
    public String getBinaryExtension () {
        return ".exe";
    }

    @Override
    public boolean verifyExecutable (Path path) {
        return true;
    }

}
