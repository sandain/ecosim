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

import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Path;

/**
 *  Performs operations specific to the Linux Operating System.
 *
 *  @author Jason M. Wood
 *  @copyright GNU General Public License
 */

public class Linux implements OperatingSystem {

    @Override
    public String getBinaryExtension () {
        return "";
    }

    @Override
    public boolean verifyExecutable (Path path) {
        String osArch = System.getProperty ("os.arch").toLowerCase ();
        int[] magicBytes = new int[] { 127, 69, 76, 70 };
        int[] x86 = new int[] { 1 };
        int[] amd64 = new int[] { 2 };
        try {
            byte[] array = Files.readAllBytes (path);
            // Verify the magic bytes of the ELF header.
            if (array.length < magicBytes.length) {
                return false;
            }
            for (int i = 0; i < magicBytes.length; i ++) {
                if (magicBytes[i] != (array[i] & 0xff)) {
                    return false;
                }
            }
            // Verify the architecture.
            int archPointer = magicBytes.length;
            if (osArch.contains ("amd64") || osArch.contains ("x86_64")) {
                // Run only 64-bit applications on 64-bit Linux.
                if (array.length < archPointer + amd64.length) {
                    return false;
                }
                for (int i = 0; i < amd64.length; i ++) {
                    if (amd64[i] != (array[archPointer + i] & 0xff)) {
                        return false;
                    }
                }
            }
            else {
                // Run only 32-bit applications on 32-bit Linux.
                if (array.length < archPointer + x86.length) {
                    return false;
                }
                for (int i = 0; i < x86.length; i ++) {
                    if (x86[i] != (array[archPointer + i] & 0xff)) {
                        return false;
                    }
                }
            }
        }
        catch (IOException e) {
            System.err.println ("Unable to verify executable!");
            e.printStackTrace ();
        }
        return true;
    }

}
