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
        String osArch = System.getProperty ("os.arch").toLowerCase ();
        int[] dosMagicBytes = new int[] { 77, 90 };
        int[] prMagicBytes = new int[] { 80, 69, 0, 0 };
        int[] x86 = new int[] { 76, 1 };
        int[] amd64 = new int[] { 100, 134 };
        try {
            byte[] array = Files.readAllBytes (path);
            // Verify the magic bytes of the DOS header.
            if (array.length < dosMagicBytes.length) {
                return false;
            }
            for (int i = 0; i < dosMagicBytes.length; i ++) {
                if (dosMagicBytes[i] != (array[i] & 0xff)) {
                    return false;
                }
            }
            // Verify the magic bytes of the PR header.
            if (array.length < 60) {
                return false;
            }
            int prPointer = ((array[60] & 0xff));
            if (array.length < prPointer + prMagicBytes.length) {
                return false;
            }
            for (int i = 0; i < prMagicBytes.length; i ++) {
                if (prMagicBytes[i] != (array[prPointer + i] & 0xff)) {
                    return false;
                }
            }
            // Verify the architecture.
            int archPointer = prPointer + prMagicBytes.length;
            if (array.length < archPointer + 2) {
                return false;
            }
            if (osArch.contains ("amd64") || osArch.contains ("x86_64")) {
                boolean x86Bool = true;
                boolean amd64Bool = true;
                // Run either 64-bit or 32-bit applications on 64-bit Windows.
                for (int i = 0; i < amd64.length; i ++) {
                    if (amd64[i] != (array[archPointer + i] & 0xff)) {
                        amd64Bool = false;
                    }
                }
                for (int i = 0; i < x86.length; i ++) {
                    if (x86[i] != (array[archPointer + i] & 0xff)) {
                        x86Bool = false;
                    }
                }
                if (amd64Bool == false && x86Bool == false) return false;
            }
            else {
                // Run only 32-bit applications on 32-bit Windows.
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
