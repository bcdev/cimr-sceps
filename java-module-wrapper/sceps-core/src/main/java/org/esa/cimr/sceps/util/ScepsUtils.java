package org.esa.cimr.sceps.util;

import esa.opensf.osfi.CLP;

import java.util.List;

/**
 * Utility class for CIMR SCEPS
 *
 * @author olafd
 */
public class ScepsUtils {

    /**
     * Checks if required OSFI command line arguments exist.
     *
     * @param clp - OSFI command line parser object.
     */
    public static void checkCommandLineArgs(CLP clp) {
        if (clp.getConfFiles().isEmpty()) {
            throw new IllegalArgumentException("Invalid arguments passed to wrapper executable. " +
                    "Make sure that at least global and local config are specified.");
        }
    }

    public static String clpInputsJava2Matlab(List<String> clpInputFiles) {
        String clpMatlabString = "";
        for (int i = 0; i < clpInputFiles.size() - 1; i++) {
            clpMatlabString = clpMatlabString.concat(clpInputFiles.get(i));
            clpMatlabString = clpMatlabString.concat(",");
        }
        clpMatlabString = clpMatlabString.concat(clpInputFiles.get(clpInputFiles.size() - 1));

        return clpMatlabString;
    }
}
