package org.esa.cimr.sceps.util;

import esa.opensf.osfi.CLP;

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
}
