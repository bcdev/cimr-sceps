package org.esa.cimr.sceps.util;

import esa.opensf.osfi.CLP;

/**
 * Utility class for CIMR SCEPS
 *
 * @author olafd
 */
public class ScepsUtils {

    /**
     * dummy
     *
     */
    public static void doSomething() {
        System.out.println();
    }

    public static void checkCommandLineArgs(CLP clp) {
        if (clp.getConfFiles().isEmpty() || clp.getInputFiles().isEmpty() || clp.getOutputFiles().isEmpty()) {
            throw new IllegalArgumentException("Invalid arguments passed to wrapper executable. " +
                    "Make sure that global and local config, inputs, and outputs are specified.");
        }
    }
}
