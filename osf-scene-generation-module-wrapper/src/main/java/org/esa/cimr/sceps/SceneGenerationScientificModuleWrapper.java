package org.esa.cimr.sceps;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;

/**
 * Class containing the wrapper for the OSF Scene Generation Scientific Moduule.
 * This class shall:
 * - be made available as executable jar (see pom) to be called from openSF
 * (set as 'Executable' in corresponding openSF module General tab).
 * - properly transform openSF arguments/parameters to ScientificModule arguments/parameters
 * - do a system call to run the SGM Scientific Module in Matlab batch mode
 */
public class SceneGenerationScientificModuleWrapper {
    // todo: implement

    /**
     * The wrapper main method.
     *
     * @param args - Program arguments
     */
    public static void main(String[] args) {
        if (args.length == 0) {
            System.out.println("No arguments given - exiting.");
        } else {
            for (String arg : args) {
                System.out.println("arg = " + arg);
            }

            // make sure that createDummyFile.m is in the same folder
            // as the executable jar containing this dummy wrapper,
            // and that the folder is on the Matlab path...
            // For example, /data/workspace/cimr-opensf/test/bin
            // on VM cimr-sceps-develop

            String devSCEPSpath = "/data/sceps/SCEPSscd/SCEPScodes";
            String dataSCEPSpath = "/data/sceps/SCEPSscd/SCEPSdata";

            // String[] commands = {"/bin/bash", "-c", "for i in {1..10000} ; do echo Hello STDERR $i 1>&2 ; echo Hello STDOUT $i; done"};
            // String[] commands = {"matlab", "-batch", "\"x1='bla', x2='blubb'\""};
            String[] commands = {
                    "matlab",
                    "-batch",
                    "\"devSCEPSpath = '" + devSCEPSpath + "', " +
                            "dataSCEPSpath = '" + dataSCEPSpath + "', " +
                            "startup_matlab_SCEPScodes(devSCEPSpath, dataSCEPSpath)\""};
            try {
                Process process = Runtime.getRuntime().exec(commands);

                InputStream inputStream = process.getInputStream();
                InputStreamReader isr = new InputStreamReader(inputStream, "UTF-8");
                BufferedReader reader = new BufferedReader(isr);
                String line;
                while ((line = reader.readLine()) != null) {
                    System.out.println(line);
                }

                reader.close();

            } catch (IOException e) {
                e.printStackTrace();
            }
        }
    }
}