package org.esa.cimr.sceps;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;

/**
 * Class containing just a main method as entry for first simple wrapper tests.
 * This class shall
 * - be made available as executable jar (see pom) to be called from openSF
 * (set as 'Executable' in corresponding openSF module General tab).
 * - properly transform openSF arguments/parameters to ScientificModule arguments/parameters
 * - do a system call to run a simple Matlab 'Scientific Module'
 */
public class SimpleScientificMatlabModuleWrapper {
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

//            String command = "matlab -batch pwd";
            // make sure that createDummyFile.m is in the same folder
            // as the executable jar containing this dummy wrapper,
            // and that the folder is on the Matlab path...
            // For example, /data/workspace/cimr-opensf/test/bin
            // on VM cimr-sceps-develop
//            String command = "matlab -batch createDummyFile";  // this works
//            String command = "matlab -batch \"createDummyFile\"";  // this does not work
//            String[] commands = {"matlab", "-batch", "createDummyFile"};  // this works
//            String[] commands = {"matlab", "-batch", "\"createDummyFile\""};   // this does not work
            String[] commands = {"matlab", "-batch", "mypath = '/home/olaf'; createDummyFile"};  // this works

            try {
                Process process = Runtime.getRuntime().exec(commands);

                BufferedReader reader = new BufferedReader(
                        new InputStreamReader(process.getInputStream()));
                String line;
                while ((line = reader.readLine()) != null) {
                    System.out.println(line);
                }

                reader.close();

            } catch (IOException e) {
                e.printStackTrace();
            }
        }

        // 1. read an xml parameters input file, such as Forward_Model_Local_Configuration_geocard2_day1.xml
        // 2. validate parameters and write into Java variables/arrays
        // 3. take these Java variables/arrays and write into xml parameters output file
        // 4. System call to ScepsSimpleScientificModule which takes this xml parameters output file as input parameter
        // 5. in ScepsSimpleScientificModule, extract again the parameters, do some fake action, and write a
        //      'scientific result' output file, e.g. as netcdf

    }
}
