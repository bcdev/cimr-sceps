package org.esa.cimr.sceps;

import java.io.BufferedReader;
import java.io.IOException;
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

//            String command = "matlab -batch pwd";
            // make sure that createDummyFile.m is in the same folder
            // as the executable jar containing this dummy wrapper,
            // and that the folder is on the Matlab path...
            // For example, /data/workspace/cimr-opensf/test/bin
            // on VM cimr-sceps-develop
            String command = "matlab -batch createDummyFile";

            try {
                Process process = Runtime.getRuntime().exec(command);

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
    }
}