package org.esa.cimr.sceps;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStreamReader;
import java.nio.charset.StandardCharsets;
import java.util.Arrays;

/**
 * Class containing the wrapper for the OSF Scene Generation Scientific Moduule.
 * This class shall:
 * - be made available as executable jar (see pom) to be called from openSF
 * (set as 'Executable' in corresponding openSF module General tab).
 * - properly transform openSF arguments/parameters to ScientificModule arguments/parameters
 * - do a system call to run the SGM Scientific Module in Matlab batch mode
 */
public class SceneGenerationModuleWrapper {
    // todo: implement

//    *** openSF Simulation --> <moduleFile1.m>, .. ,<moduleFileN.m>  :
//
//            - openSF session folder is <Workspce path>/simulations/<Simulation name>.yyyyMMddThhmmssxxxx
//            - global config file is <session folder>/GlobalConfiguration.xml, automatically copied from manual setting
//            - local config file is <session folder>/<filename>.xml, automatically copied from manual setting
//            --> these are passed as parameters to executable jar, so can we extract the openSF session folder from one of them
//            --> equivalent to D:\olaf\bc\cimr\SCEPSscd\SCEPSdata\sessions\L1_Scene_Forward_Model_Simulation_GeoCard2_Day1 (variable 'sfo') where
//                L1_Scene_Forward_Model_Simulation_GeoCard2_Day1 is the SESSION_ID defined in L1_Scene_Forward_Model_Simulation.m
//
//            - we need to let openSF know where the SCEPSscd software is located:
//                        --> suggestion: set in global configuration: <scepsScdRoot>scepsScdRoot</scepsScdRoot>
//
//                        - we need to identify the moduleFile name (without extension):
//                        --> suggestion: default could be same as local config file. If not, specify in local config: <moduleName>mymodule</moduleName>
//
//                Wrapper needs to do:
//                        - main method receives as input parameters:
//                        # global config file (see above)
//                      # local config file (see above)
//                      # input files 1..N
//                      # output files 1..M

    /**
     * The wrapper main method.
     *
     * @param args - Program arguments. These are:
     *             # full path of global config file
     *             # full path of local config file
     *             # full paths of input files 1..N
     *             # full paths of output files 1..M
     */
    public static void main(String[] args) {
        if (args.length < 2) {
            System.out.println("Too few arguments given - must be at least paths of " +
                    "global and local config file. Exiting.");
        } else {
            for (String arg : args) {
                System.out.println("arg = " + arg);
            }

            // set relevant paths:
            String scepsScdRoot = "/data/sceps/SCEPSscd"; // "<parsed from global config>";  TODO: add and test Java method for xml parsing
            String moduleName = "GeoInputs_Extract";      // "<default or parsed from local config>";  // here, GeoInputs_Extract or Forward_Model
            String devSCEPSpath = scepsScdRoot + File.separator + "SCEPScodes";  // make configurable in global config? If not, set as a constant.
            String dataSCEPSpath = scepsScdRoot + File.separator + "SCEPSdata";  // make configurable in global config? If not, set as a constant.
            String modulesParentName = devSCEPSpath + File.separator + "SceGenMod";  // make configurable in global config? If not, set as a constant.
            String modulePath = modulesParentName + File.separator + "Modules" + File.separator + moduleName;

            // set relevant parameters to match module name signature (see e.g. GeoInputs_Extract.m):
            String globalConfigFilename = "";  //  ... get filename without path...
            String localConfigFilename = "";  // ... get filename without path...

            String configurationParameters = globalConfigFilename + "," + localConfigFilename;
            String inputs = "";  // <openSF sessionFolder>
            String outputs = "";  // <openSF sessionFolder>

            String[] commands = {
                    "matlab",
                    "-batch",
                    "devSCEPSpath = '" + devSCEPSpath + "'; " +
                            "addpath '" + devSCEPSpath + "'; " +
                            "dataSCEPSpath = '" + dataSCEPSpath + "'; " +
                            "cd " + modulePath + "; " +
                            "addpath '" + modulePath + "'; " +
                            "moduleName(" + configurationParameters + ", " + inputs + ", " + outputs + ");"
            };

//            StringBuffer sb = new StringBuffer();
//            for (int i = 0; i < commands.length; i++) {
//                sb.append(commands[i]);
//            }
            String str = Arrays.toString(commands);
            System.out.println("Command sequence: " + str);

            try {
                Process process = Runtime.getRuntime().exec(commands);

                BufferedReader reader = new BufferedReader(new InputStreamReader(process.getInputStream(),
                        StandardCharsets.UTF_8));
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