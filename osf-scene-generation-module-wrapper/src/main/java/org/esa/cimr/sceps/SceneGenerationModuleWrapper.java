package org.esa.cimr.sceps;

import org.apache.commons.io.FilenameUtils;
import org.w3c.dom.Document;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStreamReader;
import java.nio.charset.StandardCharsets;
import java.util.Arrays;

import static org.esa.cimr.sceps.ScepsConstants.SCEPS_MODULE_NAME_CONFIG_ITEM_NAME;
import static org.esa.cimr.sceps.ScepsConstants.SCEPS_SCD_ROOT_CONFIG_ITEM_NAME;

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
     * @param args - Program arguments. Two single strings with comma separated entries. These are:
     *             # first string: <full path of global config file>,<full path of global config file>
     *             # second string: <full paths of input files 1..N>,<full paths of output files 1..M>
     *             # Optionally, a 'simulation' mode can be set: last arg is 'simulation=true'
     */
    public static void main(String[] args) {
        if (args.length < 1) {
            System.out.println("Wrong number of arguments given - must be one comma separated string containing " +
                    "global and local config file, all input and output files. Exiting.");
            System.exit(1);
        } else {
            final String[] argConfigItems = args[0].split(",");

            // 'simulation' mode for testing. Omits execution of the Matlab batch command.
            final boolean simulation = args[args.length-1].equals("simulation=true");

            if (simulation) {
                for (String arg : argConfigItems) {
                    System.out.println("argItem = " + arg);
                }
            }

            // set relevant paths:
            String scepsScdRoot;
            String moduleName;
            final String globalConfigXmlPath = argConfigItems[0];
            final String localConfigXmlPath = argConfigItems[1];
            try {
                final Document globalConfigDoc = ScepsConfig.readXMLDocumentFromFile(globalConfigXmlPath);
                scepsScdRoot = ScepsConfig.getDocumentElementTextItemByName(globalConfigDoc,
                        ScepsConstants.SCEPS_CONFIG_ELEMENTS_TAG_NAME, SCEPS_SCD_ROOT_CONFIG_ITEM_NAME);
                // todo: add this to global config:
                // <parameter description="text" name="sceps_scd_root" type="STRING">/data/sceps/SCEPSscd</parameter>

                final Document localConfigDoc = ScepsConfig.readXMLDocumentFromFile(localConfigXmlPath);
                moduleName = ScepsConfig.getDocumentElementTextItemByName(localConfigDoc,
                        ScepsConstants.SCEPS_CONFIG_ELEMENTS_TAG_NAME, SCEPS_MODULE_NAME_CONFIG_ITEM_NAME);
                // todo: add this to GeoInputs_Extractlocal config:
                // <parameter description="text" name="module_name" type="STRING">GeoInputs_Extract</parameter>
                // todo: add this to Forward_Model local config:
                // <parameter description="text" name="module_name" type="STRING">Forward_Model</parameter>
                if (moduleName == null) {
                    moduleName = FilenameUtils.removeExtension((new File(localConfigXmlPath)).getName());
                    // strip extension '_Local_Configuration':
                    int index = moduleName.indexOf("_Local_Configuration");
                    moduleName = moduleName.substring(0, index);
                }

            } catch (Exception e) {
                // todo
                throw new RuntimeException(e);
            }

            String devSCEPSpath = scepsScdRoot + File.separator + ScepsConstants.SCEPS_CODES_FOLDER_NAME;
            String dataSCEPSpath = scepsScdRoot + File.separator + ScepsConstants.SCEPS_DATA_FOLDER_NAME;
            String modulesParentName = devSCEPSpath + File.separator +
                    ScepsConstants.SCENE_GENERATION_MODULE_FOLDER_NAME + File.separator +
                    ScepsConstants.SCENE_GENERATION_MODULE_MODULES_SUBFOLDER_NAME;

            // set relevant parameters to match module name signature (see e.g. GeoInputs_Extract.m):
            final File globalConfigXmlFile = new File(globalConfigXmlPath);

            String configurationParameters = globalConfigXmlPath + "," + localConfigXmlPath;
            String inputs = globalConfigXmlFile.getParent();  // everything is in the <openSF sessionFolder>
            String outputs = globalConfigXmlFile.getParent();  // same for outputs

            String[] commands = {
                    "matlab",
                    "-batch",
                    "devSCEPSpath = '" + devSCEPSpath + "'; " +
                            "addpath '" + devSCEPSpath + "'; " +
                            "global E2E_HOME; E2E_HOME = '" + dataSCEPSpath + "'; " +
                            "global LOG; LOG = Logger(); " +
                            "cd " + modulesParentName + "; " +
                            "addpath '" + modulesParentName + "'; " +
                            moduleName + "('" + configurationParameters + "','" + inputs + "','" + outputs + "');"
            };

            String str = Arrays.toString(commands);
            System.out.println("Command sequence: " + str);

            if (!simulation) {
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
}