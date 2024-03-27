package org.esa.cimr.sceps;

import esa.opensf.osfi.CLP;
import esa.opensf.osfi.Logger;
import esa.opensf.osfi.ParamReader;
import org.apache.commons.io.FilenameUtils;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStreamReader;
import java.nio.charset.StandardCharsets;
import java.util.Arrays;

import static org.esa.cimr.sceps.ScepsConstants.SCEPS_DEVALGO_L2_L2GRID_CONFIG_ITEM_NAME;
import static org.esa.cimr.sceps.ScepsConstants.SCEPS_SCD_ROOT_CONFIG_ITEM_NAME;

/**
 * Class containing the wrapper for the Devalgo L2 Seaice Concentration Moduule.
 * This class shall:
 * - be made available as executable jar (see pom) to be called from openSF
 * (set as 'Executable' in corresponding openSF module General tab).
 * - properly transform openSF arguments/parameters to ScientificModule arguments/parameters
 * - do a system call to run the Seaice Concentration Scientific Module in Python
 */
public class DevalgoL2SeaiceConcentrationModuleWrapper {

//    *** openSF Simulation --> <moduleFile1.m>, .. ,<moduleFileN.m>  :
//
//            - openSF session folder is <Workspce path>/simulations/<Simulation name>.yyyyMMddThhmmssxxxx
//            - global config file is <session folder>/GlobalConfiguration.xml, automatically copied from manual setting
//            - local config file is <session folder>/<filename>.xml, automatically copied from manual setting
//            --> these are passed as parameters to executable jar, so can we extract the openSF session folder from one of them
//
//            - we need to let openSF know where the SCEPSscd software is located:
//                        --> suggestion: set in global configuration: <scepsScdRoot>scepsScdRoot</scepsScdRoot>
//
//                        - we need to identify the moduleFile name (without extension):
//                        --> suggestion: default could be same as local config file. If not, specify in local config: <moduleName>mymodule</moduleName>
//
//                Wrapper needs to do:
//                        - main method receives as input parameters:
//                        # global config file
//                        # local config file
//                        # input files 1..N
//                        # output files 1..M

//    java -jar /data/workspace/cimr-sceps/devalgo-l2-module-module-wrapper/target/devalgo-l2-module-wrapper-1.0-SNAPSHOT.jar
//   /data/workspace/cimr-opensf/simulations/Global_Configuration.xml,
//   /data/workspace/cimr-opensf/simulations/SCEPS_DevalgoL2SeaiceConcentration_test/Devalgo_L2_Seaice_Concentration_Local_Configuration.xml
//   /data/workspace/cimr-opensf/simulations/SCEPS_DevalgoL2SeaiceConcentration_test/Devalgo_L2_Seaice_Concentration_Output

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
            Logger.error("Wrong number of arguments given - must be one comma separated string containing " +
                    "global and local config file, all input and output files. Exiting.");
            System.exit(1);
        } else {
            // 'simulation' mode for testing. Omits execution of the Python batch command.
            final boolean simulation = args[args.length - 1].equals("simulation=true");

            // get program arguments from OSFI command line parser:
            final CLP clp = new CLP(args);

            final String globalConfigXmlPath = clp.getConfFiles().get(0);
            final String localConfigXmlPath = clp.getConfFiles().get(1);
            final String inputL1bPath = clp.getInputFiles().get(0);
            final String outputL2Dir = clp.getOutputFiles().get(0);

            Logger.info("globalConfigXmlPath: " + globalConfigXmlPath);
            Logger.info("localConfigXmlPath: " + localConfigXmlPath);
            Logger.info("inputL1bPath: " + inputL1bPath);
            Logger.info("outputL2Dir: " + outputL2Dir);

            // set relevant paths:
            String scepsScdRoot;
            String moduleName;
            String pythonScriptName;
            String l2Grid;
            try {
                // this was added to global config:
                // <parameter description="text" name="sceps_scd_root" type="STRING">/data/sceps/SCEPSscd</parameter>
                final ParamReader globalParamReader = new ParamReader(globalConfigXmlPath);
                scepsScdRoot = globalParamReader.getParameter(SCEPS_SCD_ROOT_CONFIG_ITEM_NAME).getStringValue();

                moduleName = FilenameUtils.removeExtension((new File(localConfigXmlPath)).getName());
                // strip extension '_Local_Configuration':
                if (moduleName.contains("_Local_Configuration")) {
                    int index = moduleName.indexOf("_Local_Configuration");
                    moduleName = moduleName.substring(0, index);
                }
                pythonScriptName = moduleName.toLowerCase() + ".py";

                final ParamReader localParamReader = new ParamReader(localConfigXmlPath);
                l2Grid = localParamReader.getParameter(SCEPS_DEVALGO_L2_L2GRID_CONFIG_ITEM_NAME).getStringValue();
            } catch (Exception e) {
                // todo
                throw new RuntimeException(e);
            }

            String devSCEPSpath = scepsScdRoot + File.separator + ScepsConstants.SCEPS_CODES_FOLDER_NAME;
            String modulesParentName = devSCEPSpath + File.separator +
                    ScepsConstants.DEVALGO_L2_MODULE_FOLDER_NAME + File.separator +
                    ScepsConstants.DEVALGO_L2_MODULE_MODULES_SUBFOLDER_NAME;

            String l2GridString = "";
            if (l2Grid != null && !l2Grid.isEmpty()) {
                l2GridString = " -g " + l2Grid;
            }

            final String pyArgs =
                    pythonScriptName + " -i " + inputL1bPath + " -o " + outputL2Dir + l2GridString;

            // test sequence of remote commands on Windows:
            /*
            final String chdirCmd = "cd " + modulesParentName + " && ";
            final String mkdirCmd = "mkdir " + outputL2Dir + " && ";
            final String condaCmd = "D:\\olaf\\miniconda3_py39\\condabin\\conda.bat activate base && ";
            final String pyCmd = "D:\\olaf\\miniconda3_py39\\python.exe ";

           String[] commands = {
                    "cmd",
                    "/c",
//                    "start /b dir && ping localhost && echo end"
//                    "start cmd.exe /K \"dir && ping localhost && echo end\""  // keeps cmd open
//                    "start cmd.exe /C" +
//                            "\"dir && ping localhost && echo end\""  // closes cmd after commands
//                    "start cmd.exe /K " + chdirCmd + mkdirCmd + pyCmd + pyArgs
                    chdirCmd + "start cmd.exe /C " + mkdirCmd + condaCmd + pyCmd + pyArgs
//                    "start cmd.exe /C " + chdirCmd + mkdirCmd + pyCmd + " --version;"
            };
            */

            // build commands for remote execution:
            final String chdirCmdSh = "cd " + modulesParentName + " ; ";
            final String mkdirCmdSh = "mkdir -p " + outputL2Dir + " ; ";
            final String pyCmdSh = "python ";

            final String[] commands = {
                    "/bin/sh",
                    "-c",
                    chdirCmdSh + mkdirCmdSh + pyCmdSh + pyArgs + ";"
            };
            Logger.info("Command sequence: " + Arrays.toString(commands));

            if (!simulation) {
                try {
                    Process process = Runtime.getRuntime().exec(commands);
                    process.waitFor();

                    BufferedReader reader = new BufferedReader(new InputStreamReader(process.getInputStream(),
                            StandardCharsets.UTF_8));

                    String line;
                    while ((line = reader.readLine()) != null) {
                        Logger.info(line);
                    }
                    Logger.info("Java Runtime.getRuntime(): exitValue = " + process.exitValue());
                    reader.close();

                    if (process.exitValue() != 0) {
                        throw new RuntimeException("Devalgo L2 Seaice Conc Python process terminated with an error");
                    }
                } catch (IOException | InterruptedException e) {
                    Logger.error("Executing Devalgo L2 openSF simulation from Java wrapper failed: " + e.getMessage());
                }
            }
        }
    }
}