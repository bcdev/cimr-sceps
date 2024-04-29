package org.esa.cimr.sceps;

import esa.opensf.osfi.CLP;
import esa.opensf.osfi.Logger;
import esa.opensf.osfi.ParamReader;
import org.apache.commons.io.FilenameUtils;
import org.esa.cimr.sceps.util.ScepsUtils;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStreamReader;
import java.nio.charset.StandardCharsets;
import java.util.Arrays;

import static org.esa.cimr.sceps.ScepsConstants.*;

/**
 * Class containing the wrapper for the OSF Scene Generation Scientific Moduule.
 * This class shall:
 * - be made available as executable jar (see pom) to be called from openSF
 * (set as 'Executable' in corresponding openSF module General tab).
 * - properly transform openSF arguments/parameters to ScientificModule arguments/parameters
 * - do a system call to run the underlying Scientific Module in Matlab batch mode
 *
 * In the CIMR SCEPS context, this wrapper will usually be called twice, subsequently for the two
 * scientific modules GeoInputs_Extract and 'Forward_Model'.
 */
public class SceneGenerationModuleWrapper {

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
//                        # global config file
//                        # local config file
//                        # input files 1..N
//                        # output files 1..M

//    java -jar /data/workspace/cimr-sceps/osf-scene-generation-module-wrapper/target/osf-scene-generation-module-wrapper-1.0-SNAPSHOT.jar
//   /data/workspace/cimr-opensf/simulations/SCEPS_SceGen_test.20231229T192846d191/Global_Configuration.xml,
//   /data/workspace/cimr-opensf/simulations/SCEPS_SceGen_test.20231229T192846d191/GeoInputs_Extract_Local_Configuration_geocard2_day1.xml
//   /data/workspace/cimr-opensf/simulations/SCEPS_SceGen_test.20231229T192846d191/GeoInputs_Extract_Output

    /**
     * The wrapper main method.
     *
     * @param args - Program arguments. Two single strings with comma separated entries. These are:
     *             # first string: <full path of global config file>,<full path of global config file>
     *             # second string: <full paths of input files 1..N>,<full paths of output files 1..M>
     *             # Optionally, a 'simulation' mode can be set: last arg is 'simulation=true'
     */
    public static void main(String[] args) throws IOException {
        if (args.length < 1) {
            Logger.error("Wrong number of arguments given - must be one comma separated string containing " +
                    "global and local config file, all input and output files. Exiting.");
            System.exit(1);
        } else {
            // 'simulation' mode for testing. Omits execution of the Matlab batch command.
            final boolean simulation = args[args.length - 1].equals("simulation=true");

            // get program arguments from OSFI command line parser:
            final CLP clp = new CLP(args);

            // check if args are complete and correctly parsed.
            // Must contain global and local config, inputs and outputs are optional
            ScepsUtils.checkCommandLineArgs(clp);

            // We only parse config files but ignore inputs and outputs.
            // They are set explicitly below, following the needs of the two Matlab modules.
            final String globalConfigXmlPath = clp.getConfFiles().get(0);
            final String localConfigXmlPath = clp.getConfFiles().get(1);

            Logger.info("globalConfigXmlPath: " + globalConfigXmlPath);
            Logger.info("localConfigXmlPath: " + localConfigXmlPath);

            // set relevant paths:
            String scepsScdRoot;
            String moduleName;
            String sceneType;
            String sceneDate;
            try {
                // this was added to global config:
                // <parameter description="text" name="sceps_scd_root" type="STRING">/data/sceps/SCEPSscd</parameter>
                final ParamReader globalParamReader = new ParamReader(globalConfigXmlPath);
                scepsScdRoot = globalParamReader.getParameter("E2E_HOME").getStringValue();

                moduleName = FilenameUtils.removeExtension((new File(localConfigXmlPath)).getName());
                // strip extension '_Local_Configuration':
                if (moduleName.contains("_Local_Configuration")) {
                    int index = moduleName.indexOf("_Local_Configuration");
                    moduleName = moduleName.substring(0, index);
                }

                // we need SCENE_TYPE and SCENE_DATE as global variables from GeoInputs_Extract config:
                // It's in GeoInputs_Extract config only, thus this was added to Forward_Model local config:
                final ParamReader localParamReader = new ParamReader(localConfigXmlPath);
                sceneType = localParamReader.getParameter(SCEPS_SCENE_TYPE_CONFIG_ITEM_NAME).getStringValue();
                sceneDate = localParamReader.getParameter(SCEPS_SCENE_DATE_CONFIG_ITEM_NAME).getStringValue();
            } catch (Exception e) {
                // todo
                throw new RuntimeException(e);
            }

            final String devSCEPSpath = scepsScdRoot + File.separator + ScepsConstants.SCEPS_CODES_FOLDER_NAME;
            final String dataSCEPSpath = scepsScdRoot + File.separator + ScepsConstants.SCEPS_DATA_FOLDER_NAME;
            final String modulesParentPath = devSCEPSpath + File.separator +
                    ScepsConstants.SCENE_GENERATION_MODULE_FOLDER_NAME + File.separator +
                    ScepsConstants.MODULES_SUBFOLDER_NAME;
            final String subModulesParentPath = devSCEPSpath + File.separator +
                    ScepsConstants.SCENE_GENERATION_MODULE_FOLDER_NAME + File.separator +
                    ScepsConstants.MODULE_SUBMODULES_SUBFOLDER_NAME;

            // set relevant parameters to match module name signature (see e.g. GeoInputs_Extract.m):
            final String configurationParameters = globalConfigXmlPath + "," + localConfigXmlPath;
            final File globalConfigXmlFile = new File(globalConfigXmlPath);
            final String geoinputSimulationOutputFolder = globalConfigXmlFile.getParent() + File.separator +
                    GEOINPUTS_EXTRACT_MODULE_NAME + "_Output";
            final String geoinputSimulationMatlabGlobal = geoinputSimulationOutputFolder + File.separator +
                    GEOINPUTS_EXTRACT_MODULE_NAME;

            final String outputs = globalConfigXmlFile.getParent() + File.separator + moduleName + "_Output";
            String inputs;
            if (moduleName.equals(ScepsConstants.GEOINPUTS_EXTRACT_MODULE_NAME)) {
                inputs = ScepsUtils.clpInputsJava2Matlab(clp.getInputFiles());
            } else if (moduleName.equals(FORWARD_MODEL_MODULE_NAME)) {
                inputs = geoinputSimulationOutputFolder;
            } else {
                throw new IOException("Module name " + moduleName + " not known.");
            }
            Logger.info("moduleName: " + moduleName);
            Logger.info("inputs: " + inputs);
            Logger.info("outputs: " + outputs);

            final String scepsPathCmdSh = "devSCEPSpath = '" + devSCEPSpath + "'; ";
            final String addpath1CmdSh =
                    "addpath '" + devSCEPSpath + File.separator + SCEPS_CODES_GENERAL_SUBMODULES_FOLDER_NAME + "'; ";
            final String addpath2CmdSh =
                    "addpath '" + devSCEPSpath + File.separator + SCEPS_CODES_OSFI_MATLAB_FOLDER_NAME + "'; ";
            final String addpath3CmdSh = "addpath '" + modulesParentPath + "'; ";
            final String addpath4CmdSh = "addpath '" + subModulesParentPath + "'; ";
            final String chdirCmdSh = "cd " + modulesParentPath + "; ";
            final String matlabGlobalVarsString =
                    "global E2E_HOME; E2E_HOME = '" + scepsScdRoot + "'; " +
                            "global SCENE_TYPE; SCENE_TYPE = '" + sceneType + "'; " +
                            "global SCENE_DATE; SCENE_DATE = '" + sceneDate + "'; " +
                            "global GEOINPUT_SIMULATION; GEOINPUT_SIMULATION = '" + geoinputSimulationMatlabGlobal + "'; " +
                            "global LOG; LOG = Logger(); ";

            final String[] commands = {
                    "matlab",
                    "-batch",
                    scepsPathCmdSh +
                            addpath1CmdSh + addpath2CmdSh + addpath3CmdSh + addpath4CmdSh +
                            chdirCmdSh + matlabGlobalVarsString +
                            moduleName + "('" + configurationParameters + "','" + inputs + "','" + outputs + "');"
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
                        throw new RuntimeException("SceGen Matlab process terminated with an error");
                    }

                } catch (IOException | InterruptedException e) {
                    throw new RuntimeException("Executing SceGen openSF simulation from Java wrapper failed: " + e.getMessage());
                }
            }
        }
    }
}