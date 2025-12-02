package org.esa.cimr.sceps;

import esa.opensf.osfi.CLP;
import esa.opensf.osfi.Logger;
import esa.opensf.osfi.ParamReader;
import esa.opensf.osfi.xmlutils.XmlParseException;
import org.apache.commons.io.FilenameUtils;
import org.esa.cimr.sceps.util.ScepsUtils;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStreamReader;
import java.nio.charset.StandardCharsets;

import static org.esa.cimr.sceps.ScepsConstants.*;

/**
 * Class containing the wrapper for the OSF Scene Generation Scientific Moduule.
 * This class shall:
 * - be made available as executable jar (see pom) to be called from openSF
 * (set as 'Executable' in corresponding openSF module General tab).
 * - properly transform openSF arguments/parameters to ScientificModule arguments/parameters
 * - do a system call to run the underlying Scientific Module in Matlab batch mode
 * <p>
 * In the CIMR SCEPS context, this wrapper will usually be called twice, subsequently for the two
 * scientific modules GeoInputs_Extract and 'Forward_Model'.
 */
public class SceneGenerationModuleMrtWrapper {

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
        System.out.println("Wrapped command: " + execute(args));
    }

    static String execute(String[] args) throws IOException {
        // default mode
        return execute(args,  null);
    }

    static String execute(String[] args, String e2eHome) throws IOException {
        if (args.length < 1) {
            Logger.error("Wrong number of arguments given - must be one comma separated string containing " +
                    "global and local config file, all input and output files. Exiting.");
            System.exit(1);
        } else {
            // for testing, E2E_HOME can pe passed
            final boolean simulation = e2eHome != null;

            // get program arguments from OSFI command line parser:
            final CLP clp = new CLP(args);

            // check if args are complete and correctly parsed.
            // Must contain global and local config, inputs and outputs are optional
            ScepsUtils.checkCommandLineArgs(clp);

            // We only parse config files but ignore inputs and outputs.
            // They are set explicitly below, following the needs of the two Matlab modules.
            final String globalConfigXmlPath = clp.getConfFiles().get(0);
            final String localConfigXmlPath = clp.getConfFiles().get(1);

            if (!simulation) {
                final ParamReader globalParamReader;
                try {
                    globalParamReader = new ParamReader(globalConfigXmlPath);
                    e2eHome = globalParamReader.getParameter("E2E_HOME").getStringValue();
                } catch (XmlParseException e) {
                    throw new RuntimeException(e);
                }
            }

            String moduleName = FilenameUtils.removeExtension((new File(localConfigXmlPath)).getName());
            // strip extension '_Local_Configuration':
            if (moduleName.contains("_Local_Configuration")) {
                int index = moduleName.indexOf("_Local_Configuration");
                moduleName = moduleName.substring(0, index);
            }
            final String modulePath = e2eHome + "/" + SCEPS_CODES_FOLDER_NAME + "/SceGenMod/Modules/" + moduleName;
            if (!simulation && !(new File(modulePath).exists())) {
                throw new IOException("Wrong structure of SCEPS package:  Module path: " + modulePath +
                        ". Must be: E2E_HOME/" + SCEPS_CODES_FOLDER_NAME +
                        "/SceGenMod/Modules/<module name>");
            }

            // set relevant parameters to match module name signature (see e.g. GeoInputs_Extract.m):
            final String configurationParameters = globalConfigXmlPath + "," + localConfigXmlPath;
            final String inputs = ScepsUtils.clpInputsOutputsJava2Matlab(clp.getInputFiles());
            final String outputs = ScepsUtils.clpInputsOutputsJava2Matlab(clp.getOutputFiles());

            final String[] commands = {
                    modulePath,
                    configurationParameters,
                    inputs,
                    outputs
            };

            final String commandsString = modulePath + " " + configurationParameters + " " +
                    inputs + " " + outputs;
            Logger.info("Wrapped command: " + commandsString);

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
                        throw new RuntimeException("SceGen Matlab Runtime process terminated with an error");
                    }

                } catch (IOException | InterruptedException e) {
                    throw new RuntimeException("Executing SceGen openSF Matlab Runtime " +
                            "simulation from Java wrapper failed: " + e.getMessage());
                }
            }
            return commandsString;
        }
        return null;
    }
}