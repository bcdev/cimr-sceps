package org.esa.cimr.sceps;

import org.junit.Test;
import org.w3c.dom.Document;
import org.w3c.dom.Element;
import org.w3c.dom.Node;
import org.w3c.dom.NodeList;

import java.io.IOException;
import java.util.Objects;

import static org.junit.Assert.*;

public class SceneGenerationModuleMrtWrapperTest {

    @Test
    public void testGetWrappedCommand_GeoInputs_Extract() {

        final String E2E_HOME = "/MY/E2E/HOME";
        final String mrtExecPath = E2E_HOME + "/SCEPScodes/SceGenMod/Modules/GeoInputs_Extract";
        final String globalConfigPath = E2E_HOME + "/executions/SCEPS_SceGenSim.last/Global_Configuration.xml";
        final String localConfigPath = E2E_HOME + "/executions/SCEPS_SceGenSim.last/GeoInputs_Extract_Local_Configuration_bla.xml";
        final String firstInputFile = E2E_HOME + "/executions/SCEPS_SceGenSim.last/first.nc";
        final String secondInputFile = E2E_HOME + "/executions/SCEPS_SceGenSim.last/second.nc";
        final String thirdInputFile = E2E_HOME + "/executions/SCEPS_SceGenSim.last/third.nc";
        final String outputDir = E2E_HOME + "/executions/SCEPS_SceGenSim.last/GeoInputs_Extract_Output";

        final String[] args = new String[]{
                "--global",
                globalConfigPath,
                "--local",
                localConfigPath,
                "--input",
                firstInputFile,
                "--input",
                secondInputFile,
                "--input",
                thirdInputFile,
                "--output",
                outputDir
        };

        try {
            final String expectedCommandString =
                    mrtExecPath + " " +
                    globalConfigPath + "," + localConfigPath + " " +
                    firstInputFile + "," + secondInputFile + "," + thirdInputFile + " " +
                    outputDir;
            final String commandString = SceneGenerationModuleMrtWrapper.execute(args, E2E_HOME);

            assertEquals(commandString, expectedCommandString);
        } catch (IOException e) {
            fail(e.getMessage());
        }
    }

    @Test
    public void testGetWrappedCommand_Forward_Model() {

        final String E2E_HOME = "/MY/E2E/HOME";
        final String mrtExecPath = E2E_HOME + "/SCEPScodes/SceGenMod/Modules/Forward_Model";
        final String globalConfigPath = E2E_HOME + "/executions/SCEPS_SceGenSim.last/Global_Configuration.xml";
        final String localConfigPath = E2E_HOME + "/executions/SCEPS_SceGenSim.last/Forward_Model_Local_Configuration_blubb.xml";
        final String firstInputFile = E2E_HOME + "/executions/SCEPS_SceGenSim.last/GeoInputs_Extract_Output";
        final String outputDir = E2E_HOME + "/executions/SCEPS_SceGenSim.last/Forward_Model_Output";

        final String[] args = new String[]{
                "--global",
                globalConfigPath,
                "--local",
                localConfigPath,
                "--input",
                firstInputFile,
                "--output",
                outputDir
        };

        try {
            final String expectedCommandString =
                    mrtExecPath + " " +
                            globalConfigPath + "," + localConfigPath + " " +
                            firstInputFile + " " +
                            outputDir;
            final String commandString = SceneGenerationModuleMrtWrapper.execute(args, E2E_HOME);

            assertEquals(commandString, expectedCommandString);
        } catch (IOException e) {
            fail(e.getMessage());
        }
    }
}