package org.esa.cimr.sceps;

import org.junit.Test;

import java.io.IOException;

import static org.esa.cimr.sceps.ScepsConstants.SCEPS_CODES_FOLDER_NAME;
import static org.junit.Assert.assertEquals;
import static org.junit.Assert.fail;

public class ObservingSystemSimulationModuleMrtWrapperTest {

    @Test
    public void testGetWrappedCommand_Sensor_Apply_Antenna() {

        final String E2E_HOME = "/MY/E2E/HOME";
        final String mrtExecPath = E2E_HOME + "/" + SCEPS_CODES_FOLDER_NAME + "/RunTimeEnv/Sensor_Apply_Antenna";
        final String globalConfigPath = E2E_HOME + "/executions/SCEPS_ObsSysSim.last/Global_Configuration.xml";
        final String localConfigPath = E2E_HOME + "/executions/SCEPS_ObsSysSim.last/Sensor_Apply_Antenna_Local_Configuration_bla.xml";
        final String firstInputFile = E2E_HOME + "/executions/SCEPS_ObsSysSim.last/first.nc";
        final String secondInputFile = E2E_HOME + "/executions/SCEPS_ObsSysSim.last/second.nc";
        final String thirdInputFile = E2E_HOME + "/executions/SCEPS_ObsSysSim.last/third.nc";
        final String outputDir = E2E_HOME + "/executions/SCEPS_ObsSysSim.last/Sensor_Apply_Antenna_Output";

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
            final String commandString = ObservingSystemSimulationModuleMrtWrapper.execute(args, E2E_HOME);

            assertEquals(commandString, expectedCommandString);
        } catch (IOException e) {
            fail(e.getMessage());
        }
    }


}