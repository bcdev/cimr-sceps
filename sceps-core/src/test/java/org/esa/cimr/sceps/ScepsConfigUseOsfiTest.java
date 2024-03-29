package org.esa.cimr.sceps;

import esa.opensf.osfi.CLP;
import esa.opensf.osfi.Logger;
import esa.opensf.osfi.ParamReader;
import esa.opensf.osfi.Parameter;
import esa.opensf.osfi.xmlutils.XmlParseException;
import org.junit.After;
import org.junit.Before;
import org.junit.Test;

import java.io.FileNotFoundException;
import java.util.List;
import java.util.Map;
import java.util.Objects;

import static org.junit.Assert.*;

public class ScepsConfigUseOsfiTest {

    private String globalConfigXmlPath;
    private String localConfigXmlPath;
    private String inputFile1;
    private String inputFile2;
    private String outputFile1;
    private String outputFile2;
    private String inputFiles;
    private String outputFiles;
    private String configFiles;

    @Before
    public void setUp() {
        globalConfigXmlPath =
                Objects.requireNonNull(getClass().getResource("Global_Configuration.xml")).getPath();
        localConfigXmlPath =
                Objects.requireNonNull(getClass().getResource("Some_Local_Configuration.xml")).getPath();
        configFiles = globalConfigXmlPath + "," + localConfigXmlPath;

        inputFile1 = "/path/to/input1.txt";
        inputFile2 = "/path/to/input2.csv";
        inputFiles = inputFile1 + "," + inputFile2;

        outputFile1 = "/path/to/output1.nc";
        outputFile2 = "/path/to/output2.png";
        outputFiles = outputFile1 + "," + outputFile2;
    }

    @After
    public void tearDown() {
        // nothing to do
    }

    @Test
    public void testCLP() {
        Logger.info("Starting testCLP...");

        final String[] args = new String[]{configFiles, inputFiles, outputFiles};
        CLP clp = new CLP(args);

        final List<String> clpConfFiles = clp.getConfFiles();
        final List<String> clpInputFiles = clp.getInputFiles();
        final List<String> clpOutputFiles = clp.getOutputFiles();

        assertNotNull(clpConfFiles);
        assertNotNull(clpInputFiles);
        assertNotNull(clpOutputFiles);

        assertEquals(2, clpConfFiles.size());
        assertEquals(2, clpInputFiles.size());
        assertEquals(2, clpOutputFiles.size());

        assertEquals(globalConfigXmlPath, clpConfFiles.get(0));
        assertEquals(localConfigXmlPath, clpConfFiles.get(1));

        assertEquals(inputFile1, clp.getInputFiles().get(0));
        assertEquals(inputFile2, clp.getInputFiles().get(1));
        assertEquals(outputFile1, clp.getOutputFiles().get(0));
        assertEquals(outputFile2, clp.getOutputFiles().get(1));

        Logger.info("Finished testCLP.");
    }

    @Test
    public void testReadGlobalConfigFile() {
        Logger.info("Starting testReadGlobalConfigFile...");

        try {
            final ParamReader globalParamReader = new ParamReader(globalConfigXmlPath);
            Map<String, Parameter> globalParams = globalParamReader.getAllParameters();
            assertNotNull(globalParams);
            assertEquals(4, globalParams.size());

            final Parameter geodataVersionParam = globalParams.get("geodata_version");
            assertNotNull(geodataVersionParam);
            assertEquals("v1", geodataVersionParam.getStringValue());
            assertEquals("text", geodataVersionParam.getDescription());
            assertEquals("STRING", geodataVersionParam.getType());

            final Parameter softwareVersionParam = globalParams.get("software_version");
            assertNotNull(softwareVersionParam);
            assertEquals("v1", softwareVersionParam.getStringValue());
            assertEquals("text", softwareVersionParam.getDescription());
            assertEquals("STRING", softwareVersionParam.getType());

            final Parameter workersNumberParam = globalParams.get("workers_number");
            assertNotNull(workersNumberParam);
            assertEquals(0, workersNumberParam.getIntValue());
            assertEquals("text", workersNumberParam.getDescription());
            assertEquals("INTEGER", workersNumberParam.getType());

            final Parameter scepsScdRootParam = globalParams.get("sceps_scd_root");
            assertNotNull(scepsScdRootParam);
            assertEquals("/data/sceps/SCEPSscd", scepsScdRootParam.getStringValue());
            assertEquals("text", scepsScdRootParam.getDescription());
            assertEquals("STRING", scepsScdRootParam.getType());

        } catch (FileNotFoundException | XmlParseException e) {
            Logger.error("testReadGlobalConfigFile failed!");
            fail(e.getMessage());
        }

        Logger.info("Finished testReadGlobalConfigFile.");
    }

    @Test
    public void testReadLocalConfigFile() {
        Logger.info("Starting testReadLocalConfigFile...");
        try {
            final ParamReader localParamReader = new ParamReader(localConfigXmlPath);
            Map<String, Parameter> localParams = localParamReader.getAllParameters();
            assertNotNull(localParams);
            assertEquals(14, localParams.size());

            final Parameter frequencyParam = localParams.get("frequencies");
            assertNotNull(frequencyParam);
            assertEquals(5, frequencyParam.getDims().get(0).intValue());
            assertEquals("ARRAY", frequencyParam.getComplexType());
            assertEquals("FLOAT", frequencyParam.getElementType().name());
            assertEquals("", frequencyParam.getDescription());
            final double[] expectedFrequencies = new double[]{1.410, 6.925, 10.650, 18.700, 36.500};
            assertArrayEquals(expectedFrequencies, frequencyParam.getRootNode().getTreeDouble().getData(), 1.E-6);

            final Parameter zenithAngleParam = localParams.get("zenith_angle");
            assertNotNull(zenithAngleParam);
            assertEquals(2, zenithAngleParam.getDims().size());
            assertEquals(1, zenithAngleParam.getDims().get(0).intValue());
            assertEquals(5, zenithAngleParam.getDims().get(1).intValue());
            assertEquals("FLOAT", zenithAngleParam.getElementType().name());
            assertEquals("freqsxscans", zenithAngleParam.getDescription());
            final double[] expectedZenithAngles = new double[]{40, 55, 55, 55, 55};
            assertArrayEquals(expectedZenithAngles, zenithAngleParam.getRootNode().getTreeDouble().getData(), 1.E-6);

        } catch (FileNotFoundException | XmlParseException e) {
            Logger.error("testReadLocalConfigFile failed!");
            fail(e.getMessage());
        }

        Logger.info("Finished testReadLocalConfigFile.");
    }
}