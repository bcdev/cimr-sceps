package org.esa.cimr.sceps;

import org.junit.Ignore;
import org.junit.Test;
import org.w3c.dom.Document;
import org.w3c.dom.Element;
import org.w3c.dom.Node;
import org.w3c.dom.NodeList;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.nio.charset.StandardCharsets;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Objects;

import static org.junit.Assert.*;

public class DevalgoL2SeaiceConcentrationModuleConfigTest {

    @Test
    public void testReadDevalgoL2SeaiceConcentrationLocalConfigFile() {

        final String localalConfigXmlPath =
                Objects.requireNonNull(getClass().
                        getResource("Devalgo_L2_Seaice_Concentration_Local_Configuration.xml")).getPath();
        Document document = null;
        try {
            document = ScepsConfig.readXMLDocumentFromFile(localalConfigXmlPath);
        } catch (Exception e) {
            fail("Cannot parse '" + localalConfigXmlPath + "'");
        }

        //Verify XML Content

        Element root = ScepsConfig.getDocumentRootElement(document);
        assertNotNull(root);
        assertEquals("Devalgo_L2_Seaice_Concentration_Local_Configuration", root.getNodeName());
        assertNotNull(root.getAttributes());
        assertEquals(1, root.getAttributes().getLength());
        assertEquals("0.1", root.getAttribute("version"));

        NodeList nList = ScepsConfig.getDocumentElementsByTagName(document, "parameter");
        assertNotNull(nList);
        assertEquals(1, nList.getLength());

        for (int i = 0; i < nList.getLength(); i++) {
            Node node = nList.item(i);
            assertNotNull(node);
            assertEquals("parameter", node.getNodeName());
            assertEquals(Node.ELEMENT_NODE, node.getNodeType());
        }

        Element elem = (Element) nList.item(0);
        assertEquals(3, elem.getAttributes().getLength());
        assertEquals("", elem.getAttribute("description"));
        assertEquals("l2_grid", elem.getAttribute("name"));
        assertEquals("STRING", elem.getAttribute("type"));
        assertNotNull(elem.getFirstChild());
        assertNotNull(elem.getFirstChild().getTextContent());
        assertEquals("testgrid", elem.getFirstChild().getTextContent());

        // todo: further fill config parameters if needed and continue testing here
    }

    @Test
    @Ignore
    public void testGetPythonSitePackagesPath_main_install_dir_linux() {
        //test Windows with python.exe in main install dir:
        // before enabling the test, make sure this is an existing Python installation:
        String pyInstallLinmac = "/data/software/miniconda3_py310/bin/";
        String pyExec = pyInstallLinmac + "python";

        final String[] commands = {
                pyExec,
                "-c",
                "\"import sysconfig; print(sysconfig.get_path('platlib')); print('bla')\""
        };

        System.out.println("Command sequence: " + Arrays.toString(commands));

        try {
            ProcessBuilder processBuilder = new ProcessBuilder().command(commands);
            final Process process1 = processBuilder.start();
            process1.waitFor();
            BufferedReader reader3 = new BufferedReader(new InputStreamReader(process1.getInputStream(),
                    StandardCharsets.UTF_8));

            String line3;
            List lines3= new ArrayList();

            while ((line3 = reader3.readLine()) != null) {
                System.out.printf(line3);
                // store output in a list of lines. Anyway, the result should be just one line.
                lines3.add(line3);
            }
            reader3.close();
            System.out.println("lines3.size() = " + lines3.size());
            final String[] linesArr3 = (String[]) lines3.toArray(new String[0]);
            System.out.println("linesArr3.length = " + linesArr3.length);


            Process process = Runtime.getRuntime().exec(commands);
            process.waitFor();

            BufferedReader reader = new BufferedReader(new InputStreamReader(process.getInputStream(),
                    StandardCharsets.UTF_8));

            String line;
            List lines = new ArrayList();

            while ((line = reader.readLine()) != null) {
                System.out.printf(line);
                // store output in a list of lines. Anyway, the result should be just one line.
                lines.add(line);
            }
            reader.close();

            if (process.exitValue() != 0) {
                throw new RuntimeException("External Python call terminated with an error.");
            }

            System.out.println("lines.size() = " + lines.size());
            final String[] linesArr = (String[]) lines.toArray(new String[0]);
            System.out.println("linesArr.length = " + linesArr.length);

            System.out.println("process.exitValue = " + process.exitValue());

            System.out.println("result: " + Paths.get(linesArr[0]));

        } catch (IOException | InterruptedException e) {
            throw new RuntimeException("Determination of Python site-packages path failed: " + e.getMessage());
        }
    }

}