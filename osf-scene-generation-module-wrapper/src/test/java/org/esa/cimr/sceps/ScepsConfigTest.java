package org.esa.cimr.sceps;

import org.esa.cimr.sceps.ScepsConfig;
import org.junit.Test;
import org.w3c.dom.Document;
import org.w3c.dom.Element;
import org.w3c.dom.Node;
import org.w3c.dom.NodeList;
import org.xml.sax.SAXException;

import javax.xml.parsers.DocumentBuilder;
import javax.xml.parsers.DocumentBuilderFactory;
import javax.xml.parsers.ParserConfigurationException;

import java.io.File;
import java.io.IOException;

import static org.junit.Assert.*;

public class ScepsConfigTest {

    @Test
    public void testReadXmlGlobalConfigFile() {

        final String globalConfigXmlPath = getClass().getResource("Global_Configuration.xml").getPath();
        Document document = null;
        try {
            document = ScepsConfig.readXMLDocumentFromFile(globalConfigXmlPath);
        } catch (Exception e) {
            fail("Cannot parse '" + globalConfigXmlPath + "'");
        }

        //Verify XML Content

        Element root = document.getDocumentElement();
        assertNotNull(root);
        assertEquals("Global_Configuration", root.getNodeName());
        assertNotNull(root.getAttributes());
        assertEquals(1, root.getAttributes().getLength());
        assertEquals("01.00.00", root.getAttribute("version"));

        NodeList nList = document.getElementsByTagName("parameter");
        assertNotNull(nList);
        assertEquals(3, nList.getLength());

        for (int i = 0; i < nList.getLength(); i++) {
            Node node = nList.item(i);
            assertNotNull(node);
            assertEquals("parameter", node.getNodeName());
            assertEquals(Node.ELEMENT_NODE, node.getNodeType());
        }

        Element elem = (Element) nList.item(0);
        assertEquals(3, elem.getAttributes().getLength());
        assertEquals("text", elem.getAttribute("description"));
        assertEquals("geodata_version", elem.getAttribute("name"));
        assertEquals("STRING", elem.getAttribute("type"));
        assertNotNull(elem.getFirstChild());
        assertNotNull(elem.getFirstChild().getTextContent());
        assertEquals("v1", elem.getFirstChild().getTextContent());

        elem = (Element) nList.item(1);
        assertEquals(3, elem.getAttributes().getLength());
        assertEquals("text", elem.getAttribute("description"));
        assertEquals("software_version", elem.getAttribute("name"));
        assertEquals("STRING", elem.getAttribute("type"));
        assertNotNull(elem.getFirstChild());
        assertNotNull(elem.getFirstChild().getTextContent());
        assertEquals("v1", elem.getFirstChild().getTextContent());

        elem = (Element) nList.item(2);
        assertEquals(3, elem.getAttributes().getLength());
        assertEquals("text", elem.getAttribute("description"));
        assertEquals("workers_number", elem.getAttribute("name"));
        assertEquals("INTEGER", elem.getAttribute("type"));
        assertNotNull(elem.getFirstChild());
        assertNotNull(elem.getFirstChild().getTextContent());
        assertEquals("0", elem.getFirstChild().getTextContent());
    }

    @Test
    public void testReadXmlLocalConfigFile() {
//        <?xml version="1.0" encoding="UTF-8"?>
//        <Forward_Model_Local_Configuration version="04.15.33">
//           <parameter name="frequencies" description="" type="ARRAY" elementType="FLOAT" dims="5">1.410 6.925 10.650 18.700 36.500</parameter>
//           <parameter name="bands" description="" type="STRING" dims="1">L C X KuKa</parameter>
//           <parameter name="v_polarization" description="" type="ARRAY" elementType="INTEGER" dims="5">1 1 1 1 1</parameter>
//           <parameter name="h_polarization" description="" type="ARRAY" elementType="INTEGER" dims="5">1 1 1 1 1</parameter>
//           <parameter name="atmosphere_adding" description="" type="ARRAY" elementType="INTEGER" dims="5">1 1 1 1 1</parameter>
//           <parameter name="land_only_coast" description="" type="ARRAY" elementType="INTEGER" dims="1">0</parameter>
//           <parameter name="zenith_angle" description="freqsxscans" type="FLOAT" dims="1 5" >40 55 55 55 55</parameter>
//           <parameter name="azimuth_angle" description="freqsxscans" type="FLOAT" dims="1 5" >0 0 0 0 0  </parameter>
//           <parameter name="sensor_height" description="" type="FLOAT">853</parameter>
//           <parameter name="data_thinning" description="" type="INTEGER">-2</parameter>
//           <parameter name="emis_ice_smoothing" description="" type="FLOAT">1</parameter>
//           <parameter name="emis_land_external" type="FLOAT" dims="1 1" >-999 </parameter>
//        </Forward_Model_Local_Configuration>

        final String fwModelTestString =
                "<Forward_Model_Local_Configuration version=\"04.15.33\">\n" +
                        "<parameter name=\"frequencies\" description=\"\" type=\"ARRAY\" elementType=\"FLOAT\" dims=\"5\">1.410 6.925 10.650 18.700 36.500</parameter>\n" +
                        "<parameter name=\"bands\" description=\"\" type=\"STRING\" dims=\"1\">L C X KuKa</parameter>\n" +
                        "<parameter name=\"v_polarization\" description=\"\" type=\"ARRAY\" elementType=\"INTEGER\" dims=\"5\">1 1 1 1 1</parameter>\n" +
                        "<parameter name=\"h_polarization\" description=\"\" type=\"ARRAY\" elementType=\"INTEGER\" dims=\"5\">1 1 1 1 1</parameter>\n" +
                        "<parameter name=\"atmosphere_adding\" description=\"\" type=\"ARRAY\" elementType=\"INTEGER\" dims=\"5\">1 1 1 1 1</parameter>\n" +
                        "<parameter name=\"land_only_coast\" description=\"\" type=\"ARRAY\" elementType=\"INTEGER\" dims=\"1\">0</parameter>\n" +
                        "<parameter name=\"zenith_angle\" description=\"freqsxscans\" type=\"FLOAT\" dims=\"1 5\" >40 55 55 55 55</parameter>\n" +
                        "<parameter name=\"azimuth_angle\" description=\"freqsxscans\" type=\"FLOAT\" dims=\"1 5\" >0 0 0 0 0  </parameter>\n" +
                        "<parameter name=\"sensor_height\" description=\"\" type=\"FLOAT\">853</parameter>\n" +
                        "<parameter name=\"data_thinning\" description=\"\" type=\"INTEGER\">-2</parameter>\n" +
                        "<parameter name=\"emis_ice_smoothing\" description=\"\" type=\"FLOAT\">1</parameter>\n" +
                        "<parameter name=\"emis_land_external\" type=\"FLOAT\" dims=\"1 1\" >-999 </parameter>\n" +
                        "</Forward_Model_Local_Configuration>";

        // todo
        assertTrue(true);
    }
}