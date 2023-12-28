package org.esa.cimr.sceps;

import org.junit.Test;
import org.w3c.dom.Document;
import org.w3c.dom.Element;
import org.w3c.dom.Node;
import org.w3c.dom.NodeList;

import java.util.Objects;

import static org.junit.Assert.*;

public class ScepsConfigTest {

    @Test
    public void testReadGlobalConfigFile() {

        final String globalConfigXmlPath =
                Objects.requireNonNull(getClass().getResource("Global_Configuration.xml")).getPath();
        Document document = null;
        try {
            document = ScepsConfig.readXMLDocumentFromFile(globalConfigXmlPath);
        } catch (Exception e) {
            fail("Cannot parse '" + globalConfigXmlPath + "'");
        }

        //Verify XML Content

        Element root = ScepsConfig.getDocumentRootElement(document);
        assertNotNull(root);
        assertEquals("Global_Configuration", root.getNodeName());
        assertNotNull(root.getAttributes());
        assertEquals(1, root.getAttributes().getLength());
        assertEquals("01.00.00", root.getAttribute("version"));

        NodeList nList = ScepsConfig.getDocumentElementsByTagName(document, "parameter");
        assertNotNull(nList);
        assertEquals(4, nList.getLength());

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

        assertEquals("v1", ScepsConfig.getDocumentElementTextItemByName(document, "parameter", "geodata_version"));

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

        elem = (Element) nList.item(3);
        assertEquals(3, elem.getAttributes().getLength());
        assertEquals("text", elem.getAttribute("description"));
        assertEquals("sceps_scd_root", elem.getAttribute("name"));
        assertEquals("STRING", elem.getAttribute("type"));
        assertNotNull(elem.getFirstChild());
        assertNotNull(elem.getFirstChild().getTextContent());
        assertEquals("/data/sceps/SCEPSscd", elem.getFirstChild().getTextContent());
    }

    @Test
    public void testReadGeoInputsLocalConfigFile() {

        final String globalConfigXmlPath =
                Objects.requireNonNull(getClass().
                        getResource("GeoInputs_Extract_Local_Configuration_geocard2_day1.xml")).getPath();
        Document document = null;
        try {
            document = ScepsConfig.readXMLDocumentFromFile(globalConfigXmlPath);
        } catch (Exception e) {
            fail("Cannot parse '" + globalConfigXmlPath + "'");
        }

        //Verify XML Content

        Element root = ScepsConfig.getDocumentRootElement(document);
        assertNotNull(root);
        assertEquals("GeoInputs_Extract_Local_Configuration", root.getNodeName());
        assertNotNull(root.getAttributes());
        assertEquals(1, root.getAttributes().getLength());
        assertEquals("0.1", root.getAttribute("version"));

        NodeList nList = ScepsConfig.getDocumentElementsByTagName(document, "parameter");
        assertNotNull(nList);
        assertEquals(6, nList.getLength());

        for (int i = 0; i < nList.getLength(); i++) {
            Node node = nList.item(i);
            assertNotNull(node);
            assertEquals("parameter", node.getNodeName());
            assertEquals(Node.ELEMENT_NODE, node.getNodeType());
        }

        Element elem = (Element) nList.item(0);
        assertEquals(3, elem.getAttributes().getLength());
        assertEquals("", elem.getAttribute("description"));
        assertEquals("scene_type", elem.getAttribute("name"));
        assertEquals("STRING", elem.getAttribute("type"));
        assertNotNull(elem.getFirstChild());
        assertNotNull(elem.getFirstChild().getTextContent());
        assertEquals("geo-nemo", elem.getFirstChild().getTextContent());

        elem = (Element) nList.item(2);
        assertEquals(3, elem.getAttributes().getLength());
        assertEquals("", elem.getAttribute("description"));
        assertEquals("data_folder", elem.getAttribute("name"));
        assertEquals("STRING", elem.getAttribute("type"));
        assertNotNull(elem.getFirstChild());
        assertNotNull(elem.getFirstChild().getTextContent());
        assertEquals("E2E_HOME/InputData/GeoInputData/GeoCardScenes/SCEPS_ts2", elem.getFirstChild().getTextContent());

        elem = (Element) nList.item(4);
        assertEquals(5, elem.getAttributes().getLength());
        assertEquals("", elem.getAttribute("description"));
        assertEquals("latitude_filter", elem.getAttribute("name"));
        assertEquals("ARRAY", elem.getAttribute("type"));
        assertEquals("FLOAT", elem.getAttribute("elementType"));
        assertEquals("2", elem.getAttribute("dims"));
        assertEquals(1, elem.getChildNodes().getLength());
        assertEquals("-90 90", elem.getChildNodes().item(0).getTextContent());
        final String[] latitudes = elem.getChildNodes().item(0).getTextContent().split(" ");
        assertEquals(2, latitudes.length);
        assertEquals(-90, Integer.parseInt(latitudes[0]));
        assertEquals(90, Integer.parseInt(latitudes[1]));
    }

    @Test
    public void testReadForwardModelLocalConfigFile() {

        final String globalConfigXmlPath =
                Objects.requireNonNull(getClass().getResource("Forward_Model_Local_Configuration_geocard2_day1.xml")).getPath();
        Document document = null;
        try {
            document = ScepsConfig.readXMLDocumentFromFile(globalConfigXmlPath);
        } catch (Exception e) {
            fail("Cannot parse '" + globalConfigXmlPath + "'");
        }

        //Verify XML Content

        Element root = ScepsConfig.getDocumentRootElement(document);
        assertNotNull(root);
        assertEquals("Forward_Model_Local_Configuration", root.getNodeName());
        assertNotNull(root.getAttributes());
        assertEquals(1, root.getAttributes().getLength());
        assertEquals("04.15.33", root.getAttribute("version"));

        NodeList nList = ScepsConfig.getDocumentElementsByTagName(document, "parameter");
        assertNotNull(nList);
        assertEquals(12, nList.getLength());

        for (int i = 0; i < nList.getLength(); i++) {
            Node node = nList.item(i);
            assertNotNull(node);
            assertEquals("parameter", node.getNodeName());
            assertEquals(Node.ELEMENT_NODE, node.getNodeType());
        }

        Element elem = (Element) nList.item(0);
        assertEquals(5, elem.getAttributes().getLength());
        assertEquals("", elem.getAttribute("description"));
        assertEquals("frequencies", elem.getAttribute("name"));
        assertEquals("ARRAY", elem.getAttribute("type"));
        assertEquals("FLOAT", elem.getAttribute("elementType"));
        assertEquals("5", elem.getAttribute("dims"));
        assertNotNull(elem.getChildNodes());
        assertNotNull(elem.getChildNodes().item(0).getTextContent());
        assertEquals("1.410 6.925 10.650 18.700 36.500", elem.getChildNodes().item(0).getTextContent());
        final String[] frequencies = elem.getChildNodes().item(0).getTextContent().split(" ");
        assertEquals(5, frequencies.length);
        assertEquals(1.410, Double.parseDouble(frequencies[0]), Double.MIN_VALUE);
        assertEquals(6.925, Double.parseDouble(frequencies[1]), Double.MIN_VALUE);
        assertEquals(10.650, Double.parseDouble(frequencies[2]), Double.MIN_VALUE);
        assertEquals(18.700, Double.parseDouble(frequencies[3]), Double.MIN_VALUE);
        assertEquals(36.500, Double.parseDouble(frequencies[4]), Double.MIN_VALUE);

        elem = (Element) nList.item(1);
        assertEquals(4, elem.getAttributes().getLength());
        assertEquals("", elem.getAttribute("description"));
        assertEquals("bands", elem.getAttribute("name"));
        assertEquals("STRING", elem.getAttribute("type"));
        assertEquals("1", elem.getAttribute("dims"));
        assertNotNull(elem.getChildNodes());
        assertEquals("L C X KuKa", elem.getChildNodes().item(0).getTextContent());
        final String[] bands = elem.getChildNodes().item(0).getTextContent().split(" ");
        assertEquals(4, bands.length);
        assertEquals("L", bands[0]);
        assertEquals("C", bands[1]);
        assertEquals("X", bands[2]);
        assertEquals("KuKa", bands[3]);

        elem = (Element) nList.item(2);
        assertEquals(5, elem.getAttributes().getLength());
        assertEquals("", elem.getAttribute("description"));
        assertEquals("v_polarization", elem.getAttribute("name"));
        assertEquals("ARRAY", elem.getAttribute("type"));
        assertEquals("INTEGER", elem.getAttribute("elementType"));
        assertEquals("5", elem.getAttribute("dims"));
        assertNotNull(elem.getChildNodes());
        assertEquals("1 1 1 1 1", elem.getChildNodes().item(0).getTextContent());
        final String[] vpol = elem.getChildNodes().item(0).getTextContent().split(" ");
        assertEquals(5, vpol.length);
        assertEquals(1, Integer.parseInt(vpol[0]));
        assertEquals(1, Integer.parseInt(vpol[1]));
        assertEquals(1, Integer.parseInt(vpol[2]));
        assertEquals(1, Integer.parseInt(vpol[3]));
        assertEquals(1, Integer.parseInt(vpol[4]));

        elem = (Element) nList.item(5);
        assertEquals(5, elem.getAttributes().getLength());
        assertEquals("", elem.getAttribute("description"));
        assertEquals("land_only_coast", elem.getAttribute("name"));
        assertEquals("ARRAY", elem.getAttribute("type"));
        assertEquals("INTEGER", elem.getAttribute("elementType"));
        assertEquals("1", elem.getAttribute("dims"));
        assertNotNull(elem.getChildNodes());
        final String loc = elem.getChildNodes().item(0).getTextContent();
        assertEquals(0, Integer.parseInt(loc));

        elem = (Element) nList.item(6);
        assertEquals(4, elem.getAttributes().getLength());
        assertEquals("freqsxscans", elem.getAttribute("description"));
        assertEquals("zenith_angle", elem.getAttribute("name"));
        assertEquals("FLOAT", elem.getAttribute("type"));
        assertEquals("1 5", elem.getAttribute("dims"));
        assertNotNull(elem.getChildNodes());
        assertNotNull(elem.getChildNodes().item(0).getTextContent());
        assertEquals("40 55 55 55 55", elem.getChildNodes().item(0).getTextContent());
        final String[] za = elem.getChildNodes().item(0).getTextContent().split(" ");
        assertEquals(5, frequencies.length);
        assertEquals(40, Double.parseDouble(za[0]), Double.MIN_VALUE);
        assertEquals(55, Double.parseDouble(za[1]), Double.MIN_VALUE);
        assertEquals(55, Double.parseDouble(za[2]), Double.MIN_VALUE);
        assertEquals(55, Double.parseDouble(za[3]), Double.MIN_VALUE);
        assertEquals(55, Double.parseDouble(za[4]), Double.MIN_VALUE);

        elem = (Element) nList.item(11);
        assertEquals(3, elem.getAttributes().getLength());
        assertEquals("emis_land_external", elem.getAttribute("name"));
        assertEquals("FLOAT", elem.getAttribute("type"));
        assertEquals("1 1", elem.getAttribute("dims"));
        assertNotNull(elem.getChildNodes());
        assertNotNull(elem.getChildNodes().item(0).getTextContent());
        final String ele = elem.getChildNodes().item(0).getTextContent();
        assertEquals("-999 ", ele);
        assertEquals(-999, Double.parseDouble(ele), Double.MIN_VALUE);

    }
}