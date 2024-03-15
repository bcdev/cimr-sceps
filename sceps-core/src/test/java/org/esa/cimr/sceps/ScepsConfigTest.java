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

}