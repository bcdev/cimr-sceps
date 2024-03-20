package org.esa.cimr.sceps;

import org.junit.Test;
import org.w3c.dom.Document;
import org.w3c.dom.Element;
import org.w3c.dom.Node;
import org.w3c.dom.NodeList;

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
}