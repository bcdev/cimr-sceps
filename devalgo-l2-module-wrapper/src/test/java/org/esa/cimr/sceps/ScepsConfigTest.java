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
    public void testReadDevalgoL2SeaiceConcentrationLocalConfigFile() {

        final String localalConfigXmlPath =
                Objects.requireNonNull(getClass().
                        getResource("Devalgo_L2_Seaice_Concentration_Local_Configuration_template.xml")).getPath();
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

        // todo: fill config parameters and continue testing here
    }
}