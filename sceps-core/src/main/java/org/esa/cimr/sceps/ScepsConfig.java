package org.esa.cimr.sceps;

import org.w3c.dom.Document;
import org.w3c.dom.Element;
import org.w3c.dom.Node;
import org.w3c.dom.NodeList;

import javax.xml.parsers.DocumentBuilder;
import javax.xml.parsers.DocumentBuilderFactory;
import java.io.File;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;

/**
 * Utility class for CIMR SCEPS
 *
 * @author olafd
 */
public class ScepsConfig {

    /**
     * dummy
     *
     */
    public static void doSomething() {
        System.out.println();
    }

    public static Document readXMLDocumentFromFile(String fileNameWithPath) throws Exception {

        //Get Document Builder
        DocumentBuilderFactory factory = DocumentBuilderFactory.newInstance();
        DocumentBuilder builder = factory.newDocumentBuilder();

        //Build Document
        Document document = builder.parse(new File(fileNameWithPath));

        //Normalize the XML Structure; It's just too important !!
        document.getDocumentElement().normalize();

        return document;
    }

    public static NodeList getDocumentElementsByTagName(Document document, String tagName) {
        return document.getElementsByTagName(tagName);
    }

    public static Element getDocumentRootElement(Document document) {
        return document.getDocumentElement();
    }

    public static String getDocumentElementTextItemByName(Document document, String tagName, String elementName) {
        NodeList nList = getDocumentElementsByTagName(document, tagName);

        for (int i = 0; i < nList.getLength(); i++) {
            Node node = nList.item(i);

            Element elem = (Element) node;
            if (elem.getFirstChild() != null && elem.getAttribute("name").equals(elementName)) {
                return elem.getFirstChild().getTextContent();
            }
        }
        return null;
    }

}
