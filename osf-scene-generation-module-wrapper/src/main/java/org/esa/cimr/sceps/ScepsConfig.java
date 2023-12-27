package org.esa.cimr.sceps;

import org.w3c.dom.Document;

import javax.xml.parsers.DocumentBuilder;
import javax.xml.parsers.DocumentBuilderFactory;
import java.io.File;

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
}
