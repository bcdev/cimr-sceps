package org.esa.cimr.sceps.util;

import org.jdom2.Document;
import org.jdom2.Element;
import org.jdom2.JDOMException;
import org.jdom2.input.SAXBuilder;

import java.io.IOException;
import java.io.InputStream;

public class JDomHelper {
    public static void addElement(String tag, long value, Element element) {
        addElement(tag, String.valueOf(value), element);
    }

    public static void addElement(String tag, float value, Element element) {
        addElement(tag, String.valueOf(value), element);
    }

    public static void addElement(String tag, double value, Element element) {
        addElement(tag, String.valueOf(value), element);
    }

    public static void addElement(String tag, boolean value, Element element) {
        addElement(tag, String.valueOf(value), element);
    }

    public static void addElement(String tag, String content, Element element) {
        element.addContent(createElement(tag, content));
    }

    public static Element createElement(String tag, String content) {
        Element element = new Element(tag);
        element.addContent(content);
        return element;
    }

    public static Document parse(InputStream stream) throws IOException, JDOMException {
        SAXBuilder builder = new SAXBuilder();
        builder.setExpandEntities(false);
        builder.setIgnoringElementContentWhitespace(false);
        builder.setValidation(false);
        return builder.build(stream);
    }
}
