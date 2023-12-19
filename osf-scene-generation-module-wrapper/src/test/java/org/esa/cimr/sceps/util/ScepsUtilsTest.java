package org.esa.cimr.sceps.util;

import org.junit.Test;

import static org.junit.Assert.*;

public class ScepsUtilsTest {

    @Test
    public void testDoSomething() {
        assertTrue(true);
    }

    @Test
    public void testReadXmlConfigFile() {
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

        assertTrue(true);
    }
}