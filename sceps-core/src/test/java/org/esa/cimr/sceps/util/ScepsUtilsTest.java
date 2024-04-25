package org.esa.cimr.sceps.util;

import org.junit.Test;

import java.util.ArrayList;
import java.util.List;

import static org.junit.Assert.*;

public class ScepsUtilsTest {

    @Test
    public void testDoSomething() {
        assertTrue(true);
    }

    @Test
    public void testClpInputsJava2Matlab() {
        List<String> javaList = new ArrayList<>();
        javaList.add("abc");
        javaList.add("def");
        javaList.add("ghi");

        final String matlabInputString = ScepsUtils.clpInputsJava2Matlab(javaList);
        assertNotNull(matlabInputString);
        final String expected = "abc,def,ghi";
        assertEquals(expected, matlabInputString);
    }
}