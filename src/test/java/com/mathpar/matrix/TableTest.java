package com.mathpar.matrix;

import com.mathpar.number.NumberZ;
import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;

import static com.mathpar.matrix.Table.fromString;
import static org.junit.Assert.*;

public class TableTest {
    public TableTest() {
    }

    @BeforeClass
    public static void setUpClass() {
    }

    @AfterClass
    public static void tearDownClass() {
    }

    @Before
    public void setUp() {
    }

    @After
    public void tearDown() {
    }

    /**
     * Test of fromString method, of class Table.
     */
    @Test
    public void testFromStringHeader() throws Exception {
        NumberZ z = NumberZ.MINUS_ONE;
        Table table = fromString("t, с	U, В\n"
                + "0.1	-0.02\n"
                + "0.3	-0.02\n"
                + "0.5	-0.02\n"
                + "0.7	-0.02\n"
                + "0.9	0");
        assertArrayEquals(new String[] {"t (с)", "U (В)"}, table.axesSign);
    }
}
