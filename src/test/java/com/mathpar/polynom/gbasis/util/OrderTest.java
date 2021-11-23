package com.mathpar.polynom.gbasis.util;

import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;
import com.mathpar.polynom.Term;
import static org.junit.Assert.*;

/**
 *
 * @author ivan
 */
public class OrderTest {

    public OrderTest() {
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

    @Test
    public void equal() {
        Term t1 = new Term(1, 0, 2); // z^2x
        Term t2 = new Term(1, 0, 2); // z^2x
        assertEquals(0, Order.REVLEX.compare(t1, t2));
        assertEquals(0, Order.LEX.compare(t1, t2));
        assertEquals(0, Order.DEGREVLEX.compare(t1, t2));
    }

    @Test
    public void revlex() {
        Term t1 = new Term(1, 0, 2); // z^2x
        Term t2 = new Term(0, 0, 3); // z^3
        assertTrue(Order.REVLEX.compare(t1, t2) < 0);
        assertTrue(Order.REVLEX.compare(t2, t1) > 0);
    }

    @Test
    public void lex() {
        Term t1 = new Term(1, 0, 2); // xz^2
        Term t2 = new Term(0, 0, 4); // z^4
        // now x > y > z, xz^2 > z^4
        assertTrue(Order.LEX.compare(t1, t2) > 0);
        assertTrue(Order.LEX.compare(t2, t1) < 0);
    }

    @Test
    public void degrevlex() {
        Term t1 = new Term(1, 0, 2); // z^2x
        Term t2 = new Term(3, 0, 1); // zx^3
        // deg(t2) > deg(t1)
        assertTrue(Order.DEGREVLEX.compare(t1, t2) < 0);
        assertTrue(Order.DEGREVLEX.compare(t2, t1) > 0);
        t1 = new Term(1, 2, 2); // z^2y^2x
        t2 = new Term(2, 1, 2); // z^2yx^2
        // deg(t1) == deg(t2), but t1 > (revlex) t2
        assertTrue(Order.DEGREVLEX.compare(t1, t2) > 0);
        assertTrue(Order.DEGREVLEX.compare(t2, t1) < 0);
    }
}
