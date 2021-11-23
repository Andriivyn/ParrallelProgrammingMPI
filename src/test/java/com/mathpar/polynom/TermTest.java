package com.mathpar.polynom;

import com.mathpar.number.Ring;
import org.junit.Before;
import org.junit.Test;

import static org.junit.Assert.*;

/**
 * @author ivan
 */
public class TermTest {

    private Ring rZxyz;

    public TermTest() {
    }

    @Before
    public void initialize() {
        rZxyz = new Ring("Z[x, y, z]");
    }

    @Test
    public void size() {
        Term t = new Term(1, 1, 1);
        assertEquals(3, t.size());
    }

    @Test
    public void multiplyByTerm() {
        Term t = new Term(1, 1, 1);
        assertEquals(new Term(1, 2, 1), t.multiply(new Term(0, 1)));
    }

    @Test
    public void multiplyByPolynom() {
        Term t = new Term(1, 1, 1);
        assertTrue(new Polynom("x^2yz + xy^2z + xyz^2", rZxyz).equals(
                t.multiply(new Polynom("x + y + z", rZxyz)), rZxyz));
    }

    @Test
    public void div() {
        Term t = new Term(1, 1, 1);
        assertEquals(new Term(1, 0, 1), t.div(new Term(0, 1)));
    }

    @Test
    public void trim() {
        Term t = new Term(1, 1, 0);
        assertEquals(new Term(1, 1), t.trim());
    }
}
