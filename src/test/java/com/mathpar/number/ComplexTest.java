package com.mathpar.number;

import org.junit.Test;
import static org.junit.Assert.*;

public class ComplexTest {
    @Test
    public void testAbsSquareIsReal() {
        Ring r = new Ring("R64[]");
        r.MachineEpsilonR = new NumberR("0.001");
        Complex c = new Complex(new NumberR("1.61"), new NumberR("-0.78"));
        Element absSquare = c.absSquare(r);

        assertTrue(absSquare instanceof NumberR);
        assertTrue(new NumberR("3.2005").equals(absSquare, r));
    }
}
