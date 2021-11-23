package com.mathpar.number;

import static org.junit.Assert.*;
import org.junit.Test;

public class NumberZpTest {

    private static final Ring r10 = new Ring("Zp[]");

    static {
        r10.setMOD(new NumberZ("10"));
    }

    @Test
    public void constructZeroModuloModHasSignum0AndEmptyMagnitude() {
        Element n10 = new NumberZp("10", r10);
        assertEquals(0, n10.signum());
        assertArrayEquals(new int[]{}, ((NumberZp) n10).mag);
    }

    @Test
    public void constructPositiveHasSignum1() {
        Element n15 = new NumberZp("15", r10);
        assertEquals(1, n15.signum());
        assertArrayEquals(new int[]{5}, ((NumberZp) n15).mag);
    }

    @Test
    public void constructNegativeHasSignumMinus1AndPositiveMagnitude() {
        Element n_15 = new NumberZp("-15", r10);
        assertEquals(-1, n_15.signum());
        assertArrayEquals(new int[]{5}, ((NumberZp) n_15).mag);
    }

    @Test
    public void addLowerTypeElementProducesZp() {
        Element n4 = new NumberZp("19", r10).add(new NumberZ("5"), r10);
        assertTrue(n4 instanceof NumberZ);
        assertArrayEquals(new int[]{14}, ((NumberZ) n4).mag);
    }

    @Test
    public void addSameTypeElementProducesZp() {
        Element n1 = new NumberZp("19", r10).add(new NumberZp("12", r10), r10);
        assertTrue(n1 instanceof NumberZp);
        assertArrayEquals(new int[]{1}, ((NumberZp) n1).mag);
    }

    @Test
    public void addHigherTypeElementProducesHigherType() {
        Element n11 = new NumberZp("19", r10).add(new NumberR("10.2"), r10);
        assertTrue(n11 instanceof NumberR);
    }
}
