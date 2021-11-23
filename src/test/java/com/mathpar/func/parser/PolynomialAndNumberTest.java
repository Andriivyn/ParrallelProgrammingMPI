package com.mathpar.func.parser;

import com.mathpar.func.F;
import com.mathpar.number.*;
import org.junit.Test;
import com.mathpar.polynom.Polynom;

import static org.junit.Assert.*;

public class PolynomialAndNumberTest {
    private static Polynom checkPolAtRoot(Element f) {
        assertNotNull(f);
        assertTrue(f instanceof F);
        F f_ = (F) f;
        assertEquals(F.ID, f_.name);
        assertTrue(f_.X.length == 1);
        assertTrue(f_.X[0] instanceof Polynom);
        return (Polynom) f_.X[0];
    }

    private static void checkNumber(Ring r, String s, Element el) {
        F f = Parser.getF(s, r);
        Polynom pol = checkPolAtRoot(f);
        assertTrue(el.equals(pol.coeffs[0], r));
    }

    @Test
    public void polNumbersAtDifferentRings() {
        checkNumber(new Ring("R64[]"), "2", new NumberR64("2"));
        checkNumber(new Ring("R[]"), "2", new NumberR("2"));
        checkNumber(new Ring("Z[]"), "2", new NumberZ("2"));
        checkNumber(new Ring("R[]Z[]"), "2", new NumberZ("2"));
        checkNumber(new Ring("R[]Z64[]Z[]"), "2", new NumberZ64("2"));
        checkNumber(new Ring("R[]"), "2.5", new NumberR("2.5"));
        checkNumber(new Ring("Z[]"), "2.5", new NumberZ("2"));
        checkNumber(new Ring("R[]Z[]"), "2.5", new NumberR("2.5"));
        checkNumber(new Ring("Z[]R[]"), "2.5", new NumberZ("2"));
        checkNumber(new Ring("Q[]"), "2", new NumberZ(2));
        checkNumber(new Ring("Q[]"), "1/2", new Fraction(new NumberZ(1), new NumberZ(2)));
        checkNumber(new Ring("Q[]"), "2.67", new NumberZ(2));
        checkNumber(new Ring("C64[]"), "2",
                new Complex(new NumberR64("2"), new Ring("C64[]")));
        checkNumber(new Ring("C[]"), "2",
                new Complex(new NumberR("2"), new Ring("C[]")));
        checkNumber(new Ring("ZMaxMin[]"), "2",
                new NumberZMaxMin("2", new Ring("ZMaxMin[]")));
        checkNumber(new Ring("ZMaxPlus[]"), "2",
                new NumberZMaxPlus("2", new Ring("ZMaxPlus[]")));
        checkNumber(new Ring("ZMinMax[]"), "2",
                new NumberZMinMax("2", new Ring("ZMinMax[]")));
        checkNumber(new Ring("ZMinPlus[]"), "2",
                new NumberZMinPlus("2", new Ring("ZMinPlus[]")));
        checkNumber(new Ring("RMaxPlus[]"), "2",
                new NumberRMaxPlus("2", new Ring("RMaxPlus[]")));
        checkNumber(new Ring("RMinPlus[]"), "2",
                new NumberRMinPlus("2", new Ring("RMinPlus[]")));
        checkNumber(new Ring("RMaxMin[]"), "2",
                new NumberRMaxMin("2", new Ring("RMaxMin[]")));
        checkNumber(new Ring("RMinMax[]"), "2",
                new NumberRMinMax("2", new Ring("RMinMax[]")));
        checkNumber(new Ring("RMaxMult[]"), "2",
                new NumberRMaxMult("2", new Ring("RMaxMult[]")));
        checkNumber(new Ring("RMinMult[]"), "2",
                new NumberRMinMult("2", new Ring("RMinMult[]")));
        checkNumber(new Ring("R64MaxPlus[]"), "2",
                new NumberR64MaxPlus("2", new Ring("R64MaxPlus[]")));
        checkNumber(new Ring("R64MinPlus[]"), "2",
                new NumberR64MinPlus("2", new Ring("R64MinPlus[]")));
        checkNumber(new Ring("R64MaxMin[]"), "2",
                new NumberR64MinPlus("2", new Ring("R64MaxMin[]")));
        checkNumber(new Ring("R64MinMax[]"), "2",
                new NumberR64MinMax("2", new Ring("R64MinMax[]")));
        checkNumber(new Ring("R64MaxMult[]"), "2",
                new NumberR64MaxMult("2", new Ring("R64MaxMult[]")));
        checkNumber(new Ring("R64MinMult[]"), "2",
                new NumberR64MinMult("2", new Ring("R64MinMult[]")));
    }

    @Test
    public void polNegativeNumbersAtDifferentRings() {
        checkNumber(new Ring("R64[]"), "-2", new NumberR64("-2"));
        checkNumber(new Ring("R[]"), "-2", new NumberR("-2"));
        checkNumber(new Ring("Z[]"), "-2", new NumberZ("-2"));
        checkNumber(new Ring("R[]Z[]"), "-2", new NumberZ("-2"));
        checkNumber(new Ring("R[]Z64[]Z[]"), "-2", new NumberZ64("-2"));
        checkNumber(new Ring("R[]"), "-2.5", new NumberR("-2.5"));
        checkNumber(new Ring("Z[]"), "-2.5", new NumberZ("-2"));
        checkNumber(new Ring("R[]Z[]"), "-2.5", new NumberR("-2.5"));
        checkNumber(new Ring("Z[]R[]"), "-2.5", new NumberZ("-2"));
        checkNumber(new Ring("Q[]"), "-2", new NumberZ(-2));
        checkNumber(new Ring("Q[]"), "-1/2", new Fraction(new NumberZ(-1), new NumberZ(2)));
        checkNumber(new Ring("Q[]"), "-2.67", new NumberZ(-2));
        checkNumber(new Ring("C64[]"), "-2",
                new Complex(new NumberR64("-2"), new Ring("C64[]")));
        checkNumber(new Ring("C[]"), "-2",
                new Complex(new NumberR("-2"), new Ring("C[]")));
        checkNumber(new Ring("ZMaxMin[]"), "-2",
                new NumberZMaxMin("-2", new Ring("ZMaxMin[]")));
        checkNumber(new Ring("ZMaxPlus[]"), "-2",
                new NumberZMaxPlus("-2", new Ring("ZMaxPlus[]")));
        checkNumber(new Ring("ZMinMax[]"), "-2",
                new NumberZMinMax("-2", new Ring("ZMinMax[]")));
        checkNumber(new Ring("ZMinPlus[]"), "-2",
                new NumberZMinPlus("-2", new Ring("ZMinPlus[]")));
    }

    @Test
    public void polFractionsInCq() {
        // + (F.ADD)
        //      * (F.MULTIPLY)
        //          1/2 (Polynom)
        //          \i (Polynom)
        //          x (Polynom)
        //      * (F.MULTIPLY)
        //          3/4 (Polynom)
        //          x (Polynom)
        Ring r = new Ring("CQ[x]");
        String s = "1/2 \\i x + 3/4 x";
        F f = Parser.getF(s, r);
        assertEquals(F.ADD, f.name);
        assertEquals(2, f.X.length);
        F mult1 = (F) f.X[0];
        F mult2 = (F) f.X[1];
        assertEquals(F.MULTIPLY, mult1.name);
        assertEquals(3, mult1.X.length);
        assertTrue(((Polynom) mult1.X[0]).coeffs[0] instanceof Fraction); // 1/2 is Fraction
        assertEquals(F.MULTIPLY, mult2.name);
        assertEquals(2, mult2.X.length);
        assertTrue(((Polynom) mult2.X[0]).coeffs[0] instanceof Fraction); // 3/4 is Fraction
    }

    @Test
    public void polFractionInShorExpr() {
        Ring r = new Ring("Q[]");
        String s = "\\value(1/2)";
        F f = Parser.getF(s, r);

        Polynom pol = (Polynom) f.X[0];
        assertEquals(1, pol.coeffs.length);
        assertTrue(pol.coeffs[0].equals(new Fraction(new NumberZ("1"), new NumberZ("2")), r));
    }

    @Test
    public void polOnlyOneFloatNumber() {
        Ring r = new Ring("R[x]");
        String s = "1.1";
        F f = Parser.getF(s, r);

        Polynom pol = checkPolAtRoot(f);
        assertEquals(0, pol.powers.length);
        assertTrue(new NumberR("1.1").equals(pol.coeffs[0], r));
    }

    @Test
    public void polSingleVar() {
        String s = "y";
        Ring r = new Ring("R64[x, y]");
        F f = Parser.getF(s, r);

        Polynom pol = checkPolAtRoot(f);
        assertArrayEquals(new int[] {0, 1}, pol.powers);
    }

    @Test
    public void polWithGreekVariable() {
        String s = "\\alpha^2";
        Ring r = new Ring("R64[\\alpha]");
        F f = Parser.getF(s, r);

        Polynom pol = checkPolAtRoot(f);
        assertArrayEquals(new int[] {2}, pol.powers);
    }

    @Test
    public void polUnusedHighestVarsShouldBeTruncated() {
        String s = "x";
        Ring r = new Ring("R64[x, y, z]");
        F f = Parser.getF(s, r);

        Polynom pol = checkPolAtRoot(f);
        assertArrayEquals(new int[] {1}, pol.powers);
    }

    @Test
    public void polShouldNotExpand() {
        // ^ (F.POW)
        //      +
        //          x^2 (Polynom)
        //          1.0 (Polynom)
        //      2.0 (Polynom)
        String s = "(x^2 + 1)^2";
        Ring r = new Ring("R64[x]");
        F f = Parser.getF(s, r);

        assertEquals(F.POW, f.name);
        assertEquals(2, f.X.length);
        assertTrue(f.X[0] instanceof F);
        F plus = (F) f.X[0];
        assertEquals(F.ADD, plus.name);
        assertEquals(2, plus.X.length);
        assertTrue(plus.X[0] instanceof Polynom);
        assertTrue(plus.X[1] instanceof Polynom);
        assertTrue(f.X[1] instanceof Polynom);
    }

    @Test
    public void polSimpleInQ() {
        // * (F.MULTIPLY)
        //      1/2 (Polynom)
        //      y^2 (Polynom)
        String s = "1/2y^2";
        Ring r = new Ring("Q[x, y]");
        F f = Parser.getF(s, r);
        assertEquals(F.MULTIPLY, f.name);
        assertEquals(2, f.X.length);
        Polynom oneHalf = (Polynom) f.X[0];
        Polynom ySquare = (Polynom) f.X[1];
        assertTrue(new Fraction(new NumberZ(1), new NumberZ(2)).
                equals(oneHalf.coeffs[0], r));
        assertArrayEquals(new int[] {0, 2}, ySquare.powers);
    }

    @Test
    public void polComplexSimple() {
        // + (F.ADD)
        //      * (F.MULTIPLY)
        //          5.2 (Polynom)
        //          \i (Polynom)
        //          y (Polynom)
        //      1 (Polynom)
        String s = "5.2\\iy + 1";
        Ring r = new Ring("C64[x, y]");
        F f = Parser.getF(s, r);

        assertEquals(F.ADD, f.name);
        assertEquals(2, f.X.length);
        assertTrue(f.X[0] instanceof F);
        F mult = (F) f.X[0];
        assertEquals(F.MULTIPLY, mult.name);
        assertEquals(3, mult.X.length);
        assertTrue(f.X[1] instanceof Polynom);
    }

    @Test
    public void polComplex() {
        // + (F.ADD)
        //      \i (Polynom)
        //      * (F.MULTIPLY)
        //          y^2 (Polynom)
        //          \i (Polynom)
        //          x (Polynom)
        String s = "\\i + y^2\\ix";
        Ring r = new Ring("C64[x, y]");
        F f = Parser.getF(s, r);
        assertTrue(f instanceof F);
        assertEquals(F.ADD, f.name);
        assertEquals(2, f.X.length);
        assertTrue(f.X[0] instanceof Polynom);
        Polynom i = (Polynom) f.X[0];
        assertTrue(((Complex) i.coeffs[0]).equals(new Complex(0.0, 1.0), r));
        assertTrue(f.X[1] instanceof F);
        F mult = (F) f.X[1];
        assertEquals(F.MULTIPLY, mult.name);
        assertEquals(3, mult.X.length);
        assertTrue(mult.X[0] instanceof Polynom);
        assertTrue(mult.X[1] instanceof Polynom);
        assertTrue(mult.X[2] instanceof Polynom);
    }

    @Test
    public void polBig() {
        //  / (F.DIVIDE)
        //      ^ (F.POW)
        //          * (F.MULTIPLY)
        //              2.0 (Polynom)
        //              y (Polynom)
        //          3.0 (Polynom)
        //      ^ (F.POW)
        //          * (F.MULTIPLY)
        //              z^2 (Polynom)
        //              y   (Polynom)
        //              x^2 (Polynom)
        //          2.00 (Polynom)
        String s = "(2 y)^3 / (z^2yx^2)^2";
        Ring r = new Ring("R64[x, y, z]");
        try {
            F f = Parser.getF(s, r);

            assertEquals(F.DIVIDE, f.name);
            assertEquals(2, f.X.length);
            assertTrue(f.X[0] instanceof F);
            F pow1 = (F) f.X[0];
            assertEquals(F.POW, pow1.name);
            assertTrue(pow1.X[0] instanceof F);
            F mult1 = (F) pow1.X[0];
            assertEquals(F.MULTIPLY, mult1.name);
            assertEquals(2, mult1.X.length);
            assertTrue(pow1.X[1] instanceof Polynom);
            Polynom p12 = (Polynom) pow1.X[1];
            assertEquals(0, p12.powers.length);
            assertTrue(new NumberR64("3.0").equals(p12.coeffs[0], r));
            assertTrue(f.X[1] instanceof F);
            F pow2 = (F) f.X[1];
            assertTrue(pow2.X[0] instanceof F);
            F mult2 = (F) pow2.X[0];
            assertEquals(F.MULTIPLY, mult2.name);
            assertEquals(3, mult2.X.length);
            assertTrue(pow2.X[1] instanceof Polynom);
            Polynom p22 = (Polynom) pow2.X[1];
            assertEquals(0, p22.powers.length);
            assertTrue(new NumberR64("2.0").equals(p22.coeffs[0], r));
        } catch (ParserException pe) {
            pe.getCause().printStackTrace(System.out);
        }
    }

    @Test
    public void polWithWeirdVariables() {
        // + (F.ADD)
        //      C$_1 (Polynom)
        //      C$_0^2 (Polynom)
        String s = "C$_1 + C$_0^2";
        Ring r = new Ring("Z[x, y, C$_0, C$_1]");

        F f = Parser.getF(s, r);
        assertTrue(f instanceof F);
        assertEquals(F.ADD, f.name);
        assertTrue(f.X[0] instanceof Polynom);
        Polynom p0 = (Polynom) f.X[0];
        assertTrue(f.X[1] instanceof Polynom);
        Polynom p1 = (Polynom) f.X[1];
        assertArrayEquals(new int[] {0, 0, 0, 1}, p0.powers);
        assertArrayEquals(new int[] {0, 0, 2}, p1.powers);
        assertTrue(NumberZ.ONE.equals(p0.coeffs[0], r));
        assertTrue(NumberZ.ONE.equals(p1.coeffs[0], r));
    }

    @Test
    public void polWithWeirdVariablesAndUnderscore() {
        String s = "$c_1";
        Ring r = new Ring("Z[$c_1]");

        F f = Parser.getF(s, r);
        Polynom pol = checkPolAtRoot(f);
        assertArrayEquals(new int[] {1}, pol.powers);
        assertTrue(NumberZ.ONE.equals(pol.coeffs[0], r));
    }

    @Test
    public void polInZpPowerGreaterThanP() {
        // + (F.ADD)
        //      x^5 (Polynom) // 5 stays as it power
        //      1 (Polynom)   // 13 = 1 (mod 3)
        String s = "x^5 + 13";
        Ring r = new Ring("Zp[x]");
        r.setMOD(new NumberZ("3"));

        F f = Parser.getF(s, r);
        assertEquals(F.ADD, f.name);
        assertEquals(2, f.X.length);
        Polynom pol1 = (Polynom) f.X[0], pol2 = (Polynom) f.X[1];
        assertArrayEquals(new int[] {5}, pol1.powers);
        assertTrue(new NumberZp("1", r).equals(pol2.coeffs[0], r));
    }

    @Test
    public void polNumberDualMinusToUnary() {
        // - (F.SUBTRACT)
        //      1 (Polynom)
        //      -1 (Polynom)
        String s = "1 - -1";
        Ring r = new Ring("Z[x]");

        F f = Parser.getF(s, r);
        assertEquals(F.SUBTRACT, f.name);
        assertEquals(2, f.X.length);
        Polynom one = (Polynom) f.X[0], minusOne = (Polynom) f.X[1];
        assertTrue(!one.isNegative());
        assertTrue(minusOne.isNegative());
    }

    @Test
    public void polMonomialShouldContainsOnlyOneMultiply() {
        String s = " x^1^3 y^2^3^4 z^3 t^4 x^5 y^6";
        Ring r = new Ring("Z[x, y, z, t]");

        F f = Parser.getF(s, r);

        assertEquals(F.MULTIPLY, f.name);
        assertEquals(6, f.X.length);
    }
}
