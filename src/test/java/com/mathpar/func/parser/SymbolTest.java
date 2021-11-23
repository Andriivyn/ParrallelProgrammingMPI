package com.mathpar.func.parser;

import com.mathpar.func.F;
import com.mathpar.func.Fname;
import com.mathpar.number.Element;
import com.mathpar.number.NumberR;
import com.mathpar.number.NumberR64;
import com.mathpar.number.Ring;
import com.mathpar.polynom.Polynom;
import org.junit.Test;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import static org.junit.Assert.*;

public class SymbolTest {
    private static void checkFnameString(String expected, String input) {
        // F.ID
        //      input (Fname)
        F f = Parser.getF(input, Ring.ringZxyz);

        assertNotNull(f);
        assertEquals(F.ID, f.name);
        assertEquals(1, f.X.length);
        assertTrue(f.X[0] instanceof Fname);
        assertEquals(expected, ((Fname) f.X[0]).name);
    }

    private static void checkFnameString(String input) {
        // F.ID
        //      input (Fname)
        checkFnameString(input, input);
    }

    @Test
    public void symbSimple() {
        checkFnameString("a");
    }

    @Test
    public void symbLongNameWithUnderscores() {
        checkFnameString("some_long_name");
    }

// XXX: can't have such names
//    @Ignore("To be reviewed")
//    @Test
//    public void symbLongNameWithUnderscoresAndPolynomialWithoutSpaces() {
//        checkFnameString("f_xxx");
//    }

    @Test
    public void symbLongNameWithoutUnderscoresAndPolynomial() {
        checkFnameString("fxxx");
    }

    @Test
    public void symbPolynomWithNumber() {
        checkFnameString("xy1");
    }

    @Test
    public void symbUnderScoreWithRingVars() {
        checkFnameString("x_1");
    }

    @Test
    public void symbMixedNameSymbolWithPolWithoutSpaces() {
        checkFnameString("ax");
    }

    @Test
    public void symbMixedSymbolWithLongPolynomialWithoutSpaces() {
        checkFnameString("xyabc");
    }

    @Test
    public void symbMixedNameSymbolWithPolWithSpaces() {
        // * (F.MULTIPLY)
        //      F.ID
        //          a (Fname)
        //      x (Polynom)

        String s = "a x";
        Ring r = new Ring("Z[x]");
        F f = Parser.getF(s, r);

        assertNotNull(f);
        assertEquals(F.MULTIPLY, f.name);
        assertEquals(2, f.X.length);
        F a = (F) f.X[0];
        assertEquals(F.ID, a.name);
        assertEquals("a", ((Fname) a.X[0]).name);
        assertTrue(f.X[1] instanceof Polynom);
    }

    @Test
    public void symbMixedNamePolWithSymbolWithSpaces() {
        // * (F.MULTIPLY)
        //      x (Polynom)
        //      F.ID
        //          a (Fname)

        String s = "x a";
        Ring r = new Ring("Z[x]");
        F f = Parser.getF(s, r);

        assertNotNull(f);
        assertEquals(F.MULTIPLY, f.name);
        assertEquals(2, f.X.length);
        assertTrue(f.X[0] instanceof Polynom);
        F a = (F) f.X[1];
        assertEquals(F.ID, a.name);
        assertEquals("a", ((Fname) a.X[0]).name);
    }

    @Test
    public void symbMixedNamePolWithSymbolWithoutSpaces() {
        checkFnameString("xzabyc");
    }

    @Test
    public void symbMixedNameStartingWithNumber() {
        // * (F.MULTIPLY)
        //      2 (Polynom)
        //      F.ID
        //          "a" (Fname)
        String s = "2a";
        Ring r = new Ring("Z[x]");
        F f = Parser.getF(s, r);

        assertNotNull("Failed to parse", f);
        assertEquals(F.MULTIPLY, f.name);
        assertEquals(2, f.X.length);
        assertTrue(f.X[0] instanceof Polynom);
        Fname a = (Fname) ((F) f.X[1]).X[0];
        assertEquals("a", a.name);
    }

    @Test
    public void symbMixedNameWithSpaces() {
        // * (F.MULTIPLY)
        //      F.ID
        //          xza (Fname)
        //      F.ID
        //          byc (Fname)
        String s = "xza byc";
        Ring r = new Ring("Z[x, y, z]");
        F f = Parser.getF(s, r);

        assertNotNull(f);
        assertEquals(F.MULTIPLY, f.name);
        assertEquals(2, f.X.length);
        assertEquals("xza", ((Fname) ((F) f.X[0]).X[0]).name);
        assertEquals("byc", ((Fname) ((F) f.X[1]).X[0]).name);
    }

    @Test
    public void symbSpacesAndBackslash() {
        // * (F.MULTIPLY)
        //      F.ID
        //          kg (Fname)
        //      F.ID
        //          \degreeC (Fname)
        String s = "kg \\degreeC";
        Ring r = Ring.ringR64xyzt;
        F f = Parser.getF(s, r);

        assertNotNull(f);
        assertEquals(F.MULTIPLY, f.name);
        assertEquals(2, f.X.length);
        assertEquals("kg", ((Fname) ((F) f.X[0]).X[0]).name);
        assertEquals("\\degreeC", ((Fname) ((F) f.X[1]).X[0]).name);
    }

    @Test
    public void symbHboxEmpty() {
        // hbox (F.HBOX)
        //      (F.ID)
        //          "" (Fname)
        Ring r = new Ring("Z[x]");
        F f = Parser.getF("\\hbox{}", r);

        assertNotNull(f);
        assertEquals(F.HBOX, f.name);
        assertEquals(1, f.X.length);
        F f1 = (F) f.X[0];
        assertEquals(F.ID, f1.name);
        assertEquals("", ((Fname) f1.X[0]).name);
    }

    @Test
    public void symbHbox() {
        // hbox (F.HBOX)
        //      (F.ID)
        //          some words (Fname)
        F f = Parser.getF("\\hbox{some words}", Ring.ringZxyz);

        assertNotNull(f);
        assertEquals(F.HBOX, f.name);
        assertEquals(1, f.X.length);
        F f1 = (F) f.X[0];
        assertEquals(F.ID, f1.name);
        assertEquals("some words", ((Fname) f1.X[0]).name);
    }

    @Test
    public void symbHboxWithSpecialSymbols() {
        // \hbox (F.HBOX)
        //      (F.ID)
        //          Test test c_1 = 1: (Fname)
        String s = "\\hbox{Test test c_1 = 1:}";
        Ring r = new Ring("Z[x, y, z]");
        F f = Parser.getF(s, r);

        assertNotNull(f);
        assertEquals(F.HBOX, f.name);
        assertEquals(1, f.X.length);
        F fId = (F) f.X[0];
        assertEquals(F.ID, fId.name);
        assertEquals(1, fId.X.length);
        assertEquals("Test test c_1 = 1:", ((Fname) fId.X[0]).name);
    }

    @Test
    public void symbGreekDoesnotConflictWithFunctions() {
        checkFnameString("\\delta");
    }

    @Test
    public void symbPolSumWithSymbol() {
        // + (F.ADD)
        //      x (Polynom)
        //      y (Polynom)
        //      F.ID
        //          b (Fname)
        //      F.ID
        //          c (Fname)
        String s = "x+y+b+c";
        Ring r = new Ring("Z[x,y,z]");
        F f = Parser.getF(s, r);

        assertNotNull(f);
        assertEquals(F.ADD, f.name);
        assertEquals(4, f.X.length);
        assertTrue(f.X[0] instanceof Polynom);
        assertTrue(f.X[1] instanceof Polynom);
        assertTrue(f.X[2] instanceof F);
        assertEquals(F.ID, ((F) f.X[2]).name);
        assertTrue(f.X[3] instanceof F);
        assertEquals(F.ID, ((F) f.X[3]).name);
    }

    @Test
    public void symbLongNameWithUnderscorsWithNulls() {
        checkFnameString("f_1_{,j,,l,}", "f_1_{, j, , l, }");
    }

    @Test
    public void symbLongNameWithUnderscorsWithNullsWithNonCommuteList() {
        // F.ID
        //      f_1 (Fname)
        //          null
        //          F.VECTORS
        //              null
        //              j
        //              null
        //              k
        //              null
        String s = "f_1_{, j, , l, }";
        Ring r = new Ring("R64[x]");
        ArrayList<Fname> nonCommut = new ArrayList<>();
        nonCommut.add(new Fname("f_1"));
        F f = Parser.getF(s, r, Collections.<String>emptyList(), nonCommut);

        assertNotNull(f);
        assertEquals(F.ID, f.name);
        assertEquals(1, f.X.length);
        assertTrue(f.X[0] instanceof Fname);
        Fname fname = (Fname) f.X[0];
        assertEquals("f_1", fname.name);
        assertNull(fname.X);
        Element[] expectedIndices = new Element[]{
                null, new F(new Fname("j")), null, new F(new Fname("l")), null
        };
        assertArrayEquals(expectedIndices, fname.lowerIndices());
    }

    @Test
    public void symbNameWithArgumentsNullQuestion() {
        // Spaces are gone here.
        checkFnameString("A_{?,j,?,l,?}", "A_{?, j, ?, l, ?}");
    }

    @Test
    public void symbNameWithArgumentsNullQuestionWithNonCommuteList() {
        // F.ID
        //      A (Fname)
        //          null
        //          F.VECTORS
        //              null
        //              j
        //              null
        //              l
        //              null
        String s = "A_{?, j, ?, l, ?}";
        Ring r = new Ring("R64[x]");
        ArrayList<Fname> nonCommut = new ArrayList<>();
        nonCommut.add(new Fname("A"));
        F f = Parser.getF(s, r, Collections.<String>emptyList(), nonCommut);

        assertNotNull(f);
        assertEquals(F.ID, f.name);
        assertEquals(1, f.X.length);
        assertTrue(f.X[0] instanceof Fname);
        Fname fname = (Fname) f.X[0];
        assertEquals("A", fname.name);
        assertNull(fname.X);
        Element[] expectedIndices = new Element[]{
                null, new F(new Fname("j")), null, new F(new Fname("l")), null
        };
        assertArrayEquals(expectedIndices, fname.lowerIndices());
    }

    @Test
    public void symbUnaryMinus() {
        // F.MULTIPLY
        //      -1 (R64)
        //      F.SUBTRACT
        //          F.MULTIPLY
        //              -1 (R64)
        //              F.ID
        //                  A (Fname)
        //          F.ID
        //              B (Fname)
        String s = "-(-A-B)";
        Ring r = new Ring("R64[x]");
        F f = Parser.getF(s, r);

        assertNotNull(f);
        assertEquals(F.MULTIPLY, f.name);
        assertTrue(new Polynom(new NumberR64("-1")).equals(f.X[0], r));
        F f1 = (F) f.X[1];
        assertEquals(F.SUBTRACT, f1.name);
        F f2 = (F) f1.X[0];
        assertEquals(F.MULTIPLY, f2.name);
        assertTrue(new Polynom(new NumberR64("-1")).equals(f2.X[0], r));
        F f4 = (F) f2.X[1];
        assertEquals(F.ID, f4.name);
        Fname fname1 = (Fname) f4.X[0];
        assertEquals("A", fname1.name);
        F f3 = (F) f1.X[1];
        assertEquals(F.ID, f3.name);
        Fname fname2 = (Fname) f3.X[0];
        assertEquals("B", fname2.name);
    }

    @Test
    public void symbStringLiteral() {
        // F.ID
        //      Hello, world. (Fname)
        String s = "'Hello, world.'";
        Ring r = new Ring("R64[x]");
        F f = Parser.getF(s, r);

        assertEquals(F.ID, f.name);
        Fname fname = (Fname) f.X[0];
        assertEquals("Hello, world.", fname.name);
    }

    @Test
    public void symbDependentOnVars() {
        // F.ID
        //      a (Fname)
        //          [0, 1]
        String s = "a(x, y)";
        Ring r = new Ring("R[x, y, z]");
        F f = Parser.getF(s, r);

        assertEquals(F.ID, f.name);
        Fname fname = (Fname) f.X[0];
        assertEquals("a", fname.name);
        assertNull(fname.X);
        Element[] expectedIndices = new Element[]{
                Polynom.polynom_zero(NumberR.ONE), new Polynom(new NumberR("1"))
        };
        assertArrayEquals(expectedIndices, fname.varsIndices());
     }

    @Test
    public void symbWithUnderscores() {
        checkFnameString("ab_{b,c}", "ab_{b,c}");
    }

    @Test
    public void symbWithBracesSimple() {
        checkFnameString("a", "(a)");
    }

    @Test
    public void symbWithBracesWithIndices() {
        checkFnameString("b_{k+1}", "(b_{k+1})");
    }
}
