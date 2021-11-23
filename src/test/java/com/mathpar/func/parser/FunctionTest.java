package com.mathpar.func.parser;

import com.mathpar.func.F;
import com.mathpar.func.Fname;
import java.util.Arrays;
import java.util.Collections;

import com.mathpar.number.Element;
import com.mathpar.number.NumberR64;
import com.mathpar.number.NumberZ;
import com.mathpar.number.Ring;
import org.junit.Test;
import com.mathpar.polynom.Polynom;

import static org.junit.Assert.*;

public class FunctionTest {
    public static void assertFuncNonId(Element f) {
        assertNotNull(f);
        assertTrue(f instanceof F);
        assertTrue(((F) f).name > F.ID);
    }

    @Test
    public void funcUserDefinedFunction() {
        // myFunc (F.MAX_F_NUMB + indexOf(myFunc1) [1000 + 0])
        //      a (Fname)
        //      b (Fname)
        Ring r = new Ring("R[x]");
        F f = Parser.getF("\\myFunc1(a, b)", r, Arrays.asList("myFunc1", "myFunc2"), Collections.<Fname>emptyList());
        assertFuncNonId(f);
        assertEquals(F.MAX_F_NUMB + 0, f.name);
        assertEquals(2, f.X.length);
        F f1 = (F) f.X[0];
        F f2 = (F) f.X[1];
        assertEquals("$2a", ((Fname) f1.X[0]).name);
        assertEquals("$2b", ((Fname) f2.X[0]).name);
    }

    @Test
    public void funcUserDefinedFunctionDoesnotConflictWithBuiltin() {
        Ring r = Ring.ringR64xyzt;
        F f = Parser.getF("\\print()", r, Arrays.asList("p"), Collections.<Fname>emptyList());
        assertFuncNonId(f);
        assertEquals(F.PRINT, f.name);
    }

    @Test
    public void funcBoolOperationsSimple() {
        Ring r = new Ring("R[x]");
        F f;
        f = Parser.getF("A == B", r);
        assertFuncNonId(f);
        assertEquals(F.B_EQ, f.name);
        assertEquals(2, f.X.length);
        f = Parser.getF("A < B", r);
        assertFuncNonId(f);
        assertEquals(F.B_LESS, f.name);
        assertEquals(2, f.X.length);
        f = Parser.getF("A > B", r);
        assertFuncNonId(f);
        assertEquals(F.B_GT, f.name);
        assertEquals(2, f.X.length);
        f = Parser.getF("A \\& B", r);
        assertFuncNonId(f);
        assertEquals(F.B_AND, f.name);
        assertEquals(2, f.X.length);
        f = Parser.getF("A \\le B", r);
        assertFuncNonId(f);
        assertEquals(F.B_LE, f.name);
        assertEquals(2, f.X.length);
        f = Parser.getF("A \\ge B", r);
        assertFuncNonId(f);
        assertEquals(F.B_GE, f.name);
        assertEquals(2, f.X.length);
        f = Parser.getF("A \\ne B", r);
        assertFuncNonId(f);
        assertEquals(F.B_NE, f.name);
        assertEquals(2, f.X.length);
        f = Parser.getF("A \\lor B", r);
        assertFuncNonId(f);
        assertEquals(F.B_OR, f.name);
        assertEquals(2, f.X.length);
        f = Parser.getF("\\neg B", r);
        assertFuncNonId(f);
        assertEquals(F.B_NOT, f.name);
        assertEquals(1, f.X.length);
    }

    @Test
    public void funcBooleanHasCorrectPriority() {
        // > (F.B_GT)
        //      x^2 (Polynom)
        //      x^3 (Polynom)
        String s = "x^2 > x^3";
        Ring r = new Ring("R64[x]");
        F f = Parser.getF(s, r);
        assertFuncNonId(f);
        assertEquals(F.B_GT, f.name);
        assertEquals(2, f.X.length);
    }

    @Test
    public void funcBooleanHasCorrectPriority2() {
        // < (F.B_LT)
        //      + (F.ADD)
        //          a (Fname)
        //          b (Fname)
        //      + (F.ADD)
        //          x (Polynom)
        //          y (Polynom)
        String s = "a + b < x + y";
        Ring r = new Ring("Z[x, y]");
        F f = Parser.getF(s, r);
        assertFuncNonId(f);
        assertEquals(2, f.X.length);
        assertTrue(f.X[0] instanceof F);
        assertTrue(f.X[1] instanceof F);
        assertEquals((long)F.B_LESS, (long)f.name);
    }

    @Test
    public void funcBooleanHasCorrectPriority3() {
        // solve (F.SOLVE)
        //      > (F.B_GT)
        //          + (F.ADD)
        //              x (Polynom)
        //              1 (Polynom)
        //          5 (Polynom)
        String s = "\\solve(x+1 > 5)";
        Ring r = new Ring("R64[x]");
        F f = Parser.getF(s, r);
        assertFuncNonId(f);
        assertEquals(F.SOLVE, f.name);
        assertEquals(1, f.X.length);
        F fGt = (F) f.X[0];
        assertEquals(F.B_GT, fGt.name);
        assertEquals(2, fGt.X.length);
        assertTrue(fGt.X[0] instanceof F);
        assertTrue(fGt.X[1] instanceof Polynom);
        Polynom five = (Polynom) fGt.X[1];
        assertArrayEquals(new int[] {}, five.powers);
    }

    @Test
    public void funcPostfix() {
        // ! (F.FACTORIAL)
        //      5.00 (Polynom)
        String s = "5!";
        Ring r = new Ring("R64[x]");
        F f = Parser.getF(s, r);
        assertFuncNonId(f);
        assertEquals(F.FACTORIAL, f.name);
        Polynom p1 = (Polynom) f.X[0];
        assertEquals(0, p1.powers.length);
        assertTrue(new NumberR64("5.0").equals(p1.coeffs[0], r));
    }

    @Test
    public void funcFactorialClash() {
        // * (F.MULTIPLY)
        //      ! (F.FACTORIAL)
        //          5 (Polynom)
        //      6 (Polynom)
        String s = "5!6";
        Ring r = new Ring("Z[x]");
        F f = Parser.getF(s, r);
        assertFuncNonId(f);
        assertEquals(F.MULTIPLY, f.name);
        F fact = (F) f.X[0];
        Polynom five = (Polynom) fact.X[0];
        assertEquals(F.FACTORIAL, fact.name);
        assertTrue(new NumberZ("5").equals(five.coeffs[0], r));
        Polynom six = (Polynom) f.X[1];
        assertTrue(new NumberZ("6").equals(six.coeffs[0], r));
    }

    @Test
    public void funcSumToInfinityBecomesSeries() {
        // series (F.SERIES)
        //      + (F.ADD)
        //          (F.ID)
        //              i (Fname)
        //          + (F.ADD)
        //              cos (F.COS)
        //                  x (Polynom)
        //              (F.ID)
        //                  i (Fname)
        //      i
        //      8
        String s = "\\sum_{i = 8}^{\\infty} (\\cos(x) + i)";
        Ring r = new Ring("R64[x]");
        F f = Parser.getF(s, r);
        assertFuncNonId(f);
        assertEquals(F.SERIES, f.name);
        assertEquals(3, f.X.length);
        F f1 = (F) f.X[0];
        assertEquals(F.ADD, f1.name);
        F f11 = (F) f1.X[0];
        assertEquals(F.COS, f11.name);
        Polynom p111 = (Polynom) f11.X[0];
        assertArrayEquals(new int[] {1}, p111.powers);
        assertTrue(new NumberR64("1.0").equals(p111.coeffs[0], r));
        F f12 = (F) f1.X[1];
        assertEquals("i", ((Fname) f12.X[0]).name);
        F f2 = (F) f.X[1];
        assertEquals("i", ((Fname) f2.X[0]).name);
        Polynom p3 = (Polynom) f.X[2];
        assertEquals(0, p3.powers.length);
        assertTrue(new NumberR64("8.0").equals(p3.coeffs[0], r));
    }

    @Test
    public void funcSum() {
        // sum (F.SUM)
        //      + (F.ADD)
        //          cos(x) + ...)
        //      (F.ID)
        //          i (Fname)
        //      a
        //      b
        String s = "\\sum_{i = a}^{b} (\\cos(x) + x)";
        Ring r = new Ring("R64[x]");
        F f = Parser.getF(s, r);
        assertFuncNonId(f);
        assertEquals(F.SUM, f.name);
        assertEquals(4, f.X.length);
        F f1 = (F) f.X[0];
        assertEquals(F.ADD, f1.name);
        F f11 = (F) f1.X[0];
        assertEquals(F.COS, f11.name);
        assertTrue(f11.X[0].equals(f1.X[1], r));
        Polynom p12 = (Polynom) f1.X[1];
        assertArrayEquals(new int[] {1}, p12.powers);
        assertTrue(new NumberR64("1").equals(p12.coeffs[0], r));
        F f2 = (F) f.X[1];
        assertEquals("i", ((Fname) f2.X[0]).name);
        F f3 = (F) f.X[2];
        assertEquals("a", ((Fname) f3.X[0]).name);
        F f4 = (F) f.X[3];
        assertEquals("b", ((Fname) f4.X[0]).name);
    }

    @Test
    public void funcSequence() {
        // SEQUENCE
        //      +
        //          *
        //              ^
        //                  x
        //                  4
        //          4
        //      i
        //      0
        //      (F.ID)
        //          \\infty (Fname)
        String s = "\\sequence_{i = 0}^{\\infty} (3x^4 + 4)";
        Ring r = new Ring("R64[x]");
        F f = Parser.getF(s, r);
        assertFuncNonId(f);
        assertEquals(F.SEQUENCE, f.name);
        assertEquals(4, f.X.length);
        F f2 = (F) f.X[1];
        assertEquals("i", ((Fname) f2.X[0]).name);
        assertEquals(0, ((Polynom) f.X[2]).powers.length);
        assertEquals(0, ((Polynom) f.X[2]).coeffs.length);
        F f4 = (F) f.X[3];
        assertEquals("\\infty", ((Fname) f4.X[0]).name);
    }

    @Test
    public void funcWithoutArgs() {
        // (F.TIME)
        String s = "\\time()";
        Ring r = new Ring("R[x]");
        F f = Parser.getF(s, r);

        assertFuncNonId(f);
        assertEquals(F.TIME, f.name);
        assertEquals(0, f.X.length);
    }

    @Test
    public void funcVector() {
        // (F.VECTORS)
        //      1 (Polynom)
        //      2 (Polynom)
        //      3 (Polynom)
        String s = "[1, 2, 3]";
        Ring r = new Ring("Z[x]");
        F f = Parser.getF(s, r);
        assertEquals(F.VECTORS, f.name);
        Polynom p1 = (Polynom) f.X[0];
        assertTrue(new NumberZ("1").equals(p1.coeffs[0], r));
        Polynom p2 = (Polynom) f.X[1];
        assertTrue(new NumberZ("2").equals(p2.coeffs[0], r));
        Polynom p3 = (Polynom) f.X[2];
        assertTrue(new NumberZ("3").equals(p3.coeffs[0], r));
    }

    @Test
    public void funcVectorWithNulls() {
        // (F.VECTORS)
        //      null
        //      null
        //      3 (Polynom)
        //      null
        String s = "[, , 3, ]";
        Ring r = new Ring("Z[x]");
        F f = Parser.getF(s, r);
        assertEquals(F.VECTORS, f.name);
        assertEquals(4, f.X.length);
        assertNull(f.X[0]);
        assertNull(f.X[1]);
        assertTrue(((Polynom) f.X[2]).coeffs[0].equals(new NumberZ("3")));
        assertNull(f.X[3]);
    }

    @Test
    public void funcLogLatexNotation() {
        // log (F.LOG)
        //      + (F.ADD)
        //          x (Polynom)
        //          sin (F.SIN)
        //              x (Polynom)
        //      * (F.MULTIPLY)2x (Polynom)
        String s = "\\log_{2x} (x + \\sin(x))";
        Ring r = new Ring("Z[x]");
        F f = Parser.getF(s, r);
        assertEquals(F.LOG, f.name);
        F fPlus = (F) f.X[1];
        assertEquals(F.ADD, fPlus.name);
    }

    @Test
    public void funcLog() {
        // log (F.LOG)
        //      + (F.ADD)
        //          x (Polynom)
        //          sin (F.SIN)
        //              x (Polynom)
        //      * (F.MULTIPLY)
        //          2 (Polynom)
        //          x (Polynom)
        String s = "\\log(2x, x + \\sin(x))";
        Ring r = new Ring("Z[x]");
        F f = Parser.getF(s, r);
        assertEquals(F.LOG, f.name);
        F fPlus = (F) f.X[1];
        assertEquals(F.ADD, fPlus.name);
    }

    @Test
    public void funcSingleFunction() {
        // sin (F.SIN)
        //      x (Polynom)
        String s = "\\sin(x)";
        Ring r = new Ring("R64[x]");
        F f = Parser.getF(s, r);
        assertFuncNonId(f);
        Polynom p1 = (Polynom) f.X[0];
        assertArrayEquals(new int[] {1}, p1.powers);
        assertTrue(r.numberONE().equals(p1.coeffs[0], r));
    }

    @Test
    public void funcD() {
        // D (F.D)
        //      y (Polynom)
        //      x^2 (Polynom)
        String s = "\\D(y, x^2)";
        Ring r = new Ring("Z[x, y]");
        F f = Parser.getF(s, r);
        assertEquals(F.D, f.name);
        assertEquals(2, f.X.length);
    }

    @Test
    public void funcD01() {
        // D (F.D)
        //      f (Fname)
        //      x^2 (Polynom)
        String s = "\\D_{x^2} (f)";
        Ring r = new Ring("Z[x, y]");
        F f = Parser.getF(s, r);
        assertEquals(F.D, f.name);
        assertEquals(2, f.X.length);
        assertNotNull("Func part must be not null", f.X[0]);
        assertNotNull("Var part must be not null", f.X[1]);
    }

    @Test
    public void funcD02() {
        // D (F.D)
        //      null
        //      x^2 (Polynom)
        String s = "\\D_{x^2}";
        Ring r = new Ring("Z[x, y]");
        F f = Parser.getF(s, r);
        assertEquals(F.D, f.name);
        assertEquals(2, f.X.length);
        assertNull("Func part must be null", f.X[0]);
        assertNotNull("Var part must be not null", f.X[1]);
    }

    @Test
    public void funcD03() {
        // D (F.D)
        //      f (Fname)
        //      x (Polynom)
        String s = "\\D_x (f)";
        Ring r = new Ring("Z[x, y]");
        F f = Parser.getF(s, r);
        assertEquals(F.D, f.name);
        assertEquals(2, f.X.length);
        assertNotNull("Func part must be not null", f.X[0]);
        assertNotNull("Var part must be not null", f.X[1]);
    }

    @Test
    public void funcD04() {
        // D (F.D)
        //      null
        //      x (Polynom)
        String s = "\\D_x";
        Ring r = new Ring("Z[x, y]");
        F f = Parser.getF(s, r);
        assertEquals(F.D, f.name);
        assertEquals(2, f.X.length);
        assertNull("Func part must be null", f.X[0]);
        assertNotNull("Var part must be not null", f.X[1]);
    }

    @Test
    public void funcSystLDE() {
        // systLDE (F.SYSTLDE)
        //      + (F.ADD)
        //          d (F.d)
        //              x (Fname)
        //              t (Polynom)
        //          x (Fname)
        //      0 (Polynom)
        //      + (F.ADD)
        //          d (F.d)
        //              (F.ID)
        //                  x (Fname)
        //              t (Polynom)
        //          (F.ID)
        //              x (Fname)
        //      0 (Polynom)
        String s = "\\systLDE(\\d(x, t) + x = 0, \\d(x, t) + x = 0)";
        Ring r = new Ring("R64[t]");
        F f = Parser.getF(s, r);
        assertFuncNonId(f);
        assertEquals(F.SYSTLDE, f.name);
        assertEquals(4, f.X.length);
        F fPlus = (F) f.X[0];
        assertEquals(F.ADD, fPlus.name);
        F d = (F) fPlus.X[0];
        assertEquals(F.d, d.name);
        F x = (F) d.X[0];
        Polynom t = (Polynom) d.X[1];
        assertEquals(F.ID, x.name);
        assertTrue(r.numberONE.equals(t.coeffs[0], r));
        x = (F) fPlus.X[1];
        assertEquals(F.ID, x.name);
        Polynom zero = (Polynom) f.X[1];
        assertEquals(0, zero.powers.length);
    }

    @Test
    public void funcIntegralLimitedSingleArg() {
        // int (F.INT)
        //      sin (F.SIN)
        //          x (Polynom)
        //      x (Polynom)
        //      1 (Polynom)
        //      2 (Polynom)
        String[] inp = {
            "\\int_{1}^{2}\\sin(x)dx",
            "\\int_{1}^{2} \\sin(x)dx",
            "\\int_{1}^{2} \\sin(x) dx",
            "\\int_{1}^{2} \\sin(x) d x",
            "\\int_{1}^{2}(\\sin(x))dx",
            "\\int_{1}^{2} (\\sin(x))dx",
            "\\int_{1}^{2} (\\sin(x)) dx",
            "\\int_{1}^{2} (\\sin(x)) d x"
        };
        for (String s : inp) {
            Ring r = new Ring("R64[x]");
            F f = Parser.getF(s, r);
            assertFuncNonId(f);
            assertEquals(F.INT, f.name);
            assertEquals(4, f.X.length);
            assertFuncNonId(f.X[0]);
            F sinX = (F) f.X[0];
            assertEquals(F.SIN, sinX.name);
        }
    }

    @Test
    public void funcIntegralSimple() {
        // int (F.INT)
        //      sin (F.SIN)
        //          x (Polynom)
        //      x (Polynom)
        String[] inp = { // Эти строки дают одинаковый результат
            "\\int\\sin(x)dx",
            "\\int \\sin(x)dx",
            "\\int \\sin(x) dx",
            "\\int \\sin(x) d x",
            "\\int(\\sin(x))dx",
            "\\int (\\sin(x))dx",
            "\\int (\\sin(x)) dx",
            "\\int (\\sin(x)) d x"
        };
        Ring r = new Ring("R64[x]");
        for (String s : inp) {
            F f = Parser.getF(s, r);
            assertFuncNonId(f);
            assertEquals(s, F.INT, f.name);
            assertEquals(2, f.X.length);
            assertFuncNonId(f.X[0]);
            F sinX = (F) f.X[0];
            assertEquals(F.SIN, sinX.name);
        }
    }

    @Test
    public void funcInitCond() {
        // initCond (F.INITCOND)
        //      d (F.D)
        //          F.ID
        //              x (Fname)
        //          t (Polynom)
        //          0 (Polynom)
        //          0 (Polynom)
        //      1 (Polynom)
        //      d (F.D)
        //          F.ID
        //              y (Fname)
        //          t (Polynom)
        //          0 (Polynom)
        //          0 (Polynom)
        //      2 (Polynom)
        //      d (F.D)
        //          F.ID
        //              z (Fname)
        //          t (Polynom)
        //          0 (Polynom)
        //          0 (Polynom)
        //      3 (Polynom)
        String s = "\\initCond(\\d(x,t,0,0)=1,\\d(y,t,0,0)=2,\\d(z,t,0,0)=3)";
        Ring r = new Ring("R64[t]");
        F f = Parser.getF(s, r);

        Polynom zeroPol = new Polynom(new int[] {}, new Element[] {});
        Polynom tPol = new Polynom(new int[] {1}, new Element[] {NumberR64.ONE});

        assertFuncNonId(f);
        assertEquals(F.INITCOND, f.name);
        assertEquals(6, f.X.length);
        F fd1 = (F) f.X[0];
        assertEquals(4, fd1.X.length);
        Fname x = (Fname) ((F) fd1.X[0]).X[0];
        assertEquals("x", x.name);
        assertTrue(fd1.X[1].equals(tPol, r));
        assertTrue(fd1.X[2].equals(zeroPol, r));
        assertTrue(fd1.X[3].equals(zeroPol, r));
        Polynom fp1 = (Polynom) f.X[1];
        assertArrayEquals(new int[] {}, fp1.powers);
        assertTrue(NumberR64.ONE.equals(fp1.coeffs[0], r));
        F fd2 = (F) f.X[2];
        assertEquals(4, fd2.X.length);
        Fname y = (Fname) ((F) fd2.X[0]).X[0];
        assertEquals("y", y.name);
        assertTrue(fd2.X[1].equals(tPol, r));
        assertTrue(fd2.X[2].equals(zeroPol, r));
        assertTrue(fd2.X[3].equals(zeroPol, r));
        Polynom fp2 = (Polynom) f.X[3];
        assertArrayEquals(new int[] {}, fp2.powers);
        assertTrue(new NumberR64("2").equals(fp2.coeffs[0], r));
        F fd3 = (F) f.X[4];
        assertEquals(4, fd3.X.length);
        Fname z = (Fname) ((F) fd3.X[0]).X[0];
        assertEquals("z", z.name);
        assertTrue(fd3.X[1].equals(tPol, r));
        assertTrue(fd3.X[2].equals(zeroPol, r));
        assertTrue(fd3.X[3].equals(zeroPol, r));
        Polynom fp3 = (Polynom) f.X[5];
        assertArrayEquals(new int[] {}, fp3.powers);
        assertTrue(new NumberR64("3").equals(fp3.coeffs[0], r));
    }

    @Test
    public void funcInitCondWithIcWithNull() {
        String s = "\\initCond(\\ic(f,x,0,0,y,,0)=1,\\ic(f,x,0,1,y,,0)=-y^2/2)";
        Ring r = new Ring("R64[x,y]");
        F f = Parser.getF(s, r);
        assertFuncNonId(f);
        assertEquals(F.INITCOND, f.name);
        assertEquals(4, f.X.length);
        F fIc = (F) f.X[0];
        assertEquals(7, fIc.X.length);
        assertEquals(F.IC, fIc.name);
    }

    @Test
    public void funcSolveWithEqualsShouldBecomeToSolveLEQ() {
        String s = "\\solve(x + 1 = 0)";
        Ring r = new Ring("R64[x]");
        F f = Parser.getF(s, r);
        assertFuncNonId(f);
        assertEquals(F.SOLVEEQ, f.name);
        assertEquals(2, f.X.length);
    }

    @Test
    public void funcSolveWithEqualsAndVectorsShouldBecomeToSolveLEQ() {
        String s = "\\solve([x+3y-1 = 0, y+5x+7 = 0], [x, y])";
        Ring r = new Ring("R64[x]");
        F f = Parser.getF(s, r);
        assertFuncNonId(f);
        assertEquals(F.SOLVEEQ, f.name);
        assertEquals(2, f.X.length);
        F vector = (F) f.X[0];
        assertEquals(F.VECTORS, vector.name);
        assertEquals(4, vector.X.length);
        F vector1 = (F) f.X[1];
        assertEquals(F.VECTORS, vector1.name);
        assertEquals(2, vector1.X.length);
    }

    @Test
    public void funcMod() {
        String s = "\\mod(a, b)";
        Ring r = new Ring("Z[x]");
        F f = Parser.getF(s, r);
        assertFuncNonId(f);
        assertEquals(F.MOD, f.name);
        assertEquals(2, f.X.length);
    }

    @Test
    public void funcFloor() {
        String s = "\\floor(a)";
        Ring r = new Ring("Z[x]");
        F f = Parser.getF(s, r);
        assertFuncNonId(f);
        assertEquals(F.FLOOR, f.name);
        assertEquals(1, f.X.length);
    }

    @Test
    public void funcCeil() {
        String s = "\\ceil(a)";
        Ring r = new Ring("Z[x]");
        F f = Parser.getF(s, r);
        assertFuncNonId(f);
        assertEquals(F.CEIL, f.name);
        assertEquals(1, f.X.length);
    }

    @Test
    public void funcRound() {
        String s = "\\round(a)";
        Ring r = new Ring("Z[x]");
        F f = Parser.getF(s, r);
        assertFuncNonId(f);
        assertEquals(F.ROUND, f.name);
        assertEquals(1, f.X.length);
    }

    @Test
    public void funcSolveTrig() {
        String s = "\\solveTrig(\\sin(x))";
        Ring r = new Ring("R64[x]");
        F f = Parser.getF(s, r);
        assertFuncNonId(f);
        assertEquals(F.SOLVETRIG, f.name);
        assertEquals(1, f.X.length);
    }

    @Test
    public void funcSolveTrigWithEqual() {
        String s = "\\solveTrig(\\sin(x) = 0)";
        Ring r = new Ring("R64[x]");
        F f = Parser.getF(s, r);
        assertFuncNonId(f);
        assertEquals(F.SOLVETRIG, f.name);
        assertEquals(1, f.X.length);
    }

    @Test
    public void funcSolveTrigWithVar() {
        String s = "\\solveTrig(\\sin(x), x)";
        Ring r = new Ring("R64[x]");
        F f = Parser.getF(s, r);
        assertFuncNonId(f);
        assertEquals(F.SOLVETRIG, f.name);
        assertEquals(2, f.X.length);
    }

    @Test
    public void funcSolveTrigWithEqualWithVar() {
        String s = "\\solveTrig(\\sin(x) = 0, x)";
        Ring r = new Ring("R64[x]");
        F f = Parser.getF(s, r);
        assertFuncNonId(f);
        assertEquals(F.SOLVETRIG, f.name);
        assertEquals(2, f.X.length);
        assertEquals(F.SUBTRACT, ((F) f.X[0]).name);
    }

    @Test
    public void funcSolveTrigWithEqualWithVarWithIndex() {
        String s = "\\solveTrig(\\sin(x) = 0, x, k)";
        Ring r = new Ring("R64[x]");
        F f = Parser.getF(s, r);
        assertFuncNonId(f);
        assertEquals(F.SOLVETRIG, f.name);
        assertEquals(3, f.X.length);
        assertEquals(F.SUBTRACT, ((F) f.X[0]).name);
    }

    @Test
    public void funcSolveTrigWithVarWithIndex() {
        String s = "\\solveTrig(\\sin(x), x, k)";
        Ring r = new Ring("R64[x]");
        F f = Parser.getF(s, r);
        assertFuncNonId(f);
        assertEquals(F.SOLVETRIG, f.name);
        assertEquals(3, f.X.length);
    }

    @Test
    public void funcGetStatusDoesntConflictWithGeBooleanOperation() {
        String s = "\\getStatus(0)";
        Ring r = new Ring("Z[x]");
        F f = Parser.getF(s, r);
        assertFuncNonId(f);
        assertEquals(F.GETSTATUS, f.name);
        assertEquals(1, f.X.length);
    }

    @Test
    public void funcCupInfix() {
        String s = "a \\cup b";
        Ring r = new Ring("Z[x]");
        F f = Parser.getF(s, r);
        assertFuncNonId(f);
        assertEquals(F.CUP, f.name);
        assertEquals(2, f.X.length);
        assertTrue(f.X[0].equals(new F(new Fname("a"))));
        assertTrue(f.X[1].equals(new F(new Fname("b"))));
    }

    @Test
    public void funcCapInfix() {
        String s = "a \\cap b";
        Ring r = new Ring("Z[x]");
        F f = Parser.getF(s, r);
        assertFuncNonId(f);
        assertEquals(F.CAP, f.name);
        assertEquals(2, f.X.length);
        assertTrue(f.X[0].equals(new F(new Fname("a"))));
        assertTrue(f.X[1].equals(new F(new Fname("b"))));
    }

    @Test
    public void funcMultManyArgs() {
        String s = "a b c";
        Ring r = new Ring("Z[x]");
        F f = Parser.getF(s, r);
        assertFuncNonId(f);
        assertEquals(F.MULTIPLY, f.name);
        assertEquals(3, f.X.length);
    }

    @Test
    public void funcMultManyArgsWithParens() {
        String s = "a (b c d) e";
        Ring r = new Ring("Z[x]");
        F f = Parser.getF(s, r);
        assertFuncNonId(f);
        assertEquals(F.MULTIPLY, f.name);
        assertEquals(3, f.X.length);
        assertFuncNonId(f.X[1]);
        F insideParens = (F) f.X[1];
        assertEquals(F.MULTIPLY, insideParens.name);
        assertEquals(3, insideParens.X.length);
    }

    @Test
    public void funcAddManyArgs() {
        String s = "a + b + c";
        Ring r = new Ring("Z[x]");
        F f = Parser.getF(s, r);
        assertFuncNonId(f);
        assertEquals(F.ADD, f.name);
        assertEquals(3, f.X.length);
    }

    @Test
    public void funcAddManyArgsWithParens() {
        String s = "a + (b + c + d) + e";
        Ring r = new Ring("Z[x]");
        F f = Parser.getF(s, r);
        assertFuncNonId(f);
        assertEquals(F.ADD, f.name);
        assertEquals(3, f.X.length);
        assertFuncNonId(f.X[1]);
        F insideParens = (F) f.X[1];
        assertEquals(F.ADD, insideParens.name);
        assertEquals(3, insideParens.X.length);
    }
}
