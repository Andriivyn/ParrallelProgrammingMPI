package com.mathpar.func.parser;

import com.mathpar.func.F;
import com.mathpar.func.Fname;
import com.mathpar.number.NumberZ;
import com.mathpar.number.Ring;
import com.mathpar.number.SubsetR;
import org.junit.BeforeClass;
import org.junit.Test;
import com.mathpar.polynom.Polynom;

import static com.mathpar.func.parser.FunctionTest.assertFuncNonId;
import static org.junit.Assert.*;

/**
 * "Special" functions tests.
 *
 * @author ivan
 */
public class FunctionSpecTest {
    private static Ring rZx;

    @BeforeClass
    public static void setUpClass() throws Exception {
        rZx = new Ring("Z[x]");
    }

    private static void testMatrix1234(F f) {
        F f1 = (F) f.X[0];
        assertEquals(F.VECTORS, f1.name);
        assertEquals(2, f1.X.length);
        F f11 = (F) f1.X[0];
        assertEquals(F.VECTORS, f11.name);
        assertEquals(2, f11.X.length);
        Polynom p111 = (Polynom) f11.X[0];
        assertEquals(0, p111.powers.length);
        assertTrue(new NumberZ("1").equals(p111.coeffs[0], rZx));
        Polynom p112 = (Polynom) f11.X[1];
        assertEquals(0, p112.powers.length);
        assertTrue(new NumberZ("2").equals(p112.coeffs[0], rZx));
        F f12 = (F) f1.X[1];
        assertEquals(F.VECTORS, f12.name);
        assertEquals(2, f12.X.length);
        Polynom p121 = (Polynom) f12.X[0];
        assertEquals(0, p121.powers.length);
        assertTrue(new NumberZ("3").equals(p121.coeffs[0], rZx));
        Polynom p122 = (Polynom) f12.X[1];
        assertEquals(0, p122.powers.length);
        assertTrue(new NumberZ("4").equals(p122.coeffs[0], rZx));
    }

    private static void testFnameA(F f) {
        F f1 = (F) f.X[0];
        assertEquals(F.ID, f1.name);
        assertEquals("A", ((Fname) f1.X[0]).name);
    }

    @Test
    public void funcSpecTransposeName() {
        // transpose (F.TRANSPOSE)
        //      (F.ID)
        //          A (Fname)
        F f = Parser.getF("A^T", rZx);

        assertFuncNonId(f);
        assertEquals(F.TRANSPOSE, f.name);
        testFnameA(f);
    }

    @Test
    public void funcSpecTransposeMatrixLiteral() {
        // transpose (F.TRANSPOSE)
        //      (F.VECTORS)
        //          (F.VECTORS)
        //              1 (Polynom)
        //              2 (Polynom)
        //          (F.VECTORS)
        //              3 (Polynom)
        //              4 (Polynom)
        F f = Parser.getF("[[1, 2], [3, 4]]^{T}", rZx);

        assertFuncNonId(f);
        assertEquals(F.TRANSPOSE, f.name);
        testMatrix1234(f);
    }

    @Test
    public void funcSpecConjugateMatrixName() {
        // conjugate (F.CONJUGATE)
        //      (F.ID)
        //          A (Fname)
        F f = Parser.getF("A^{\\ast}", rZx);

        assertFuncNonId(f);
        assertEquals(F.CONJUGATE, f.name);
        testFnameA(f);
    }

    @Test
    public void funcSpecConjugateMatrixLiteral() {
        // conjugate (F.CONJUGATE)
        //      (F.ID)
        //          A (Fname)
        F f = Parser.getF("[[1, 2], [3, 4]]^{\\ast}", rZx);

        assertFuncNonId(f);
        assertEquals(F.CONJUGATE, f.name);
        testMatrix1234(f);
    }

    @Test
    public void funcSpecConjugatePolynomial() {
        // conjugate (F.CONJUGATE)
        //      x (Polynom)
        F f = Parser.getF("x^{\\ast}", rZx);

        assertFuncNonId(f);
        assertEquals(F.CONJUGATE, f.name);
        Polynom p = (Polynom) f.X[0];
        assertArrayEquals(new int[] {1}, p.powers);
        assertTrue(rZx.numberONE().equals(p.coeffs[0], rZx));
    }

    @Test
    public void funcSpecAdjointMatrixName() {
        // adjoint (F.ADJOINT)
        //      (F.ID)
        //          A (Fname)
        F f = Parser.getF("A^\\star", rZx);

        assertFuncNonId(f);
        assertEquals(F.ADJOINT, f.name);
        testFnameA(f);
    }

    @Test
    public void funcSpecAdjointMatrixLiteral() {
        // adjoint (F.ADJOINT)
        //      (F.VECTORS)
        //          (F.VECTORS)
        //              1 (Polynom)
        //              2 (Polynom)
        //          (F.VECTORS)
        //              3 (Polynom)
        //              4 (Polynom)
        F f = Parser.getF("[[1, 2], [3, 4]]^{\\star}", rZx);

        assertFuncNonId(f);
        assertEquals(F.ADJOINT, f.name);
        testMatrix1234(f);
    }

    @Test
    public void funcSpecGenInverseMatrixName() {
        // genInverse (F.GENINVERSE)
        //      (F.ID)
        //          A (Fname)
        F f = Parser.getF("A^{+}", rZx);

        assertFuncNonId(f);
        assertEquals(F.GENINVERSE, f.name);
        testFnameA(f);
    }

    @Test
    public void funcSpecGenInversMatrixLiteral() {
        // genInverse (F.GENINVERSE)
        //      (F.VECTORS)
        //          (F.VECTORS)
        //              1 (Polynom)
        //              2 (Polynom)
        //          (F.VECTORS)
        //              3 (Polynom)
        //              4 (Polynom)
        F f = Parser.getF("[[1, 2], [3, 4]]^{+}", rZx);

        assertFuncNonId(f);
        assertEquals(F.GENINVERSE, f.name);
        testMatrix1234(f);
    }

    @Test
    public void funcSpecInverseMatrixName() {
        // inverse (F.INVERSE)
        //      (F.ID)
        //          A (Fname)
        F f = Parser.getF("A^{-1}", rZx);

        assertFuncNonId(f);
        assertEquals(F.POW, f.name);
        testFnameA(f);
    }

    @Test
    public void funcSpecInverseMatrixLiteral() {
        // inverse (F.INVERSE)
        //      (F.VECTORS)
        //          1 (Polynom)
        //          2 (Polynom)
        //      (F.VECTORS)
        //          3 (Polynom)
        //          4 (Polynom)
        F f = Parser.getF("[[1, 2], [3, 4]]^{-1}", rZx);

        assertFuncNonId(f);
        assertEquals(F.POW, f.name);
        testMatrix1234(f);
    }

    @Test
    public void funcSpecClosureName() {
        // closure (F.CLOSURE)
        //      (F.ID)
        //          A (Fname)
        F f = Parser.getF("A^\\times", rZx);

        assertFuncNonId(f);
        assertEquals(F.CLOSURE, f.name);
        testFnameA(f);
    }

    @Test
    public void funcSpecClosureMatrixLiteral() {
        // closure (F.CLOSURE)
        //      (F.VECTORS)
        //          (F.VECTORS)
        //              1 (Polynom)
        //              2 (Polynom)
        //          (F.VECTORS)
        //              3 (Polynom)
        //              4 (Polynom)
        F f = Parser.getF("[[1, 2], [3, 4]]^\\times", rZx);

        assertFuncNonId(f);
        assertEquals(F.CLOSURE, f.name);
        testMatrix1234(f);
    }

    @Test(timeout = 500)
    public void funcSpecLim() {
        // lim (F.LIM)
        //      / (F.DIVIDE)
        //          1 (Polynom)
        //          x (Polynom)
        //      x (Polynom)
        //      (F.ID)
        //          \infty (Fname)
        String s = "\\lim_{x \\to \\infty} (1 / x)";
        Ring r = new Ring("R[x, y]");
        F f = Parser.getF(s, r);
        assertFuncNonId(f);

        assertEquals(F.LIM, f.name);
        assertEquals(3, f.X.length);
        F f1 = (F) f.X[0];
        assertEquals(F.DIVIDE, f1.name);
        Polynom p11 = (Polynom) f1.X[0];
        assertEquals(0, p11.powers.length);
        assertTrue(r.numberONE().equals(p11.coeffs[0], r));
        Polynom p12 = (Polynom) f1.X[1];
        assertArrayEquals(new int[] {1}, p12.powers);
        assertTrue(r.numberONE().equals(p12.coeffs[0], r));
        Polynom p2 = (Polynom) f.X[1];
        assertArrayEquals(new int[] {1}, p2.powers);
        assertTrue(r.numberONE().equals(p2.coeffs[0], r));
        F f3 = (F) f.X[2];
        assertEquals("\\infty", ((Fname) f3.X[0]).name);

        // lim (F.LIM)
        //      x (Polynom)
        //      x (Polynom)
        //      (F.ID)
        //          \infty (Fname)
        s = "\\lim_{x \\to \\infty} x";
        f = Parser.getF(s, r);
        assertFuncNonId(f);

        assertEquals(F.LIM, f.name);
        assertEquals(3, f.X.length);
        Polynom p1 = (Polynom) f.X[0];
        assertArrayEquals(new int[] {1}, p1.powers);
        assertTrue(r.numberONE().equals(p1.coeffs[0], r));
        p2 = (Polynom) f.X[1];
        assertArrayEquals(new int[] {1}, p2.powers);
        assertTrue(r.numberONE().equals(p2.coeffs[0], r));
        f3 = (F) f.X[2];
        assertEquals("\\infty", ((Fname) f3.X[0]).name);

        // lim (F.LIM)
        //      sin (F.SIN)
        //          x (Polynom)
        //      x (Polynom)
        //      (F.ID)
        //          \infty (Fname)
        s = "\\lim_{x \\to \\infty} \\sin(x)";
        f = Parser.getF(s, r);
        assertFuncNonId(f);

        assertEquals(F.LIM, f.name);
        assertEquals(3, f.X.length);
        f1 = (F) f.X[0];
        assertEquals(F.SIN, f1.name);
        p11 = (Polynom) f1.X[0];
        assertArrayEquals(new int[] {1}, p11.powers);
        assertTrue(r.numberONE().equals(p11.coeffs[0], r));
        p2 = (Polynom) f.X[1];
        assertArrayEquals(new int[] {1}, p2.powers);
        assertTrue(r.numberONE().equals(p2.coeffs[0], r));
        f3 = (F) f.X[2];
        assertEquals("\\infty", ((Fname) f3.X[0]).name);
    }

    @Test
    public void funcSpecSetBasic() {
        String s = "\\set((1, 5), (6, 7), {7})";
        Ring r = new Ring("Z[x]");

        F f = Parser.getF(s, r);

        assertEquals(F.ID, f.name);
        assertEquals(1, f.X.length);
        assertTrue(f.X[0] instanceof SubsetR);
    }

    @Test
    public void funcSpecSetBasicWithSymbols() {
        String s = "\\set((1, 2) \\cup [3, 4] \\cup {5} \\cup [6, 7) \\cup (7, 8])";
        Ring r = new Ring("Z[x]");

        F f = Parser.getF(s, r);

        assertEquals(F.ID, f.name);
        assertEquals(1, f.X.length);
        assertTrue(f.X[0] instanceof SubsetR);
    }

    @Test
    public void funcSpecSetSingle() {
        String s = "\\(1, 2\\]";
        Ring r = new Ring("Z[x]");

        F f = Parser.getF(s, r);

        assertEquals(F.ID, f.name);
        assertEquals(1, f.X.length);
        assertTrue(f.X[0] instanceof SubsetR);
    }

    @Test
    public void funcSpecVectorSet() {
        // F.VECTOR_SET
        // F.VECTORS
        //      1
        //      Element.LESS
        //      a
        //      Element.LESS
        //      2
        // F.VECTORS
        //      0
        //      Element.LESS_OR_EQUAL
        //      b
        //      Element.LESS
        //      c
        String s = "\\{1 < a < 2, 0 \\le b < c\\}";
        Ring r = new Ring("Z[x]");

        try {
            F f = Parser.getF(s, r);
            assertEquals(F.VECTOR_SET, f.name);
            assertEquals(2, f.X.length); // Должно быть два вектора
            // Первый вектор
            assertTrue(f.X[0] instanceof F);
            F vec0 = (F) f.X[0];
            assertEquals(F.VECTORS, vec0.name);
            assertEquals(5, vec0.X.length);
            // Второй вектор
            assertTrue(f.X[1] instanceof F);
            F vec1 = (F) f.X[1];
            assertEquals(F.VECTORS, vec0.name);
            assertEquals(5, vec1.X.length);
            assertEquals("c", vec1.X[4].toString(r)); // Проверка, что "c" не слиплось со слешем.
        } catch (ParserException e) {
            e.getCause().printStackTrace();
        }
    }
}
