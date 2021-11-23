package com.mathpar.polynom.gbasis.util;

import com.mathpar.matrix.MatrixS;
import com.mathpar.number.Element;
import com.mathpar.number.Fraction;
import com.mathpar.number.Ring;
import com.mathpar.number.VectorS;
import org.junit.Test;
import com.mathpar.web.exceptions.MathparException;

import static org.junit.Assert.*;
import static com.mathpar.polynom.gbasis.util.Modular.*;

public class ModularTest {
    @Test
    public void rationalReconstructionSimple() {
        assertEquals(new Fraction(1, 2), rationalReconstruction(7, 13));
        assertEquals(new Fraction(2, 3), rationalReconstruction(34, 100));
        assertEquals(new Fraction(-1, 2), rationalReconstruction(8, 17));
        assertEquals(new Fraction(7, 4), rationalReconstruction(-215, 289));
        assertEquals(new Fraction(45, 97), rationalReconstruction(204977, 292393));
        assertEquals(new Fraction(1, 4), rationalReconstruction(-72, 289));
    }

    @Test(expected = MathparException.class)
    public void rationalReconstruction400mod1000DoesntExist() {
        rationalReconstruction(400, 1000);
    }

    @Test
    public void solveLaePadic1() {
        Ring ringZ = new Ring("Z[]");
        MatrixS a = new MatrixS(new int[][] {
            {3, 1},
            {2, 2}}, ringZ);
        VectorS b = new VectorS(new long[] {-1, -3}, ringZ);

        VectorS x = solveLaePadic(a, b, 17);
        VectorS expected = new VectorS(new Element[] {new Fraction(1, 4), new Fraction(-7, 4)});
        assertTrue(expected.equals(x, ringZ));
    }

    /**
     * Example from
     * http://zealint.ru/education/how-to-solve-integer-linear-system.html .
     */
    @Test
    public void solveLaePadic2() {
        Ring ringZ = new Ring("Z[]");
        MatrixS a = new MatrixS(new int[][] {
            {-10, 12},
            {16, -27}}, ringZ);
        VectorS b = new VectorS(new long[] {-4, 8}, ringZ);

        VectorS x = solveLaePadic(a, b, 5);
        VectorS expected = new VectorS(new Element[] {new Fraction(2, 13), new Fraction(-8, 39)});
        assertTrue(expected.equals(x, ringZ));
    }

    @Test
    public void solveLaePadic3() {
        Ring r = Ring.ringZxyz;
        MatrixS a = new MatrixS(new int[][] {
            {2, 1, 3},
            {2, 6, 8},
            {6, 8, 18}}, r);
        VectorS b = new VectorS(new long[] {1, 3, 5}, r);

        VectorS x = solveLaePadic(a, b, 7);
        VectorS expected = new VectorS(new Element[] {
            new Fraction(3, 10),
            new Fraction(2, 5),
            new Fraction(0)});
        assertTrue(expected.equals(x, r));
    }

    @Test(expected = MathparException.class)
    public void solveLaePadic3Modulo5IsSingular() {
        Ring r = Ring.ringZxyz;
        MatrixS a = new MatrixS(new int[][] {
            {2, 1, 3},
            {2, 6, 8},
            {6, 8, 18}}, r);
        VectorS b = new VectorS(new long[] {1, 3, 5}, r);

        VectorS x = solveLaePadic(a, b, 5);
    }
}
