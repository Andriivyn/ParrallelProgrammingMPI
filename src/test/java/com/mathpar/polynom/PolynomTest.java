package com.mathpar.polynom;

import com.mathpar.func.F;
import com.mathpar.number.Element;
import com.mathpar.number.NumberZ;
import com.mathpar.number.Ring;
import org.junit.Test;

import static org.hamcrest.MatcherAssert.assertThat;
import static org.junit.Assert.*;

public class PolynomTest {
    
        @Test(timeout = 5000)
    public void factotInR64() {
        Ring r = new Ring("R64[x]");
        Polynom p = new Polynom("(1/9)x^4+ (1/9)x^2", r);
        Element fp = p.factor(r);
        System.out.println(fp);
        Polynom p1 = new Polynom("((1/3)x)", r);
        Polynom p2 = new Polynom("x^2+1", r);
        int[] pow = new int[]{2, 1};
        Polynom[] mm = new Polynom[]{p1, p2};
        FactorPol expected = new FactorPol(pow, mm);
        Element exp = expected.toPolynomOrFraction(r);
        Element res = ((FactorPol)fp).toPolynomOrFraction(r);
        res = res.subtract(exp, r);
        assertTrue(res.isZero(r));
    }
    
    
    @Test(timeout = 5000)
    public void factorInCFirstVar() {
        Ring r = new Ring("Z[x, y, z]");
        Polynom p = new Polynom("x^2 + 1", r);

        FactorPol res = p.factorOfPol_inC(r);

        assertEquals(2, res.multin.length);
        assertArrayEquals(new int[]{1, 1}, res.powers);
        // x - i
        assertArrayEquals(new int[]{1, 0}, res.multin[0].powers);
        // x + i
        assertArrayEquals(new int[]{1, 0}, res.multin[1].powers);
    }

    @Test(timeout = 5000)
    public void factorInCNotFirstVar() {
        Ring r = new Ring("Z[x, y, z]");
        Polynom p = new Polynom("y^2 + 1", r);

        FactorPol res = p.factorOfPol_inC(r);

        assertEquals(2, res.multin.length);
        assertArrayEquals(new int[]{1, 1}, res.powers);
        // y - i
        assertArrayEquals(new int[]{0, 1, 0, 0}, res.multin[0].powers);
        // y + i
        assertArrayEquals(new int[]{0, 1, 0, 0}, res.multin[1].powers);

        res = new Polynom("z^2 + 1", r).factorOfPol_inC(r);
        assertEquals(2, res.multin.length);
        assertArrayEquals(new int[]{1, 1}, res.powers);
        // z - i
        assertArrayEquals(new int[]{0, 0, 1, 0, 0, 0}, res.multin[0].powers);
        // z + i
        assertArrayEquals(new int[]{0, 0, 1, 0, 0, 0}, res.multin[1].powers);
    }

    @Test(timeout = 5000)
    public void factorInCNonfactorizablePolynomialStaysTheSame() {
        Ring r = new Ring("Z[x, y, z]");
        Polynom p = new Polynom("x + y", r);
        FactorPol res = p.factorOfPol_inC(r);
        assertEquals(1, res.multin.length);
        assertArrayEquals(new int[]{1}, res.powers);
        // p lives in FactorPol without changes, but its coefficients go from Z to C
        assertTrue(p.subtract(res.multin[0], r).isZero(r));
    }

    @Test(timeout = 5000)
    public void testFactorInCExample002() {
        Ring r = new Ring("Z[x, y, z]");
        Polynom p = new Polynom("-y^4 - 2y^2", r);

        FactorPol res = p.factorOfPol_inC(r);

        assertEquals(4, res.multin.length);
        assertArrayEquals(new int[]{1, 2, 1, 1}, res.powers);
    }

    @Test
    public void homogenize() {
        Ring r = new Ring("Z[x, y, z]");

        Polynom p = new Polynom("y^3x + yx^6 + yx + x^20 + 1", r);
        Polynom expected = new Polynom("z^3yx^16 + zy^6x^13 + zyx^18 + y^20 + x^20", r);
        assertTrue(expected.equals(p.homogenize(2), r));

        Polynom p2 = new Polynom("y^2 + y x^2", r);
        Polynom expected2 = new Polynom("z^2x + zy^2", r);
        assertTrue(expected2.equals(p2.homogenize(2), r));

        Polynom p3 = new Polynom("x^2 + 1", r);
        assertArrayEquals(new int[]{0, 2, 0, 2, 0, 0}, p3.homogenize(2).powers);
    }

    @Test
    public void dehomogenize() {
        Ring r = new Ring("Z[x, y, z]");

        Polynom p = new Polynom("z^3yx^16 + zy^6x^13 + zyx^18 + y^20 + x^20", r);
        Polynom expected = new Polynom("y^3x + yx^6 + yx + x^20 + 1", r);
        assertTrue(expected.equals(p.dehomogenize(2), r));

        Polynom p2 = new Polynom("x", r);
        assertTrue(new Polynom("1", r).equals(p2.dehomogenize(2)));
    }

    @Test
    public void rootOf() {
        Ring r = new Ring("Z[x, y]");

        Polynom x = new Polynom("x", r);
        Polynom p = new Polynom("x^2", r);
        assertThat("sqrt[2](x^2) = x", p.rootOf(2, r).equals(new FactorPol(x), r));
        assertThat("sqrt[3](x^2) = sqrt[3](x^2)", p.rootOf(3, r)
                .equals(new F(F.ROOTOF, new FactorPol(new int[]{2}, new Polynom[]{x}),
                        new NumberZ("3")), r));

        p = new Polynom("(x+1)^2(x+5)^4", r);
        assertThat("sqrt[2]((x+1)^2(x+5)^4) = (x+1)(x+5)^2", p.rootOf(2, r)
                .equals(new FactorPol(new int[]{1, 2}, new Polynom[]{
                        new Polynom("x+1", r), new Polynom("x+5", r)
                }), r));
        assertThat("sqrt[3]((x+1)^2(x+5)^4) = (x+5)*sqrt[3]((x+1)^2(x+5))",
                p.rootOf(3, r).equals(
                        new F(
                                F.MULTIPLY, new FactorPol(new Polynom("x+5", r)),
                                new F(
                                        F.ROOTOF,
                                        new FactorPol(new int[]{2, 1}, new Polynom[]{
                                                new Polynom("x+1", r),
                                                new Polynom("x+5", r)
                                        }), new NumberZ("3"))), r));
    }

    @Test
    public void factotInZ() {
        Ring r = new Ring("Z[x]");
        Polynom p = new Polynom("1-x-285*x^2-411*x^3+18027*x^4+20689*x^5-472275*x^6-271027*x^7+6149853*x^8+"
                + "471319*x^9-42303393*x^10+10402780*x^11+157353820*x^12-58545372*x^13-335484428*x^14+"
                + "123321948*x^15+429447820*x^16-123321948*x^17-335484428*x^18+58545372*x^19+"
                + "157353820*x^20-10402780*x^21-42303393*x^22-471319*x^23+"
                + "6149853*x^24+271027*x^25-472275*x^26-20689*x^27+18027*x^28+411*x^29-285*x^30+x^31+x^32", r);
        FactorPol expected = new FactorPol(p, r);
        FactorPol res = p.factorOfPol_inQ(true, r);
        assertTrue(expected.equals(p.FactorPol_SquareFreeOneVar(r), r));
    }

    @Test
    public void factotInC() {
        Ring r = new Ring("C64[x]");
        Polynom p = new Polynom("x^4+ x^2", r);
        FactorPol fp = p.factorOfPol_inC(r);
        System.out.println(fp);
        Polynom p1 = new Polynom("(x)", r);
        Polynom p2 = new Polynom("x-\\i", r);
        Polynom p3 = new Polynom("(x+\\i", r);
        int[] pow = new int[]{2, 1, 1};
        Polynom[] mm = new Polynom[]{p1, p2, p3};
        FactorPol expected = new FactorPol(pow, mm);
        Element exp = expected.toPolynomOrFraction(r);
        Element res = fp.toPolynomOrFraction(r);
        res = res.subtract(exp, r);
        assertTrue(res.isZero(r));
    }
    
  
    
    
    
//        @Test(timeout = 5000)
//    public void factorInCFirstVar3() {
//        Ring r = new Ring("Z[x, y, z]");
//        Polynom p = new Polynom("x^3 + 1", r);
//
//        FactorPol res = p.factorOfPol_inC(r);
//
//        assertEquals(2, res.multin.length);
//        assertArrayEquals(new int[]{1, 1,1}, res.powers);
//        // x - i
//        assertArrayEquals(new int[]{1, 0}, res.multin[0].powers);
//        // x + i
//        assertArrayEquals(new int[]{ 1, 0}, res.multin[1].powers);
//        assertArrayEquals(new int[]{ 1, 0}, res.multin[2].powers);
//    }
    
}
