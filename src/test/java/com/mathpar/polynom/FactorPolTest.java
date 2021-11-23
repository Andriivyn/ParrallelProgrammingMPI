package com.mathpar.polynom;

import com.mathpar.func.F;
import com.mathpar.number.Element;
import com.mathpar.number.NumberZ;
import com.mathpar.number.Ring;
import static org.junit.Assert.*;
import org.junit.Test;

public class FactorPolTest {
    private static String msg(Element fp, Element res, Element expected, Ring r) {
        return String.format("D(%s) = %s%nExpected: %s", fp.toString(r), res.toString(r), expected.toString(r));
    }

    @Test
    public void d() {
        Ring r = new Ring("Z[x,y]");
        Element fp, expected, res;
//-----%<------
        fp = new FactorPol(new int[] {1}, new Polynom[] {
            new Polynom(new NumberZ("42"))
        });
        expected = NumberZ.ZERO;
        res = fp.D(0, r);
        assertTrue(msg(fp, res, expected, r), res.equals(expected, r));
//-----%<------
        fp = new FactorPol(new int[] {1, 1}, new Polynom[] {
            new Polynom("x", r),
            new Polynom("y", r)
        });
        expected = new FactorPol(new int[] {1, 1}, new Polynom[] {
            new Polynom("y", r)
        });
        res = fp.D(0, r);
        assertTrue(msg(fp, res, expected, r), res.equals(expected, r));
//-----%<------
        expected = new FactorPol(new int[] {1, 1}, new Polynom[] {
            new Polynom("x", r)
        });
        res = fp.D(1, r);
        assertTrue(msg(fp, res, expected, r), res.equals(expected, r));
//-----%<------
        fp = new FactorPol(new int[] {1, 2, 3}, new Polynom[] {
            new Polynom("42", r),
            new Polynom("5x^3 + 2y^4 + 5y", r),
            new Polynom("x + 1", r)
        });
        expected = new F(F.ADD, new Element[] {
            new FactorPol(new int[] {1, 3, 1}, new Polynom[] {
                new Polynom("84", r),
                new Polynom("x+1", r),
                new Polynom("15x^2", r)
            }),
            new FactorPol(new int[] {1, 2}, new Polynom[] {
                new Polynom("126", r),
                new Polynom("2y^4+5y+5x^3", r)
            })
        });
        res = fp.D(0, r);
        assertTrue(msg(fp, res, expected, r), res.equals(expected, r));
//-----%<------
        fp = new FactorPol(new int[] {2, 3}, new Polynom[] {
            new Polynom("x+1", r),
            new Polynom("x+2", r)
        });
        expected = new F(F.ADD, new Element[] {
            new FactorPol(new int[] {1, 1, 3}, new Polynom[] {
                new Polynom("2", r),
                new Polynom("x+1", r),
                new Polynom("x+2", r)
            }),
            new FactorPol(new int[] {1, 2, 2}, new Polynom[] {
                new Polynom("3", r),
                new Polynom("x+1", r),
                new Polynom("x+2", r)
            })
        });
        res = fp.D(0, r);
        assertTrue(msg(fp, res, expected, r), res.equals(expected, r));
    }


}
