package com.mathpar.polynom.gbasis.algorithms;

import com.mathpar.polynom.gbasis.util.Utils;
import com.mathpar.number.Element;
import com.mathpar.number.Ring;
import org.junit.BeforeClass;
import org.junit.Test;

import static org.junit.Assert.assertTrue;

public class FaugereNewTest {
    @BeforeClass
    public static void beforeClass() {
    }

    private void doTest(Ring r, String[] ideal, String[] exp) {
        assertTrue(Element.equals(r,
                Utils.polList(r, exp),
                new FaugereNew(r, ideal).gbasis()));
    }

    @Test
    public void katsura3() {
        doTest(new Ring("Z[z, y, x]"), new String[] {
            "x + 2*y + 2*z - 1",
            "x^2 - x + 2*y^2 + 2*z^2",
            "2*x*y + 2*y*z - y"
        }, new String[] {
            "7x-420z^3+158z^2+8z-7",
            "7y+210z^3-79z^2+3z",
            "84z^4-40z^3+z^2+z"
        });
    }

    @Test
    public void katsura4() {
        /*
         * sage (Ring is Q):
         *
         * a - 53230079232/1971025*d^7 + 10415423232/1971025*d^6
         * + 9146536848/1971025*d^5 - 2158574456/1971025*d^4
         * - 838935856/5913075*d^3 + 275119624/5913075*d^2 + 4884038/5913075*d - 1,
         *
         * b - 97197721632/1971025*d^7 + 73975630752/1971025*d^6 -
         * 12121915032/1971025*d^5 - 2760941496/1971025*d^4 +
         * 814792828/1971025*d^3 - 1678512/1971025*d^2 - 9158924/1971025*d,
         *
         * c + 123812761248/1971025*d^7 - 79183342368/1971025*d^6 +
         * 7548646608/1971025*d^5 + 3840228724/1971025*d^4 -
         * 2024910556/5913075*d^3 - 132524276/5913075*d^2 + 30947828/5913075*d,
         *
         * d^8 - 8/11*d^7 + 4/33*d^6 + 131/5346*d^5 - 70/8019*d^4 + 1/3564*d^3 +
         * 5/42768*d^2 - 1/128304*d
         */
        doTest(new Ring("Z[d, c, b, a]"), new String[] {
            "a + 2*b + 2*c + 2*d - 1",
            "a^2 - a + 2*b^2 + 2*c^2 + 2*d^2",
            "2*a*b + 2*b*c - b + 2*c*d",
            "2*a*c + b^2 + 2*b*d - c"
        }, new String[] {
            "5913075a-159690237696d^7+31246269696d^6+27439610544d^5"
            + "-6475723368d^4-838935856d^3+275119624d^2+4884038d-5913075",
            "1971025b-97197721632d^7+73975630752d^6-12121915032d^5"
            + "-2760941496d^4+814792828d^3-1678512d^2-9158924d",
            "5913075c+371438283744d^7-237550027104d^6+22645939824d^5"
            + "+11520686172d^4-2024910556d^3-132524276d^2+30947828d",
            "128304d^8-93312d^7+15552d^6+3144d^5-1120d^4+36d^3+15d^2-d"
        });
    }

    @Test
    public void cyclic4() {
        /*
         * sage:
         *
         * R.<a, b, c, d> = PolynomialRing(ZZ, order='lex')
         *
         * I = R.ideal(a + b + c + d, a b + a d + b d + c d,
         * a b c + a b d + a c d + b c d, a b c d - 1) I.groebner_basis()
         *
         * [-a - b - d^5 - d, b^2 + b*d^5 + b*d + d^2, -c + d^5, -d^8 + 1]
         */
        doTest(new Ring("Z[d, c, b, a]"), new String[] {
            "a*b*c*d - 1",
            "a*b*c + a*b*d + a*c*d + b*c*d",
            "a*b + a*d + b*d + c*d",
            "a + b + c + d"
        }, new String[] {
            "a + b + d^5 + d",
            "b^2 + d^5 b + d b + d^2",
            "c - d^5",
            "d^8 - 1"
        });
    }

    @Test
    public void ex1() {
        /*
         * sage:
         *
         * R.<z, y, x> = PolynomialRing(ZZ, order='lex')
         * I = R.ideal(y^2 + 2 y x^2, y x + 2 x^3 - 1) I.groebner_basis()
         *
         * [y, -2*x^3 + 1]
         */
        doTest(new Ring("Z[x, y, z]"), new String[] {
            "y^2 + 2 y x^2",
            "y x + 2 x^3 - 1"
        }, new String[] {
            "y",
            "2x^3-1"
        });
    }

    @Test
    public void ex2() {
        doTest(new Ring("Z[x, y]"), new String[] {
            "x^2 - y^2",
            "x^2 + y"
        }, new String[] {
            "y + x^2",
            "x^4 - x^2"
        });
    }

    @Test
    public void ex3() {
        doTest(new Ring("Z[z, y, x]"), new String[] {
            "y - x^2",
            "z - x^3"
        }, new String[] {
            "x^2 - y",
            "y x - z",
            "z x - y^2",
            "y^3 - z^2"
        });
    }

    @Test
    public void mathematica1() {
        doTest(new Ring("Z[x, y, z]"), new String[] {
            "y^2 - 2 x^2",
            "x y - 3"
        }, new String[] {
            "3 y - 2 x^3",
            "2 x^4 - 9"
        });
    }

    @Test
    public void mathematica2() {
        doTest(new Ring("Z[x, y, z]"), new String[] {
            "x + y",
            "x^2 - 1",
            "y^2 - 2 x"
        }, new String[] {
            "1"
        });
    }

    @Test
    public void cox1() {
        doTest(new Ring("Z[w, z, y, x]"), new String[] {
            "3 x - 6 y - 2 z",
            "2 x - 4 y + 4 w",
            "x - 2 y - z - w"
        }, new String[] {
            "x - 2 y + 2 w",
            "z + 3 w"
        });
    }

    @Test
    public void cox2() {
        doTest(new Ring("Z[z, y, x]"), new String[] {
            "x^2 + y^2 + z^2 - 1",
            "x^2 + z^2 - y",
            "x - z"
        }, new String[] {
            "x-z",
            "y-2z^2",
            "4z^4+2z^2-1"
        });
    }

    @Test
    public void sage1() {
        doTest(new Ring("Z[z, y, x]"), new String[] {
            "x^6 + 4 x^4 y^2 + 4 x^2 y^4",
            "x^2 y^2"
        }, new String[] {
            "x^6",
            "y^2 x^2"
        });
    }

    @Test
    public void sage2() {
        /*
         * sage:
         *
         * R.<x, y, z> = PolynomialRing(ZZ, order='lex')
         * I = R.ideal(x^2 y - z, y^2 z - x, z^2 x - y) I.groebner_basis()
         *
         * [-x + z^9, y - z^11, z^15 - z]
         */
        doTest(new Ring("Z[z, y, x]"), new String[] {
            "x^2 y - z",
            "y^2 z - x",
            "z^2 x - y"
        }, new String[] {
            "x - z^9",
            "y - z^11",
            "z^15 - z"
        });
    }

    @Test
    public void sage3() {
        /*
         * sage:
         *
         * R.<x,y,z> = PolynomialRing(ZZ, order='lex')
         * I = R.ideal(z x + y^3, z + y^3, z+x y) I.groebner_basis()
         *
         * [x*y + z, x*z - z, y^3 + z, -y*z - z^2, -z^3 + z]
         */
        doTest(new Ring("Z[z, y, x]"), new String[] {
            "z x + y^3",
            "z + y^3",
            "z + x * y"
        }, new String[] {
            "z x - z",
            "y^3 + z",
            "y x + z",
            "z y + z^2",
            "z^3 - z"
        });
    }

    @Test
    public void sage001() {
        /*
         * sage
         *
         * [-4*x*y^2*z + y*z, -x*z + 21*y^2]
         *
         * [-x*z + 21*y^2, -420*y^4 + 5*y*z]
         */
        doTest(new Ring("Z[z, y, x]"), new String[] {
            "-4*x*y^2*z + y*z",
            "-x*z + 21*y^2"
        }, new String[] {
            "zx-21y^2",
            "84y^4-zy"
        });
    }

    @Test
    public void sage002() {
        /*
         * sage [x*y^2 - x*y*z^2, -z] [x*y^2, -z]
         */
        doTest(new Ring("Z[z, y, x]"), new String[] {
            "x*y^2 - x*y*z^2",
            "-z"
        }, new String[] {
            "y^2x",
            "z"
        });
    }

    @Test
    public void sage003() {
        /*
         * sage [-x^3, -16*x*z + y^2]
         * [-x^3, x^2*y^2, -x*y^4, -16*x*z + y^2, y^6]
         */
        doTest(new Ring("Z[z, y, x]"), new String[] {
            "-x^3",
            "-16*x*z + y^2"
        }, new String[] {
            "x^3",
            "16zx-y^2",
            "y^2x^2",
            "y^4x",
            "y^6"
        });
    }

    @Test
    public void sage004() {
        /*
         * sage [3*x^2*z + x*y, -2*z^2 - 7*z]
         * [-21*x^2*y + 2*x*y^2, 3*x^2*z + x*y, -2*x*y*z - 7*x*y, -2*z^2 - 7*z]
         */
        doTest(new Ring("Z[z, y, x]"), new String[] {
            "3*x^2*z + x*y",
            "-2*z^2 - 7*z"
        }, new String[] {
            "3zx^2 + yx",
            "2z^2 + 7z",
            "2zyx + 7yx",
            "21yx^2 - 2y^2x"
        });
    }

    @Test
    public void sage005() {
        /*
         * sage [x^2*z^2 - x, -70*x*y - 2*y]
         * [x^2*z^2 - x, -70*x*y - 2*y, -2*y*z^2 - 70*y]
         */
        doTest(new Ring("Z[z, y, x]"), new String[] {
            "x^2*z^2 - x",
            "-70*x*y - 2*y"
        }, new String[] {
            "z^2x^2-x",
            "35yx+y",
            "z^2y+35y"
        });
    }

    @Test
    public void sage006() {
        /*
         * sage [-x*y^2*z - 2*x*z^2 + 3*y*z^2, x] [x, 3*y*z^2]
         */
        doTest(new Ring("Z[z, y, x]"), new String[] {
            "-x*y^2*z - 2*x*z^2 + 3*y*z^2",
            "x"
        }, new String[] {
            "x",
            "z^2y"
        });
    }

    @Test
    public void sage007() {
        /*
         * sage [-x*y^3 + 6*x*y^2*z + z^3, y^2] [y^2, z^3]
         */
        doTest(new Ring("Z[z, y, x]"), new String[] {
            "-x*y^3 + 6*x*y^2*z + z^3",
            "y^2"
        }, new String[] {
            "y^2",
            "z^3"
        });
    }

    @Test
    public void sage008() {
        /*
         * sage [-10*x^2 + x, -2*z - 37] [-10*x^2 + x, -2*z - 37]
         */
        doTest(new Ring("Z[z, y, x]"), new String[] {
            "-10*x^2 + x",
            "-2*z - 37"
        }, new String[] {
            "10x^2-x",
            "2z+37"
        });
    }

    @Test
    public void sage009() {
        /*
         * sage [-3*y^4 - y^3*z, y^2 - z] [y^2 - z, -y*z^2 - 3*z^2, z^3 - 9*z^2]
         */
        doTest(new Ring("Z[z, y, x]"), new String[] {
            "-3*y^4 - y^3*z",
            "y^2 - z"
        }, new String[] {
            "y^2-z",
            "z^2y+3z^2",
            "z^3-9z^2"
        });
    }

    @Test
    public void sage010() {
        /*
         * sage [-x^2*z^2 + 5*x*y^2*z, 5*x*z - 3*z^2]
         * [-3*x*z^3 + 15*y^2*z^2, 5*x*z - 3*z^2, 75*y^2*z^2 - 9*z^4]
         */
        doTest(new Ring("Z[z, y, x]"), new String[] {
            "-x^2*z^2 + 5*x*y^2*z",
            "5*x*z - 3*z^2"
        }, new String[] {
            "5zx - 3z^2",
            "25z^2y^2 - 3z^4"
        });
    }

    @Test
    public void sage011() {
        /* Sage:
         * Ideal (1/4*y*x^2 + 1/5*y - 1/2*x^3, 11/2*y^2*x - 7/2*y^2 - 2*y*x)
         * of Multivariate Polynomial Ring in y, x over Rational Field
         [y - 9625/5832*x^6 + 19975/5832*x^5 - 275/729*x^4 - 5/2*x^3, x^7 - 9/11*x^6 - 8/55*x^4]
         */
        doTest(new Ring("Z[x, y]"), new String[] {
            "5*y*x^2 + 4*y - 10*x^3", "11*y^2*x - 7*y^2 - 4*y*x"
        }, new String[] {
            "5832y-9625x^6+19975x^5-2200x^4-14580x^3",
            "55x^7-45x^6-8x^4"
        });
    }

    @Test
    public void cocoa1() {
        doTest(new Ring("Z[z, y, x]"), new String[] {
            "x-y+z-2",
            "3x-z+6",
            "x+y-1"
        }, new String[] {
            "5x+3",
            "5y-8",
            "5z-21"
        });
    }

    @Test
    public void ourHelp01() {
        doTest(new Ring("Z[x, y, z]"), new String[] {
            "x^4y^3+2xy^2+3x+1",
            "x^3y^2+x^2",
            "x^4y+z^2+x y^4+3"
        }, new String[] {
            "z^2-x^4+3x^2-10x+9",
            "y-9x^4-3x^3-x^2-81x+27",
            "x^5+9x^2-6x+1"
        });
    }
}
