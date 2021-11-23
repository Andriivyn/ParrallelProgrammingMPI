package com.mathpar.polynom.gbasis.algorithms;

import com.mathpar.polynom.gbasis.util.Utils;
import com.mathpar.number.Element;
import com.mathpar.number.Ring;
import org.junit.Test;

import static org.junit.Assert.assertTrue;

public class BuchbergerTest {
    private void doTest(Ring r, String[] ideal, String[] exp) {
        assertTrue(Element.equals(r,
                Utils.polList(r, exp),
                new Buchberger(r, ideal).gbasis()));
    }

    @Test
    public void katsura3() {
        doTest(new Ring("Q[z, y, x]"),
                new String[] {
            "x + 2*y + 2*z - 1",
            "x^2 - x + 2*y^2 + 2*z^2",
            "2*x*y + 2*y*z - y"
        }, new String[] {
            "x-60z^3+(158/7)z^2+(8/7)z-1",
            "y+30z^3+(-79/7)z^2+(3/7)z",
            "z^4+(-10/21)z^3+(1/84)z^2+(1/84)z"
        });
    }

    @Test
    public void katsura4() {
        doTest(new Ring("Q[d, c, b, a]"),
                new String[] {
            "a + 2*b + 2*c + 2*d - 1",
            "a^2 - a + 2*b^2 + 2*c^2 + 2*d^2",
            "2*a*b + 2*b*c - b + 2*c*d",
            "2*a*c + b^2 + 2*b*d - c"
        }, new String[] {
            "a+(-53230079232/1971025)d^7+(10415423232/1971025)d^6+(9146536848/1971025)d^5+(-2158574456/1971025)d^4+(-838935856/5913075)d^3+(275119624/5913075)d^2+(4884038/5913075)d-1",
            "b+(-97197721632/1971025)d^7+(73975630752/1971025)d^6+(-12121915032/1971025)d^5+(-2760941496/1971025)d^4+(814792828/1971025)d^3+(-1678512/1971025)d^2+(-9158924/1971025)d",
            "c+(123812761248/1971025)d^7+(-79183342368/1971025)d^6+(7548646608/1971025)d^5+(3840228724/1971025)d^4+(-2024910556/5913075)d^3+(-132524276/5913075)d^2+(30947828/5913075)d",
            "d^8+(-8/11)d^7+(4/33)d^6+(131/5346)d^5+(-70/8019)d^4+(1/3564)d^3+(5/42768)d^2+(-1/128304)d"
        });
    }

    @Test
    public void cyclic4() {
        doTest(new Ring("Q[d, c, b, a]"),
                new String[] {
            "abcd - 1",
            "abc + abd + acd + bcd",
            "ab + ad + bd + cd",
            "a + b + c + d"
        }, new String[] {
            "a + b + d^5 + d",
            "b^2 + d^5b + db + d^2",
            "c - d^5",
            "d^8 - 1"
        });
    }

    @Test
    public void ex1() {
        doTest(new Ring("Q[z, y, x]"),
                new String[] {
            "x^2 + 2xy^2",
            "xy + 2y^3 -1"
        }, new String[] {
            "x",
            "y^3+(-1/2)"
        });
    }

    @Test
    public void ex2() {
        doTest(new Ring("Q[x, y]"),
                new String[] {
            "x^2 - y^2",
            "x^2 + y"
        }, new String[] {
            "y+x^2",
            "x^4-x^2"
        });
    }

    @Test
    public void ex3() {
        doTest(new Ring("Q[z, y, x]"),
                new String[] {
            "y - x^2",
            "z - x^3"
        }, new String[] {
            "x^2-y",
            "yx-z",
            "zx-y^2",
            "y^3-z^2"
        });
    }

    @Test
    public void mathematica1() {
        doTest(new Ring("Q[x, y, z]"),
                new String[] {
            "y^2 - 2x^2",
            "xy - 3"
        }, new String[] {
            "y-2/3x^3",
            "x^4-9/2"
        });
    }

    @Test
    public void mathematica2() {
        doTest(new Ring("Q[x, y, z]"),
                new String[] {
            "x + y",
            "x^2 - 1",
            "y^2 - 2x"
        }, new String[] {
            "1"
        });
    }

    @Test
    public void cox1() {
        doTest(new Ring("Q[w, z, y, x]"),
                new String[] {
            "3x - 6y - 2z",
            "2x - 4y + 4w",
            "x - 2y - z - w"
        }, new String[] {
            "x-2y+2w",
            "z+3w"
        });
    }

    @Test
    public void cox2() {
        doTest(new Ring("Q[z, y, x]"),
                new String[] {
            "x^2 + y^2 + z^2 - 1",
            "x^2 + z^2 - y",
            "x - z"
        },
                new String[] {
            "x-z",
            "y-2z^2",
            "z^4+(1/2)z^2+(-1/4)"
        });
    }

    @Test
    public void sage1() {
        doTest(new Ring("Q[z, y, x]"),
                new String[] {
            "x^6 + 4*x^4*y^2 + 4*x^2*y^4",
            "x^2*y^2"
        },
                new String[] {
            "x^6",
            "y^2x^2"
        });
    }

    @Test
    public void sage2() {
        doTest(new Ring("Q[z, y, x]"),
                new String[] {
            "x^2*y - z",
            "y^2*z - x",
            "z^2*x - y"
        },
                new String[] {
            "x-z^9",
            "y-z^11",
            "z^15-z"
        });
    }

    @Test
    public void sage3() {
        /*
         * sage:
         *
         * R.<x,y,z> = PolynomialRing(ZZ, order='lex') I =
         * R.ideal(z*x+y^3,z+y^3,z+x*y) I.groebner_basis() [x*y + z, x*z - z,
         * y^3 + z, -y*z - z^2, -z^3 + z]
         */
        doTest(new Ring("Q[z, y, x]"),
                new String[] {
            "z*x+y^3",
            "z+y^3",
            "z+x*y"
        },
                new String[] {
            "zx-z",
            "y^3+z",
            "yx+z",
            "zy+z^2",
            "z^3-z"
        });
    }

    @Test
    public void sage001() {
        /*
         * sage [-4*x*y^2*z + y*z, -x*z + 21*y^2] [-x*z + 21*y^2, -420*y^4 +
         * 5*y*z]
         */
        doTest(new Ring("Q[z, y, x]"),
                new String[] {
            "-4*x*y^2*z + y*z",
            "-x*z + 21*y^2"
        },
                new String[] {
            "zx-21y^2",
            "y^4+(-1/84)zy"
        });
    }

    @Test
    public void sage002() {
        /*
         * sage [x*y^2 - x*y*z^2, -z] [x*y^2, -z]
         */
        doTest(new Ring("Q[z, y, x]"),
                new String[] {
            "x*y^2 - x*y*z^2",
            "-z"
        },
                new String[] {
            "y^2x",
            "z"
        });
    }

    @Test
    public void sage003() {
        /*
         * sage [-x^3, -16*x*z + y^2] [-x^3, x^2*y^2, -x*y^4, -16*x*z + y^2,
         * y^6]
         */
        doTest(new Ring("Q[z, y, x]"),
                new String[] {
            "-x^3",
            "-16*x*z + y^2"
        },
                new String[] {
            "x^3",
            "zx+(-1/16)y^2",
            "y^2x^2",
            "y^4x",
            "y^6"
        });
    }

    @Test
    public void sage004() {
        /*
         * sage [3*x^2*z + x*y, -2*z^2 - 7*z] [-21*x^2*y + 2*x*y^2, 3*x^2*z +
         * x*y, -2*x*y*z - 7*x*y, -2*z^2 - 7*z]
         */
        doTest(new Ring("Q[z, y, x]"),
                new String[] {
            "3*x^2*z + x*y",
            "-2*z^2 - 7*z"
        },
                new String[] {
            "zx^2+(1/3)yx",
            "z^2+(7/2)z",
            "zyx+(7/2)yx",
            "yx^2+(-2/21)y^2x"
        });
    }

    @Test
    public void sage005() {
        /*
         * sage [x^2*z^2 - x, -70*x*y - 2*y] [x^2*z^2 - x, -70*x*y - 2*y,
         * -2*y*z^2 - 70*y]
         */
        doTest(new Ring("Q[z, y, x]"),
                new String[] {
            "x^2*z^2 - x",
            "-70*x*y - 2*y"
        },
                new String[] {
            "z^2x^2-x",
            "yx+(1/35)y",
            "z^2y+35y"
        });
    }

    @Test
    public void sage006() {
        /*
         * sage [-x*y^2*z - 2*x*z^2 + 3*y*z^2, x] [x, 3*y*z^2]
         */
        doTest(new Ring("Q[z, y, x]"),
                new String[] {
            "-x*y^2*z - 2*x*z^2 + 3*y*z^2",
            "x"
        },
                new String[] {
            "x",
            "z^2y"
        });
    }

    @Test
    public void sage007() {
        /*
         * sage [-x*y^3 + 6*x*y^2*z + z^3, y^2] [y^2, z^3]
         */
        doTest(new Ring("Q[z, y, x]"),
                new String[] {
            "-x*y^3 + 6*x*y^2*z + z^3",
            "y^2"
        },
                new String[] {
            "y^2",
            "z^3"
        });
    }

    @Test
    public void sage008() {
        /*
         * sage [-10*x^2 + x, -2*z - 37] [-10*x^2 + x, -2*z - 37]
         */
        doTest(new Ring("Q[z, y, x]"),
                new String[] {
            "-10*x^2 + x",
            "-2*z - 37"
        },
                new String[] {
            "x^2+(-1/10)x",
            "z+(37/2)"
        });
    }

    @Test
    public void sage009() {
        /*
         * sage [-3*y^4 - y^3*z, y^2 - z] [y^2 - z, -y*z^2 - 3*z^2, z^3 - 9*z^2]
         */
        doTest(new Ring("Q[z, y, x]"),
                new String[] {
            "-3*y^4 - y^3*z",
            "y^2 - z"
        },
                new String[] {
            "y^2-z",
            "z^2y+3z^2",
            "z^3-9z^2"
        });
    }

    @Test
    public void sage010() {
        /*
         * sage [-x^2*z^2 + 5*x*y^2*z, 5*x*z - 3*z^2] [-3*x*z^3 + 15*y^2*z^2,
         * 5*x*z - 3*z^2, 75*y^2*z^2 - 9*z^4]
         */
        doTest(new Ring("Q[z, y, x]"),
                new String[] {
            "-x^2*z^2 + 5*x*y^2*z",
            "5*x*z - 3*z^2"
        },
                new String[] {
            "zx+(-3/5)z^2",
            "z^2y^2+(-3/25)z^4"
        });
    }

    @Test
    public void cocoa1() {
        doTest(new Ring("Q[z, y, x]"),
                new String[] {
            "x-y+z-2",
            "3x-z+6",
            "x+y-1"
        },
                new String[] {
            "x+(3/5)",
            "y+(-8/5)",
            "z+(-21/5)"
        });
    }
}
