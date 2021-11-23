package com.mathpar.polynom.gbasis.datastructures;

import com.mathpar.polynom.gbasis.util.Utils;
import java.util.Arrays;
import java.util.List;
import com.mathpar.number.Element;
import com.mathpar.matrix.MatrixS;
import com.mathpar.number.NumberZ;
import com.mathpar.number.Ring;
import static org.junit.Assert.*;
import org.junit.Test;
import com.mathpar.polynom.Polynom;
import com.mathpar.polynom.Term;

public class PolyListDefaultTest {
    private static final Ring rZxyz = new Ring("Z[x, y, z]"); // z > y > x
    private static final List<Polynom> pols1 = Utils.polList(rZxyz,
            "2z + yx^2",
            "z^3 - 3y^2",
            "-10z");
    private static final TermsContainer terms1 = new TermsList(rZxyz,
            Arrays.asList(
            new Term(0, 0, 3), // z^3   #0
            new Term(0, 0, 1), // z     #1
            new Term(0, 2), //    y^2   #2
            new Term(2, 1) //     yx^2  #3
            ), true);
    private static final Element[][] rows1 = new Element[][] {
        new Element[] {
            new NumberZ("2"),
            new NumberZ("1")
        },
        new Element[] {
            new NumberZ("1"),
            new NumberZ("-3")
        },
        new Element[] {
            new NumberZ("-10")
        }
    };
    private static final MatrixS matr1 = new MatrixS(rows1, new int[][] {
        new int[] {1, 3},
        new int[] {0, 2},
        new int[] {1}
    });

    @Test
    public void constructFromListOfPolynomials() {
        PolyListDefault plist = new PolyListDefault(rZxyz, pols1);
        assertEquals(terms1, plist.terms());
        assertArrayEquals(rows1, plist.rows());
    }

    @Test
    public void toMatrixS() {
        PolyListDefault plist = new PolyListDefault(rZxyz, pols1);
        MatrixS m1 = plist.toMatrixS();
        assertTrue(matr1.equals(m1, rZxyz));
    }

    @Test
    public void getPolynomial() {
        PolyListDefault plist = new PolyListDefault(rZxyz, pols1);
        assertTrue(pols1.get(0).equals(plist.get(0), rZxyz));
    }

    @Test
    public void toPolynomialList() {
        PolyListDefault plist = new PolyListDefault(rZxyz, pols1);
        List<Polynom> res = plist.toList();
        assertTrue(Element.equals(rZxyz, pols1, res));
    }

    @Test
    public void addUpdatesTermsListAndTermsIndexes() {
        PolyListDefault plist = new PolyListDefault(rZxyz, pols1);
        // This adds new terms:
        // z^2      #1
        // y        #5
        plist.add(new Polynom("5z^2 + 3y^2 + 10y", rZxyz));

        // Terms
        assertEquals(new TermsList(rZxyz, Arrays.asList(
                new Term(0, 0, 3), // z^3 #0
                new Term(0, 0, 2), // z^2 #1
                new Term(0, 0, 1), // z #2
                new Term(0, 2), // y^2 #3
                new Term(2, 1), // yx^2 #4
                new Term(0, 1) // y #5
                ), true), plist.terms());
        // Matrix
        assertTrue(new MatrixS(new Element[][] {
            new Element[] {new NumberZ(2), new NumberZ(1)},
            new Element[] {new NumberZ(1), new NumberZ(-3)},
            new Element[] {new NumberZ(-10)},
            new Element[] {new NumberZ(5), new NumberZ(3), new NumberZ(10)}
        }, new int[][] {
            new int[] {2, 4},
            new int[] {0, 3},
            new int[] {2},
            new int[] {1, 3, 5}
        }).equals(plist.toMatrixS(), rZxyz));
    }

    @Test
    public void addAllUpdatesTermsListAndTermsIndexes() {
        PolyListDefault plist = new PolyListDefault(rZxyz, pols1);
        plist.addAll(Utils.polList(rZxyz,
                "5z^2 + 3y^2 + 10y",
                "z^5 + z^4 + z^3 + y^2 + x^2 + x"));

        // Terms
        assertEquals(new TermsList(rZxyz, Arrays.asList(
                new Term(0, 0, 5), // z^5 #0
                new Term(0, 0, 4), // z^4 #1
                new Term(0, 0, 3), // z^3 #2
                new Term(0, 0, 2), // z^2 #3
                new Term(0, 0, 1), // z #4
                new Term(0, 2), // y^2 #5
                new Term(2, 1), // yx^2 #6
                new Term(0, 1), // y #7
                new Term(2), // x^2 #8
                new Term(1) // x #9
                ), true), plist.terms());
        // Matrix
        assertTrue(new MatrixS(new Element[][] {
            new Element[] {new NumberZ(2), new NumberZ(1)},
            new Element[] {new NumberZ(1), new NumberZ(-3)},
            new Element[] {new NumberZ(-10)},
            new Element[] {
                new NumberZ(5),
                new NumberZ(3),
                new NumberZ(10)
            },
            new Element[] {
                new NumberZ(1),
                new NumberZ(1),
                new NumberZ(1),
                new NumberZ(1),
                new NumberZ(1),
                new NumberZ(1)
            }
        }, new int[][] {
            new int[] {4, 6},
            new int[] {2, 5},
            new int[] {4},
            new int[] {3, 5, 7},
            new int[] {0, 1, 2, 5, 8, 9}
        }).equals(plist.toMatrixS(), rZxyz));
    }

    @Test
    public void multiplyPolynomialOnTerm() {
        PolyListDefault plist = new PolyListDefault(rZxyz, pols1);
        plist.multiply(0, new Term(0, 0, 2) // z^2
                );

        // Terms
        assertEquals(new TermsList(rZxyz, Arrays.asList(
                new Term(0, 0, 3), // z^3       #0
                new Term(2, 1, 2), // z^2yx^2   #1
                new Term(0, 0, 1), // z         #2
                new Term(0, 2), // y^2          #3
                new Term(2, 1) // yx^2         #4
                ), true), plist.terms());

        // Matrix
        assertTrue(new MatrixS(rows1,
                new int[][] {
            new int[] {0, 1},
            new int[] {0, 3},
            new int[] {2}
        }).equals(plist.toMatrixS(), rZxyz));
    }

    @Test
    public void ht() {
        PolyListDefault plist = new PolyListDefault(rZxyz, pols1);
        // ht of z^3 - 3y^2 == z^3
        assertEquals(new Term(0, 0, 3), plist.ht(1));
    }

    @Test
    public void lcmHT() {
        PolyListDefault plist = new PolyListDefault(rZxyz, pols1);
        // lcm of z and z^3 == z^3
        assertEquals(new Term(0, 0, 3), plist.lcmHT(0, 1));
    }

    @Test
    public void copyPols() {
        PolyListDefault src = new PolyListDefault(rZxyz, pols1);
        PolyListDefault dest = new PolyListDefault(rZxyz);

        PolyListDefault.appendPols(src, dest, 0, 2);
        assertEquals(new TermsList(rZxyz, Arrays.asList(
                new Term(0, 0, 1), // z #0
                new Term(2, 1) // yx^2 #1
                ), true), dest.terms());
        assertTrue(new MatrixS(new Element[][] {
            new Element[] {new NumberZ("2"), new NumberZ("1")},
            new Element[] {new NumberZ("-10")}
        }, new int[][] {
            new int[] {0, 1},
            new int[] {0}
        }).equals(dest.toMatrixS(), rZxyz));
    }

    @Test
    public void addLR() {
        PolyListDefault plist = new PolyListDefault(rZxyz, pols1);
        PolyListDefault lr = new PolyListDefault(rZxyz);
        // Compute pair of 2z + yx^2 and z^3 - 3y^2
        // L = 2z^3 + z^2yx^2; R = z^3 - 3y^2
        PolyListDefault.computeAndAddLRPair(plist, lr, 0, 1);
        assertEquals(new TermsList(rZxyz, Arrays.asList(
                new Term(0, 0, 3), // z^3 # 0
                new Term(2, 1, 2), // z^2yx^2 # 1
                new Term(0, 0, 1), // z # 2
                new Term(0, 2), // y^2 # 3
                new Term(2, 1) // yx^2 # 4
                ), true), lr.terms());
        assertTrue(new MatrixS(new Element[][] {
            new Element[] {new NumberZ("2"), NumberZ.ONE},
            new Element[] {NumberZ.ONE, new NumberZ("-3")}
        }, new int[][] {
            new int[] {0, 1},
            new int[] {0, 3}
        }).equals(lr.toMatrixS(), rZxyz));
    }

    @Test
    public void copyToAnotherAndRemovePols() {
        PolyListDefault plist = new PolyListDefault(rZxyz, pols1);
        PolyListDefault dest = new PolyListDefault(rZxyz);
        plist.copyToAndRemovePols(dest, 0, 1);

        assertTrue(Element.equals(rZxyz,
                Arrays.asList(pols1.get(2)), plist.toList()));
    }
}
