package com.mathpar.polynom.gbasis.algorithms;

import com.mathpar.number.*;
import com.mathpar.polynom.Polynom;
import org.junit.*;

import java.util.*;

import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

/**
 * Contains few examples to test polynomial systems solving with Groebner bases.
 */
public class SNLESolverTest {
    public SNLESolverTest() {
    }

    @BeforeClass
    public static void setUpClass() {
    }

    @AfterClass
    public static void tearDownClass() {
    }

    @Before
    public void setUp() {
    }

    @After
    public void tearDown() {
    }
  
    /**
     * Checks if given unordered lists of Maps which hold solutions of
     * polynomial system are equal. Position of corresponding roots may differ
     * due to random nature of polynomial's roots calculation.
     *
     * @param ring ring
     * @param exp  expected solutions
     * @param got  got solutions
     * @return {@code true} if lists are equal, {@code false} otherwise
     */
    private static boolean fuzzyEqual(Ring ring, List<SortedMap<Polynom, Element>> exp,
                                      List<SortedMap<Polynom, Element>> got) {
        int sz = exp.size();
        int equalMapsCnt = 0;
        if (sz != got.size()) {
            return false;
        }
        for (int i = 0; i < sz; i++) {
            SortedMap<Polynom, Element> currRoots = exp.get(i);
            for (int j = 0; j < sz; j++) {
                SortedMap<Polynom, Element> gotRoots = got.get(j);
                if (mapsEqual(ring, currRoots, gotRoots)) {
                    equalMapsCnt++;
                    break;
                }
            }
        }
        return equalMapsCnt == sz;
    }

    /**
     * Checks if given Maps are equal. Uses Element-specific equals() method
     * (with Ring) if needed.
     *
     * @param ring   ring
     * @param first  first map
     * @param second second map
     * @return {@code true} if maps are equal {@code false} otherwise
     */
    private static boolean mapsEqual(Ring ring, Map<Polynom, Element> first,
                                     SortedMap<Polynom, Element> second) {
        if (first.size() != second.size()) {
            return false;
        }
        Set<Polynom> expKeys = first.keySet();
        for (Iterator<Polynom> itExp = expKeys.iterator(),
                     itGot = second.keySet().iterator();
             itExp.hasNext(); ) {
            Polynom keyExp = itExp.next();
            Polynom keyGot = itGot.next();
            Element valExp = first.get(keyExp);
            Element valGot = second.get(keyGot);
            boolean keysEqual = keyExp.equals(keyGot, ring);
            boolean valsEqual = valExp.equals(valGot, ring);
            if (valGot instanceof Complex && valExp instanceof Complex) {
                Complex cValExp = (Complex) valExp;
                Complex cValGot = (Complex) valGot;
                boolean reEqual = cValExp.re.equals(cValGot.re, ring);
                boolean imEqual = cValExp.im.equals(cValGot.im, ring);
                valsEqual = reEqual && imEqual;
            }
            if (!(keysEqual && valsEqual)) {
                return false;
            }
        }
        return true;
    }

    // Test testing methods - hell, yeah!
    @Test
    public void elementsMapEquals() {
        Ring ring = new Ring("R[x]");
        SortedMap<Polynom, Element> one = new TreeMap<>();
        one.put(new Polynom("x", ring), new NumberR("1.0"));
        one.put(new Polynom("x^2", ring), new NumberR("2.0"));
        SortedMap<Polynom, Element> two = new TreeMap<>();
        two.put(new Polynom("x", ring), new NumberR("1.0"));
        two.put(new Polynom("x^2", ring), new NumberR("2.0"));
        assertTrue(mapsEqual(ring, one, two));
    }

    @Test
    public void elementsMapNotEquals() {
        Ring ring = new Ring("R[x]");
        SortedMap<Polynom, Element> one = new TreeMap<>();
        one.put(new Polynom("x", ring), new NumberR("1.0"));
        one.put(new Polynom("x^2", ring), new NumberR("2.0"));
        SortedMap<Polynom, Element> two = new TreeMap<>();
        two.put(new Polynom("x", ring), new NumberR("1.0"));
        two.put(new Polynom("x^2", ring), new NumberR("3.0"));
        assertFalse(mapsEqual(ring, one, two));
    }

    @Test
    public void example01() {
        /*
         * (x, y)
         *
         * (0, 1); (0.8, -0.6)
         */
        final Ring r = new Ring("R[y, x]");
        List<SortedMap<Polynom, Element>> res = new SNLESolver(r,
                "x^2 + y^2 - 1",
                "(x - 2)^2 + (y - 1)^2 - 4").solve();

        List<SortedMap<Polynom, Element>> expecting = new ArrayList<>();
        SortedMap<Polynom, Element> s1 = new TreeMap<>();
        s1.put(new Polynom("x", r), new Complex(new NumberR("0.0"), NumberR.ZERO));
        s1.put(new Polynom("y", r), new Complex(new NumberR("1.0"), NumberR.ZERO));
        SortedMap<Polynom, Element> s2 = new TreeMap<>();
        s2.put(new Polynom("x", r), new Complex(new NumberR("0.8"), NumberR.ZERO));
        s2.put(new Polynom("y", r), new Complex(new NumberR("-0.6"), NumberR.ZERO));
        expecting.add(s1);
        expecting.add(s2);

        assertTrue(fuzzyEqual(r, expecting, res));
    }

    @Test
    public void example02() {
        /*
         * (x, y)
         *
         * (-1.89005, 0.653986); (0.653986, -1.89005);
         * (1.61803 - -0.78615i, 1.61803 + -0.78615i); (1.61803 + -0.78615i, 1.61803 - -0.78615i)
         */
        Ring r = new Ring("R[y, x]");
        List<SortedMap<Polynom, Element>> res = new SNLESolver(r,
                "x^2 + y^2 - 4",
                "(x - 1) (y - 1) - 1").solve();

        List<SortedMap<Polynom, Element>> expecting = new ArrayList<>();
        SortedMap<Polynom, Element> s1 = new TreeMap<>();
        s1.put(new Polynom("x", r), new Complex(new NumberR("-1.89005"), NumberR.ZERO));
        s1.put(new Polynom("y", r), new Complex(new NumberR("0.653986"), NumberR.ZERO));
        SortedMap<Polynom, Element> s2 = new TreeMap<>();
        s2.put(new Polynom("x", r), new Complex(new NumberR("0.653986"), NumberR.ZERO));
        s2.put(new Polynom("y", r), new Complex(new NumberR("-1.89005"), NumberR.ZERO));
        SortedMap<Polynom, Element> s3 = new TreeMap<>();
        s3.put(new Polynom("x", r), new Complex(new NumberR("1.61803"), new NumberR("-0.78615")));
        s3.put(new Polynom("y", r), new Complex(new NumberR("1.61803"), new NumberR("0.78615")));
        SortedMap<Polynom, Element> s4 = new TreeMap<>();
        s4.put(new Polynom("x", r), new Complex(new NumberR("1.61803"), new NumberR("0.78615")));
        s4.put(new Polynom("y", r), new Complex(new NumberR("1.61803"), new NumberR("-0.78615")));
        expecting.add(s1);
        expecting.add(s2);
        expecting.add(s3);
        expecting.add(s4);

        r.MachineEpsilonR = new NumberR("0.001");
        assertTrue(fuzzyEqual(r, expecting, res));
    }

    @Test
    public void example03() {
        /*
         * (a, b, c)
         *
         * 1) (1, -0.5 - 0.866025i, -0.5 + 0.866025i);
         * 2) (1, -0.5 + 0.866025i, -0.5 - 0.866025i);
         * 3) (-0.5 - 0.866025i, 1, -0.5 + 0.866025i);
         * 4) (-0.5 - 0.866025i, -0.5 + 0.866025i, 1);
         * 5) (-0.5 + 0.866025i, 1, -0.5 - 0.866025i);
         * 6) (-0.5 + 0.866025i, -0.5 - 0.866025i, 1);
         */
        Ring r = new Ring("R[a, b, c]");
        List<SortedMap<Polynom, Element>> res = new SNLESolver(r,
                "a + b + c",
                "a b + a c + b c",
                "a b c - 1").solve();

        List<SortedMap<Polynom, Element>> expecting = new ArrayList<>();
        SortedMap<Polynom, Element> s1 = new TreeMap<>();
        s1.put(new Polynom("a", r), new Complex(new NumberR("1.0"), NumberR.ZERO));
        s1.put(new Polynom("b", r), new Complex(new NumberR("-0.5"), new NumberR("-0.866025")));
        s1.put(new Polynom("c", r), new Complex(new NumberR("-0.5"), new NumberR("0.866025")));
        SortedMap<Polynom, Element> s2 = new TreeMap<>();
        s2.put(new Polynom("a", r), new Complex(new NumberR("1.0"), NumberR.ZERO));
        s2.put(new Polynom("b", r), new Complex(new NumberR("-0.5"), new NumberR("0.866025")));
        s2.put(new Polynom("c", r), new Complex(new NumberR("-0.5"), new NumberR("-0.866025")));
        SortedMap<Polynom, Element> s3 = new TreeMap<>();
        s3.put(new Polynom("a", r), new Complex(new NumberR("-0.5"), new NumberR("-0.866025")));
        s3.put(new Polynom("b", r), new Complex(new NumberR("1.0"), new NumberR("0.0")));
        s3.put(new Polynom("c", r), new Complex(new NumberR("-0.5"), new NumberR("0.866025")));
        SortedMap<Polynom, Element> s4 = new TreeMap<>();
        s4.put(new Polynom("a", r), new Complex(new NumberR("-0.5"), new NumberR("-0.866025")));
        s4.put(new Polynom("b", r), new Complex(new NumberR("-0.5"), new NumberR("0.866025")));
        s4.put(new Polynom("c", r), new Complex(new NumberR("1.0"), new NumberR("0.0")));
        SortedMap<Polynom, Element> s5 = new TreeMap<>();
        s5.put(new Polynom("a", r), new Complex(new NumberR("-0.5"), new NumberR("0.866025")));
        s5.put(new Polynom("b", r), new Complex(new NumberR("1.0"), new NumberR("0.0")));
        s5.put(new Polynom("c", r), new Complex(new NumberR("-0.5"), new NumberR("-0.866025")));
        SortedMap<Polynom, Element> s6 = new TreeMap<>();
        s6.put(new Polynom("a", r), new Complex(new NumberR("-0.5"), new NumberR("0.866025")));
        s6.put(new Polynom("b", r), new Complex(new NumberR("-0.5"), new NumberR("-0.866025")));
        s6.put(new Polynom("c", r), new Complex(new NumberR("1.0"), new NumberR("0.0")));
        expecting.add(s1);
        expecting.add(s2);
        expecting.add(s3);
        expecting.add(s4);
        expecting.add(s5);
        expecting.add(s6);

        r.FLOATPOS = 6;
        r.MachineEpsilonR = new NumberR("1e-3");
        assertTrue(fuzzyEqual(r, expecting, res));
    }

    @Test
    public void example04() {
        /*
         * (x, y, z)
         *
         * (0, -1, -1); (0, 1, -1); (3, -2, -4); (3, 2, -4); (8, -3, -9); (8, 3, -9)
         */
        Ring r = new Ring("R[z, y, x]");
        List<SortedMap<Polynom, Element>> res = new SNLESolver(r,
                "x^2 + x z + y^2 - 1",
                "x^3 - 11 x^2 + 23 x + y^2 - 1",
                "x + z + 1").solve();

        List<SortedMap<Polynom, Element>> expecting = new ArrayList<>();
        SortedMap<Polynom, Element> s1 = new TreeMap<>();
        s1.put(new Polynom("x", r), new Complex(new NumberR("0.0"), NumberR.ZERO));
        s1.put(new Polynom("y", r), new Complex(new NumberR("-1.0"), NumberR.ZERO));
        s1.put(new Polynom("z", r), new Complex(new NumberR("-1.0"), NumberR.ZERO));
        SortedMap<Polynom, Element> s2 = new TreeMap<>();
        s2.put(new Polynom("x", r), new Complex(new NumberR("0.0"), NumberR.ZERO));
        s2.put(new Polynom("y", r), new Complex(new NumberR("1.0"), NumberR.ZERO));
        s2.put(new Polynom("z", r), new Complex(new NumberR("-1.0"), NumberR.ZERO));
        SortedMap<Polynom, Element> s3 = new TreeMap<>();
        s3.put(new Polynom("x", r), new Complex(new NumberR("3.0"), NumberR.ZERO));
        s3.put(new Polynom("y", r), new Complex(new NumberR("-2.0"), NumberR.ZERO));
        s3.put(new Polynom("z", r), new Complex(new NumberR("-4.0"), NumberR.ZERO));
        SortedMap<Polynom, Element> s4 = new TreeMap<>();
        s4.put(new Polynom("x", r), new Complex(new NumberR("3.0"), NumberR.ZERO));
        s4.put(new Polynom("y", r), new Complex(new NumberR("2.0"), NumberR.ZERO));
        s4.put(new Polynom("z", r), new Complex(new NumberR("-4.0"), NumberR.ZERO));
        SortedMap<Polynom, Element> s5 = new TreeMap<>();
        s5.put(new Polynom("x", r), new Complex(new NumberR("8.0"), NumberR.ZERO));
        s5.put(new Polynom("y", r), new Complex(new NumberR("-3.0"), NumberR.ZERO));
        s5.put(new Polynom("z", r), new Complex(new NumberR("-9.0"), NumberR.ZERO));
        SortedMap<Polynom, Element> s6 = new TreeMap<>();
        s6.put(new Polynom("x", r), new Complex(new NumberR("8.0"), NumberR.ZERO));
        s6.put(new Polynom("y", r), new Complex(new NumberR("3.0"), NumberR.ZERO));
        s6.put(new Polynom("z", r), new Complex(new NumberR("-9.0"), NumberR.ZERO));
        expecting.add(s1);
        expecting.add(s2);
        expecting.add(s3);
        expecting.add(s4);
        expecting.add(s5);
        expecting.add(s6);

        r.MachineEpsilonR = new NumberR("0.001");
        assertTrue(fuzzyEqual(r, expecting, res));
    }

    @Test
    public void example05() {
        /*
         * (x, y)
         *
         * (0.828427, -0.643594); (0.828427, 0.643594); (-4.82843, 0.0 -1.55377i); (-4.82843, 0.0 + 1.55377i)
         */
        Ring r = new Ring("R[y, x]");
        List<SortedMap<Polynom, Element>> res = new SNLESolver(r,
                "x^2 + 4 x - 4",
                "2 y^2 - x").solve();

        List<SortedMap<Polynom, Element>> expecting = new ArrayList<>();
        SortedMap<Polynom, Element> s1 = new TreeMap<>();
        s1.put(new Polynom("x", r), new Complex(new NumberR("0.828427"), NumberR.ZERO));
        s1.put(new Polynom("y", r), new Complex(new NumberR("-0.643594"), NumberR.ZERO));
        SortedMap<Polynom, Element> s2 = new TreeMap<>();
        s2.put(new Polynom("x", r), new Complex(new NumberR("0.828427"), NumberR.ZERO));
        s2.put(new Polynom("y", r), new Complex(new NumberR("0.643594"), NumberR.ZERO));
        SortedMap<Polynom, Element> s3 = new TreeMap<>();
        s3.put(new Polynom("x", r), new Complex(new NumberR("-4.82843"), NumberR.ZERO));
        s3.put(new Polynom("y", r), new Complex(NumberR.ZERO, new NumberR("-1.55377")));
        SortedMap<Polynom, Element> s4 = new TreeMap<>();
        s4.put(new Polynom("x", r), new Complex(new NumberR("-4.82843"), NumberR.ZERO));
        s4.put(new Polynom("y", r), new Complex(NumberR.ZERO, new NumberR("1.55377")));
        expecting.add(s1);
        expecting.add(s2);
        expecting.add(s3);
        expecting.add(s4);

        r.MachineEpsilonR = new NumberR("0.001");
        assertTrue(fuzzyEqual(r, expecting, res));
    }

    @Test
    public void example06() {
        /*
         * (c, s)
         *
         * (-0.328074, -0.944652); (0.328074, 0.944652); (-0.837408, -0.546579);
         * (0.837408, 0.546579); (-0.882809, 0.469733); (0.882809, -0.469733)
         */
        Ring r = new Ring("R[c, s]");
        List<SortedMap<Polynom, Element>> res = new SNLESolver(r,
                "s^3 + 4 c^3 - 3c",
                "s^2 + c^2 - 1").solve();

        List<SortedMap<Polynom, Element>> expecting = new ArrayList<>();
        SortedMap<Polynom, Element> s1 = new TreeMap<>();
        s1.put(new Polynom("c", r), new Complex(new NumberR("-0.328074"), NumberR.ZERO));
        s1.put(new Polynom("s", r), new Complex(new NumberR("-0.944652"), NumberR.ZERO));
        SortedMap<Polynom, Element> s2 = new TreeMap<>();
        s2.put(new Polynom("c", r), new Complex(new NumberR("0.328074"), NumberR.ZERO));
        s2.put(new Polynom("s", r), new Complex(new NumberR("0.944652"), NumberR.ZERO));
        SortedMap<Polynom, Element> s3 = new TreeMap<>();
        s3.put(new Polynom("c", r), new Complex(new NumberR("-0.837408"), NumberR.ZERO));
        s3.put(new Polynom("s", r), new Complex(new NumberR("-0.546579"), NumberR.ZERO));
        SortedMap<Polynom, Element> s4 = new TreeMap<>();
        s4.put(new Polynom("c", r), new Complex(new NumberR("0.837408"), NumberR.ZERO));
        s4.put(new Polynom("s", r), new Complex(new NumberR("0.546579"), NumberR.ZERO));
        SortedMap<Polynom, Element> s5 = new TreeMap<>();
        s5.put(new Polynom("c", r), new Complex(new NumberR("-0.882809"), NumberR.ZERO));
        s5.put(new Polynom("s", r), new Complex(new NumberR("0.469733"), NumberR.ZERO));
        SortedMap<Polynom, Element> s6 = new TreeMap<>();
        s6.put(new Polynom("c", r), new Complex(new NumberR("0.882809"), NumberR.ZERO));
        s6.put(new Polynom("s", r), new Complex(new NumberR("-0.469733"), NumberR.ZERO));
        expecting.add(s1);
        expecting.add(s2);
        expecting.add(s3);
        expecting.add(s4);
        expecting.add(s5);
        expecting.add(s6);

        r.MachineEpsilonR = new NumberR("0.001");
        assertTrue(fuzzyEqual(r, expecting, res));
    }

    @Test
    public void example07() {
        /*
         * (x, y)
         * (-1.24962, 1.56155); (1.24962, 1.56155);
         * (-1.60049 i, -2.56155); (1.60049 i, -2.56155)
         */
        Ring r = new Ring("R[x, y]");
        List<SortedMap<Polynom, Element>> res = new SNLESolver(r,
                "x^2 + y^2 - 4",
                "y - x^2").solve();

        List<SortedMap<Polynom, Element>> expecting = new ArrayList<>();
        SortedMap<Polynom, Element> s1 = new TreeMap<>();
        s1.put(new Polynom("x", r), new Complex(new NumberR("-1.24962"), NumberR.ZERO));
        s1.put(new Polynom("y", r), new Complex(new NumberR("1.56155"), NumberR.ZERO));
        SortedMap<Polynom, Element> s2 = new TreeMap<>();
        s2.put(new Polynom("x", r), new Complex(new NumberR("1.24962"), NumberR.ZERO));
        s2.put(new Polynom("y", r), new Complex(new NumberR("1.56155"), NumberR.ZERO));
        SortedMap<Polynom, Element> s3 = new TreeMap<>();
        s3.put(new Polynom("x", r), new Complex(NumberR.ZERO, new NumberR("-1.60049")));
        s3.put(new Polynom("y", r), new Complex(new NumberR("-2.56155"), NumberR.ZERO));
        SortedMap<Polynom, Element> s4 = new TreeMap<>();
        s4.put(new Polynom("x", r), new Complex(NumberR.ZERO, new NumberR("1.60049")));
        s4.put(new Polynom("y", r), new Complex(new NumberR("-2.56155"), NumberR.ZERO));
        expecting.add(s1);
        expecting.add(s2);
        expecting.add(s3);
        expecting.add(s4);

        r.MachineEpsilonR = new NumberR("0.001");
        assertTrue(fuzzyEqual(r, expecting, res));
    }

    @Test
    public void example08() {
        Ring r = new Ring("R64[t,x,e]");
        List<SortedMap<Polynom, Element>> res = new SNLESolver(r,
                "t^2+1-x", "3t-e", "1-e t").solve();
        List<SortedMap<Polynom, Element>> expecting = new ArrayList<>();
        SortedMap<Polynom, Element> s1 = new TreeMap<>();
        s1.put(new Polynom("x", r), new Complex(new Fraction(new NumberZ("4"), new NumberZ("3")), NumberR64.ZERO));
        s1.put(new Polynom("t", r), new Complex(new NumberR64("-0.57735"), NumberR64.ZERO));
        s1.put(new Polynom("e", r), new Complex(new NumberR64("-1.73205"), NumberR64.ZERO));
        SortedMap<Polynom, Element> s2 = new TreeMap<>();
        s2.put(new Polynom("x", r), new Complex(new Fraction(new NumberZ("4"), new NumberZ("3")), NumberR64.ZERO));
        s2.put(new Polynom("t", r), new Complex(new NumberR64("0.57735"), NumberR64.ZERO));
        s2.put(new Polynom("e", r), new Complex(new NumberR64("1.73205"), NumberR64.ZERO));
        expecting.add(s1);
        expecting.add(s2);

        r.FLOATPOS = 6;
        r.MachineEpsilonR = new NumberR("1e-3");
        r.MachineEpsilonR64 = new NumberR64("1e-3");
        assertTrue(fuzzyEqual(r, expecting, res));
    }

    @Ignore("Not ready.")
    @Test
    public void should_return_the_rest_of_the_system_if_there_is_no_onevar_polynomials() {
        // sage006 test from GB tests.
        Ring r = new Ring("R[z, y, x]");
        List<SortedMap<Polynom, Element>> res = new SNLESolver(r,
                "z*x+y^3",
                "z+y^3",
                "z+x*y"
        ).solve();

        // TODO: finish test.
        System.out.println(res);
    }
}
