package com.mathpar.polynom.gbasis.datastructures;

import java.util.Arrays;
import java.util.List;
import com.mathpar.number.Ring;
import static org.junit.Assert.*;
import org.junit.Before;
import org.junit.Test;
import com.mathpar.polynom.Term;

/**
 * @author Ivan Borisov
 */
public class TermsListTest {

    private static final Ring rZxy = new Ring("Z[x, y]"); // y > x
    private static final List<Term> termsForAdd = Arrays.asList(
            new Term(5, 2), // y^2x^5
            new Term(3, 4), // y^4x^3
            new Term(2, 6) //  y^6x^2
            );
    private TermsList termsList;

    public TermsListTest() {
    }

    @Before
    public void setUp() {
        termsList = new TermsList(rZxy, termsForAdd, true);
    }

    private static void assertAddResult(int[][] expected, int[][] actual) {
        assertEquals(2, expected.length);
        assertEquals(2, actual.length);
        assertArrayEquals(expected[0], actual[0]);
        assertArrayEquals(expected[1], actual[1]);
    }

    private static int[] arr(int... arr) {
        return arr;
    }

    private static int[][] addRes(int[] arr1, int[] arr2) {
        return new int[][]{arr1, arr2};
    }

    @Test
    public void testAddLowestTermGoesToEnd() {
        Term termToAdd = new Term(10, 1); // yx^10 goes to pos 3

        int[][] res = termsList.add(termToAdd);

        assertEquals(termsList.size(), 4);
        assertEquals(termToAdd, termsList.get(3));
        assertAddResult(addRes(arr(0, 1, 2), arr(3)), res);
    }

    @Test
    public void testAddMiddleTermGoesToCorrectPosition() {
        Term termToAdd = new Term(0, 4); // y^4 goes to pos 2

        int[][] res = termsList.add(termToAdd);

        assertEquals(termsList.size(), 4);
        assertEquals(termToAdd, termsList.get(2));
        assertAddResult(addRes(arr(0, 1, 3), arr(2)), res);
    }

    @Test
    public void testAddHighestTermGoesToStart() {
        Term termToAdd = new Term(20, 10); // y^10x^20 goes to pos 0

        int[][] res = termsList.add(termToAdd);

        assertEquals(termsList.size(), 4);
        assertEquals(termToAdd, termsList.get(0));
        assertAddResult(addRes(arr(1, 2, 3), arr(0)), res);
    }

    @Test
    public void testAddExistingTermShouldNotBeAdded() {
        Term termToAdd = new Term(3, 4); // y^4x^3 DON'T add it

        int[][] res = termsList.add(termToAdd);

        assertEquals(termsList.size(), 3);
        assertAddResult(addRes(arr(0, 1, 2), arr(-1)), res);
    }

    @Test
    public void testAddAllLowMidHighTermsGoToCorrectPositionsAndCorrectIndexesAreReturned() {
        List<Term> termsToAdd = Arrays.asList(
                new Term(10, 1), // yx^10       goes to pos 5
                new Term(0, 4), //  y^4         goes to pos 3
                new Term(20, 10) // y^10x^20    goes to pos 0
                );
        // Result should be [(y^10x^20), y^6x^2, y^4x^3, (y^4), y^2x^5, (yx^10)]
        int[][] res = termsList.addAll(termsToAdd, true);

        assertEquals(termsList.size(), 6);
        assertAddResult(addRes(arr(1, 2, 4), arr(0, 3, 5)), res);
    }

    @Test
    public void testAddAllExistingTermsShouldNotBeAdded() {
        List<Term> termsToAdd = Arrays.asList(
                new Term(2, 6), // y^6x^2
                new Term(3, 4) // y^4x^3
                );

        int[][] res = termsList.addAll(termsToAdd, true);

        assertEquals(termsList.size(), 3);
        // Identity map + no terms were added.
        assertAddResult(addRes(arr(0, 1, 2), arr(-1)), res);
    }
}
