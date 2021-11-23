/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package com.mathpar.matrix;

import com.mathpar.number.Element;
import com.mathpar.number.Ring;
import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;

import static org.junit.Assert.*;

/**
 *
 * @author gennadi
 */
public class BigMatrixSTest {

    public BigMatrixSTest() {
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
     * Test of add method, of class BigMatrixS.
     */
    @Test
    public void testAdd() {
        Ring ring = new Ring("Z[x]");
        int[][] mat1 = {{1, 3}, {2, 4}};
        MatrixD[][] MD2x2 = new MatrixD[2][2];
        for (int i = 0; i < 2; i++) {
            for (int j = 0; j < 2; j++) {
                MD2x2[i][j] = new MatrixD(mat1, ring);
            }
        }
        BigMatrixS BMS = new BigMatrixS(MD2x2, ring);
        BigMatrixS sum = BMS.add(BMS, ring);
        System.out.println("BMS+BMS= " + sum);
    }

    @Test
    public void testSubtr() {
        Ring ring = new Ring("Z[x]");
        int[][] mat1 = {{2, 3}, {5, 4}};
        MatrixD[][] MD2x2 = new MatrixD[2][2];
        for (int i = 0; i < 2; i++) {
            for (int j = 0; j < 2; j++) {
                MD2x2[i][j] = new MatrixD(mat1, ring);
            }
        }
        BigMatrixS BMS = new BigMatrixS(MD2x2, ring);
        BigMatrixS sub = BMS.subtract(BMS, ring);
        System.out.println("BMS-BMS= " + sub);
    }

    @Test
    public void testMultByNumber() {
        Ring ring = new Ring("Z[x]");
        Element e = ring.numberMINUS_ONE;
        System.out.println("e= " + e);
        int[][] mat1 = {{2, 3}, {5, 4}};
        int[][] mat1Neg = {{-2, -3}, {-5, -4}};
        MatrixD[][] MD2x2 = new MatrixD[2][2];
        MatrixD[][] MD2x2Neg = new MatrixD[2][2];
        for (int i = 0; i < 2; i++) {
            for (int j = 0; j < 2; j++) {
                MD2x2[i][j] = new MatrixD(mat1, ring);
                MD2x2Neg[i][j] = new MatrixD(mat1Neg, ring);
            }
        }
        BigMatrixS BMS = new BigMatrixS(MD2x2, ring);
        BigMatrixS multN = BMS.multiplyByNumber(e, ring);
        BigMatrixS BMSNeg = new BigMatrixS(MD2x2Neg, ring);
        System.out.println("BMS*e= " + multN + "\n BMSNeg=" + BMSNeg + "\ndifference=" + BMSNeg.subtract(multN, ring));
        assertTrue("Somthing wrong in BigMatrixS multiplyByNumber",
                (BMSNeg.subtract(multN, ring)).isZero(ring));
    }
 
    /**
     * Test of add method, of class BigMatrixS.
     */
    /**
     * Test of subtract method, of class BigMatrixS.
     */


    /**
     * Test of toMatrixD method, of class BigMatrixS.
     */
    @Test
    public void testToMatrixD() {
        Ring ring = new Ring("Z[x]");
        int[][] mat1 = {{1, 3}, {2, 4}};
        MatrixD[][] MD2x2 = new MatrixD[2][2];
        for (int i = 0; i < 2; i++) {
            for (int j = 0; j < 2; j++) {
                MD2x2[i][j] = new MatrixD(mat1, ring);}}
        BigMatrixS instance = new BigMatrixS(MD2x2, ring);
        MatrixD expResult = new MatrixD(new int[][]{{1, 3, 1, 3}, {2, 4, 2, 4},{1, 3, 1, 3}, {2, 4, 2, 4}}, ring);
        MatrixD result = instance.toMatrixD(ring);
        assertTrue(expResult.subtract(result, ring).isZero(ring));
    }

    /**
     * Test of toString method, of class BigMatrixS.
     */
    @Test
    public void testToString_Ring() {
        System.out.println("toString");
        Ring ring = new Ring("Z[x]");
        int[][] mat1 = {{1, 3}, {2, 4}};
        MatrixD[][] MD2x2 = new MatrixD[2][2];
        for (int i = 0; i < 2; i++) {
            for (int j = 0; j < 2; j++) {
                MD2x2[i][j] = new MatrixD(mat1, ring); }}
        BigMatrixS bm = new BigMatrixS(MD2x2, ring);
        String expResult = "[[1, 3, 1, 3],\n[2, 4, 2, 4],\n[1, 3, 1, 3],\n[2, 4, 2, 4]]";
        String result = bm.toString(ring);
        assertEquals(expResult, result);
    }
  
}
