/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package com.mathpar.matrix.LDU;

import com.mathpar.matrix.MatrixS;
import com.mathpar.number.NFunctionZ32;
import com.mathpar.number.Ring;
import java.io.IOException;
import java.util.logging.Level;
import java.util.logging.Logger;
import org.junit.After;
import org.junit.AfterClass;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;

import static org.junit.Assert.*;

/**
 *
 * @author ridkeim
 */
public class ETDTest {
    
    public ETDTest() {
    }
    
    @BeforeClass
    public static void setUpClass() {
        try {
            System.out.println("generating primes in setUpClass");
            NFunctionZ32.doFileOfIntPrimes();
        } catch (IOException ex) {
            Logger.getLogger(ETDTest.class.getName()).log(Level.SEVERE, null, ex);
        }
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
     * Test of ETDmodLDU method, of class ETD.
     */
    @Test
    public void testETDmodLDU() {
        System.out.println("ETDmodLDU");
        Ring ring = Ring.ringZxyz;
        MatrixS T = new MatrixS(
                new int[][]{
                    {422, 753,  890,  0,   0,   0,    0,    -721},
                    {0,   0,    0,    699, 0,   -652, -245, 829 },
                    {0,   -388, 0,    0,   0,   0,    0,    0   },
                    {873, 0,    0,    0,   0,   -937, 230,  -825},
                    {0,   0,    -212, 0,   0,   -562, 0,    0   },
                    {95,  -303, 0,    0,   628, 0,    0,    -310},
                    {0,   98,   0,    0,   0,   0,    75,   -376},
                    {0,   235,  -26,  797, 660, -803, 0,    -8  },
                },ring);
        MatrixS[] result = ETD.ETDmodLDU(T);
        MatrixS L,D,U;
        L = result[0];
        D = result[1];
        U = result[2];
        MatrixS expRes = L.multiply(D, ring).multiply(U, ring);
        MatrixS zeroMatrix = T.subtract(expRes, ring);
        assertTrue(zeroMatrix.isZero(ring));        
    }

    /**
     * Test of ETDmodWDK method, of class ETD.
     */
    @Test
    public void testETDmodWDK() {
        System.out.println("ETDmodWDK");
        Ring ring = Ring.ringZxyz;
        MatrixS T = new MatrixS(
                new int[][]{
                    {422, 753,  890,  0,   0,   0,    0,    -721},
                    {0,   0,    0,    699, 0,   -652, -245, 829 },
                    {0,   -388, 0,    0,   0,   0,    0,    0   },
                    {873, 0,    0,    0,   0,   -937, 230,  -825},
                    {0,   0,    -212, 0,   0,   -562, 0,    0   },
                    {95,  -303, 0,    0,   628, 0,    0,    -310},
                    {0,   98,   0,    0,   0,   0,    75,   -376},
                    {0,   235,  -26,  797, 660, -803, 0,    -8  },
                },ring);
        MatrixS[] result = ETD.ETDmodWDK(T);
        MatrixS W,D,K;
        W = result[0];
        D = result[1];
        K = result[2];
        MatrixS expRes = W.multiply(D, ring).multiply(K, ring);
        MatrixS zeroMatrix = T.inverse(ring).subtract(expRes, ring);
        assertTrue(zeroMatrix.isZero(ring));        
    }

    /**
     * Test of ETDmodPLDUQWDK method, of class ETD.
     */
    @Test
    public void testETDmodPLDUQWDK() {
        System.out.println("ETDmodPLDUQWDK");
        Ring ring = Ring.ringZxyz;
        MatrixS T = new MatrixS(
                new int[][]{
                    {422, 753,  890,  0,   0,   0,    0,    -721},
                    {0,   0,    0,    699, 0,   -652, -245, 829 },
                    {0,   -388, 0,    0,   0,   0,    0,    0   },
                    {873, 0,    0,    0,   0,   -937, 230,  -825},
                    {0,   0,    -212, 0,   0,   -562, 0,    0   },
                    {95,  -303, 0,    0,   628, 0,    0,    -310},
                    {0,   98,   0,    0,   0,   0,    75,   -376},
                    {0,   235,  -26,  797, 660, -803, 0,    -8  },
                },ring);
        MatrixS[] result = ETD.ETDmodPLDUQWDK(T);
        MatrixS P,L,D,U,Q,W,K,res;
        P = result[0];
        L = result[1];
        D = result[2];
        U = result[3];
        Q = result[4];
        W = result[5];
        K = result[6];
        MatrixS expRes = P.multiply(L, ring).multiply(D, ring).multiply(U, ring).multiply(Q, ring);
        res = T.subtract(expRes, ring);
        assertTrue(res.isZero(ring));
        MatrixS expInv = Q.transpose().multiply(W, ring).multiply(D, ring).multiply(K, ring).multiply(P.transpose(), ring);
        res = T.inverse(ring).subtract(expInv, ring);
        assertTrue(res.isZero(ring));
    }

    /**
     * Test of ETDmodLDUWDK method, of class ETD.
     */
    @Test
    public void testETDmodLDUWDK() {
        System.out.println("ETDmodLDUWDK");
        Ring ring = Ring.ringZxyz;
        MatrixS T = new MatrixS(
                new int[][]{
                    {422, 753,  890,  0,   0,   0,    0,    -721},
                    {0,   0,    0,    699, 0,   -652, -245, 829 },
                    {0,   -388, 0,    0,   0,   0,    0,    0   },
                    {873, 0,    0,    0,  0,   -937, 230,  -825},
                    {0,   0,    -212, 0,   0,   -562, 0,    0   },
                    {95,  -303, 0,    0,   628, 0,    0,    -310},
                    {0,   98,   0,    0,   0,   0,    75,   -376},
                    {0,   235,  -26,  797, 660, -803, 0,    -8  },
                },ring);
        MatrixS[] result = ETD.ETDmodLDUWDK(T);
        MatrixS L,D,U,W,Dd,K,res;
        L = result[0];
        D = result[1];
        U = result[2];
        W = result[3];
        Dd = result[4];
        K = result[5];
        MatrixS expRes = L.multiply(D, ring).multiply(U, ring);
        res = T.subtract(expRes, ring);
        assertTrue(res.isZero(ring));
        MatrixS expInv = W.multiply(Dd, ring).multiply(K, ring);
        res = T.inverse(ring).subtract(expInv, ring);
        assertTrue(res.isZero(ring));
    }

     /**
     * Test of ETDLDU method, of class ETD.
     */
    @Test
    public void testETDLDU() {
        System.out.println("ETDLDU");
        Ring ring = Ring.ringZxyz;
        MatrixS T = new MatrixS(new int[][]{{-1,0,0,0},{-9,0,0,-20},{-19,0,1,-5},{-5,-13,0,-19}},ring);
        MatrixS[] result = ETD.ETDLDU(T);
        MatrixS L,D,U;
        L = result[0];
        D = result[1];
        U = result[2];
        MatrixS expRes = L.multiply(D, ring).multiply(U, ring);
        MatrixS zeroMatrix = T.subtract(expRes, ring);
        assertTrue(zeroMatrix.isZero(ring));        
    }

    /**
     * Test of ETDWDK method, of class ETD.
     */
    @Test
    public void testETDWDK() {
        System.out.println("ETDWDK");
        Ring ring = Ring.ringZxyz;
        MatrixS T = new MatrixS(new int[][]{{-1,0,0,0},{-9,0,0,-20},{-19,0,1,-5},{-5,-13,0,-19}},ring);
        MatrixS[] result = ETD.ETDWDK(T);
        MatrixS W,D,K;
        W = result[0];
        D = result[1];
        K = result[2];
        MatrixS expRes = W.multiply(D, ring).multiply(K, ring);
        MatrixS zeroMatrix = T.inverse(ring).subtract(expRes, ring);
        assertTrue(zeroMatrix.isZero(ring));
    }

    /**
     * Test of ETDPLDUQWDK method, of class ETD.
     */
    @Test
    public void testETDPLDUQWDK() {
        System.out.println("ETDPLDUQWDK");
        Ring ring = Ring.ringZxyz;
        MatrixS T = new MatrixS(new int[][]{{-1,0,0,0},{-9,0,0,-20},{-19,0,1,-5},{-5,-13,0,-19}},ring);
        MatrixS[] result = ETD.ETDPLDUQWDK(T);
        MatrixS P,L,D,U,Q,W,K,res;
        P = result[0];
        L = result[1];
        D = result[2];
        U = result[3];
        Q = result[4];
        W = result[5];
        K = result[6];
        MatrixS expRes = P.multiply(L, ring).multiply(D, ring).multiply(U, ring).multiply(Q, ring);
        res = T.subtract(expRes, ring);
        assertTrue(res.isZero(ring));
        MatrixS expInv = Q.transpose().multiply(W, ring).multiply(D, ring).multiply(K, ring).multiply(P.transpose(), ring);
        res = T.inverse(ring).subtract(expInv, ring);
        assertTrue(res.isZero(ring));
    }

    /**
     * Test of ETDLDUWDK method, of class ETD.
     */
    @Test
    public void testETDLDUWDK() {
        System.out.println("ETDLDUWDK");
        Ring ring = Ring.ringZxyz;
        MatrixS T = new MatrixS(new int[][]{{-1,0,0,0},{-9,0,0,-20},{-19,0,1,-5},{-5,-13,0,-19}},ring);
        MatrixS[] result = ETD.ETDLDUWDK(T);
        MatrixS L,D,U,W,Dd,K,res;
        L = result[0];
        D = result[1];
        U = result[2];
        W = result[3];
        Dd = result[4];
        K = result[5];
        MatrixS expRes = L.multiply(D, ring).multiply(U, ring);
        res = T.subtract(expRes, ring);
        assertTrue(res.isZero(ring));
        MatrixS expInv = W.multiply(Dd, ring).multiply(K, ring);
        res = T.inverse(ring).subtract(expInv, ring);
        assertTrue(res.isZero(ring));
    }

    /**
     * Test of ETDmodLDUpar method, of class ETD.
     */
//    @Test
//    public void testETDmodLDUpar() {
//        System.out.println("ETDmodLDUpar");
//        MatrixS T = null;
//        MatrixS[] expResult = null;
//        MatrixS[] result = ETD.ETDmodLDUpar(T);
//        assertArrayEquals(expResult, result);
//        // TODO review the generated test code and remove the default call to fail.
//        fail("The test case is a prototype.");
//    }

    /**
     * Test of ETDmodWDKpar method, of class ETD.
     */
//    @Test
//    public void testETDmodWDKpar() {
//        System.out.println("ETDmodWDKpar");
//        MatrixS T = null;
//        MatrixS[] expResult = null;
//        MatrixS[] result = ETD.ETDmodWDKpar(T);
//        assertArrayEquals(expResult, result);
//        // TODO review the generated test code and remove the default call to fail.
//        fail("The test case is a prototype.");
//    }

    /**
     * Test of ETDmodPLDUQWDKpar method, of class ETD.
     */
//    @Test
//    public void testETDmodPLDUQWDKpar() {
//        System.out.println("ETDmodPLDUQWDKpar");
//        MatrixS T = null;
//        MatrixS[] expResult = null;
//        MatrixS[] result = ETD.ETDmodPLDUQWDKpar(T);
//        assertArrayEquals(expResult, result);
//        // TODO review the generated test code and remove the default call to fail.
//        fail("The test case is a prototype.");
//    }

    /**
     * Test of ETDmodLDUWDKpar method, of class ETD.
     */
//    @Test
//    public void testETDmodLDUWDKpar() {
//        System.out.println("ETDmodLDUWDKpar");
//        MatrixS T = null;
//        MatrixS[] expResult = null;
//        MatrixS[] result = ETD.ETDmodLDUWDKpar(T);
//        assertArrayEquals(expResult, result);
//        // TODO review the generated test code and remove the default call to fail.
//        fail("The test case is a prototype.");
//    }
    
   
    
}
