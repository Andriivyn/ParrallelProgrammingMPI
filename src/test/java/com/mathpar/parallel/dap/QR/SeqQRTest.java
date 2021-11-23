package com.mathpar.parallel.dap.QR;

import com.mathpar.matrix.MatrixS;
import com.mathpar.number.Ring;
import org.junit.Test;

import java.util.Random;
import java.util.logging.Logger;

import static org.junit.Assert.*;

public class SeqQRTest {
    private static Logger LOGGER = Logger.getLogger(SeqQRTest.class.getName());

    private static int size = 256;
    private static Random random = new Random();
    private static Ring ring = new Ring("R64[]");

    private static int density = 10;

    @Test
    public void SeqBlockQRTest(){
        MatrixS A = matrix(size);

        long t1 = System.currentTimeMillis();
        MatrixS[] QR = SeqBlockQR.compute(A, ring);
        long t2 = System.currentTimeMillis();

        MatrixS resultCheck = QR[0].multiply(QR[1], ring);
//        System.out.println("Sequenced QR: "+ (t2-t1)+"ms");
//        System.out.println("Q = "+QR[0]);
//        System.out.println("R = "+QR[1]);
        MatrixS diff = A.subtract(resultCheck, ring);
        System.out.println("Diff = "+diff);

//        assertTrue("Orthogonal", isOrthogonal(QR[0], ring));
        boolean check = A.equals(resultCheck, ring);
        assertTrue(check);
    }

    @Test
    public void SeqQRTest(){
        MatrixS A = matrix(size);

        long t1 = System.currentTimeMillis();
        MatrixS[] QR = SeqBlockQR.sequentialQR(A, ring);
        long t2 = System.currentTimeMillis();

        MatrixS resultCheck = QR[0].multiply(QR[1], ring);
        System.out.println("Sequenced QR: "+ (t2-t1)+"ms");
//        System.out.println("Q = "+QR[0]);
//        System.out.println("R = "+QR[1]);
        MatrixS diff = A.subtract(resultCheck, ring);
        System.out.println("Diff = "+diff);

//        assertTrue("Orthogonal",isOrthogonal(QR[0], ring));
        assertTrue(A.equals(resultCheck, ring));
    }



//    @Test
//    public void SeqBlockQPTest(){
//
//        MatrixS A = matrix(size);
//        MatrixS B = matrix(size);
//
//        makeR(B, ring);
//        MatrixS expected = MatrixS.embedBlocksOfColumn(A, B, ring);
//
//        MatrixS[] QP = SeqBlockQR.computeQP(A, B, ring);
//        MatrixS resultCheck = QP[0].multiply(expected, ring);
//
//        MatrixS diff = resultCheck.subtract(expected, ring);
//        System.out.println("Result QP diff" + diff);
//
//        assertTrue("Block Q is not orthogonal", isOrthogonal(QP[0], ring));
//        assertTrue("Wrong result of computation of Block QP", QP[1].equals(resultCheck, ring));
//    }

//    @Test
//    public void SeqQPTest(){
//
//        MatrixS A = matrix(size);
//        MatrixS B = matrix(size);
//
//        makeR(B, ring);
//        MatrixS expected = MatrixS.embedBlocksOfColumn(A, B, ring);
//
//        MatrixS[] seqQP = SeqBlockQR.sequentialQP(A, B, ring);
////        MatrixS P2 = MatrixS.embedBlocksOfColumn(seqQP[1], MatrixS.zeroMatrix(size), ring);
//        MatrixS resultCheck2 = seqQP[0].multiply(expected, ring);
//
//        assertEquals(expected.size, resultCheck2.size);
//        assertTrue("Seq Q is not orthogonal", isOrthogonal(seqQP[0], ring));
//        assertTrue("Wrong result of computation of Seq QP",seqQP[1].equals(resultCheck2, ring));
//    }

    @Test
    public void leafQRTest(){

        MatrixS A = matrix(2);
        MatrixS[] QR = SeqBlockQR.computeLeaf(A.copy(), ring);
        MatrixS resultCheck = QR[0].multiply(QR[1], ring);

        assertTrue("Orthogonal",isOrthogonal(QR[0], ring));
        boolean check = A.equals(resultCheck, ring);
        assertTrue(check);

    }

    @Test
    public void leafQRTest2(){

        MatrixS A = matrix(2);
        MatrixS[] QR = SeqBlockQR.computeLeaf(A, ring);
        MatrixS resultCheck = QR[0].transpose().multiply(A, ring);

        assertTrue("Orthogonal",isOrthogonal(QR[0], ring));
        boolean check = QR[1].equals(resultCheck, ring);
        assertTrue(check);
    }

    @Test
    public void leafQPTest(){
        MatrixS A = matrix(1);
        MatrixS B = matrix(1);

        B.putElement(ring.numberZERO(), 1, 0);

        MatrixS[] QP = SeqBlockQR.computeQPLeaf(A, B, ring);

        MatrixS expected = MatrixS.embedBlocksOfColumn(A, B, ring);
        MatrixS resultCheck = QP[0].multiply(expected, ring);

        assertTrue("Orthogonal", isOrthogonal( QP[0], ring));
        assertTrue(QP[1].equals(resultCheck, ring));

    }

    private MatrixS matrix(int size){
        return new MatrixS(size, size, density, new int[]{5}, random, ring.numberONE(), ring);
    }


    private boolean isOrthogonal(MatrixS matrix, Ring ring){
        MatrixS transposed = matrix.transpose();
//        MatrixS inverse = matrix.inverse(ring);
//        if(inverse == null) return false;
//        return transposed.equals(inverse, ring);

        MatrixS result = transposed.multiply(matrix, ring);
        boolean isOne  = result.isOne(ring);
        if(!isOne){
            System.out.println("M x MT"+result);
        }

        return isOne;
    }
}
