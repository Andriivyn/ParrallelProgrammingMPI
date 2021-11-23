package com.mathpar.parallel.dap.QR;

import com.mathpar.matrix.MatrixS;
import com.mathpar.number.Element;
import com.mathpar.number.Ring;
import org.junit.Test;

import java.util.Random;
import java.util.logging.Logger;

import static org.junit.Assert.*;

public class SeqQRTest2 {

    private static Logger LOGGER = Logger.getLogger(SeqQRTest2.class.getName());

    private static int size = 1024;
    private static int density = 1;
    private static Random random = new Random(seed());

    private static Ring ring = new Ring("R64[]");


    private static long seed(){
        long s = System.currentTimeMillis();
        System.out.println("Seed: "+s);
        return s;
    }

    @Test
    public void SeqQRTest(){
        MatrixS A = matrix(size);

        long t1 = System.currentTimeMillis();
        MatrixS[] QR = SeqQR.compute(A, ring);
        long t2 = System.currentTimeMillis();

        MatrixS resultCheck = QR[0].multiplySorted(QR[1], ring);
        System.out.println("Sequenced QR: "+ (t2-t1)+"ms");

        MatrixS diff = A.subtract(resultCheck, ring);
        Element maxElement = diff.max(ring);
        ring.setFLOATPOS(100);
        System.out.println("max element "+ maxElement.toString(ring));

        boolean isCorrect = A.equals(resultCheck, ring);
        assertTrue(isCorrect);
    }



    @Test
    public void SeqBlockQRTest(){
        MatrixS A = matrix(size);

        long t1 = System.currentTimeMillis();
        MatrixS[] QR = SeqBlockQR.compute(A, ring);
        long t2 = System.currentTimeMillis();

        MatrixS resultCheck = QR[0].multiply(QR[1], ring);
        System.out.println("Sequenced Block QR: "+ (t2-t1)+"ms");
        MatrixS diff = A.subtract(resultCheck, ring);

        Element maxElement = diff.max(ring);
        ring.setFLOATPOS(100);
        System.out.println("max element "+ maxElement.toString(ring));

        assertTrue(A.equals(resultCheck, ring));
    }

//    @Test
    public void SeqUnoptimizedQRTest(){
        MatrixS A = matrix(size);

        long t1 = System.currentTimeMillis();
        MatrixS[] QR = SeqBlockQR.sequentialQR(A, ring);
        long t2 = System.currentTimeMillis();

        MatrixS resultCheck = QR[0].multiply(QR[1], ring);
        System.out.println("Sequenced  unoptimized QR: "+ (t2-t1)+"ms");

        MatrixS diff = A.subtract(resultCheck, ring);
        Element maxElement = diff.max(ring);
        ring.setFLOATPOS(100);
        System.out.println("max element "+ maxElement.toString(ring));

        assertTrue(A.equals(resultCheck, ring));
    }

    private MatrixS matrix(int size){
        return new MatrixS(size, size, density, new int[]{5}, random, ring.numberONE(), ring);
    }
}
