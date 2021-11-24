package com.mathpar.parallel.dap.multiply.multiplyVar;

import com.mathpar.matrix.MatrixS;
import com.mathpar.number.Element;
import com.mathpar.number.Ring;
import com.mathpar.parallel.dap.core.Amin;
import org.junit.Test;

import java.util.Random;

public class MultiplyVarTest {

    private Ring ring = new Ring("Z[]");
    private int size = 64;
    private int density = 3000;
    private Random random = new Random();


    @Test
    public void test1(){
        MultiplyVar drop = new MultiplyVar(
                MultiplyVarConfig.Builder.startWith()
                .multiplicationBlock()
                .negateMultiplication()
                .setFirst(MultiplyVarMatrix.Builder.newBuilder().addQuarter(MultiplyVarMatrix.Quarter.IV).set().done().build())
                .setSecond(MultiplyVarMatrix.Builder.newBuilder().addCenter().transpose().set().done().build())
                .build()
        );


        MatrixS m1 = matrix( size);
        MatrixS m2 = matrix(size/2);

        MatrixS[] input = new MatrixS[12];

        input[3] = m1;
        input[4] = m2;

       /* assertFalse(drop.hasAdditionalInputData());
        assertFalse(drop.hasMainInputData());
        assertFalse(drop.hasFullInputData());
        assertFalse(drop.canHaveAdditionalData());*/

        drop._setInData(input);

        //assertFalse(drop.canHaveAdditionalData());


        Element[] result = drop.inputFunction(input, new Amin(drop, 1, 0), ring);

        drop._setInData(result);

       /* assertFalse(drop.hasAdditionalInputData());
        assertTrue(drop.hasMainInputData());
        assertTrue(drop.hasFullInputData());*/


    }



    private MatrixS matrix(int size){
        return new MatrixS(size, size, density, new int[]{5}, random, ring.numberONE(), ring);
    }
}
