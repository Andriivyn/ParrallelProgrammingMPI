package com.mathpar.parallel.dap.multiply.MatrixS;

import com.mathpar.log.MpiLogger;
import com.mathpar.matrix.MatrixD;
import com.mathpar.matrix.MatrixS;
import com.mathpar.number.Array;
import com.mathpar.number.Element;
import com.mathpar.number.Ring;
import com.mathpar.parallel.dap.core.Amin;
import com.mathpar.parallel.dap.core.Drop;

import java.util.ArrayList;

public class MatrSMult4 extends Drop {
    private final static MpiLogger LOGGER = MpiLogger.getLogger(MatrSMult4.class);
    protected static int leafSize = 4;


    private static int[][] _arcs = new int[][]{
            {1, 0, 0, 1, 4, 1, 1, 1, 2, 1, 6, 3, 2, 0, 0, 2, 5, 1, 2, 1, 2, 2, 7, 3,
                    3, 2, 0, 3, 4, 1, 3, 3, 2, 3, 6, 3, 4, 2, 0, 4, 5, 1, 4, 3, 2, 4, 7, 3},
            {5, 0, 0},
            {5, 0, 1},
            {5, 0, 2},
            {5, 0, 3},
            {}};

    public MatrSMult4() {

        //Дроп має тип 5
        type = 5;

        //унікальний номер дропа
        number = cnum++;
        arcs = _arcs;
    }

    //Розгортання аміну з дропами, відповідно до графу, для обрахунку поточного дропа.
    @Override
    public ArrayList<Drop> doAmin() {
        ArrayList<Drop> amin = new ArrayList<Drop>();

        amin.add(new MatrSMultiplyScalar());
        amin.add(new MatrSMultiplyScalar());
        amin.add(new MatrSMultiplyScalar());
        amin.add(new MatrSMultiplyScalar());

        return amin;
    }
    @Override
    public void setVars(){
        switch (key){
            case(0):
            case(1):
            case(102):
            case(105):
            case(111):
            case(122):
            case(114):
            case(110):
            case(119):
            case(124):
            case(123):
            case(125):
            case(126): {
                inputDataLength = 2;
                outputDataLength = 1;
                resultForOutFunctionLength = 4;
                break;
            }

            case(2): {
                inputDataLength = 2;
                outputDataLength = 2;
                resultForOutFunctionLength = 5;
                break;
            }

            //drop 12
            case(3):{
                inputDataLength = 6;
                outputDataLength = 3;
                resultForOutFunctionLength = 5;
                break;
            }
            //drop 5
            case(4):{
                inputDataLength = 6;
                outputDataLength = 1;
                resultForOutFunctionLength = 5;
                break;
            }
            // LDUMW drop 3 and 4
            case(103):
            case(104):{
                inputDataLength = 3;
                outputDataLength = 2;
                resultForOutFunctionLength = 5;
                break;
            }
            case(107):
            case(115):
            case(120):
            case(116):
            case(121):{
                inputDataLength = 3;
                outputDataLength = 1;
                resultForOutFunctionLength = 4;
                break;
            }
            case(109):{
                inputDataLength = 5;
                outputDataLength = 2;
                resultForOutFunctionLength = 5;
                break;
            }
            case(112): {
                inputDataLength = 5;
                outputDataLength = 4;
                resultForOutFunctionLength = 8;
                break;
            }
            case(113): {
                inputDataLength = 4;
                outputDataLength = 1;
                resultForOutFunctionLength = 4;
                break;
            }
            case(118): {
                inputDataLength = 5;
                outputDataLength = 1;
                resultForOutFunctionLength = 4;
                break;
            }
        }

        inData =  new Element[inputDataLength];
        outData =  new Element[outputDataLength];

        //LOGGER.info(inData);
    }

    //Послідовний обрахунок листових вершин
    @Override
    public void sequentialCalc(Ring ring) {
        // LOGGER.info("in sequentialCalc indata = " + inData[0] + ",  "+inData[1]);
        MatrixS A = (MatrixS) inData[0];
        MatrixS B = (MatrixS) inData[1];

        switch (key){
            case(0): outData[0] =A.multiply(B, ring); break;
            case(1): outData[0] = A.multiply(B, ring).negate(ring); break;
            case(2): {

                MatrixS b = ((MatrixS) inData[0]).transpose();
                outData[1] = b;
                MatrixS bbT =  b.multiply((MatrixS) inData[0], ring);
                outData[0] = ((MatrixS)inData[1]).subtract(bbT, ring);

                break;
            }
            case(3):{

                MatrixS M22_2 = A.multiply(B, ring).divideByNumber(inData[4].
                        multiply(inData[4], ring), ring);

                Element M22_3 = inData[5].multiply(M22_2, ring);//temporary(use multiplyLeftI)
                Element ds = inData[2].multiply(inData[3], ring).divide(inData[4], ring);

                outData[0]  = M22_2;
                outData[1]  =  ds;
                outData[2] = M22_3;

                break;
            }
            case(4):{
                MatrixS ET = ((MatrixS) inData[4]).transpose();
                Element d0 = inData[5];
                MatrixS Md = ((MatrixS) inData[2]).multiplyByNumber(inData[3],ring);
                MatrixS M22_1 = A.multiply(ET.multiply(B, ring),ring).negate(ring).
                        divideByNumber(d0, ring).add(Md,ring);
                outData[0] = M22_1;
                break;
            }
        }
    }

    @Override
    //Вхідна функція дропа, розбиває вхідні дані на блоки.
    public MatrixS[] inputFunction(Element[] input, Amin amin, Ring ring) {
//LOGGER.info(input[0]);
        MatrixS[] res = new MatrixS[8];
        MatrixS ms = (MatrixS) input[0];
        MatrixS ms1 = (MatrixS) input[1];


        switch (key) {
            case(2):{
                ms =  ((MatrixS) input[0]).transpose();
                ms1 = (MatrixS) input[0];
                amin.resultForOutFunction[4] = ms;
                break;
            }
            case(4):{
                MatrixS E11T = ((MatrixS) input[4]).transpose();
                ms1 = E11T.multiply(ms1, ring);
                break;
            }
        }
        Array.concatTwoArrays(ms.split(), ms1.split(), res);
        return res;

    }

    @Override
    public void independentCalc(Ring ring, Amin amin)
    {
        switch (key){
            case(0):
            case(1):
            case(2): break;
            case(3):{
                Element ds = inData[2].multiply(inData[3], ring).divide(inData[4], ring);
                amin.resultForOutFunction[4] = ds;
                break;
            }
            case(4):{
                Element md = inData[2].multiply(inData[3], ring);
                amin.resultForOutFunction[4] = md;
                break;
            }
        }
        return;
    }

    //Вихідна функція дропа, яка збирає блоки в результат
    @Override
    public Element[] outputFunction(Element[] input, com.mathpar.number.Ring ring) {
        MatrixS[] resmat = new MatrixS[4];
        for (int i = 0; i < 4; i++) {
            resmat[i] = (MatrixS) input[i];
        }

        Element[] res = new Element[outputDataLength];
        switch (key){
            case(0):
                res = new MatrixS[]{MatrixS.join(resmat)};
                break;
            case(1):
                res = new MatrixS[]{MatrixS.join(resmat).negate(ring)};
                break;
            case(2):{
                res = new Element[]{inData[1].add(MatrixS.join(resmat).negate(ring),ring), input[4]};
                break;
            }
            case(3):{

                MatrixS M22_2 = MatrixS.join(resmat).divideByNumber(inData[4].
                        multiply(inData[4], ring), ring);

                Element M22_3 = inData[5].multiply(M22_2, ring);//temporary(use multiplyLeftI)
                res = new Element[]{M22_2, input[4], M22_3};

                break;
            }
            case(4):{

                res = new MatrixS[]{MatrixS.join(resmat).negate(ring).
                        divideByNumber(inData[5], ring).add((MatrixS)input[4],ring)};
                break;
            }
        }

        return res;
    }

    //Перевіряє чи є дроп листовим
    @Override
    public boolean isItLeaf() {
        MatrixS ms = (MatrixS) inData[0];
        return (ms.size <= leafSize);
    }

    @Override
    public void setLeafSize(int dataSize) {
        leafSize = dataSize;
    }
}
