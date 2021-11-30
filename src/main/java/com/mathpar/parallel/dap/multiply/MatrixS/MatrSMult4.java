package com.mathpar.parallel.dap.multiply.MatrixS;

import com.mathpar.log.MpiLogger;
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

    /**
     * key that starts on 77 (M in ascii) is a mldtsv custom key,
     * the next 2 digits form the number of the step
     */
    @Override
    public void setVars(){
        switch (key){
            // a*b (mldtsv step 8, 13, 20)
            case(0):
            // -a*b (mldtsv step 9)
            case(1): {
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

            case(7702):
            case(7707):
            case(7710):
            case(7718):
            case(7724): {
                inputDataLength = 3;
                outputDataLength = 1;
                resultForOutFunctionLength = 4;
                break;
            }

            case(7703): {
                inputDataLength = 4;
                outputDataLength = 2;
                resultForOutFunctionLength = 4;
                break;
            }

            case(7705): {
                inputDataLength = 6;
                outputDataLength = 1;
                resultForOutFunctionLength = 5;
                break;
            }

            case(7711):
            case(7721):
            case(7725): {
                inputDataLength = 4;
                outputDataLength = 1;
                resultForOutFunctionLength = 4;
                break;
            }

            case(7712): {
                inputDataLength = 7;
                outputDataLength = 3;
                resultForOutFunctionLength = 5;
                break;
            }

            case(7714):
            case(7717):
            case(7722):
            case(7723): {
                inputDataLength = 6;
                outputDataLength = 1;
                resultForOutFunctionLength = 4;
                break;
            }

            case(7716): {
                inputDataLength = 7;
                outputDataLength = 1;
                resultForOutFunctionLength = 5;
                break;
            }

            case(7719): {
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
