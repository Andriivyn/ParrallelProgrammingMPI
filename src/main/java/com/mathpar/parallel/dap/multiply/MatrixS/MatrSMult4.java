package com.mathpar.parallel.dap.multiply.MatrixS;

import com.mathpar.log.MpiLogger;
import com.mathpar.matrix.MatrixD;
import com.mathpar.matrix.MatrixS;
import com.mathpar.number.Array;
import com.mathpar.number.Element;
import com.mathpar.number.Ring;
import com.mathpar.parallel.dap.core.Amin;
import com.mathpar.parallel.dap.core.Drop;
import com.mathpar.parallel.dap.ldumw.LdumwDto;

import java.util.ArrayList;
import java.util.Objects;

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
    public void setVars() {
        switch (key) {
            case (0):
            case (1):
            case (102):
            case (105):
            case (111):
            case (122):
            case (114):
            case (110):
            case (119):
            case (124):
            case (123):
            case (125):
            case (126): {
                inputDataLength = 2;
                outputDataLength = 1;
                resultForOutFunctionLength = 4;
                break;
            }

            case (2): {
                inputDataLength = 2;
                outputDataLength = 2;
                resultForOutFunctionLength = 5;
                break;
            }

            //drop 12
            case (3): {
                inputDataLength = 6;
                outputDataLength = 3;
                resultForOutFunctionLength = 5;
                break;
            }
            //drop 5
            case (4): {
                inputDataLength = 6;
                outputDataLength = 1;
                resultForOutFunctionLength = 5;
                break;
            }
            // LDUMW drop 3 and 4
            case (103):
            case (104): {
                inputDataLength = 3;
                outputDataLength = 2;
                resultForOutFunctionLength = 5;
                break;
            }
            case (107):
            case (115):
            case (120):
            case (116):
            case (121): {
                inputDataLength = 3;
                outputDataLength = 1;
                resultForOutFunctionLength = 4;
                break;
            }
            case (109): {
                inputDataLength = 5;
                outputDataLength = 2;
                resultForOutFunctionLength = 5;
                break;
            }
            case (112): {
                inputDataLength = 5;
                outputDataLength = 4;
                resultForOutFunctionLength = 8;
                break;
            }
            case (113): {
                inputDataLength = 4;
                outputDataLength = 1;
                resultForOutFunctionLength = 4;
                break;
            }
            case (118): {
                inputDataLength = 5;
                outputDataLength = 1;
                resultForOutFunctionLength = 4;
                break;
            }
        }

        inData = new Element[inputDataLength];
        outData = new Element[outputDataLength];

        //LOGGER.info(inData);
    }

    //Послідовний обрахунок листових вершин
    @Override
    public void sequentialCalc(Ring ring) {
        // LOGGER.info("in sequentialCalc indata = " + inData[0] + ",  "+inData[1]);

        MatrixS A = null;
        MatrixS B = null;
        if (inData[0] instanceof MatrixS && inData[1] instanceof MatrixS) {
            A = (MatrixS) inData[0];
            B = (MatrixS) inData[1];
        }

        switch (key) {
            case (0):
                outData[0] = A.multiply(B, ring);
                break;
            case (1):
                outData[0] = A.multiply(B, ring).negate(ring);
                break;
            case (2): {

                MatrixS b = ((MatrixS) inData[0]).transpose();
                outData[1] = b;
                MatrixS bbT = b.multiply((MatrixS) inData[0], ring);
                outData[0] = ((MatrixS) inData[1]).subtract(bbT, ring);

                break;
            }
            case (3): {

                MatrixS M22_2 = A.multiply(B, ring).divideByNumber(inData[4].
                        multiply(inData[4], ring), ring);

                Element M22_3 = inData[5].multiply(M22_2, ring);//temporary(use multiplyLeftI)
                Element ds = inData[2].multiply(inData[3], ring).divide(inData[4], ring);

                outData[0] = M22_2;
                outData[1] = ds;
                outData[2] = M22_3;

                break;
            }
            case (4): {
                MatrixS ET = ((MatrixS) inData[4]).transpose();
                Element d0 = inData[5];
                MatrixS Md = ((MatrixS) inData[2]).multiplyByNumber(inData[3], ring);
                MatrixS M22_1 = A.multiply(ET.multiply(B, ring), ring).negate(ring).
                        divideByNumber(d0, ring).add(Md, ring);
                outData[0] = M22_1;
                break;
            }
            case (102): {
                LdumwDto F11 = ((LdumwDto) inData[0]);
                Element A12 = inData[1];
                Element X_U2 = F11.getJ().multiply(F11.getM(), ring).multiply(A12, ring).divide(F11.getA_n(), ring);

                outData[0] = X_U2;
                break;
            }
            case (103): {
                LdumwDto F11 = ((LdumwDto) inData[0]);
                Element A12 = inData[1];
                Element a = inData[2];

                Element A12_0 = F11.getM().multiply(A12, ring);
                Element A12_2 = F11.getD().multiply(A12_0, ring).divide(a, ring);

                outData[0] = A12_0;
                outData[1] = A12_2;
                break;
            }
            case (104): {
                LdumwDto F11 = ((LdumwDto) inData[0]);
                Element A21 = inData[1];
                Element a = inData[2];

                Element A21_0 = A21.multiply(F11.getW(), ring);
                Element A21_2 = A21_0.multiply(F11.getD(), ring).divide(a, ring);

                outData[0] = A21_0;
                outData[1] = A21_2;
                break;
            }
            case (105): {
                Element A21 = inData[0];
                LdumwDto F11 = ((LdumwDto) inData[1]);

                Element X_L3 = A21.multiply(F11.getW(), ring).multiply(F11.getI(), ring).divide(F11.getA_n(), ring);

                outData[0] = X_L3;
                break;
            }
            case (107): {
                LdumwDto F11 = ((LdumwDto) inData[0]);
                Element A21_0 = inData[1];
                Element A12_0 = inData[2];

                Element A21_1 = F11.getA_n().multiply(A21_0, ring).multiply(F11.getDhat(), ring);
                Element A12_1 = F11.getA_n().multiply(F11.getDhat(), ring).multiply(A12_0, ring);

                MatrixS D11PLUS = F11.getD().transpose();
                Element A22_0 = A21_1.multiply(D11PLUS, ring).multiply(A12_1, ring);

                outData[0] = A22_0;
                break;
            }
            case (109): {
                LdumwDto F11 = ((LdumwDto) inData[0]);
                Element A22 = inData[1];
                Element A22_0 = inData[2];
                LdumwDto F21 = ((LdumwDto) inData[3]);
                Element a = inData[4];

                Element ak2 = F11.getA_n().multiply(F11.getA_n(), ring);
                Element aak2 = a.multiply(ak2, ring);
                Element A22_1left = aak2.multiply(A22, ring).subtract(A22_0, ring);
                Element A22_1 = A22_1left.divide(a.multiply(F11.getA_n(), ring), ring);

                Element X_A22_2 = F21.getDbar().multiply(F21.getM(), ring).multiply(A22_1, ring);

                outData[0] = A22_1;
                outData[1] = X_A22_2;
                break;
            }
            case (110): {
                LdumwDto F11 = ((LdumwDto) inData[0]);
                LdumwDto F21 = ((LdumwDto) inData[1]);

                Element UU = F21.getU().multiply(F11.getU(), ring);

                outData[0] = UU;
                break;
            }
            case (111): {
                LdumwDto F11 = ((LdumwDto) inData[0]);
                LdumwDto F21 = ((LdumwDto) inData[1]);

                Element U1_m1 = F11.getW()
                        .multiply(F11.getDhat(), ring)
                        .multiply(F21.getW(), ring)
                        .multiply(F21.getDhat(), ring);

                outData[0] = U1_m1;
                break;
            }
            case (112): {
                Element X_A22_2 = inData[0];
                LdumwDto F12 = ((LdumwDto) inData[1]);
                LdumwDto F21 = ((LdumwDto) inData[2]);
                LdumwDto F11 = ((LdumwDto) inData[3]);
                Element a = inData[4];

                Element A22_2 = X_A22_2.multiply(F12.getW(), ring).multiply(F12.getDbar(), ring);
                Element lambda = F21.getA_n().divide(F11.getA_n(), ring);
                Element as = lambda.multiply(F12.getA_n(), ring);

                Element I12lambdaM2 = F12.getI().divideByNumbertoFraction(lambda, ring).add(F12.getIbar(), ring);
                Element invD12hat = I12lambdaM2.multiply(F12.getDhat(), ring);

                Element ak2 = F11.getA_n().multiply(F11.getA_n(), ring);
                Element aak2 = ak2.multiply(a, ring);
                Element A22_3 = A22_2.divide(aak2, ring);

                outData[0] = lambda;
                outData[1] = as;
                outData[2] = A22_3;
                outData[3] = invD12hat;
                break;
            }
            case (113): {
                LdumwDto F21 = ((LdumwDto) inData[0]);
                Element A22_1 = inData[1];
                Element X_U2 = inData[2];
                Element a = inData[3];

                Element U2 = F21.getJ().multiply(F21.getM(), ring)
                        .multiply(A22_1, ring)
                        .divide(F21.getA_n().multiply(a, ring), ring)
                        .add(X_U2, ring);

                outData[0] = U2;
                break;
            }
            case (114): {
                LdumwDto F21 = ((LdumwDto) inData[0]);
                Element A22_1 = inData[1];

                Element Y_L3 = F21.getDbar().multiply(F21.getM(), ring).multiply(A22_1, ring);

                outData[0] = Y_L3;
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
            case (2): {
                ms = ((MatrixS) input[0]).transpose();
                ms1 = (MatrixS) input[0];
                amin.resultForOutFunction[4] = ms;
                break;
            }
            case (4): {
                MatrixS E11T = ((MatrixS) input[4]).transpose();
                ms1 = E11T.multiply(ms1, ring);
                break;
            }
        }
        Array.concatTwoArrays(ms.split(), ms1.split(), res);
        return res;

    }

    @Override
    public void independentCalc(Ring ring, Amin amin) {
        switch (key) {
            case (0):
            case (1):
            case (2):
                break;
            case (3): {
                Element ds = inData[2].multiply(inData[3], ring).divide(inData[4], ring);
                amin.resultForOutFunction[4] = ds;
                break;
            }
            case (4): {
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
        switch (key) {
            case (0):
                res = new MatrixS[]{MatrixS.join(resmat)};
                break;
            case (1):
                res = new MatrixS[]{MatrixS.join(resmat).negate(ring)};
                break;
            case (2): {
                res = new Element[]{inData[1].add(MatrixS.join(resmat).negate(ring), ring), input[4]};
                break;
            }
            case (3): {

                MatrixS M22_2 = MatrixS.join(resmat).divideByNumber(inData[4].
                        multiply(inData[4], ring), ring);

                Element M22_3 = inData[5].multiply(M22_2, ring);//temporary(use multiplyLeftI)
                res = new Element[]{M22_2, input[4], M22_3};

                break;
            }
            case (4): {

                res = new MatrixS[]{MatrixS.join(resmat).negate(ring).
                        divideByNumber(inData[5], ring).add((MatrixS) input[4], ring)};
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
