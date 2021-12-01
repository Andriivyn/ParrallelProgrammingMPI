package com.mathpar.parallel.dap.multiply.MatrixS;

import com.mathpar.log.MpiLogger;
import com.mathpar.matrix.AdjMatrixS;
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
            // a*b
            case(0):
            // -a*b
            case(1):
            case(7708):
            case(7709):
            case(7713):
            case(7720): {
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


        switch (key){
            case(0):
            case(7713): {
                MatrixS A = (MatrixS) inData[0];
                MatrixS B = (MatrixS) inData[1];
                outData[0] =A.multiply(B, ring);
                break;
            }
            case(1): {
                MatrixS A = (MatrixS) inData[0];
                MatrixS B = (MatrixS) inData[1];
                outData[0] = A.multiply(B, ring).negate(ring);
                break;
            }
            case(2): {
                MatrixS b = ((MatrixS) inData[0]).transpose();
                outData[1] = b;
                MatrixS bbT =  b.multiply((MatrixS) inData[0], ring);
                outData[0] = ((MatrixS)inData[1]).subtract(bbT, ring);
                break;
            }
            // todo ask if is needed to change recursive methods to non-recursive
            case(7702):
            case(7718): {
                MatrixS A = (MatrixS) inData[0];
                MatrixS B = (MatrixS) inData[1];
                Element d = inData[2];
                outData[0] = A.multiplyDivRecursive(B, d.negate(ring), ring);
                break;
            }

            case(7703): {
                // todo be aware it can possibly not work where is used AdjMatrix
                AdjMatrixS m11 = (AdjMatrixS) inData[0];
                MatrixS M12 = (MatrixS) inData[1];
                Element d0 = inData[2];
                Element finalN = inData[3];
                MatrixS M12_1 = m11.A.multiplyDivRecursive(M12, d0, ring);
                MatrixS M12_2 = M12_1.multiplyLeftI(Array.involution(m11.Ei, (int) finalN.value));
                outData[0] = M12_1;
                outData[1] = M12_2;
                break;
            }

            case(7705): {
                MatrixS M22 = (MatrixS) inData[0];
                Element d11 = inData[1];
                MatrixS M21 = (MatrixS) inData[2];
                MatrixS M12_1 = (MatrixS) inData[3];
                AdjMatrixS m11 = (AdjMatrixS) inData[4];
                Element d0 = inData[5];
                outData[0] = ((M22.multiplyByNumber(d11, ring))
                        .subtract(M21.multiplyRecursive(M12_1.multiplyLeftE(m11.Ej, m11.Ei), ring), ring))
                        .divideByNumber(d0, ring);
                break;
            }

            case(7707): {
                MatrixS A = ((AdjMatrixS) inData[0]).S;
                MatrixS B = (MatrixS) inData[1];
                Element d = inData[2];
                outData[0] = A.multiplyDivRecursive(B, d.negate(ring), ring);
                break;
            }

            case(7708):
            case(7720): {
                MatrixS A1 = ((AdjMatrixS) inData[0]).A;
                MatrixS A2 = ((AdjMatrixS) inData[1]).A;
                outData[0] = A1.multiply(A2, ring);
                break;
            }

            case(7709): {
                MatrixS A21 = ((AdjMatrixS) inData[0]).A;
                MatrixS M22_1 = (MatrixS) inData[1];
                outData[0] = A21.multiply(M22_1, ring).negate(ring);
                break;
            }

            case(7710): {
                AdjMatrixS m11 = (AdjMatrixS) inData[0];
                AdjMatrixS m21 = (AdjMatrixS) inData[1];
                Element d11 = inData[2];
                outData[0] = m11.S.multiplyDivRecursive(m21.A.multiplyLeftE(m21.Ej, m21.Ei), d11, ring);
                break;
            }

            case(7711):
            case(7721): {
                MatrixS M22_1 = (MatrixS) inData[0];
                MatrixS A1 = (MatrixS) inData[1];
                AdjMatrixS m12 = (AdjMatrixS) inData[2];
                Element d11 = inData[3];
                outData[0] = M22_1.multiplyDivRecursive(A1.multiplyLeftE(m12.Ej, m12.Ei), d11, ring);
                break;
            }

            case(7712): {
                MatrixS B = (MatrixS) inData[0];
                MatrixS y12 = (MatrixS) inData[1];
                Element d11 = inData[2];
                Element d21 = inData[3];
                Element d12 = inData[4];
                AdjMatrixS m21 = (AdjMatrixS) inData[5];
                Element finalN = inData[6];
                Element d11_2 = d11.multiply(d11, ring);
                MatrixS M22_2 = B.multiplyDivRecursive(y12, d11_2.negate(ring), ring);
                Element ds = d12.multiply(d21, ring).divide(d11, ring);
                MatrixS M22_3 = M22_2.multiplyLeftI(Array.involution(m21.Ei, (int) finalN.value));
                outData[0] = M22_2;
                outData[1] = ds;
                outData[2] = M22_3;
                break;
            }

            case(7714): {
                MatrixS M21 = (MatrixS) inData[0];
                AdjMatrixS m11 = (AdjMatrixS) inData[1];
                Element d0 = inData[2];
                Element d12 = inData[3];
                MatrixS K2 = (MatrixS) inData[4];
                Element d11 = inData[5];
                outData[0] = (M21.multiplyDivMulRecursive(m11.A.multiplyLeftE(m11.Ej, m11.Ei), d0, d12, ring).add(K2, ring))
                        .divideByNumber(d11.negate(ring), ring);
                break;
            }

            case(7716): {
                MatrixS Q1 = (MatrixS) inData[0];
                MatrixS M12_1 = (MatrixS) inData[1];
                AdjMatrixS m11 = (AdjMatrixS) inData[2];
                Element d21 = inData[3];
                Element d11 = inData[4];
                MatrixS y12 = (MatrixS) inData[5];
                AdjMatrixS m12 = (AdjMatrixS) inData[6];
                outData[0] = (
                ((Q1.subtract((M12_1.multiplyLeftI(m11.Ei).multiplyByNumber(d21, ring)), ring))
                    .divideByNumber(d11, ring).multiplyRecursive(y12, ring))
                        .add((m12.S).multiplyByNumber(d21, ring), ring)
                )
                        .divideByNumber(d11, ring);
                break;
            }

            case(7717): {
                MatrixS A1 = (MatrixS) inData[0];
                MatrixS M12_1 = (MatrixS) inData[1];
                AdjMatrixS m11 = (AdjMatrixS) inData[2];
                AdjMatrixS m12 = (AdjMatrixS) inData[3];
                Element d11 = inData[4];
                Element d22 = inData[5];
                outData[0] = (A1.subtract((M12_1.multiplyLeftI(m11.Ei)).
                        multiplyDivRecursive(A1.multiplyLeftE(m12.Ej, m12.Ei), d11, ring), ring)
                ).divideMultiply(d11, d22, ring);
                break;
            }

            case(7719): {
                MatrixS M22_2 = (MatrixS) inData[0];
                AdjMatrixS m21 = (AdjMatrixS) inData[1];
                MatrixS y22 = (MatrixS) inData[2];
                Element ds = inData[3];
                AdjMatrixS m22 = (AdjMatrixS) inData[4];
                outData[0] = ((M22_2.multiplyLeftI(m21.Ei))
                        .multiplyDivRecursive(y22, ds.negate(ring), ring)).add(m22.S, ring);
                break;
            }

            case(7722): {
                MatrixS A2 = (MatrixS) inData[0];
                MatrixS M22_2 = (MatrixS) inData[1];
                AdjMatrixS m21 = (AdjMatrixS) inData[2];
                AdjMatrixS m22 = (AdjMatrixS) inData[3];
                Element ds = inData[4];
                Element d21 = inData[5];
                outData[0] = (A2.subtract((M22_2.multiplyLeftI(m21.Ei)).
                        multiplyDivRecursive(A2.multiplyLeftE(m22.Ej, m22.Ei), ds, ring), ring)
                ).divideByNumber(d21, ring);
                break;
            }

            case(7723): {
                AdjMatrixS m11 = (AdjMatrixS) inData[0];
                AdjMatrixS m21 = (AdjMatrixS) inData[1];
                Element d11 = inData[2];
                Element d22 = inData[3];
                MatrixS K1 = (MatrixS) inData[4];
                Element d21 = inData[5];
                outData[0] = (m11.S.multiplyDivMulRecursive(m21.A.multiplyLeftE(m21.Ej, m21.Ei), d11, d22, ring).add(K1, ring))
                        .divideByNumber(d21.negate(ring), ring);
                break;
            }

            case(7724): {
                MatrixS P = (MatrixS) inData[0];
                MatrixS G = (MatrixS) inData[1];
                Element d12 = inData[2];
                outData[0] = P.multiplyDivRecursive(G, d12, ring);
                break;
            }

            case(7725): {
                MatrixS L = (MatrixS) inData[0];
                MatrixS F = (MatrixS) inData[1];
                MatrixS G = (MatrixS) inData[2];
                Element d12 = inData[3];
                outData[0] = (L.add(F.multiplyRecursive(G, ring), ring)).divideByNumber(d12, ring);
                break;
            }
        }
    }

    @Override
    //Вхідна функція дропа, розбиває вхідні дані на блоки.
    public MatrixS[] inputFunction(Element[] input, Amin amin, Ring ring) {
        //LOGGER.info(input[0]);
        MatrixS[] res = new MatrixS[8];
        MatrixS v1;
        MatrixS v2;

        switch (key) {
            case(0):
            case(1):
            case(7702):
            case(7712):
            case(7713):
            case(7718):
            case(7724):
            default: {
                v1 = (MatrixS) input[0];
                v2 = (MatrixS) input[1];
                break;
            }
            case(2):{
                v1 =  ((MatrixS) input[0]).transpose();
                v2 = (MatrixS) input[0];
                amin.resultForOutFunction[4] = v1;
                break;
            }
            case(7703):{
                v1 = ((AdjMatrixS) input[0]).A;
                v2 = (MatrixS) input[1];
                break;
            }
            case(7705): {
                v1 = (MatrixS) inData[2];
                MatrixS M12_1 = (MatrixS) inData[3];
                AdjMatrixS m11 = (AdjMatrixS) inData[4];
                v2 = M12_1.multiplyLeftE(m11.Ej, m11.Ei);
                break;
            }
            case(7707): {
                v1 = ((AdjMatrixS) input[0]).S;
                v2 = (MatrixS) inData[1];
                break;
            }
            case(7708):
            case(7720): {
                v1 = ((AdjMatrixS) input[0]).A;
                v2 = ((AdjMatrixS) input[1]).A;
                break;
            }

            case(7709): {
                v1 = ((AdjMatrixS) input[0]).A;
                v2 = (MatrixS) inData[1];
                break;
            }
            case(7710): {
                MatrixS S11 = ((AdjMatrixS) input[0]).S;
                AdjMatrixS m21 = (AdjMatrixS) inData[1];
                v1 = S11;
                v2 = m21.A.multiplyLeftE(m21.Ej, m21.Ei);
                break;
            }
            case(7711):
            case(7721): {
                MatrixS M22_1 = (MatrixS) inData[0];
                MatrixS A1 = (MatrixS) inData[1];
                AdjMatrixS m12 = (AdjMatrixS) inData[2];
                v1 = M22_1;
                v2 = A1.multiplyLeftE(m12.Ej, m12.Ei);
                break;
            }
            case(7714): {
                v1 = (MatrixS) inData[0];
                AdjMatrixS m11 = (AdjMatrixS) inData[1];
                v2 = m11.A.multiplyLeftE(m11.Ej, m11.Ei);
                break;
            }
            case(7716): {
                MatrixS Q1 = (MatrixS) inData[0];
                MatrixS M12_1 = (MatrixS) inData[1];
                AdjMatrixS m11 = (AdjMatrixS) inData[2];
                Element d21 = inData[3];
                Element d11 = inData[4];
                v1 = (Q1.subtract((M12_1.multiplyLeftI(m11.Ei).multiplyByNumber(d21, ring)), ring))
                        .divideByNumber(d11, ring);
                v2 = (MatrixS) inData[5];
                break;
            }
            case(7717):
            case(7722): {
                MatrixS A1 = (MatrixS) inData[0];
                MatrixS M12_1 = (MatrixS) inData[1];
                AdjMatrixS m11 = (AdjMatrixS) inData[2];
                AdjMatrixS m12 = (AdjMatrixS) inData[3];
                v1 = M12_1.multiplyLeftI(m11.Ei);
                v2 = A1.multiplyLeftE(m12.Ej, m12.Ei);
                break;
            }
            case(7719): {
                MatrixS M22_2 = (MatrixS) inData[0];
                AdjMatrixS m21 = (AdjMatrixS) inData[1];
                MatrixS y22 = (MatrixS) inData[2];
                v1 = M22_2.multiplyLeftI(m21.Ei);
                v2 = y22;
                break;
            }
            case(7723): {
                v1 = ((AdjMatrixS) input[0]).S;
                AdjMatrixS m21 = (AdjMatrixS) inData[1];
                v2 = m21.A.multiplyLeftE(m21.Ej, m21.Ei);
                break;
            }
            case(7725): {
                v1 = (MatrixS) inData[1];
                v2 = (MatrixS) inData[2];
            }
        }

        Array.concatTwoArrays(v1.split(), v2.split(), res);
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
