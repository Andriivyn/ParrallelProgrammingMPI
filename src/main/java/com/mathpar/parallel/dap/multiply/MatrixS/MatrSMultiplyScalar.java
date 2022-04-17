package com.mathpar.parallel.dap.multiply.MatrixS;

import com.mathpar.log.MpiLogger;
import com.mathpar.matrix.MatrixS;
import com.mathpar.number.Element;
import com.mathpar.number.Ring;
import com.mathpar.parallel.dap.core.Amin;
import com.mathpar.parallel.dap.core.Drop;

import java.util.ArrayList;

public class MatrSMultiplyScalar extends Drop {
    private final static MpiLogger LOGGER = MpiLogger.getLogger(MatrSMultiplyScalar.class);
    protected static int leafSize = 4;

    private static int[][] _arcs = new int[][]{
            //Зв'язки від вхідної функції до всіх інших дропів
            {1, 0, 0, 1, 1, 1, 2, 2, 0, 2, 3, 1},
            {3, 0, 0},
            {3, 0, 1},
            {}};

    public MatrSMultiplyScalar() {
        inData =  new Element[4];
        outData =  new Element[1];
        //Дроп має тип 7
        type = 7;
        //кількість блоків, для формування результату
        resultForOutFunctionLength = 2;
        inputDataLength = 4;
        //унікальний номер дропа
        number = cnum++;
        arcs = _arcs;
    }

   /* @Override
    public void setVars(){
       return;
    }
*/
    //Розгортання аміну з дропами, відповідно до графу, для обрахунку поточного дропа.
    @Override
    public ArrayList<Drop> doAmin() {
        ArrayList<Drop> amin = new ArrayList<Drop>();

        amin.add(new MatrSMult4());
        amin.add(new MatrSMult4());

        return amin;
    }

    //Послідовний обрахунок листових вершин
    @Override
    public void sequentialCalc(Ring ring) {
        // LOGGER.info("in sequentialCalc indata = " + inData[0] + ",  "+inData[1]);
        MatrixS A = (MatrixS) inData[0];
        MatrixS B = (MatrixS) inData[1];
        MatrixS C = (MatrixS) inData[2];
        MatrixS D = (MatrixS) inData[3];
        MatrixS R = A.multiply(B, ring).add(C.multiply(D, ring), ring);

        outData[0] = R;
    }

    @Override
    //Вхідна функція дропа, розбиває вхідні дані на блоки.
    public MatrixS[] inputFunction(Element[] input, Amin amin, Ring ring) {
        MatrixS[] res = {(MatrixS) input[0], (MatrixS) input[1], (MatrixS) input[2], (MatrixS) input[3]};
        return res;
    }

    //Вихідна функція дропа, яка збирає блоки в результат
    @Override
    public MatrixS[] outputFunction(Element[] input, com.mathpar.number.Ring ring) {

        MatrixS A = (MatrixS) input[0];
        MatrixS B = (MatrixS) input[1];
//        LOGGER.info("A: " + A.getElement(0, 0, ring));
//        LOGGER.info("B: " + B);
        MatrixS[] res = new MatrixS[]{A.add(B, ring)};
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
