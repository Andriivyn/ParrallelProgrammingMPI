package com.mathpar.parallel.dap.adjmatrix.MatrixS.Tests;

import com.mathpar.log.MpiLogger;
import com.mathpar.matrix.AdjMatrixS;
import com.mathpar.matrix.MatrixS;
import com.mathpar.number.Element;
import com.mathpar.number.Ring;
import com.mathpar.parallel.dap.core.DispThread;
import com.mathpar.parallel.dap.multiply.MatrixS.Tests.MatrSMult4Test;
import com.mathpar.parallel.dap.test.DAPTest;
import mpi.MPIException;
import org.javatuples.Pair;

import java.io.IOException;
import java.util.Random;
import java.util.logging.Logger;

public class MatrSAdjMatrixTest extends DAPTest {
    private final static MpiLogger LOGGER = MpiLogger.getLogger(MatrSMult4Test.class);
    protected MatrSAdjMatrixTest() {
        super("MatrSAdjMatrix", 7701, 0);
        //ring = new Ring("Z[]");
    }

    @Override
    protected Element[] initData(int size, int density, int maxBits, Ring ring) {
        MatrixS M = matrix(size, density, maxBits, ring);
        return new Element[]{M, ring.numberONE};
    }

    @Override
    protected Pair<Boolean, Element> checkResult(DispThread dispThread, String[] args, Element[] initData, Element[] resultData, Ring ring) {
        MatrixS initM = (MatrixS) initData[0];
        LOGGER.info("Input matrix: " + initM);
        AdjMatrixS resAdj = (AdjMatrixS) resultData[0];
        //AdjMatrixS resAdj = new AdjMatrixS(initM, ring.numberONE,  ring);
        LOGGER.info("Adj matrix = " + resAdj.A);
        LOGGER.info("Output matrix det = " + resAdj.Det);
        MatrixS divided = resAdj.A.divideByNumbertoFraction(resAdj.Det, ring);
        MatrixS rr = initM.multiply(divided, ring);
        LOGGER.info("Output: " + rr);
        boolean succeed = rr.isOne(ring);
        return new Pair<>(succeed, null);
    }

    @Override
    protected int dispRunsOnOtherProc() {
        return 0;
    }

    @Override
    protected MatrixS matrix(int size, int density, int maxBits, Ring ring){
//        int [][]mat = {{6, -5, 8, 4}, {9,7,5,2}, {7,5,3,7}, {-4,8,-8,-3}};
//        MatrixS matrix = new MatrixS(mat, ring);
        MatrixS matrix = new MatrixS(size, size, density, new int[]{maxBits}, new Random(),ring.numberONE(), ring);
        // LOGGER.trace("bef matrix = " + matrix);
//        for (int i = 0; i < size; i++) {
//            for (int j = 0; j < size; j++) {
//                if(!matrix.getElement(i,j, ring).equals(ring.numberZERO))
//                    matrix.putElement( ring.numberONE.divide(matrix.getElement(i,j, ring),  ring), i, j);
//            }
//        }
        //LOGGER.info("matrix = " + matrix);
        return matrix;
    }

    @Override
    protected Element[] sequentialExecute(Element[] data, Ring ring) {
        MatrixS m = (MatrixS) data[0];
        Element d0 = data[1];
        AdjMatrixS adjM = new AdjMatrixS(m, d0,  ring);
        Element resD = adjM.Det;
        MatrixS y = adjM.S.ES_min_dI(resD, adjM.Ei, adjM.Ej, ring);
        return new Element[] {adjM, y, resD};
    }

    public static void main(String[] args) throws InterruptedException, ClassNotFoundException, MPIException, IOException {
        MatrSAdjMatrixTest test = new MatrSAdjMatrixTest();
        test.runTests(args);
    }
}
