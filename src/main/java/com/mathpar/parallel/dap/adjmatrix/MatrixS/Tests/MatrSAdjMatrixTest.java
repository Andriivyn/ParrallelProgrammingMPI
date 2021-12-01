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

public class MatrSAdjMatrixTest extends DAPTest {
    private final static MpiLogger LOGGER = MpiLogger.getLogger(MatrSMult4Test.class);
    protected MatrSAdjMatrixTest() {
        super("MatrSAdjMatrix", 7701, 0);
    }

    @Override
    protected Element[] initData(int size, int density, int maxBits, Ring ring) {
        MatrixS M = matrix(size, density, maxBits, ring);
        return new Element[]{M, M.det(ring)};
    }

    @Override
    protected Pair<Boolean, Element> checkResult(DispThread dispThread, String[] args, Element[] initData, Element[] resultData, Ring ring) {
        AdjMatrixS resAdj = (AdjMatrixS) resultData[0];
        MatrixS initM = (MatrixS) initData[0];
        Element d0 = initData[1];
        AdjMatrixS toCompare = new AdjMatrixS(initM, d0, ring);
        boolean succeed = resAdj.A.equals(toCompare.A);
        return new Pair<>(succeed, resAdj.Det);
    }

    @Override
    protected MatrixS matrix(int size, int density, int maxBits, Ring ring){
        MatrixS matrix = new MatrixS(size, size, density, new int[]{maxBits}, new Random(),ring.numberONE(), ring);
        // LOGGER.trace("bef matrix = " + matrix);
        for (int i = 0; i < size; i++) {
            for (int j = 0; j < size; j++) {
                if(!matrix.getElement(i,j, ring).equals(ring.numberZERO))
                    matrix.putElement( ring.numberONE.divide(matrix.getElement(i,j, ring),  ring), i, j);
            }
        }
        //LOGGER.info("matrix = " + matrix);
        return matrix;
    }

    public static void main(String[] args) throws InterruptedException, ClassNotFoundException, MPIException, IOException {
        MatrSAdjMatrixTest test = new MatrSAdjMatrixTest();
        test.runTests(args);
    }
}
