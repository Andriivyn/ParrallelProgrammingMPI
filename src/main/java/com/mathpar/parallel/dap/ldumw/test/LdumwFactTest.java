package com.mathpar.parallel.dap.ldumw.test;

import com.mathpar.matrix.LDUMW;
import com.mathpar.matrix.MatrixS;
import com.mathpar.number.Element;
import com.mathpar.number.Ring;
import com.mathpar.parallel.dap.core.DispThread;
import com.mathpar.parallel.dap.ldumw.LdumwDto;
import com.mathpar.parallel.dap.ldumw.LdumwFact;
import com.mathpar.parallel.dap.test.DAPTest;
import mpi.MPIException;
import org.javatuples.Pair;

import java.io.IOException;

public class LdumwFactTest extends DAPTest {

    protected LdumwFactTest() {
        super("LdumwFactTest", 23, 0);
        //ring = new Ring("Z[]");
    }

    @Override
    protected Element[] initData(int size, int density, int maxBits, Ring ring) {
        MatrixS A = matrix(size, density, 5, ring);
        return new Element[]{A, ring.numberONE};
    }

    @Override
    protected Pair<Boolean, Element> checkResult(DispThread dispThread, String[] args, Element[] initData, Element[] resultData, Ring ring) {
        LdumwDto ldumwDto = (LdumwDto) resultData[0];
        ldumwDto.setD(LdumwFact.invForD(ldumwDto.D(), ring));

        MatrixS A = (MatrixS) initData[0];
        Element a = initData[1];

        LdumwDto ldumwDtoSequential = LDUMW.LDUMW(A, a, ring);

        if (ldumwDto.equals(ldumwDtoSequential)) {
            return new Pair<>(true, null);
        }

        return new Pair<>(false, null);
    }

    public static void main(String[] args) throws InterruptedException, ClassNotFoundException, MPIException, IOException {
        LdumwFactTest test = new LdumwFactTest();
        test.runTests(args);
    }

    @Override
    protected int dispRunsOnOtherProc() {
        return 0;
    }

    @Override
    protected Element[] sequentialExecute(Element[] data, Ring ring) {
        MatrixS A = (MatrixS) data[0];
        Element a = data[1];

        LdumwDto FF = LDUMW.LDUMW(A, a, ring);

        return new Element[] {FF, FF.A_n()};
    }
}
