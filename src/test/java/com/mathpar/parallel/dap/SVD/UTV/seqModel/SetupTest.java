package com.mathpar.parallel.dap.SVD.UTV.seqModel;

import com.mathpar.number.Ring;
import com.mathpar.parallel.dap.SVD.UTV.seqModel.setup.NSetup;
import org.junit.After;
import org.junit.Before;
import org.junit.Test;

import static org.junit.Assert.*;

public class SetupTest {
    Ring ring;
    @Before
    public void setUp() throws Exception {
        ring =  new Ring("R64[]");
    }
    @Test
    public void NSetupTest(){
        int N = 6;

        NSetup setup = new NSetup(N, ring);

        for (int rank = 0; rank < N; rank++) {
            System.out.println(String.format("Rank=%d virtual ranks=%s", rank, setup.getVirtualProcessorRanks(rank)));
        }
    }
}
