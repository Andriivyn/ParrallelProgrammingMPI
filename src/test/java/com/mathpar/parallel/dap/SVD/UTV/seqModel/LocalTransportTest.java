package com.mathpar.parallel.dap.SVD.UTV.seqModel;

import com.mathpar.parallel.dap.SVD.UTV.seqModel.transport.LocalTransport;
import com.mathpar.parallel.dap.SVD.UTV.seqModel.transport.Transport;
import org.junit.Before;
import org.junit.Test;

import static org.junit.Assert.*;

public class LocalTransportTest {
    LocalTransport transport;
    @Before
    public void setUp() throws Exception {
        transport = new LocalTransport(4);
    }

    @Test
    public void receiveObject() {

        Object msg1 = transport.receiveObject(1, 2, Transport.Tag.ANY);

        assertNull(msg1);

        transport.sendObject(1, 1, 2, Transport.Tag.ANY);

        Object msg2 = transport.receiveObject(1, 2, Transport.Tag.ANY);

        assertNotNull(msg2);
        assertEquals(msg2, 1);

    }

    @Test
    public void receiveObjectArray() {

        Object[] msg1 = transport.receiveObjectArray(1, 2, Transport.Tag.ANY);

        assertNull(msg1);

        transport.sendObjectArray(new Integer[]{1,2,3}, 1, 2, Transport.Tag.ANY);

        Object[] msg2 = transport.receiveObjectArray(1, 2, Transport.Tag.ANY);
        Object[] msg3 = transport.receiveObjectArray(1, 2, Transport.Tag.ANY);

        assertNull(msg3);
        assertNotNull(msg2);
        assertArrayEquals(msg2, new Integer[]{1,2,3});


    }

    @Test
    public void sendObjectArray() {

        transport.sendObjectArray(new Integer[]{1,1,1}, 1, 2, Transport.Tag.ANY);
        transport.sendObjectArray(new Integer[]{2,2,2}, 1, 2, Transport.Tag.ANY);
        transport.sendObjectArray(new Integer[]{3,3,3}, 1, 2, Transport.Tag.ANY);

        Object[] msg1 = transport.receiveObjectArray(1, 2, Transport.Tag.ANY);
        Object[] msg2 = transport.receiveObjectArray(1, 2, Transport.Tag.ANY);
        Object[] msg3 = transport.receiveObjectArray(1, 2, Transport.Tag.ANY);

        assertNotNull(msg1);
        assertNotNull(msg2);
        assertNotNull(msg3);
        assertArrayEquals(msg1, new Integer[]{1,1,1});
        assertArrayEquals(msg2, new Integer[]{2,2,2});
        assertArrayEquals(msg3, new Integer[]{3,3,3});

    }
}