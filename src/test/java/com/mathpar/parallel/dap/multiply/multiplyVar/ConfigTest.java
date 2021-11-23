package com.mathpar.parallel.dap.multiply.multiplyVar;

import com.mathpar.number.Ring;
import org.junit.Test;

import static org.junit.Assert.*;

public class ConfigTest {

    private static Ring ring = new Ring("Z[]");

    @Test
    public void createNegatedMultiplication(){
        MultiplyVarConfig config = MultiplyVarConfig.Builder.startWith()
                                        .multiplicationBlock()
                                        .negateMultiplication()
                                        .setFirst(MultiplyVarMatrix.Builder.newBuilder().addFull().set().done().build())
                                        .setSecond(MultiplyVarMatrix.Builder.newBuilder().addFull().transpose().copyOf(0).done().build())
                                        .build();

        byte[] configData = config.getData();
        byte[] answer = new byte[14];
        answer[0] = 0x0f;
        answer[4] = 0x11;
        answer[12] = 1;

        assertFalse(config.isAdditionBlockNegated());
        assertTrue(config.isMultiplicationBlockNegated());
        assertFalse(config.hasAdditionalData());
        assertArrayEquals(answer, configData);

    }

    @Test
    public void createQuarteredMultiplicationWithNegatedCenterAddition(){
        MultiplyVarConfig config = MultiplyVarConfig.Builder.startWith()
                .multiplicationBlock()
                .negateMultiplication()
                .setFirst(MultiplyVarMatrix.Builder.newBuilder()
                        .addQuarter(MultiplyVarMatrix.Quarter.I).set().done()
                        .addQuarter(MultiplyVarMatrix.Quarter.III).transpose().set().done().build())
                .setSecond(MultiplyVarMatrix.Builder.newBuilder().addFull().transpose().copyOf(0).done().build())
                .additionBlock()
                .negate()
                .set(MultiplyVarMatrix.Builder.newBuilder().fillDiagonal(MultiplyVarMatrix.FillElement.ZERO).addCenter().copyOf(2).done().build())
                .build();

        byte[] configData = config.getData();
        byte[] answer = new byte[14];
        answer[0] = (byte) 0b10001111;
        answer[1] = (byte) 0b10000000;
        answer[2] = (byte) 0b10011111;
        answer[3] = (byte) 0b10000000;
        answer[4] = 0b00010001;
        answer[8] = 0b01000011;
        answer[9] = 0b01000000;
        answer[10] = 0b01000000;
        answer[11] = 0b01000000;
        answer[12] = 1;
        answer[13] = 1;

        assertTrue(config.isSet(0));
        assertFalse(config.isSet(1));
        assertTrue(config.isSet(2));
        assertFalse(config.isSet(3));

        assertTrue(config.isSet(4));
        assertFalse(config.isSet(5));
        assertFalse(config.isSet(6));
        assertFalse(config.isSet(7));

        assertTrue(config.isSet(8));
        assertFalse(config.isSet(9));
        assertFalse(config.isSet(10));
        assertFalse(config.isSet(11));

        assertEquals(MultiplyVarConfig.INPUT, config.dataSource(0));
        assertEquals(MultiplyVarConfig.NONE, config.dataSource(1));
        assertEquals(MultiplyVarConfig.INPUT, config.dataSource(2));
        assertEquals(MultiplyVarConfig.NONE, config.dataSource(3));

        assertEquals(0, config.dataSource(4));
        assertEquals(MultiplyVarConfig.NONE, config.dataSource(5));
        assertEquals(MultiplyVarConfig.NONE, config.dataSource(6));
        assertEquals(MultiplyVarConfig.NONE, config.dataSource(7));

        assertEquals(2, config.dataSource(8));
        assertEquals(MultiplyVarConfig.NONE, config.dataSource(9));
        assertEquals(MultiplyVarConfig.NONE, config.dataSource(10));
        assertEquals(MultiplyVarConfig.NONE, config.dataSource(11));

        assertFalse(config.isTransposed(0));
        assertFalse(config.isTransposed(1));
        assertTrue(config.isTransposed(2));
        assertFalse(config.isTransposed(3));

        assertTrue(config.isTransposed(4));
        assertFalse(config.isTransposed(5));
        assertFalse(config.isTransposed(6));
        assertFalse(config.isTransposed(7));

        assertFalse(config.isTransposed(8));
        assertFalse(config.isTransposed(9));
        assertFalse(config.isTransposed(10));
        assertFalse(config.isTransposed(11));


        assertTrue(config.isAdditionBlockNegated());
        assertTrue(config.isMultiplicationBlockNegated());

        assertEquals(MultiplyVarMatrix.Position.QUARTER, config.getPosition(0));
        assertEquals(MultiplyVarMatrix.Position.QUARTER, config.getPosition(1));
        assertEquals(MultiplyVarMatrix.Position.QUARTER, config.getPosition(2));
        assertEquals(MultiplyVarMatrix.Position.QUARTER, config.getPosition(3));

        assertEquals(MultiplyVarMatrix.Position.FULL, config.getPosition(4));
        assertEquals(MultiplyVarMatrix.Position.FULL, config.getPosition(5));
        assertEquals(MultiplyVarMatrix.Position.FULL, config.getPosition(6));
        assertEquals(MultiplyVarMatrix.Position.FULL, config.getPosition(7));

        assertEquals(MultiplyVarMatrix.Position.CENTER, config.getPosition(8));
        assertEquals(MultiplyVarMatrix.Position.CENTER, config.getPosition(9));
        assertEquals(MultiplyVarMatrix.Position.CENTER, config.getPosition(10));
        assertEquals(MultiplyVarMatrix.Position.CENTER, config.getPosition(11));

        assertEquals(ring.numberZERO(), config.getDiagonalElement(0, ring));
        assertEquals(ring.numberZERO(), config.getDiagonalElement(1, ring));
        assertEquals(ring.numberZERO(), config.getDiagonalElement(2, ring));
        assertEquals(ring.numberZERO(), config.getDiagonalElement(3, ring));

        assertEquals(ring.numberZERO(), config.getDiagonalElement(4, ring));
        assertEquals(ring.numberZERO(), config.getDiagonalElement(5, ring));
        assertEquals(ring.numberZERO(), config.getDiagonalElement(6, ring));
        assertEquals(ring.numberZERO(), config.getDiagonalElement(7, ring));

        assertEquals(ring.numberZERO(), config.getDiagonalElement(8, ring));
        assertEquals(ring.numberZERO(), config.getDiagonalElement(9, ring));
        assertEquals(ring.numberZERO(), config.getDiagonalElement(10, ring));
        assertEquals(ring.numberZERO(), config.getDiagonalElement(11, ring));

        assertTrue(config.hasAdditionalData());

        assertArrayEquals(answer, configData);

    }
}
