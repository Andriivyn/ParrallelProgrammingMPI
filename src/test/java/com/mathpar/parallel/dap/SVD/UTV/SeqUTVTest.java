package com.mathpar.parallel.dap.SVD.UTV;

import com.mathpar.matrix.MatrixS;
import com.mathpar.number.Element;
import com.mathpar.number.Ring;
import org.junit.Before;
import org.junit.Test;

import java.util.Arrays;
import java.util.Random;

import static org.junit.Assert.*;

public class SeqUTVTest {
    private static int size = 16;
    private static Random random = new Random();
    private static Ring ring = new Ring("R64[]");

    private static int density = 2500;
    private MatrixS M;

    @Before
    public void setUp(){
        M = matrix(size);
    }

    @Test
    public void testUTV(){

        MatrixS[] UTV = SeqUTV.compute(M, ring);

        MatrixS transposedU = UTV[0].transpose();
        MatrixS transposedV = UTV[2].transpose();
        MatrixS resultCheck = transposedU.multiply(UTV[1].multiply(transposedV, ring), ring);
        System.out.println("M = "+M+";");
        System.out.println("U = "+UTV[0]+";");
        System.out.println("T = "+UTV[1]+";");
        System.out.println("V = "+UTV[2]+";");
        assertTrue(M.equals(resultCheck, ring));
        assertTrue(UTV[1].equals(UTV[0].multiply(M.multiply(UTV[2], ring), ring), ring));
    }

    @Test
    public void testEliminateColumn(){

        for (int col = 0; col < M.size - 1; col++) {
            MatrixS[] UT = SeqUTV.eliminateColumn(M, col, col+2, ring);

            Element[] column =  UT[1].getCol(col, ring);

            for(int row = column.length - 1; row > col + 1; row--){
                assertTrue(column[row].isZero(ring));
            }

            for(int row = col; row >= 0; row-- ){
                assertFalse(column[row].isZero(ring));
            }
            System.out.println("Eliminate Column (shift = "+ (col +2)+ "):\n"+ Arrays.toString(column));
        }


    }

    @Test
    public void testEliminateRow(){

        for (int rowN = 0; rowN < M.size - 1; rowN++) {

            MatrixS[] TV = SeqUTV.eliminateRow(M, rowN, rowN+2, ring);

            Element[] row =  TV[0].getRow(rowN, ring);

            for(int col = row.length - 1; col > rowN + 1; col--){
                assertTrue(row[col].isZero(ring));
            }

            System.out.println("Eliminate Row (shift = "+ (rowN +2)+ "):\n"+ Arrays.toString(row));
        }


    }

    private MatrixS matrix(int size) {
        return new MatrixS(size, size, density, new int[]{5}, random, ring.numberONE(), ring);
    }
}
