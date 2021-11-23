package com.mathpar.matrix;

import com.mathpar.number.*;
import com.mathpar.parallel.dap.QR.Utils;
import org.junit.Test;
import com.mathpar.polynom.Polynom;

import java.util.Random;

import static org.junit.Assert.assertTrue;

public class MatrixSTest {

    private static int size = 32;
    private static int density = 100;
    private static Random random = new Random();
    private static Ring ring = new Ring("R64[]");

    @Test
    public void adjointPolynomialShouldBeExact() {
        Ring r = new Ring("Z[p]");
        MatrixS m = new MatrixS(new Element[][] {
            {new Polynom("3p + 2", r), new Polynom("p", r)},
            {new Polynom("p", r), new Polynom("4p + 3", r)}
        }, r);

        MatrixS adj = m.adjoint(r);
        assertTrue("Element at pos (0, 1) should be polynomial -p, not fraction",
                adj.getElement(0, 1, r).equals(new Polynom("-p", r), r));
    }

    @Test
    public void adjointTest() {
        int[][] mi={{1,0,2,3},{2,0,0,1},{1,2,0,3},{1,3,0,3}};
        Ring r = new Ring("Z[]");         MatrixS m=new MatrixS(mi,r);
        Element det=m.det(r);
        MatrixS adj = m.adjoint(r);
      //  System.out.println("m="+adj);
      //  System.out.println("det="+det);
        int[][] res={{0,-6,6,-4},{0,0,10,-10},{-5,0,15,-10},{0,2,-12,8}};
        assertTrue("Somthing wrong in Matrix Adjoint",
        (adj.subtract((new MatrixS(res, r)),r)).isZero(r));
    }

    @Test
    public void multiplySorted(){
        MatrixS A = matrix(size);
        MatrixS B = matrix(size);
      //  System.out.println("A="+A);
      //  System.out.println("B="+B);
        long t1 = System.currentTimeMillis();
        MatrixS AB = A.multiply(B, ring);
        long t2 = System.currentTimeMillis();
        MatrixS sortedAB = A.multiplySorted(B, ring);
        long t3 = System.currentTimeMillis();
        assertTrue(AB.subtract(sortedAB, ring).isZero(ring));
        long multiplyTime = t2 - t1;
        long multiplySortedTime = t3 - t2;
        long diff = Math.abs(multiplySortedTime - multiplyTime);
        System.out.println(String.format("multiply: %dms multiplySorted: %dms. Diff=%d", multiplyTime, multiplySortedTime, diff));

        for (int row = 0; row < sortedAB.size; row++) {
            for (int col = 0; col < sortedAB.col[row].length - 1; col++) {
                assertTrue("row="+row+" col="+col+  " "+sortedAB.col[row][col]+" "+sortedAB.col[row][col+1] ,sortedAB.col[row][col] < sortedAB.col[row][col + 1]);
            }
        }
    }

    @Test
    public void multiplySorted2(){
        MatrixS[] matrices = Utils.readMatrixArrayFromFile("crash_matrix");
        MatrixS rotationMatrix = matrices[0];
        MatrixS rowMatrix = matrices[1];
        MatrixS result = matrices[2];                
        MatrixS reCheckResult = rotationMatrix.multiplySorted(rowMatrix, ring);
        MatrixS simpleMultiplication = rotationMatrix.multiply(rowMatrix, ring);
 
    //    assertTrue(result.subtract(reCheckResult, ring).isZero(ring));
         assertTrue(simpleMultiplication.subtract(reCheckResult, ring).isZero(ring));

         for (int row = 0; row < reCheckResult.size; row++) {
            for (int col = 0; col < reCheckResult.col[row].length - 1; col++) {
               assertTrue("row="+row+" col="+col+  " cols values="+reCheckResult.col[row][col]+", "+reCheckResult.col[row][col+1] ,reCheckResult.col[row][col] < reCheckResult.col[row][col + 1]);
          }
        }


    }


    private MatrixS matrix(int size){
        return new MatrixS(size, size, density, new int[]{5}, random, ring.numberONE(), ring);
    }
    
  private MatrixS matrixDouble1( ){
          double[][] aD=new double[][]{{0.84, 0.54},{-0.54, 0.84}};
          return new MatrixS(new MatrixD(aD, ring), ring);
    }
          
   private MatrixS matrixDouble2( ){        
          double[][] bD=new double[][]{{ 0.91, 0.54, 0.98, 0.51, 0.41, 0, 0.29, 0.3, 0.24,
              0.81, 0, 0.96, 0.39, 0.34, 0.39, 0.95, 0.96, 0.56, 0.05, 0, 0, 0.24, 
              0.61, 0.95, 0, 0.73, 0.4, 0.24, 0.11, 0.15, 0, 0.83}
          ,{0.59, 0.54, 0.75, 0.66, 0.4, 0.94, 0.67, 0.77, 0.3, 0.41, 0.96, 0.36, 0.61, 
              0.07, 1.05, 0.74, 1.11, 0.39, 0.86, 0.78, 0.71, 0, 0.16, 0.74, 0.97, 1.08, 
              0.49, 0.62, 0.64, 0.96, 0.37, 1.01}};
      return new MatrixS(new MatrixD(bD, ring),  ring);
    }    
     
}
 