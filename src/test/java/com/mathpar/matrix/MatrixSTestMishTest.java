/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package com.mathpar.matrix;

import com.mathpar.number.Element;
import com.mathpar.number.Ring;
import com.mathpar.polynom.Polynom;
import org.junit.Test;

import static org.junit.Assert.assertTrue;

/**
 *
 * @author gennadi
 */
public class MatrixSTestMishTest {
    
 
    /**
     * @param args the command line arguments
     */
    
     @Test
    public void  misha16() {
        Ring ring = new Ring("R[t]");
        ring.setMachineEpsilonR(50);
        ring.setDefaulRing();
        int n = 32;
        Element[][] A = new Element[n][];
        for (int i = 0; i < n; i++) {
            A[i] = new Element[n];
            for (int j = 0; j < n; j++) {
                if (j % 2 == 0) {
                    if (i == j) {
                        A[i][j] = new Polynom("t", ring);
                    } else {
                        A[i][j] = new Polynom("1", ring);
                    }
                } else {
                    if (i == j) {
                        A[i][j] = new Polynom("-t", ring);
                    } else {
                        A[i][j] = new Polynom("-1", ring);
                    }
                }
            }
        }
        //Преобразование матрицы полиномов в объект MatrixS
        MatrixS matrix = new MatrixS(A, ring);
        //Получение присоединенной матрицы
        MatrixS adjMatrix = matrix.adjoint(ring);
        System.out.println("inp="+matrix);
        System.out.println("result="+adjMatrix);
        assertTrue("Size=n:",
                (adjMatrix.size==n));
    }
}