package com.mathpar.NAUKMA.MAG21.mishakushka;

import mpi.MPI;
import mpi.MPIException;

import java.util.Arrays;

public class TestReduceScatter {
    public static void main(String[] args) throws MPIException {
        // ініціалізація MPI
        MPI.Init(args);
        // визначення номера процесора
        int myrank = MPI.COMM_WORLD.getRank();
        int n = 10;
        int[] a = new int[n];
        for (int i = 0; i < n; i++) {
            a[i] = myrank*10+i;
        }
        System.out.println("myrank = " + myrank + ": a = " + Arrays.toString(a));
        int[] q = new int[n];
        MPI.COMM_WORLD.reduceScatter(a, q, new int[]{1, 2, 3, 4},
                MPI.INT, MPI.SUM);

       // if (myrank == 0) {
            System.out.println("myrank = " + myrank + ": q = " + Arrays.toString(q));
       // }
        // завершення паралельної частини
        MPI.Finalize();
    }
}

/**
 * Output:
 * myrank = 0: a = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
 * myrank = 0: q = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
 */