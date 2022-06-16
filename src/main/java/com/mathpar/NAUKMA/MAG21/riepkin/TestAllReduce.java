package com.mathpar.NAUKMA.MAG21.riepkin;

import mpi.MPI;
import mpi.MPIException;

import java.util.Arrays;
/*
Each process creates array with 5 elements and fills it.
Sends own array to another processes.
MPI accumulates arrays from all processes and does multiplication of array elements.
*/
public class TestAllReduce {
    public static void main(String[] args) throws MPIException {
        // ініціалізація MPI
        MPI.Init(args);
        // визначення номера процесора
        int myrank = MPI.COMM_WORLD.getRank();
        // визначення числа процесорів в групі
        int np = MPI.COMM_WORLD.getSize();
        int n = 5;
        int[] a = new int[n];
        for (int i = 0; i < n; i++) {
            a[i] = myrank*10+i;
        }
        System.out.println("myrank = " + myrank + ": a = " + Arrays.toString(a));
        int[] q = new int[n];
        // accumulates results in q, multiplies elements
        MPI.COMM_WORLD.allReduce(a, q, n, MPI.INT, MPI.PROD);
        for (int j = 0; j < np; j++) {
            if (myrank == j) {
                System.out.println("myrank = " + myrank + ": q = " + Arrays.toString(q));
            }
        }
        // завершення параалельної частини
        MPI.Finalize();
    }
}
/*
Output:
root@b1a37b43d450:/DAP/src/main/java/com/mathpar/NAUKMA/riepkin# sh run_test.sh 2 TestAllReduce     
myrank = 0: a = [0, 1, 2, 3, 4]
myrank = 1: a = [10, 11, 12, 13, 14]
myrank = 0: q = [0, 11, 24, 39, 56]
myrank = 1: q = [0, 11, 24, 39, 56]
root@b1a37b43d450:/DAP/src/main/java/com/mathpar/NAUKMA/riepkin# sh run_test.sh 3 TestAllReduce
myrank = 1: a = [10, 11, 12, 13, 14]
myrank = 0: a = [0, 1, 2, 3, 4]
myrank = 2: a = [20, 21, 22, 23, 24]
myrank = 1: q = [0, 231, 528, 897, 1344]
myrank = 0: q = [0, 231, 528, 897, 1344]
myrank = 2: q = [0, 231, 528, 897, 1344]
*/
