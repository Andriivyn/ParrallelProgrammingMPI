package com.mathpar.NAUKMA.MAG21.riepkin;

import mpi.MPI;
import mpi.MPIException;

import java.util.Arrays;
/*
Creates array and fills it on all processes.
All processes send it's data.
All data is accumulated on main process.
Generates array with summ of all values on the same position.
*/
public class TestReduce {
    public static void main(String[] args) throws MPIException {
        // ініціалізація MPI
        MPI.Init(args);
        // визначення номера процесора
        int myrank = MPI.COMM_WORLD.getRank();
        int n = 5;
        int[] a = new int[n];
        for (int i = 0; i < n; i++) {
            a[i] = myrank*10+i;
        }
        System.out.println("myrank = " + myrank + ": a = " + Arrays.toString(a));
        int[] q = new int[n];
        MPI.COMM_WORLD.reduce(a, q, n, MPI.INT, MPI.SUM, 0);
        if (myrank == 0)
            System.out.println("myrank = " + myrank + ": q = " + Arrays.toString(q));
        // завершення паралельної частини
        MPI.Finalize();
    }
}
/*
Output:
root@b1a37b43d450:/DAP/src/main/java/com/mathpar/NAUKMA/riepkin# sh run_test.sh 2 TestReduce
myrank = 0: a = [0, 1, 2, 3, 4]
myrank = 1: a = [10, 11, 12, 13, 14]
myrank = 0: q = [10, 12, 14, 16, 18]
root@b1a37b43d450:/DAP/src/main/java/com/mathpar/NAUKMA/riepkin# sh run_test.sh 4 TestReduce
myrank = 1: a = [10, 11, 12, 13, 14]
myrank = 2: a = [20, 21, 22, 23, 24]
myrank = 0: a = [0, 1, 2, 3, 4]
myrank = 3: a = [30, 31, 32, 33, 34]
myrank = 0: q = [60, 64, 68, 72, 76]
*/
