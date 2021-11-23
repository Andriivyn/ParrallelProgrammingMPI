package com.mathpar.NAUKMA.MAG21.nguyen;

import mpi.MPI;
import mpi.MPIException;

import java.util.Arrays;

//Input:
//mpirun --hostfile hostfile -np 2 java -cp target/classes com/mathpar/NAUKMA/MAG21/nguyen/TestAllGatherv

//Output:
//myrank = 0: a = [0, 1, 2, 3, 4]
//myrank = 1: a = [10, 11, 12, 13, 14]
//myrank = 0: q = [10, 11, 12, 13, 14, 0, 1, 2, 3, 4]
//myrank = 1: q = [10, 11, 12, 13, 14, 0, 1, 2, 3, 4]

public class TestAllGatherv {
    public static void main(String[] args) throws MPIException,
            InterruptedException {

        // ініціалізація MPI
        MPI.Init(args);
        // визначення номера процесора
        int myrank = MPI.COMM_WORLD.getRank();
        int n = 5;
        int[] a = new int[n];
        int np = MPI.COMM_WORLD.getSize();
        for (int i = 0; i < n; i++) {
            a[i] = myrank*10+i;
        }
        System.out.println("myrank = " + myrank + ": a = " + Arrays.toString(a));
        int[] q = new int[n * np];
        MPI.COMM_WORLD.allGatherv(a, n, MPI.INT,
                q, new int[]{n, n}, new int[]{5, 0}, MPI.INT);

        System.out.println("myrank = " + myrank + ": q = " + Arrays.toString(q));
        // завершення паралельної частини
        MPI.Finalize();
    }
}
