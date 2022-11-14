package com.mathpar.NAUKMA.anvyn;

import mpi.MPI;
import mpi.MPIException;

import java.util.Arrays;

public class Lab_4_6 {
    public static void main(String[] args) throws MPIException {
        MPI.Init(args);
        int myRank = MPI.COMM_WORLD.getRank();
        int n = 5;
        int[] a = new int[n];

        for (int i = 0; i < n; i++) {
            a[i] = myRank*10+i;
        }

        System.out.println("myRank = " + myRank + ": a = " + Arrays.toString(a));

        int[] q = new int[n];
        MPI.COMM_WORLD.reduce(a, q, a.length, MPI.SHORT_INT, MPI.MINLOC, 2);

        if (myRank == 2) {
            //мiнiмальне значення i його iндекс.
            System.out.println("MINLOC = " + myRank + ": q = " + Arrays.toString(q));
        }

        MPI.Finalize();
    }
}
