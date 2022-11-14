package com.mathpar.NAUKMA.anvyn;

import java.util.Arrays;

import mpi.MPI;
import mpi.MPIException;

public class example {
    public static void main(String[] args)
            throws MPIException {

        MPI.Init(args);
        int myrank = MPI.COMM_WORLD.getRank();
        int np = MPI.COMM_WORLD.getSize();
        int n = 12;
        int[] a = new int[n];
        for (int i = 0; i < n; i++) a[i] = myrank * 10 + i;
        System.out.println("myrank = " + myrank + ": a = " + Arrays.toString(a));
        int[] q = new int[n * np];
        MPI.COMM_WORLD.gatherv(a, n, MPI.INT, q, new int[]{n, n},
                new int[]{0, n}, MPI.INT, np - 1);
        if (myrank == np - 1)
            System.out.println("myrank = " + myrank + ": q = " + Arrays.toString(q));
        MPI.Finalize();
    }
}
