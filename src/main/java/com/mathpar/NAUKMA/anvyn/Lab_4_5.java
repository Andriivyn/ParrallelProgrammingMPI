package com.mathpar.NAUKMA.anvyn;
import java.util.Arrays;
import java.math.BigInteger;
import java.util.Random;
import mpi.*;

public class Lab_4_5 {
    public static void main(String[] args)
            throws MPIException {
        MPI.Init(args);
        int myrank = MPI.COMM_WORLD.getRank();
        int n = 5;
        char[] a = new char[n];
        if (myrank == 3) {
            for (int i = 0; i < n; i++) {
                a[i] = (char) (myrank *22+i);
            }
            System.out.println("Sender rank = " + myrank + " : a = "
                    + Arrays.toString(a));
        }
        MPI.COMM_WORLD.bcast(a, a.length, MPI.INT, 3);
        if (myrank != 3)
            System.out.println(" Receiver rank = " + myrank + " : a = "
                    + Arrays.toString(a));
        MPI.Finalize();
    }
    }
