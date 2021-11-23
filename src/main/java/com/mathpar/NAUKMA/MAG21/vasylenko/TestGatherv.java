package com.mathpar.NAUKMA.MAG21.vasylenko;

import mpi.MPI;
import mpi.MPIException;

import java.util.Arrays;

// mpirun --hostfile /home/vasya/hostfile -np 2 java -cp /home/vasya/DAP/target/classes com/mathpar/NAUKMA/MAG21/vasylenko/TestGatherv

public class TestGatherv {
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
        int[] q = new int[n * np];
        // Процесор root збирає та зберiгає данi в
        //буферi recvbuf, розташовуючи їх у порядку зростання номерiв
        //процесорiв i вiдповiдно до змiщень, заданих у масивi displs
        MPI.COMM_WORLD.gatherv(a, n, MPI.INT, q,
                new int[]{n, n}, new int[]{5, 0},
                MPI.INT, np - 1);

        if (myrank == np - 1) {
            System.out.println("myrank = " + myrank + ": q = " + Arrays.toString(q));
        }
        // завершення паралельної частини
        MPI.Finalize();
    }
}

//        myrank = 0: a = [0, 1, 2, 3, 4]
//        myrank = 1: a = [10, 11, 12, 13, 14]
//        myrank = 1: q = [10, 11, 12, 13, 14, 0, 1, 2, 3, 4]
