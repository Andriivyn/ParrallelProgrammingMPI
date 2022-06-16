package com.mathpar.NAUKMA.MAG21.vasylenko;

import mpi.MPI;
import mpi.MPIException;

import java.util.Arrays;

// mpirun --hostfile /home/vasya/hostfile -np 8 java -cp /home/vasya/DAP/target/classes com/mathpar/NAUKMA/MAG21/vasylenko/TestBcast

/**
 * Процесор з номером 0 пересилає масив чисел
 * іншим процесорам,
 * використовуючи роздачу за бінарним деревом.
 */
public class TestBcast {
    public static void main(String[] args)
            throws MPIException {
        // ініціалізація MPI
        MPI.Init(args);
        // визначення номера процесора
        int myrank = MPI.COMM_WORLD.getRank();
        int n = 5;
        int[] a = new int[n];
        if (myrank == 0) {
            for (int i = 0; i < n; i++) {
                a[i] = myrank*10+i;
            }
            System.out.println("myrank = " + myrank + " : a = "+ Arrays.toString(a));
        }
        // передача даних від 0
        MPI.COMM_WORLD.bcast(a, a.length, MPI.INT, 0);
        if (myrank != 0) {
            System.out.println("myrank = " + myrank + " : a = "+ Arrays.toString(a));
        }
        // завершення параленьої частини
        MPI.Finalize();
    }
}

//        myrank = 0 : a = [0, 1, 2, 3, 4]
//        myrank = 1 : a = [0, 1, 2, 3, 4]
//        myrank = 2 : a = [0, 1, 2, 3, 4]
//        myrank = 4 : a = [0, 1, 2, 3, 4]
//        myrank = 5 : a = [0, 1, 2, 3, 4]
//        myrank = 7 : a = [0, 1, 2, 3, 4]
//        myrank = 6 : a = [0, 1, 2, 3, 4]
//        myrank = 3 : a = [0, 1, 2, 3, 4]
