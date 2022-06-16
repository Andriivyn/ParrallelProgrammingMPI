package com.mathpar.NAUKMA.MAG21.riepkin;

import mpi.MPI;
import mpi.MPIException;

import java.util.Arrays;
import java.util.Random;

/*
Main process creates array with 5 elements and sends it to all other processes.
Slave processes print received array.
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
        // передача даних від 0 процесора іншим
        MPI.COMM_WORLD.bcast(a, a.length, MPI.INT, 0);
        if (myrank != 0) {
            System.out.println("myrank = " + myrank + " : a = "+ Arrays.toString(a));
        }
        // завершення параленьої частини
        MPI.Finalize();
    }
}
/*
Output:
root@b1a37b43d450:/DAP/src/main/java/com/mathpar/NAUKMA/riepkin# sh run_test.sh 4 TestBcast    
myrank = 0 : a = [0, 1, 2, 3, 4]
myrank = 1 : a = [0, 1, 2, 3, 4]
myrank = 3 : a = [0, 1, 2, 3, 4]
myrank = 2 : a = [0, 1, 2, 3, 4]
root@b1a37b43d450:/DAP/src/main/java/com/mathpar/NAUKMA/riepkin# sh run_test.sh 2 TestBcast
myrank = 0 : a = [0, 1, 2, 3, 4]
myrank = 1 : a = [0, 1, 2, 3, 4]
*/
