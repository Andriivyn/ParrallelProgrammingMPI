package com.mathpar.NAUKMA.MAG21.vasylenko;

import mpi.MPI;

import java.util.Arrays;

// mpirun --hostfile /home/vasya/hostfile -np 4 java -cp /home/vasya/DAP/target/classes com/mathpar/NAUKMA/MAG21/vasylenko/TestGather

public class TestGather {
    public static void main(String[] args) throws Exception{
        // ініціалізація MPI
        MPI.Init(args);
        // визначення номера процесора
        int myrank = MPI.COMM_WORLD.getRank();
        // визначення числа процесорів в групі
        int np = MPI.COMM_WORLD.getSize();
        int n = 5;
        int[] a = new int[n];
        for(int i = 0; i < n; i++) a[i] = myrank*10+i;
        System.out.println("myrank = " + myrank + " : a = "+ Arrays.toString(a));
        int[] q = new int[n*np];
        // останній процесор збирає масив з елементів усіх масивів,
        // отриманих з усіх процесорів
        MPI.COMM_WORLD.gather(a, n, MPI.INT, q, n, MPI.INT, np-1);
        if(myrank == np-1) {
            System.out.println("myrank = " + myrank + " : q = "+ Arrays.toString(q));
        }
        // завершення параленьої частини
        MPI.Finalize();
    }
}

//        myrank = 2 : a = [20, 21, 22, 23, 24]
//        myrank = 0 : a = [0, 1, 2, 3, 4]
//        myrank = 1 : a = [10, 11, 12, 13, 14]
//        myrank = 3 : a = [30, 31, 32, 33, 34]
//        myrank = 3 : q = [0, 1, 2, 3, 4, 10, 11, 12, 13, 14, 20, 21, 22, 23, 24, 30, 31, 32, 33, 34]
