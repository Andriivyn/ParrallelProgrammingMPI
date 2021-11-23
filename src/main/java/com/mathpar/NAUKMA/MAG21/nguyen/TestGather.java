package com.mathpar.NAUKMA.MAG21.nguyen;

import mpi.MPI;

import java.util.Arrays;

//Input:
//mpirun --hostfile hostfile -np 4 java -cp target/classes com/mathpar/NAUKMA/MAG21/nguyen/TestGather

//Output:
//myrank = 1 : a = [10, 11, 12, 13, 14]
//myrank = 2 : a = [20, 21, 22, 23, 24]
//myrank = 3 : a = [30, 31, 32, 33, 34]
//myrank = 0 : a = [0, 1, 2, 3, 4]
//myrank = 3 : q = [0, 1, 2, 3, 4, 10, 11, 12, 13, 14, 20, 21, 22, 23, 24, 30, 31, 32, 33, 34]

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
        MPI.COMM_WORLD.gather(a, n, MPI.INT, q, n, MPI.INT, np-1);
        if(myrank == np-1) {
            System.out.println("myrank = " + myrank + " : q = "+ Arrays.toString(q));
        }
        // завершення параленьої частини
        MPI.Finalize();
    }
}
