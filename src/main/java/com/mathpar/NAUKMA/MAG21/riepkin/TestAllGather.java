package com.mathpar.NAUKMA.MAG21.riepkin;

import mpi.MPI;

import java.util.Arrays;
/*
On each process creates array with 2 elements and fills it with (num_of_process*10+number_of_element)
Sends this array with to elements to all other processes.
Accumulates data sent from all processes in one array and prints it out.
*/
public class TestAllGather {
    public static void main(String[] args)throws Exception {
        // ініціалізація MPI
        MPI.Init(args);
        // визначення номера процесора
        int myrank = MPI.COMM_WORLD.getRank();
        // визначення числа процесорів в групі
        int np = MPI.COMM_WORLD.getSize();
        int n = 2;
        int[] a = new int[n];
        // fills array
        for(int i = 0; i < n; i++) {
            a[i] = myrank*10+i;
        }
        System.out.println("myrank = " + myrank + " : a = "+ Arrays.toString(a));
        // create array for accumulation
        int[] q = new int[n * np];
        // send own data and accumulate everything in "q"
        MPI.COMM_WORLD.allGather(a, n, MPI.INT, q, n, MPI.INT);
        System.out.println("myrank = " + myrank + " : q = "+ Arrays.toString(q));
        // завершення паралельної частини
        MPI.Finalize();
    }
}
/*
Output:
root@b1a37b43d450:/DAP/src/main/java/com/mathpar/NAUKMA/riepkin# sh run_test.sh 4 TestAllGather         
myrank = 0 : a = [0, 1]myrank = 1 : a = [10, 11]
myrank = 3 : a = [30, 31]

myrank = 2 : a = [20, 21]
myrank = 3 : q = [0, 1, 10, 11, 20, 21, 30, 31]
myrank = 2 : q = [0, 1, 10, 11, 20, 21, 30, 31]
myrank = 1 : q = [0, 1, 10, 11, 20, 21, 30, 31]
myrank = 0 : q = [0, 1, 10, 11, 20, 21, 30, 31]
root@b1a37b43d450:/DAP/src/main/java/com/mathpar/NAUKMA/riepkin# sh run_test.sh 2 TestAllGather
myrank = 0 : a = [0, 1]
myrank = 1 : a = [10, 11]
myrank = 0 : q = [0, 1, 10, 11]
myrank = 1 : q = [0, 1, 10, 11]
*/
