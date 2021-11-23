package com.mathpar.NAUKMA.MAG21.riepkin;

import mpi.MPI;
import mpi.MPIException;

import java.util.Arrays;
/*
Creates and fills array on each process.
Sends it to all other processes.
Input:
Proc0: [A0, A1, A2, A3]
Proc2: [B0, B1, B2, B3]
ProcN: [N0, N1, N2, N3]

Output:
Proc0: [A0, B0, C0, D0]
Proc1: [A1, B1, C1, D1]
ProcN: [AN, BN, CN, DN]

N <= number of items in array, if n < number of items in array then will be zeros at the end.
*/
public class TestAllToAll {
    public static void main(String[] args) throws MPIException {
        MPI.Init(args);
        int myrank = MPI.COMM_WORLD.getRank();
        int np = MPI.COMM_WORLD.getSize();
        int n = 4;
        int[] a = new int[n];
        for (int i = 0; i < n; i++) a[i] = myrank*10+i;
        System.out.println("myrank = " + myrank + ": a = " + Arrays.toString(a));
        int[] q = new int[n];
        MPI.COMM_WORLD.allToAll(a, 1, MPI.INT, q, 1, MPI.INT);
        System.out.println("myrank = " + myrank + ": q = " + Arrays.toString(q));
        MPI.Finalize();
    }
}
/*
Output:
root@b1a37b43d450:/DAP/src/main/java/com/mathpar/NAUKMA/riepkin# sh run_test.sh 4 TestAllToAll      
myrank = 0: a = [0, 1, 2, 3]
myrank = 0: q = [0, 10, 20, 30]
myrank = 1: a = [10, 11, 12, 13]
myrank = 1: q = [1, 11, 21, 31]
myrank = 2: a = [20, 21, 22, 23]
myrank = 2: q = [2, 12, 22, 32]
myrank = 3: a = [30, 31, 32, 33]
myrank = 3: q = [3, 13, 23, 33]
root@b1a37b43d450:/DAP/src/main/java/com/mathpar/NAUKMA/riepkin# sh run_test.sh 2 TestAllToAll
myrank = 0: a = [0, 1, 2, 3]
myrank = 1: a = [10, 11, 12, 13]
myrank = 0: q = [0, 10, 0, 0]
myrank = 1: q = [1, 11, 0, 0]
*/
