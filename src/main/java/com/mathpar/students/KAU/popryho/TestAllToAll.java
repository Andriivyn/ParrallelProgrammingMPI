//mpirun --hostfile /home/popryho/hostfile -np 4 java -cp /home/popryho/IdeaProjects/dap/target/classes com/mathpar/students/KAU/popryho/TestAllToAll
package com.mathpar.students.KAU.popryho;

import mpi.MPI;
import mpi.MPIException;

import java.util.Arrays;

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