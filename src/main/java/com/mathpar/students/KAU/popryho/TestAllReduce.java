//mpirun --hostfile /home/popryho/hostfile -np 4 java -cp /home/popryho/IdeaProjects/dap/target/classes com/mathpar/students/KAU/popryho/TestAllReduce
package com.mathpar.students.KAU.popryho;

import mpi.MPI;
import mpi.MPIException;

import java.util.Arrays;

public class TestAllReduce {
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
        int[] q = new int[n];
        MPI.COMM_WORLD.allReduce(a, q, n, MPI.INT, MPI.PROD);
        for (int j = 0; j < np; j++) {
            if (myrank == j) {
                System.out.println("myrank = " + myrank + ": q = " + Arrays.toString(q));
            }
        }
        // завершення параалельної частини
        MPI.Finalize();
    }
}
