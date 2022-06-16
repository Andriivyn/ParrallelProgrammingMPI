package com.mathpar.NAUKMA.MAG21.nguyen;

import mpi.MPI;
import mpi.MPIException;

import java.nio.IntBuffer;
import java.util.Random;

//Input:
//mpirun --hostfile hostfile -np 4 java -cp target/classes com/mathpar/NAUKMA/MAG21/nguyen/TestISendAndIRecv

//Output:
//proc num = 0 array sent
//proc num = 1 array received
//proc num = 2 array received
//proc num = 3 array received

public class TestISendAndIRecv {
    public static void main(String[] args) throws MPIException {

        MPI.Init(args);

        int myrank = MPI.COMM_WORLD.getRank();

        int np = MPI.COMM_WORLD.getSize();

        int n = 5;

        IntBuffer b = MPI.newIntBuffer(n);

        if (myrank == 0) {
            for (int i = 0; i < n; i++) {
                b.put(new Random().nextInt(10));
            }
            for (int i = 1; i < np; i++) {
                MPI.COMM_WORLD.iSend(b, b.capacity(), MPI.INT, i, 3000);
            }
            System.out.println("proc num = " + myrank + " array sent");

        } else {
            MPI.COMM_WORLD.recv(b, b.capacity(), MPI.INT, 0, 3000);
            System.out.println("proc num = " + myrank + " array received");

        }
        MPI.Finalize();
    }
}
