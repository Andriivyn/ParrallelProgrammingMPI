package com.mathpar.NAUKMA.MAG21.vasylenko;

import mpi.MPI;
import mpi.MPIException;

import java.nio.IntBuffer;
import java.util.Random;

// mpirun --hostfile /home/vasya/hostfile -np 7 java -cp /home/vasya/DAP/target/classes com/mathpar/NAUKMA/MAG21/vasylenko/TestISendAndIRecv

public class TestISendAndIRecv {
    public static void main(String[] args) throws MPIException {
        //iнiцiалiзацiя MPI
        MPI.Init(args);
        //визначення номера процесора
        int myrank = MPI.COMM_WORLD.getRank();
        //визначення кiлькостi процесорiв у групi
        int np = MPI.COMM_WORLD.getSize();
        //вхiдний параметр - розмiр масиву
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

//    proc num = 0 array sent
//    proc num = 6 array received
//    proc num = 5 array received
//    proc num = 3 array received
//    proc num = 2 array received
//    proc num = 4 array received
//    proc num = 1 array received