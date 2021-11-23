package com.mathpar.NAUKMA.MAG21.vasylenko;

import mpi.MPI;
import mpi.MPIException;

import java.util.Random;

// mpirun --hostfile /home/vasya/hostfile -np 2 java -cp /home/vasya/DAP/target/classes com/mathpar/NAUKMA/MAG21/vasylenko/TestDeadlock

public class TestDeadlock {
    public static void main(String[] args) throws MPIException {
        //iнiцiалiзацiя паралельної частини
        MPI.Init(args);
        int myRank = MPI.COMM_WORLD.getRank();
        int n = 1010; // 1010 - працює, а на 1011 відбувається deadlock
        int[] a = new int[n];
        int[] b = new int[n];
        if (myRank == 0) {
            Random rnd = new Random();
            for (int i = 0; i < n; i++) {
                a[i] = rnd.nextInt() % n;
            }
            MPI.COMM_WORLD.send(a, n, MPI.INT, 1, 0);
            MPI.COMM_WORLD.recv(b, n, MPI.INT, 1, 1);
            /*
            iншi iнструкцiї програми для вузла 0
            ...
            */
        }
        if (myRank == 1) {
            Random rnd = new Random();
            for (int i = 0; i < n; i++) {
                b[i] = rnd.nextInt() % n;
            }
            MPI.COMM_WORLD.send(b, n, MPI.INT, 0, 1);
            MPI.COMM_WORLD.recv(a, n, MPI.INT, 0, 0);
            /*
            iншi iнструкцiї програми для вузла 1
            ...
            */
        }
        System.out.println("myrank = " + myRank + " success!");
        MPI.Finalize();
    }
}
//        до 1011
//        myrank = 0 success!
//        myrank = 1 success!
