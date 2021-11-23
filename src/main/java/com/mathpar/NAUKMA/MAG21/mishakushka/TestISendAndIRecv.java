package com.mathpar.NAUKMA.MAG21.mishakushka;

import mpi.MPI;
import mpi.MPIException;

import java.nio.IntBuffer;
import java.util.Random;

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

/**
 * Output:
 * a[0]= 0.28666949282560616
 * a[1]= 0.9660840218880319
 * a[2]= 0.8018557701678086
 * a[3]= 0.9141152365096255
 * a[4]= 0.6411419603689779
 * Proc num 0 ?????????? ??i??????????????????
 */