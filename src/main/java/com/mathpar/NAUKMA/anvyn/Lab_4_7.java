package com.mathpar.NAUKMA.anvyn;

import mpi.MPI;
import mpi.MPIException;

import java.util.Arrays;
import java.util.stream.IntStream;

public class Lab_4_7 {
    final static int INITIAL_VALUE = 4;
    final static int ROOT = 5;

    public static void main(String[] args) throws MPIException {
        MPI.Init(args);

        int np = MPI.COMM_WORLD.getSize();
        int myRank = MPI.COMM_WORLD.getRank();
        // масив кiлькості елементiв, отриманих вiд кожного процесора
        int[] recvCount = new int[np];
        recvCount[0] = 4;
        for (int i = 1; i < recvCount.length; i++) {
            recvCount[i] = recvCount[i - 1] * 2;
        }
        //змiщення, на якi можна розмiстити вхiднi елементи;
        int[] displs = new int[np];
        displs[0] = 0;
        for (int i = 1; i < displs.length; i++) {
            displs[i] = displs[i - 1] + recvCount[i-1];
        }
        //кiлькiсть елементiв, що пересилаються з поточного процесора процесору root;
        int sendCount = recvCount[myRank];

        //масив об’єктiв, якi пересилаються;
        int[] sendBuf = new int[sendCount];
        for (int i = 0; i < sendBuf.length; i++) {
            sendBuf[i] = (int)(Math.random() * 100);
        }
        //розмір масива, який повинен отримати процесор root.
        int recvBufLength = IntStream.of(recvCount).sum();
        //масив, в який процесор root отримає дані від решту процесорів
        int[] recvBuf = new int[recvBufLength];

        //надсилаємо з кожного процесора включно з root,  процесору root дані
        MPI.COMM_WORLD.gatherv(sendBuf, sendCount, MPI.INT, recvBuf, recvCount, displs, MPI.INT, ROOT);

        if (myRank == ROOT) {
            System.out.println("ROOT rank:" + myRank + "get buf size:= " + recvBuf.length);
        } else {
            System.out.println("myRank:=" + myRank + " sent buf size:= " + sendBuf.length);
        }

        MPI.Finalize();
    }
}
