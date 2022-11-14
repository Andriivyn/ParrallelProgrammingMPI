package com.mathpar.NAUKMA.anvyn;

import mpi.MPI;
import mpi.MPIException;

import java.util.Arrays;
import java.util.stream.IntStream;

public class Lab_4_8 {
    final static int ROOT = 3;

    public static void main(String[] args) throws MPIException {
        MPI.Init(args);
        int np = MPI.COMM_WORLD.getSize();
        int myRank = MPI.COMM_WORLD.getRank();

        //масив кiлькості елементiв, що пересилаються процесорам
        int[] sendCount = new int[np];
        sendCount[0] = 1;
        for (int i = 1; i < sendCount.length; i++) {
            sendCount[i] = sendCount[i - 1] * 3;
        }
        //масив змiщень
        int[] displs = new int[np];
        displs[0] = 0;
        for (int i = 1; i < displs.length; i++) {
            displs[i] = displs[i - 1] + sendCount[i - 1];
        }
        //розмір буферу надсилаючого процесора (сума всіх елементів  sendCount)
        int sendBufLength = IntStream.of(sendCount).sum();

        //масив об’єктiв, якi пересилаються;
        int[] sendBuf = new int[sendBufLength];
        if (myRank == ROOT) {
            for (int i = 0; i < sendBuf.length; i++) {
                sendBuf[i] = (int) (Math.random() * 100);
            }
        }
        //розмір масива, який повинен отримати процесор myRank від root.
        int recvCount = sendCount[myRank];
        //буфер приймання повiдомлення (масив об’єктiв);
        int[] recvBuf = new int[recvCount];

        //надсилаємо з кожного процесора включно з root,  процесору root дані
        MPI.COMM_WORLD.scatterv(sendBuf, sendCount, displs, MPI.INT, recvBuf, recvCount, MPI.INT, ROOT);

        if (myRank == ROOT) {
            System.out.println("ROOT rank:" + myRank + "sent buf size:= " + sendBuf.length);
        }
        System.out.println("myRank:=" + myRank + " get buf size:= " + recvBuf.length);

        MPI.Finalize();
    }
}
