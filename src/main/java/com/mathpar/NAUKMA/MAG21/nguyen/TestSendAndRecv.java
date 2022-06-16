package com.mathpar.NAUKMA.MAG21.nguyen;

import mpi.MPI;
import mpi.MPIException;

import java.util.Random;

//Input:
//mpirun --hostfile hostfile -np 2 java -cp target/classes com/mathpar/NAUKMA/MAG21/nguyen/TestSendAndRecv

//Output:
//a[0]= 0.8262410578190684
//a[1]= 0.7777073823713857
//a[2]= 0.7706600087939518
//a[3]= 0.27389418568651525
//a[4]= 0.09670167548994268
//Proc num 0 array is sent
//
//a[0]= 0.8262410578190684
//a[1]= 0.7777073823713857
//a[2]= 0.7706600087939518
//a[3]= 0.27389418568651525
//a[4]= 0.09670167548994268
//Proc num 1 array is received

public class TestSendAndRecv {
    public static void main(String[] args) throws MPIException {
        //iнiцiалiзацiя MPI
        MPI.Init(args);
        //визначення номера процесора
        int myrank = MPI.COMM_WORLD.getRank();
        //визначення кiлькостi процесорiв у групi
        int np = MPI.COMM_WORLD.getSize();
        //вхiдний параметр - розмiр масиву
        int n = 5;
        double[] a = new double[n];
        //синхронiзацiя процесорiв
        // MPI.COMM_WORLD.barrier();
        // якщо процесор з номером 0
        if (myrank == 0) {
            for (int i = 0; i < n; i++) {
                a[i] = (new Random()).nextDouble();
                System.out.println("a[" + i + "]= " + a[i]);
            }
            //передання 0-процесором елементiв усiм iншим процесорам у групi
            for (int i = 1; i < np; i++) {
                MPI.COMM_WORLD.send(a, n, MPI.DOUBLE, i, 3000);
            }
            System.out.println("Proc num " + myrank + " array is sent" + "\n");
        } else {
            //приймання i-м процесором повiдомлення вiд
            // процесора з номером 0 та тегом 3000.
            MPI.COMM_WORLD.recv(a, n, MPI.DOUBLE, 0, 3000);
            for (int i = 0; i < n; i++) {
                System.out.println("a[" + i + "]= " + a[i]);
            }
            System.out.println("Proc num " + myrank + " array is received" + "\n");
        }

        // завершення паралельної частини
        MPI.Finalize();
    }
}
