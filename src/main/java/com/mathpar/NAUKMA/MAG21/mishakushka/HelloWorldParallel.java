package com.mathpar.NAUKMA.MAG21.mishakushka;

import mpi.*;

public class HelloWorldParallel {
    public static void main(String[] args) throws MPIException {
        //iнiцiалiзацiя паралельної частини
        MPI.Init(args);
        //визначення номера процесора
        int myrank = MPI.COMM_WORLD.getRank();
        System.out.println("Proc num " + myrank + " Hello World");
        //завершення паралельної частини
        MPI.Finalize();
    }
}

/**
 * Output: Proc num 0 Hello World
 */
