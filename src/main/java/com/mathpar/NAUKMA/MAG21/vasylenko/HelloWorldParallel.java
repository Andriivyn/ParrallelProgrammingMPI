package com.mathpar.NAUKMA.MAG21.vasylenko;

import mpi.MPI;
import mpi.MPIException;

//mpirun --hostfile /home/vasya/hostfile -np 7 java -cp /home/vasya/DAP/target/classes com/mathpar/NAUKMA/MAG21/vasylenko/HelloWorldParallel

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

//Proc num 5 Hello World
//Proc num 6 Hello World
//Proc num 4 Hello World
//Proc num 3 Hello World
//Proc num 2 Hello World
//Proc num 1 Hello World
//Proc num 0 Hello World

