package com.mathpar.NAUKMA.MAG21.riepkin;

import mpi.*;
// to run type mpirun --allow-run-as-root -np <number-of-processes> java -cp /DAP/target/classes com/mathpar/NAUKMA/examples/HelloWorldParallel
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

// Output:
// Proc num <process_number> Hello World 
