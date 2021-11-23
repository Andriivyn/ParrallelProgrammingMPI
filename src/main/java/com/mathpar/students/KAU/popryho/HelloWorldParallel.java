//mpirun --hostfile /home/popryho/hostfile -np 8 java -cp /home/popryho/IdeaProjects/dap/target/classes com/mathpar/students/KAU/popryho/HelloWorldParallel
package com.mathpar.students.KAU.popryho;

import mpi.MPI;
import mpi.MPIException;

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
