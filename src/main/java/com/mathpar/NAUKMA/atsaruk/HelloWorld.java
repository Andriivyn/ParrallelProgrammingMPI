package com.mathpar.NAUKMA.atsaruk;
 
import mpi.*;
//mpirun -np 2 java -cp /home/teacher/mpi-dap/target/classes com.mathpar.NAUKMA/atsaruk/HelloWorld
public class HelloWorld{
public static void main(String[] args)
throws MPIException { 
MPI.Init(args);
int myrank = MPI.COMM_WORLD.getRank();
System.out.println("Proc num " + myrank + " Hello World");
 
MPI.Finalize();
}
}