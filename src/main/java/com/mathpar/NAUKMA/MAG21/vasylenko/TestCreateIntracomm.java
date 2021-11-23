package com.mathpar.NAUKMA.MAG21.vasylenko;

import mpi.Intracomm;
import mpi.MPI;
import mpi.MPIException;

import java.util.ArrayList;
import java.util.Arrays;

// mpirun --hostfile /home/vasya/hostfile -np 5 java -cp /home/vasya/DAP/target/classes com/mathpar/NAUKMA/MAG21/vasylenko/TestCreateIntracomm

public class TestCreateIntracomm {
    public static void main(String[] args) throws MPIException {
        MPI.Init(args);
        int myrank = MPI.COMM_WORLD.getRank();
        ArrayList s = new ArrayList();
        s.add(0);
        s.add(1);
        s.add(2);
        ArrayList s1 = new ArrayList();
        s1.add(3);
        s1.add(4);
        // отримання двух груп процесорів
        mpi.Group g = MPI.COMM_WORLD.getGroup().incl(new int[]{0, 1, 2});
        mpi.Group g2 = MPI.COMM_WORLD.getGroup().incl(new int[]{3, 4});
        // створення комунікаторів для кожної групи
        Intracomm COMM_NEW = MPI.COMM_WORLD.create(g);
        Intracomm COMM_NEW_1 = MPI.COMM_WORLD.create(g2);
        int n = 5;
        int[] a = new int[n];
        if (myrank == 0 || myrank == 3) {
            for (int i = 0; i < n; i++) {
                a[i] = myrank*10+i;
            }
            System.out.println("myrank = " + myrank + ": a = " + Arrays.toString(a));
        }
        // процесори 0 та 3 поширюють (bcast) свої масиви між процесорами відповідних
        // груп в яких знаходяться
        if (s.contains(myrank)) COMM_NEW.bcast(a, a.length, MPI.INT, 0);
        if (s1.contains(myrank)) COMM_NEW_1.bcast(a, a.length, MPI.INT, 0);
        if (myrank != 0 && myrank != 3) System.out.println("myrank = " + myrank + ": a = " + Arrays.toString(a));
        MPI.Finalize();
    }
}

//        myrank = 0: a = [0, 1, 2, 3, 4]
//        myrank = 2: a = [0, 1, 2, 3, 4]
//        myrank = 1: a = [0, 1, 2, 3, 4]
//        myrank = 3: a = [30, 31, 32, 33, 34]
//        myrank = 4: a = [30, 31, 32, 33, 34]

