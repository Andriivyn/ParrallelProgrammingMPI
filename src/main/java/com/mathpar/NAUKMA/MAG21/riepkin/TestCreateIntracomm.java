package com.mathpar.NAUKMA.MAG21.riepkin;

import mpi.Intracomm;
import mpi.MPI;
import mpi.MPIException;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Random;

/*
Could not start the program, because dont have 5 cores :(
*/

public class TestCreateIntracomm {
    public static void main(String[] args) throws MPIException {
        MPI.Init(args);
        int myrank = MPI.COMM_WORLD.getRank();
        ArrayList s = new ArrayList();
        s.add(0);
        s.add(1);
        ArrayList s1 = new ArrayList();
        s1.add(2);
        s1.add(3);
        mpi.Group g = MPI.COMM_WORLD.getGroup().incl(new int[]{0, 1});
        mpi.Group g2 = MPI.COMM_WORLD.getGroup().incl(new int[]{2, 3});
        Intracomm COMM_NEW = MPI.COMM_WORLD.create(g);
        Intracomm COMM_NEW_1 = MPI.COMM_WORLD.create(g2);
        int n = 4;
        int[] a = new int[n];
        if (myrank == 0 || myrank == 3) {
            for (int i = 0; i < n; i++) {
                a[i] = myrank*10+i;
            }
            System.out.println("myrank = " + myrank + ": a = " + Arrays.toString(a));
        }
        if (s.contains(myrank)) COMM_NEW.bcast(a, a.length, MPI.INT, 0);
        if (s1.contains(myrank)) COMM_NEW_1.bcast(a, a.length, MPI.INT, 0);
        if (myrank != 0 && myrank != 3) System.out.println("myrank = " + myrank + ": a = " + Arrays.toString(a));
        MPI.Finalize();
    }
}
