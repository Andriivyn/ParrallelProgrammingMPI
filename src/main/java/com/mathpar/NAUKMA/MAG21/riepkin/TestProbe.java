package com.mathpar.NAUKMA.MAG21.riepkin;

import mpi.MPI;
import mpi.MPIException;
import mpi.Status;
/*
Send array in blocking way to slave process with tag 1
Slave receives data in blocking manner.
*/
public class TestProbe {
    public static void main(String[] args) throws MPIException, InterruptedException {
        MPI.Init(args);
        int n = 5;
        int rank = MPI.COMM_WORLD.getRank();
        int size = MPI.COMM_WORLD.getSize();
        if (rank == 0) {
            int[] array = new int[n];
            for (int i = 1; i < size; i++) {
                MPI.COMM_WORLD.send(array, n, MPI.INT, i, 1);
                System.out.println("rank = " + rank + " відправлено до " + i);
            }
        } else {
            Status st = null;
            while (st == null) {
                st = MPI.COMM_WORLD.probe(0, 1);
            }
            System.out.println("st.getCount(MPI.INT) = " + st.getCount(MPI.INT));
            System.out.println("st.getTag() = " + st.getTag());
            System.out.println("st.getSource() = " + st.getSource());
            int[] array = new int[n];
            MPI.COMM_WORLD.recv(array, n, MPI.INT, 0, 1);
            System.out.println("rank = " + rank + " отримано");
        } MPI.Finalize();
    }
}

/*
Output:
rank = 0 ??????????? ?? 1
rank = 0 ??????????? ?? 2
rank = 0 ??????????? ?? 3st.getCount(MPI.INT) = 5
st.getTag() = 1
st.getSource() = 0
rank = 3 ????????
st.getCount(MPI.INT) = 5
st.getTag() = 1
st.getSource() = 0

rank = 1 ????????
st.getCount(MPI.INT) = 5
st.getTag() = 1
st.getSource() = 0
rank = 2 ????????
root@b1a37b43d450:/DAP/src/main/java/com/mathpar/NAUKMA/riepkin# sh run_test.sh 2 TestProbe
st.getCount(MPI.INT) = 5
st.getTag() = 1
st.getSource() = 0
rank = 0 ??????????? ?? 1
rank = 1 ????????
*/
