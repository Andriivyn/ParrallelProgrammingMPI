package com.mathpar.NAUKMA.MAG21.vasylenko;

import com.mathpar.matrix.MatrixS;
import com.mathpar.number.NumberZ;
import com.mathpar.number.Ring;
import com.mathpar.parallel.utils.MPITransport;
import mpi.MPI;
import mpi.MPIException;

import java.io.IOException;
import java.util.Random;

//mpirun --hostfile /home/vasya/hostfile -np 8 java -cp /home/vasya/DAP/target/classes com/mathpar/NAUKMA/MAG21/vasylenko/MatrixMul8

class MatrixMul8 {
    public static void main(String[] args)throws MPIException, IOException,ClassNotFoundException
    {
        Ring ring = new Ring("Z[]");
        MPI.Init(args);
        int rank = MPI.COMM_WORLD.getRank();
        if (rank == 0) {
            // програма виконується на нульовому процесорi
            int ord = 8;
            int den = 10000;
            Random rnd = new Random();
            // ord = розмiр матрицi, den = щiльнiсть
            MatrixS A = new MatrixS(ord, ord, den, new int[] {5}, rnd, NumberZ.ONE, ring);
            System.out.println("A = " + A);
            MatrixS B = new MatrixS(ord, ord, den, new int[] {5}, rnd, NumberZ.ONE, ring);
            System.out.println("B = " + B);
            MatrixS D = null;
            // розбиваємо матрицю A на 4 частини
            MatrixS[] AA = A.split();
            // розбиваємо матрицю B на 4 частини
            MatrixS[] BB = B.split();
            int tag = 0;
            //розподіляємо частини матриць між процесорами
            MPITransport.sendObjectArray(new Object[] {AA[1], BB[2]},0,2, 1, tag);
            MPITransport.sendObjectArray(new Object[] {AA[0], BB[1]},0,2, 2, tag);
            MPITransport.sendObjectArray(new Object[] {AA[1], BB[3]},0,2, 3, tag);
            MPITransport.sendObjectArray(new Object[] {AA[2], BB[0]},0,2, 4, tag);
            MPITransport.sendObjectArray(new Object[] {AA[3], BB[2]},0,2, 5, tag);
            MPITransport.sendObjectArray(new Object[] {AA[2], BB[1]},0,2, 6, tag);
            MPITransport.sendObjectArray(new Object[] {AA[3], BB[3]},0,2, 7, tag);
            MatrixS[] DD = new MatrixS[4];
            //збираємо результати
            DD[0] = (AA[0].multiply(BB[0], ring)).add((MatrixS) MPITransport.recvObject(1, 3),ring);
            DD[1] = (MatrixS) MPITransport.recvObject(2, 3);
            DD[2] = (MatrixS) MPITransport.recvObject(4, 3);
            DD[3] = (MatrixS) MPITransport.recvObject(6, 3);
            D = MatrixS.join(DD);
            System.out.println("RES = " + D.toString());
        }
        else
        {
            System.out.println("I'm processor " + rank);
            Object[] b = new Object[2];
            MPITransport.recvObjectArray(b,0,2,0, 0);
            MatrixS[] a = new MatrixS[b.length];
            for (int i = 0; i < b.length; i++)
                a[i] = (MatrixS) b[i];
            MatrixS res = a[0].multiply(a[1], ring);
            //здійснюємо додавання на парних процесорах та відправляємо результат на 0
            if (rank % 2 == 0)
            {
                MatrixS p = res.add((MatrixS) MPITransport.recvObject(rank + 1, 3), ring);
                MPITransport.sendObject(p, 0, 3);
            }
            else
            {
                MPITransport.sendObject(res, rank - 1, 3);
            }
        }
        MPI.Finalize();
    }
}
/*
I'm processor 5
I'm processor 7
I'm processor 2
I'm processor 1
I'm processor 6
I'm processor 4
I'm processor 3
A =
[[30, 24, 17, 22, 19, 24, 19, 28]
 [15, 29, 14, 10, 31, 19, 4,  28]
 [3,  8,  28, 29, 24, 18, 6,  3 ]
 [6,  3,  13, 5,  11, 11, 28, 0 ]
 [6,  14, 19, 9,  8,  6,  12, 20]
 [24, 27, 9,  6,  6,  0,  1,  22]
 [19, 17, 10, 10, 30, 11, 8,  16]
 [18, 11, 22, 26, 14, 10, 20, 31]]
B =
[[20, 24, 18, 14, 22, 23, 29, 21]
 [29, 27, 7,  18, 28, 31, 27, 24]
 [20, 10, 0,  5,  9,  3,  13, 2 ]
 [5,  17, 14, 27, 11, 23, 22, 0 ]
 [31, 14, 6,  29, 19, 4,  28, 22]
 [21, 21, 0,  6,  23, 8,  10, 24]
 [31, 12, 23, 5,  16, 18, 14, 16]
 [24, 13, 13, 0,  26, 9,  25, 29]]
RES =
[[4100, 3274, 1931, 2321, 3672, 2853, 3961, 3350]
 [3627, 2698, 1255, 2105, 3196, 2116, 3434, 3053]
 [2377, 1886, 837,  1943, 1905, 1443, 2316, 1454]
 [1932, 1161, 909,  863,  1298, 1021, 1344, 1178]
 [2177, 1507, 916,  1002, 1796, 1312, 1949, 1592]
 [2218, 1879, 1050, 1208, 2133, 1794, 2406, 1956]
 [2916, 2140, 1173, 1868, 2461, 1720, 2822, 2343]
 [3257, 2440, 1712, 1828, 2810, 2194, 3224, 2453]]

 */