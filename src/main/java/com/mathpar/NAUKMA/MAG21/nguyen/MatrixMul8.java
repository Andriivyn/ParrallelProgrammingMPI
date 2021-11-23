package com.mathpar.NAUKMA.MAG21.nguyen;

import com.mathpar.matrix.MatrixS;
import com.mathpar.number.NumberZ;
import com.mathpar.number.Ring;
import com.mathpar.parallel.utils.MPITransport;
import mpi.MPI;
import mpi.MPIException;

import java.io.IOException;
import java.util.Random;

//Input:
//mpirun --hostfile hostfile -np 8 java -cp target/classes com/mathpar/NAUKMA/MAG21/nguyen/MatrixMul8

//Output:
//I'm processor 5
//I'm processor 4
//I'm processor 1
//I'm processor 3
//I'm processor 7
//I'm processor 2
//I'm processor 6
//A =
//[[25, 10, 13, 1,  5,  4,  9,  6 ]
//[5,  8,  29, 21, 31, 3,  10, 25]
//[24, 24, 7,  21, 5,  30, 2,  10]
//[20, 29, 20, 11, 27, 4,  11, 31]
//[2,  13, 18, 24, 12, 20, 3,  18]
//[12, 20, 15, 15, 5,  3,  12, 17]
//[2,  16, 1,  24, 23, 27, 14, 30]
//[23, 20, 31, 16, 24, 24, 18, 14]]
//B =
//[[3,  8,  7,  5,  0,  11, 7,  11]
//[1,  6,  24, 28, 14, 25, 0,  20]
//[17, 28, 23, 13, 2,  27, 18, 26]
//[21, 14, 12, 4,  2,  29, 24, 3 ]
//[22, 19, 4,  0,  24, 1,  16, 9 ]
//[30, 22, 15, 21, 16, 28, 29, 26]
//[29, 21, 2,  12, 29, 22, 5,  24]
//[6,  3,  20, 0,  19, 12, 18, 24]]
//RES =
//[[854,  1028, 944,  770,  727,  1292, 782,  1325]
//[2169, 2134, 1835, 893,  1769, 2282, 2144, 2229]
//[1784, 1653, 1831, 1621, 1240, 2671, 1938, 2102]
//[1879, 1973, 2238, 1432, 2088, 2557, 1925, 2708]
//[1888, 1719, 1742, 1160, 1303, 2383, 2025, 1954]
//[1276, 1310, 1518, 1082, 1179, 2029, 1247, 1786]
//[2445, 1891, 1834, 1302, 2234, 2592, 2369, 2405]
//[2806, 2800, 2318, 1862, 2122, 3314, 2525, 3115]]


class MatrixMul8 {
    public static void main(String[] args)throws MPIException, IOException,ClassNotFoundException
    {
        Ring ring = new Ring("Z[]");
        MPI.Init(args);
        int rank = MPI.COMM_WORLD.getRank();
        if (rank == 0) {
            int ord = 8;
            int den = 10000;
            Random rnd = new Random();
            MatrixS A = new MatrixS(ord, ord, den, new int[] {5}, rnd, NumberZ.ONE, ring);
            System.out.println("A = " + A);
            MatrixS B = new MatrixS(ord, ord, den, new int[] {5}, rnd, NumberZ.ONE, ring);
            System.out.println("B = " + B);
            MatrixS D = null;
            MatrixS[] AA = A.split();
            MatrixS[] BB = B.split();
            int tag = 0;
            MPITransport.sendObjectArray(new Object[] {AA[1], BB[2]},0,2, 1, tag);
            MPITransport.sendObjectArray(new Object[] {AA[0], BB[1]},0,2, 2, tag);
            MPITransport.sendObjectArray(new Object[] {AA[1], BB[3]},0,2, 3, tag);
            MPITransport.sendObjectArray(new Object[] {AA[2], BB[0]},0,2, 4, tag);
            MPITransport.sendObjectArray(new Object[] {AA[3], BB[2]},0,2, 5, tag);
            MPITransport.sendObjectArray(new Object[] {AA[2], BB[1]},0,2, 6, tag);
            MPITransport.sendObjectArray(new Object[] {AA[3], BB[3]},0,2, 7, tag);
            MatrixS[] DD = new MatrixS[4];
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

