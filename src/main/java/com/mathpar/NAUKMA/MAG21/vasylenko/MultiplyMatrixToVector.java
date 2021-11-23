package com.mathpar.NAUKMA.MAG21.vasylenko;

import com.mathpar.NAUKMA.examples.Transport;
import com.mathpar.matrix.MatrixS;
import com.mathpar.number.Element;
import com.mathpar.number.NumberZ;
import com.mathpar.number.Ring;
import com.mathpar.number.VectorS;
import mpi.MPI;
import mpi.MPIException;

import java.io.IOException;
import java.util.Random;

//mpirun --hostfile /home/vasya/hostfile -np 4 java -cp /home/vasya/DAP/target/classes com/mathpar/NAUKMA/MAG21/vasylenko/MultiplyMatrixToVector

public class MultiplyMatrixToVector {
    public static void main(String[] args) throws MPIException,
            IOException, ClassNotFoundException {

        Ring ring = new Ring("Z[]");
        MPI.Init(args);
        int rank = MPI.COMM_WORLD.getRank();
        int size = MPI.COMM_WORLD.getSize();
        int ord = 16; // matrix size (ord х ord)
        int k = ord / size; // the number of rows for each processor
        int n = ord - k * (size - 1); // the last portion will be n (less than k) elements )
        if (rank == 0) {
            int den = 10000; // for density=100%
            Random rnd = new Random();
            MatrixS A = new MatrixS(ord, ord, den, new int[]{5},rnd, NumberZ.ONE, ring);
            System.out.println("Matrix A = " + A.toString());
            VectorS B = new VectorS(ord, den, new int[]{5},rnd, ring);
            System.out.println("Vector B = " + B.toString());
            Element[] result = new Element[ord]; // For resulting vector
            for (int j = 1; j < size; j++) {
                Element[][] Mj=new Element[k][];
                int[][] colj=new int[k][];
                for (int s = 0; s < k; s++) {
                    Mj[s]=A.M[s+j*k]; colj[s]=A.col[s+j*k];
                }
                MatrixS Aj=new MatrixS(Mj,colj);
                Transport.sendObject(Aj, j, 1);
                Transport.sendObject(B, j, 2);
            }
            Element[][] M0=new Element[n][];  int[][] col0=new int[n][];
            for (int s = 0; s < n; s++) {
                M0[s]=A.M[s+(size-1)*k]; col0[s]=A.col[s+(size-1)*k];
            }
            MatrixS A0=new MatrixS(M0,col0);
            VectorS V0=A0.multiplyByColumn(B.transpose(ring), ring);
            System.arraycopy(V0.V, 0, result, (size-1)*k, n);
            for (int t = 1; t < size; t++) {
                VectorS Vt = (VectorS)Transport.recvObject(t, t);
                System.arraycopy(Vt.V, 0, result,(t - 1) * k, k);
            }
            System.out.println("RESULT: A*B=" + new VectorS(result).toString(ring));
        } else {
            MatrixS AA = (MatrixS)Transport.recvObject(0, 1 );
            System.out.println(  " AA("+rank+") = " + AA.toString(ring));
            VectorS BB = (VectorS)Transport.recvObject(0, 2 );
            System.out.println(  " BB("+rank+") = " + BB.toString(ring));
            VectorS VV=AA.multiplyByColumn(BB.transpose(ring), ring);
            System.out.println(  " VV("+rank+") = " + VV.toString(ring));
            Transport.sendObject(VV, 0,  rank);
            System.out.println("send result from("+rank+")");
        }
        MPI.Finalize();
    }
}

/*
Matrix A =
[[9,  10, 11, 28, 29, 29, 3,  28, 7,  15, 28, 24, 16, 10, 0,  15]
 [6,  16, 3,  23, 15, 5,  22, 10, 8,  1,  8,  12, 2,  21, 4,  0 ]
 [29, 27, 8,  11, 6,  3,  3,  8,  14, 26, 17, 14, 12, 29, 26, 31]
 [21, 31, 9,  18, 25, 23, 9,  3,  1,  14, 6,  25, 6,  7,  1,  24]
 [20, 15, 21, 8,  0,  6,  9,  5,  6,  15, 21, 9,  7,  20, 25, 22]
 [2,  29, 26, 8,  17, 19, 0,  26, 6,  15, 13, 23, 6,  9,  31, 4 ]
 [5,  29, 9,  11, 3,  0,  14, 11, 6,  3,  30, 31, 25, 7,  20, 12]
 [5,  2,  16, 15, 13, 16, 17, 5,  15, 29, 22, 28, 25, 18, 22, 26]
 [8,  30, 29, 13, 7,  7,  9,  18, 7,  16, 5,  9,  27, 13, 5,  25]
 [22, 5,  1,  0,  8,  7,  21, 30, 21, 16, 25, 2,  23, 27, 30, 20]
 [22, 9,  10, 4,  5,  19, 20, 19, 0,  14, 2,  0,  4,  30, 2,  21]
 [14, 31, 6,  19, 23, 15, 4,  28, 15, 31, 14, 31, 24, 27, 14, 9 ]
 [6,  23, 6,  9,  20, 29, 5,  4,  9,  21, 13, 26, 12, 19, 3,  25]
 [4,  7,  23, 0,  20, 1,  22, 22, 21, 2,  13, 13, 21, 11, 0,  16]
 [16, 0,  9,  20, 13, 28, 22, 30, 26, 27, 29, 17, 11, 4,  2,  27]
 [16, 22, 12, 0,  22, 20, 3,  6,  0,  19, 13, 2,  5,  2,  10, 3 ]]
Vector B = [12, 27, 13, 30, 6, 10, 29, 25, 12, 27, 11, 28, 26, 8, 22, 20]
 AA(1) =
[[20, 15, 21, 8,  0,  6,  9,  5,  6,  15, 21, 9,  7,  20, 25, 22]
 [2,  29, 26, 8,  17, 19, 0,  26, 6,  15, 13, 23, 6,  9,  31, 4 ]
 [5,  29, 9,  11, 3,  0,  14, 11, 6,  3,  30, 31, 25, 7,  20, 12]
 [5,  2,  16, 15, 13, 16, 17, 5,  15, 29, 22, 28, 25, 18, 22, 26]]
 BB(1) = [12, 27, 13, 30, 6, 10, 29, 25, 12, 27, 11, 28, 26, 8, 22, 20]
 VV(1) = [3896, 4581, 4726, 5415]
send result from(1)
 AA(2) =
[[8,  30, 29, 13, 7,  7,  9,  18, 7,  16, 5,  9,  27, 13, 5,  25]
 [22, 5,  1,  0,  8,  7,  21, 30, 21, 16, 25, 2,  23, 27, 30, 20]
 [22, 9,  10, 4,  5,  19, 20, 19, 0,  14, 2,  0,  4,  30, 2,  21]
 [14, 31, 6,  19, 23, 15, 4,  28, 15, 31, 14, 31, 24, 27, 14, 9 ]]
 BB(2) = [12, 27, 13, 30, 6, 10, 29, 25, 12, 27, 11, 28, 26, 8, 22, 20]
 VV(2) = [4735, 4778, 3240, 6124]
send result from(2)
 AA(3) =
[[6,  23, 6,  9,  20, 29, 5,  4,  9,  21, 13, 26, 12, 19, 3,  25]
 [4,  7,  23, 0,  20, 1,  22, 22, 21, 2,  13, 13, 21, 11, 0,  16]
 [16, 0,  9,  20, 13, 28, 22, 30, 26, 27, 29, 17, 11, 4,  2,  27]
 [16, 22, 12, 0,  22, 20, 3,  6,  0,  19, 13, 2,  5,  2,  10, 3 ]]
 BB(3) = [12, 27, 13, 30, 6, 10, 29, 25, 12, 27, 11, 28, 26, 8, 22, 20]
 VV(3) = [4272, 3621, 5393, 2649]
send result from(3)
RESULT: A*B=[3896, 4581, 4726, 5415, 4735, 4778, 3240, 6124, 4272, 3621, 5393, 2649, 4272, 3621, 5393, 2649]

 */