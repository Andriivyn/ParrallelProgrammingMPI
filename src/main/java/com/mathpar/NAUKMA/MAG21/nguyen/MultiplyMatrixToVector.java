package com.mathpar.NAUKMA.MAG21.nguyen;

import com.mathpar.matrix.MatrixS;
import com.mathpar.number.Element;
import com.mathpar.number.NumberZ;
import com.mathpar.number.Ring;
import com.mathpar.number.VectorS;
import mpi.MPI;
import mpi.MPIException;

import java.io.IOException;
import java.util.Random;

//Input:
//mpirun --hostfile hostfile -np 4 java -cp target/classes com/mathpar/NAUKMA/MAG21/nguyen/MultiplyMatrixToVector

//Output:
//Matrix A =
//[[28, 8,  28, 2,  20, 3,  13, 23, 29, 22, 19, 2,  27, 23, 31, 1 ]
//[13, 22, 9,  8,  17, 6,  9,  13, 27, 13, 22, 21, 17, 9,  29, 31]
//[10, 3,  5,  15, 7,  24, 2,  30, 16, 5,  23, 0,  19, 18, 16, 21]
//[25, 3,  31, 20, 10, 27, 22, 21, 4,  17, 5,  30, 22, 23, 22, 16]
//[2,  8,  6,  18, 19, 18, 5,  30, 28, 4,  16, 10, 18, 26, 14, 25]
//[12, 6,  9,  3,  25, 1,  14, 3,  0,  18, 7,  31, 30, 29, 23, 29]
//[29, 16, 23, 3,  17, 10, 15, 15, 24, 12, 24, 25, 17, 5,  16, 4 ]
//[29, 24, 30, 13, 20, 16, 8,  13, 10, 7,  18, 17, 2,  19, 25, 28]
//[18, 18, 15, 26, 24, 23, 24, 9,  28, 21, 5,  3,  8,  11, 16, 14]
//[4,  16, 24, 8,  21, 1,  30, 13, 23, 9,  31, 10, 25, 16, 20, 7 ]
//[25, 7,  8,  16, 1,  18, 24, 4,  1,  22, 16, 1,  13, 5,  4,  16]
//[25, 8,  13, 13, 19, 12, 6,  25, 1,  15, 6,  4,  27, 13, 0,  20]
//[11, 6,  24, 27, 2,  18, 7,  24, 30, 8,  17, 12, 26, 10, 5,  0 ]
//[21, 21, 22, 20, 6,  8,  25, 11, 31, 19, 31, 21, 1,  1,  20, 29]
//[30, 8,  24, 20, 22, 16, 28, 29, 11, 2,  0,  26, 6,  20, 11, 23]
//[15, 18, 11, 21, 7,  14, 15, 26, 9,  10, 2,  23, 1,  5,  9,  25]]
//Vector B = [13, 5, 0, 19, 18, 4, 16, 12, 0, 29, 20, 5, 10, 24, 14, 13]
//AA(1) =
//[[2,  8,  6,  18, 19, 18, 5,  30, 28, 4,  16, 10, 18, 26, 14, 25]
//[12, 6,  9,  3,  25, 1,  14, 3,  0,  18, 7,  31, 30, 29, 23, 29]
//[29, 16, 23, 3,  17, 10, 15, 15, 24, 12, 24, 25, 17, 5,  16, 4 ]
//[29, 24, 30, 13, 20, 16, 8,  13, 10, 7,  18, 17, 2,  19, 25, 28]]
//BB(1) = [13, 5, 0, 19, 18, 4, 16, 12, 0, 29, 20, 5, 10, 24, 14, 13]
//VV(1) = [3073, 3469, 2799, 3290]
//send result from(1)
//AA(2) =
//[[18, 18, 15, 26, 24, 23, 24, 9,  28, 21, 5,  3,  8,  11, 16, 14]
//[4,  16, 24, 8,  21, 1,  30, 13, 23, 9,  31, 10, 25, 16, 20, 7 ]
//[25, 7,  8,  16, 1,  18, 24, 4,  1,  22, 16, 1,  13, 5,  4,  16]
//[25, 8,  13, 13, 19, 12, 6,  25, 1,  15, 6,  4,  27, 13, 0,  20]]
//BB(2) = [13, 5, 0, 19, 18, 4, 16, 12, 0, 29, 20, 5, 10, 24, 14, 13]
//VV(2) = [3308, 3238, 2663, 2815]
//send result from(2)
//AA(3) =
//[[11, 6,  24, 27, 2,  18, 7,  24, 30, 8,  17, 12, 26, 10, 5,  0 ]
//[21, 21, 22, 20, 6,  8,  25, 11, 31, 19, 31, 21, 1,  1,  20, 29]
//[30, 8,  24, 20, 22, 16, 28, 29, 11, 2,  0,  26, 6,  20, 11, 23]
//[15, 18, 11, 21, 7,  14, 15, 26, 9,  10, 2,  23, 1,  5,  9,  25]]
//BB(3) = [13, 5, 0, 19, 18, 4, 16, 12, 0, 29, 20, 5, 10, 24, 14, 13]
//VV(3) = [2396, 3397, 3247, 2444]
//send result from(3)
//RESULT: A*B=[3073, 3469, 2799, 3290, 3308, 3238, 2663, 2815, 2396, 3397, 3247, 2444, 2396, 3397, 3247, 2444]

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
                VectorS Vt = (VectorS) Transport.recvObject(t, t);
                System.arraycopy(Vt.V, 0, result,(t - 1) * k, k);
            }
            System.out.println("RESULT: A*B=" + new VectorS(result).toString(ring));
        } else {
            MatrixS AA = (MatrixS) Transport.recvObject(0, 1 );
            System.out.println(  " AA("+rank+") = " + AA.toString(ring));
            VectorS BB = (VectorS) Transport.recvObject(0, 2 );
            System.out.println(  " BB("+rank+") = " + BB.toString(ring));
            VectorS VV=AA.multiplyByColumn(BB.transpose(ring), ring);
            System.out.println(  " VV("+rank+") = " + VV.toString(ring));
            Transport.sendObject(VV, 0,  rank);
            System.out.println("send result from("+rank+")");
        }
        MPI.Finalize();
    }
}
