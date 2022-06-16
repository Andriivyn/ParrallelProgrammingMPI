package com.mathpar.NAUKMA.MAG21.riepkin;

import com.mathpar.matrix.MatrixS;
import com.mathpar.number.*;
import mpi.MPI;
import mpi.MPIException;

import java.io.IOException;
import java.util.Random;

/*
Generates random matrix 16x16 and vector with 16 elements

*/
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
	    // split matrixes and send parts to another processes
            for (int j = 1; j < size; j++) {
                Element[][] Mj=new Element[k][];
                int[][] colj=new int[k][];
                for (int s = 0; s < k; s++) {
                    Mj[s]=A.M[s+j*k]; colj[s]=A.col[s+j*k];
                }
                MatrixS Aj=new MatrixS(Mj,colj);
                Transport.sendObject(Aj, j, 1); // send matrix to process j with tag 1
                Transport.sendObject(B, j, 2); // send vector to process j with tag 2
            }
            Element[][] M0=new Element[n][];  int[][] col0=new int[n][];
            for (int s = 0; s < n; s++) {
                M0[s]=A.M[s+(size-1)*k]; col0[s]=A.col[s+(size-1)*k];
            }
            MatrixS A0=new MatrixS(M0,col0);
            // calculate main process part
            VectorS V0=A0.multiplyByColumn(B.transpose(ring), ring);
            System.arraycopy(V0.V, 0, result, (size-1)*k, n);
            // accumulate results from other processes
            for (int t = 1; t < size; t++) {
                VectorS Vt = (VectorS)Transport.recvObject(t, t);
                System.arraycopy(Vt.V, 0, result,(t - 1) * k, k);
            }
            System.out.println("RESULT: A*B=" + new VectorS(result).toString(ring));
        } else {
            MatrixS AA = (MatrixS)Transport.recvObject(0, 1 ); // receive matrix from main process with tag 1
            System.out.println(  " AA("+rank+") = " + AA.toString(ring));
            VectorS BB = (VectorS)Transport.recvObject(0, 2 ); // receive vector from main process with tag 2
            System.out.println(  " BB("+rank+") = " + BB.toString(ring));
            VectorS VV=AA.multiplyByColumn(BB.transpose(ring), ring); // calculate results
            System.out.println(  " VV("+rank+") = " + VV.toString(ring));
            Transport.sendObject(VV, 0,  rank); // send results back with tag equals to it's process number
            System.out.println("send result from("+rank+")");
        }
        MPI.Finalize();
    }
}
/*
Output:
root@b1a37b43d450:/DAP/src/main/java/com/mathpar/NAUKMA/riepkin# sh run_test.sh 4 MultiplyMatrixToVector
Matrix A = 
[[16, 31, 8,  31, 8,  6,  30, 10, 4,  7,  30, 5,  12, 22, 10, 7 ]
 [19, 18, 3,  30, 15, 13, 6,  31, 26, 31, 6,  14, 9,  11, 11, 10]
 [11, 18, 10, 29, 5,  23, 5,  3,  13, 18, 5,  24, 4,  25, 22, 6 ]
 [24, 11, 3,  15, 25, 30, 19, 24, 27, 11, 3,  22, 13, 25, 4,  29]
 [24, 7,  22, 6,  26, 30, 6,  1,  4,  22, 5,  28, 18, 7,  10, 18]
 [10, 19, 10, 22, 23, 17, 21, 1,  7,  29, 19, 27, 2,  30, 0,  7 ]
 [13, 24, 5,  18, 31, 19, 24, 29, 2,  11, 21, 8,  13, 25, 29, 0 ]
 [0,  2,  15, 12, 1,  30, 4,  21, 21, 11, 17, 20, 1,  19, 24, 1 ]
 [3,  27, 14, 24, 11, 14, 19, 14, 31, 9,  13, 20, 14, 7,  7,  11]
 [20, 16, 10, 9,  14, 29, 25, 10, 31, 15, 18, 26, 3,  17, 2,  2 ]
 [18, 10, 9,  17, 24, 7,  15, 13, 19, 8,  6,  25, 13, 10, 1,  9 ]
 [28, 21, 19, 4,  23, 6,  11, 25, 27, 31, 16, 5,  30, 14, 25, 3 ]
 [17, 27, 15, 20, 28, 22, 5,  17, 21, 0,  23, 15, 2,  24, 1,  18]
 [17, 9,  3,  22, 20, 18, 19, 25, 2,  21, 2,  12, 10, 4,  18, 27]
 [26, 8,  22, 26, 16, 20, 12, 19, 31, 9,  24, 14, 5,  11, 2,  31]
 [11, 4,  17, 22, 29, 23, 8,  3,  12, 8,  0,  8,  1,  9,  7,  3 ]]
Vector B = [3, 8, 25, 27, 12, 23, 28, 11, 7, 28, 1, 13, 13, 26, 25, 8]
 AA(1) = 
[[24, 7,  22, 6,  26, 30, 6,  1,  4,  22, 5,  28, 18, 7,  10, 18]
 [10, 19, 10, 22, 23, 17, 21, 1,  7,  29, 19, 27, 2,  30, 0,  7 ]
 [13, 24, 5,  18, 31, 19, 24, 29, 2,  11, 21, 8,  13, 25, 29, 0 ]
 [0,  2,  15, 12, 1,  30, 4,  21, 21, 11, 17, 20, 1,  19, 24, 1 ]]
 BB(1) = [3, 8, 25, 27, 12, 23, 28, 11, 7, 28, 1, 13, 13, 26, 25, 8]
 VV(1) = [3844, 4385, 4633, 3607]
send result from(1)
 AA(2) = 
[[3,  27, 14, 24, 11, 14, 19, 14, 31, 9,  13, 20, 14, 7,  7,  11]
 [20, 16, 10, 9,  14, 29, 25, 10, 31, 15, 18, 26, 3,  17, 2,  2 ]
 [18, 10, 9,  17, 24, 7,  15, 13, 19, 8,  6,  25, 13, 10, 1,  9 ]
 [28, 21, 19, 4,  23, 6,  11, 25, 27, 31, 16, 5,  30, 14, 25, 3 ]]
 BB(2) = [3, 8, 25, 27, 12, 23, 28, 11, 7, 28, 1, 13, 13, 26, 25, 8]
 VV(2) = [3732, 3866, 3044, 4373]
send result from(2)
 AA(3) = 
[[17, 27, 15, 20, 28, 22, 5,  17, 21, 0,  23, 15, 2,  24, 1,  18]
 [17, 9,  3,  22, 20, 18, 19, 25, 2,  21, 2,  12, 10, 4,  18, 27]
 [26, 8,  22, 26, 16, 20, 12, 19, 31, 9,  24, 14, 5,  11, 2,  31]
 [11, 4,  17, 22, 29, 23, 8,  3,  12, 8,  0,  8,  1,  9,  7,  3 ]]
 BB(3) = [3, 8, 25, 27, 12, 23, 28, 11, 7, 28, 1, 13, 13, 26, 25, 8]
 VV(3) = [3535, 3913, 3915, 3076]
send result from(3)
RESULT: A*B=[3844, 4385, 4633, 3607, 3732, 3866, 3044, 4373, 3535, 3913, 3915, 3076, 3535, 3913, 3915, 3076]
*/
