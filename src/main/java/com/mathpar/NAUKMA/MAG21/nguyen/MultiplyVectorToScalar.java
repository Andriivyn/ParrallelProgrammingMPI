package com.mathpar.NAUKMA.MAG21.nguyen;

import com.mathpar.number.*;
import mpi.MPI;
import mpi.MPIException;

import java.io.IOException;
import java.util.Random;

//Input:
//mpirun --hostfile hostfile -np 4 java -cp target/classes com/mathpar/NAUKMA/MAG21/nguyen/MultiplyVectorToScalar

//Output:
//I'm processor 1
//I'm processor 2
//I'm processor 3
//Vector B = [7, 22, 23, 12, 17, 9, 10, 12]
//rank = 1 B = [23, 12]
//rank = 2 B = [17, 9]
//rank = 3 B = [10, 12]
//send result
//send result
//send result
//B * S = [35, 110, 115, 60, 85, 45, 50, 60]

public class MultiplyVectorToScalar {
    public static void main(String[] args) throws MPIException,
            IOException, ClassNotFoundException {

        Ring ring = new Ring("Z[]");
        MPI.Init(args);
        int rank = MPI.COMM_WORLD.getRank();
        int size = MPI.COMM_WORLD.getSize();

        int ord = 8;

        Element s = NumberZ.valueOf(5);

        int k = ord / size;

        int n = ord - k * (size - 1);
        if (rank == 0) {
            int den = 10000;
            Random rnd = new Random();
            VectorS B = new VectorS(ord, den, new int[]{5},
                    rnd, ring);
            System.out.println("Vector B = " + B);

            Element[] res0 = new Element[n];
            for (int i = 0; i < n; i++) {
                res0[i] = B.V[i].multiply(s, ring);
            }

            for (int j = 1; j < size; j++) {
                Element[] v = new Element[k];
                System.arraycopy(B.V, n + (j - 1) * k, v, 0, k);
                Transport.sendObject(v, j, 100 + j);
            }

            Element[] result = new Element[ord];
            System.arraycopy(res0, 0, result, 0, n);

            for (int t = 1; t < size; t++) {
                Element[] resRank = (Element[])
                        Transport.recvObject(t, 100 + t);
                System.arraycopy(resRank, 0, result, n +
                        (t - 1) * k, resRank.length);

            }
            System.out.println("B * S = " +
                    new VectorS(result).toString(ring));
        } else {

            System.out.println("I'm processor " + rank);

            Element[] B = (Element[])
                    Transport.recvObject(0, 100 + rank);
            System.out.println("rank = " + rank +
                    " B = " + Array.toString(B));

            Element[] result = new Element[k];
            for (int j = 0; j < B.length; j++) {
                result[j] = B[j].multiply(s, ring);
            }

            Transport.sendObject(result, 0, 100 + rank);
            System.out.println("send result");
        }
        MPI.Finalize();
    }
}
