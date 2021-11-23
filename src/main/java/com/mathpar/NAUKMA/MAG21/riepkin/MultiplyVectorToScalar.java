package com.mathpar.NAUKMA.MAG21.riepkin;

import com.mathpar.number.*;
import mpi.MPI;
import mpi.MPIException;

import java.io.IOException;
import java.util.Random;

/*
Generates vector with 8 elements

*/

public class MultiplyVectorToScalar {
    public static void main(String[] args) throws MPIException,
            IOException, ClassNotFoundException {

        Ring ring = new Ring("Z[]");
        MPI.Init(args);
        int rank = MPI.COMM_WORLD.getRank();
        int size = MPI.COMM_WORLD.getSize();

        int ord = 8;

        Element s = NumberZ.valueOf(5);

        int k = ord / size; // portion of data per process

        int n = ord - k * (size - 1); // remains of data for the last process
        if (rank == 0) {
            int den = 10000;
            Random rnd = new Random();
            VectorS B = new VectorS(ord, den, new int[]{5},
                    rnd, ring);
            System.out.println("Vector B = " + B);

            // calculates own portion of data
            Element[] res0 = new Element[n];
            for (int i = 0; i < n; i++) {
                res0[i] = B.V[i].multiply(s, ring);
            }

            // send parts of data to other processes
            for (int j = 1; j < size; j++) {
                Element[] v = new Element[k];
                System.arraycopy(B.V, n + (j - 1) * k, v, 0, k);
                Transport.sendObject(v, j, 100 + j); // send part of data to process j with tag 100 + proccess num
            }

            Element[] result = new Element[ord];
            System.arraycopy(res0, 0, result, 0, n);

            // accumulate results from another processes
            for (int t = 1; t < size; t++) {
                Element[] resRank = (Element[])
                        Transport.recvObject(t, 100 + t);
                System.arraycopy(resRank, 0, result, n +
                        (t - 1) * k, resRank.length);

            }
            System.out.println("B * S = " +
                    new VectorS(result).toString(ring));
        } else {

            System.out.println("I’m processor " + rank);

            Element[] B = (Element[])
                    Transport.recvObject(0, 100 + rank); // receive data for calculation
            System.out.println("rank = " + rank +
                    " B = " + Array.toString(B));

            Element[] result = new Element[k];
            for (int j = 0; j < B.length; j++) {
                result[j] = B[j].multiply(s, ring); // calculates results
            }

            Transport.sendObject(result, 0, 100 + rank); // sends data back to the main process
            System.out.println("send result");
        }
        MPI.Finalize();
    }
}
/*
Output:
I?m processor 2
Vector B = [22, 25, 10, 3, 2, 4, 0, 5]
I?m processor 1
I?m processor 3
rank = 2 B = [2, 4]
rank = 1 B = [10, 3]
rank = 3 B = [0, 5]
send result
send result
send result
B * S = [110, 125, 50, 15, 10, 20, 0, 25]
*/
