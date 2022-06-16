package com.mathpar.NAUKMA.MAG21.vasylenko;

import com.mathpar.matrix.MatrixS;
import com.mathpar.number.NumberZ;
import com.mathpar.number.Ring;
import com.mathpar.parallel.utils.MPITransport;
import mpi.MPI;
import mpi.MPIException;

import java.io.IOException;
import java.util.Random;

//mpirun --hostfile /home/vasya/hostfile -np 4 java -cp /home/vasya/DAP/target/classes com/mathpar/NAUKMA/MAG21/vasylenko/MatrixMul4

public class MatrixMul4 {
    public static MatrixS mmultiply(MatrixS a, MatrixS b, MatrixS c, MatrixS d, Ring ring) {
        // помножимо a на b, с на d та додамо результати
        return (a.multiply(b, ring)).add(c.multiply(d, ring), ring);
    }

    public static void main(String[] args) throws MPIException, IOException, ClassNotFoundException {
        Ring ring = new Ring("Z[]");
        //iнiцiалiзацiя MPI
        MPI.Init(args);
        // отримання номера процесора
        int rank = MPI.COMM_WORLD.getRank();
        if (rank == 0) {
            // програма виконується на нульовому процесорi
            int ord = 4;
            int den = 10000;
            // представник класу випадкового генератора
            Random rnd = new Random();
            // ord = розмiр матрицi, den = щiльнiсть
            MatrixS A = new MatrixS(ord, ord, den, new int[]{5}, rnd, NumberZ.ONE, ring);
            MatrixS B = new MatrixS(ord, ord, den, new int[]{5}, rnd, NumberZ.ONE, ring);

            MatrixS[] DD = new MatrixS[4];
            MatrixS CC = null;
            // розбиваємо матрицю A на 4 частини
            MatrixS[] AA = A.split();
            // розбиваємо матрицю B на 4 частини
            MatrixS[] BB = B.split();
            // вiдправлення вiд нульового процесора масиву Object процесору 1 з iдентифiкатором tag = 1
            MPITransport.sendObjectArray(new Object[]{AA[0], BB[1], AA[1], BB[3]}, 0, 4, 1, 1);
            // вiдправлення вiд нульового процесора масиву Object процесору 2 з iдентифiкатором tag = 2
            MPITransport.sendObjectArray(new Object[]{AA[2], BB[0], AA[3], BB[2]}, 0, 4, 2, 2);
            // вiдправлення вiд нульового процесора масиву Object// процесору 3 з iдентифiкатором tag = 3
            MPITransport.sendObjectArray(new Object[]{AA[2], BB[1], AA[3], BB[3]}, 0, 4, 3, 3);
            // залишаємо один блок нульовому процесору для оброблення
            DD[0] = (AA[0].multiply(BB[0], ring)).add(AA[1].multiply(BB[2], ring), ring);
            // отримуємо результат вiд першого процесора
            DD[1] = (MatrixS) MPITransport.recvObject(1, 1);
            System.out.println("recv 1 to 0");// отримуємо результат вiд другого процесора
            DD[2] = (MatrixS) MPITransport.recvObject(2, 2);
            System.out.println("recv 2 to 0");
            // отримуємо результат вiд третього процесора
            DD[3] = (MatrixS) MPITransport.recvObject(3, 3);
            System.out.println("recv 3 to 0");
            //процедура збору матрицi з блокiв DD[i]//(i=0,...,3)
            CC = MatrixS.join(DD);
            System.out.println("RES= " + CC.toString());
        } else
        {
            // програма виконується на процесорi з номером rank
            System.out.println("I’m processor " + rank);
            // отримуємо масив Object з блоками матриць вiд нульового процесора
            Object[] n = new Object[4];
            MPITransport.recvObjectArray(n, 0, 4, 0, rank);
            MatrixS a = (MatrixS) n[0];
            MatrixS b = (MatrixS) n[1];
            MatrixS c = (MatrixS) n[2];
            MatrixS d = (MatrixS) n[3];

            // перемножуємо та складаємо блоки матриць
            MatrixS res = mmultiply(a, b, c, d, ring);
            // надсилаємо результат обчислень вiд процесора rank нульовому процесору
            System.out.println("res = " + res);
            MPITransport.sendObject(res, 0, rank);
            // повiдомлення на консоль про те, що результат буде надiслано
            System.out.println("send result");
        }
        MPI.Finalize();
    }
}

/*
I’m processor 2
I’m processor 3
I’m processor 1
res =
[[784, 1447]
 [653, 1567]]
res =
[[2268, 1748]
 [2353, 1816]]
res =
[[839, 1301]
 [433, 733 ]]
send result
send result
send result
recv 1 to 0
recv 2 to 0
recv 3 to 0
RES=
[[1307, 1326, 839, 1301]
 [727,  754,  433, 733 ]
 [2268, 1748, 784, 1447]
 [2353, 1816, 653, 1567]]
*/