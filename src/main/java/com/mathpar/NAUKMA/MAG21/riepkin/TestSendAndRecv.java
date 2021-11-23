package com.mathpar.NAUKMA.MAG21.riepkin;

import java.util.Random;
import mpi.*;

/*
Blocking send from main, blocking receive from slave.
*/
public class TestSendAndRecv {
    public static void main(String[] args) throws MPIException {
        //iнiцiалiзацiя MPI
        MPI.Init(args);
        //визначення номера процесора
        int myrank = MPI.COMM_WORLD.getRank();
        //визначення кiлькостi процесорiв у групi
        int np = MPI.COMM_WORLD.getSize();
        //вхiдний параметр - розмiр масиву
        int n = 5;
        double[] a = new double[n];
        //синхронiзацiя процесорiв
        // MPI.COMM_WORLD.barrier();
        // якщо процесор з номером 0
        if (myrank == 0) {
            for (int i = 0; i < n; i++) {
                a[i] = (new Random()).nextDouble();
                System.out.println("a[" + i + "]= " + a[i]);
            }
            //передання 0-процесором елементiв усiм iншим процесорам у групi
            for (int i = 1; i < np; i++) {
                MPI.COMM_WORLD.send(a, n, MPI.DOUBLE, i, 3000);
            }
            System.out.println("Proc num " + myrank + " масив вiдправлено" + "\n");
        } else {
            //приймання i-м процесором повiдомлення вiд
            // процесора з номером 0 та тегом 3000.
            MPI.COMM_WORLD.recv(a, n, MPI.DOUBLE, 0, 3000);
            for (int i = 0; i < n; i++) {
                System.out.println("a[" + i + "]= " + a[i]);
            }
            System.out.println("Proc num " + myrank + " масив прийнято" + "\n");
        }

        // завершення паралельної частини
        MPI.Finalize();
    }
}
/*
Output:
root@b1a37b43d450:/DAP/src/main/java/com/mathpar/NAUKMA/riepkin# sh run_test.sh 4 TestSendAndRecv
a[0]= 0.8028898223565718
a[1]= 0.8291859172759469
a[2]= 0.2351314137759727
a[3]= 0.5550165126932414
a[4]= 0.5320329459154761
Proc num 0 ????? ?i?????????

a[0]= 0.8028898223565718
a[1]= 0.8291859172759469
a[2]= 0.2351314137759727
a[3]= 0.5550165126932414
a[4]= 0.5320329459154761
a[0]= 0.8028898223565718
a[1]= 0.8291859172759469
a[2]= 0.2351314137759727
a[3]= 0.5550165126932414
a[4]= 0.5320329459154761
Proc num 3 ????? ????????

a[0]= 0.8028898223565718
a[1]= 0.8291859172759469
a[2]= 0.2351314137759727
a[3]= 0.5550165126932414
a[4]= 0.5320329459154761
Proc num 2 ????? ????????

Proc num 1 ????? ????????

root@b1a37b43d450:/DAP/src/main/java/com/mathpar/NAUKMA/riepkin# sh run_test.sh 2 TestSendAndRecv
a[0]= 0.16612561346286725
a[1]= 0.05104009278890875
a[2]= 0.909419266173667
a[3]= 0.9181670876027017
a[4]= 0.5963890669332914
Proc num 0 ????? ?i?????????

a[0]= 0.16612561346286725
a[1]= 0.05104009278890875
a[2]= 0.909419266173667
a[3]= 0.9181670876027017
a[4]= 0.5963890669332914
Proc num 1 ????? ????????
*/
