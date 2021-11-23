/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package test_web_parallel.file_testing;

import com.mathpar.matrix.file.dense.FileMatrixD;
import com.mathpar.parallel.webCluster.engine.QueryResult;
import com.mathpar.parallel.webCluster.engine.Tools;
import java.io.File;
import mpi.*;

/**
 \\uploadToCluster(test.zip);

\TOTALNODES = 1;
\PROCPERNODE = 1;
\runUploadedClass(test.zip, test.CheckFileSavingWithFMatrix);

 */
public class CheckFileSavingWithFMatrix {
    // mpirun -np 2 java -cp /home/r1d1/NetBeansProjects/mathpar/target/classes test.CheckFileSavingWithFMatrix
    public static void main(String[] args) throws Exception{
        MPI.Init(args);
        QueryResult queryRes=Tools.getDataFromClusterRootNode(args);
        
        int myRank = MPI.COMM_WORLD.getRank();
        Runtime rt = Runtime.getRuntime();
        String filePathPrefix="/tmp/testMatrixFolder"+String.valueOf(myRank)+"/";                        
        String[] delDir = {"rm", "-R", filePathPrefix};
        Process chmod = Runtime.getRuntime().exec(delDir);
        chmod.waitFor();
        System.out.println("tmp dir removed on "+myRank);
        String[] mkDir = {"mkdir", filePathPrefix};
        Process chmod1 = Runtime.getRuntime().exec(mkDir);
        chmod1.waitFor();
        System.out.println("tmp dir created on"+myRank);
        
        File fA = new File(filePathPrefix+"mA");
        int depth=0,mSize=128,nbits=3;
        FileMatrixD a = new FileMatrixD(fA, depth, mSize, mSize, nbits);
        MPI.COMM_WORLD.barrier();
        System.out.println("File writing is successfull on "+myRank);
        System.out.flush();
        
        Tools.sendFinishMessage(args);
        MPI.Finalize();
    }
}