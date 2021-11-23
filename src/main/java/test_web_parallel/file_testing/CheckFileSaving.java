/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package test_web_parallel.file_testing;


import com.mathpar.parallel.webCluster.engine.QueryResult;
import com.mathpar.parallel.webCluster.engine.Tools;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileOutputStream;
import java.io.InputStreamReader;
import java.io.ObjectOutputStream;
import java.io.Serializable;
import java.util.Random;
import mpi.MPI;

/**
 \\uploadToCluster(test_web_parallel.zip);

\TOTALNODES = 1;
\PROCPERNODE = 1;
\runUploadedClass(test_web_parallel.zip, test_web_parallel.file_testing.CheckFileSaving);

 */
class myClass implements Serializable{
    myClass(){        
        i=(new Random()).nextInt();
        d=(new Random()).nextDouble();
    }
    int i;
    double d;
}

public class CheckFileSaving {
    // mpirun -np 2 java -cp /home/r1d1/NetBeansProjects/mathpar/target/classes test_web_parallel.file_testing.CheckFileSaving
    public static void main(String[] args) throws Exception{
        MPI.Init(args);
        QueryResult queryRes=Tools.getDataFromClusterRootNode(args);
        
        int myRank = MPI.COMM_WORLD.getRank();      
        Runtime rt = Runtime.getRuntime();
        String filePathPrefix="/tmp/testFileFolder"+String.valueOf(myRank)+"/";                        
        String[] delDir = {"rm", "-R", filePathPrefix};
        Process chmod = Runtime.getRuntime().exec(delDir);
        chmod.waitFor();
        System.out.println("tmp dir removed on "+myRank);
        String[] mkDir = {"mkdir", filePathPrefix};
        Process chmod1 = Runtime.getRuntime().exec(mkDir);
        chmod1.waitFor();
        System.out.println("tmp dir created on"+myRank);
        
        if (myRank==0){
            System.out.println("mpijavac version checking...");
            String[] checkJavaCVersion = {"mpijavac","-version"};        
            Process javaCVersion = rt.exec(checkJavaCVersion);
            javaCVersion.waitFor();
            BufferedReader stdInput = new BufferedReader(new InputStreamReader(javaCVersion.getErrorStream()));       
            String s = null;        
            while ((s = stdInput.readLine()) != null) {
                System.out.println(s+" from "+myRank);
            }
            System.out.flush();
        }
        MPI.COMM_WORLD.barrier();
        
        for (int j=0; j<10; j++){
            File f = new File(filePathPrefix+"f"+String.valueOf(j)+".out");        
            FileOutputStream fout = new FileOutputStream(f);
            ObjectOutputStream oos = new ObjectOutputStream(fout);   
            System.out.println("File writing begin...");
            for (int i=0; i<1000000; i++){
                myClass tmp=new myClass();
                oos.writeObject(tmp);
            }            
            oos.close();
            fout.close();
            System.out.println("File writing is successfull on "+myRank+" for iteration numb="+j);
            System.out.flush();
        }
        
        Tools.sendFinishMessage(args);
        MPI.Finalize();
    }
}