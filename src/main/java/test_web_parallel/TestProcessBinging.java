/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package test_web_parallel;

import com.mathpar.parallel.webCluster.engine.QueryResult;
import com.mathpar.parallel.webCluster.engine.Tools;
import java.io.BufferedReader;
import java.io.InputStreamReader;
import mpi.MPI;

/*

\\uploadToCluster(test_web_parallel.zip);

\TOTALNODES = 3;
\PROCPERNODE = 4;
\runUploadedClass(test_web_parallel.zip, test_web_parallel.TestProcessBinging);

*/
public class TestProcessBinging {
    
    // mpirun -np 2 java -cp /home/r1d1/NetBeansProjects/mathpar/target/classes test_web_parallel.TestProcessBinging
    public static void main(String[] args) throws Exception {
        MPI.Init(args);
        QueryResult queryRes = Tools.getDataFromClusterRootNode(args);
        int myrank = MPI.COMM_WORLD.getRank();
        Runtime rt = Runtime.getRuntime();
        String[] commands = {"hostname"};
        Process proc = rt.exec(commands);

        BufferedReader stdInput = new BufferedReader(new InputStreamReader(proc.getInputStream()));

        String s = null,hostname="";        
        while ((s = stdInput.readLine()) != null) {
            hostname+=s;            
        }
        long stTime=System.currentTimeMillis();
        double v=1.0;
        for (int i=0; i<2000000000; i++){
            if (i%2==0){
                v/=2;
            }
            else{
                v*=2;
            }
        }
        System.out.println("runtime on host="+hostname+" and rank="+myrank+" is: "+(System.currentTimeMillis()-stTime)+"ms; v="+v);
        
        Tools.sendFinishMessage(args);
        MPI.Finalize();
    }
}
