/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

package test_web_parallel;

import com.mathpar.number.*;
import com.mathpar.parallel.webCluster.engine.*;
import java.nio.Buffer;
import java.nio.ByteBuffer;
import java.util.ArrayList;
import mpi.*;
/**
 *
 * @author r1d1
 */
/*
\TOTALNODES = 3;
\PROCPERNODE = 4;
\runUploadedClass(test_web_parallel.zip, test_web_parallel.test, 12, 14);
*/
public class test {
    public static void main(String[] args) throws MPIException {   
        int tmp;
        int z=3;
        int qq=10;
        MPI.Init(args);
        QueryResult queryRes = Tools.getDataFromClusterRootNode(args);
        int myRank = MPI.COMM_WORLD.getRank();
        if (myRank == 0) {
            Object[] ar = queryRes.getData();
            System.out.println("test...");
            for (int i = 0; i < ar.length; i++) {
                System.out.println(((Element) ar[i]).intValue());
            }
        }        
        System.out.println("hello from "+myRank);
        Tools.sendFinishMessage(args);
        MPI.Finalize();
    }
}
