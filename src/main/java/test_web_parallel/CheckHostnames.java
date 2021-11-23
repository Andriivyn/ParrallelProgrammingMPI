package test_web_parallel;
import com.mathpar.parallel.webCluster.engine.QueryResult;
import com.mathpar.parallel.webCluster.engine.Tools;
import java.io.BufferedReader;
import java.io.InputStreamReader;
import mpi.*;

/*

\\uploadToCluster(test.zip);

\TOTALNODES = 3;
\PROCPERNODE = 4;
\runUploadedClass(test_web_parallel.zip, test_web_parallel.CheckHostnames);

*/
public class CheckHostnames {
    public static void main(String[] args) throws Exception{
        MPI.Init(args);
        QueryResult queryRes=Tools.getDataFromClusterRootNode(args);
        int myrank = MPI.COMM_WORLD.getRank();
        Runtime rt = Runtime.getRuntime();
        String[] commands = {"hostname"};
        Process proc = rt.exec(commands);

        BufferedReader stdInput = new BufferedReader(new InputStreamReader(proc.getInputStream()));
       
        String s = null;
        int myRank=MPI.COMM_WORLD.getRank();
        while ((s = stdInput.readLine()) != null) {
            System.out.println(s+" from "+myRank);
        }
        Tools.sendFinishMessage(args);
        MPI.Finalize();
    }
    
}
