/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package test_web_parallel;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.OutputStream;
import java.io.OutputStreamWriter;
import java.util.ArrayList;
import java.util.logging.Level;
import java.util.logging.Logger;
import static test_web_parallel.MakeGetNodes.stringIsCorrectJobName;

/**
 *
 * @author derbist
 */
public class LaunchTaskOnCluster {
    public static String TASKNAME;
    public static String JOBFOLDER;
    public static String CLASSPATH;
    public static int TOTALNODES;
    public static int PROCPERNODE;
    public static int WALLTIME;

    public LaunchTaskOnCluster(String filePath) throws FileNotFoundException, IOException {
        File fileSettings = new File(filePath);
        BufferedReader brSet = new BufferedReader(new InputStreamReader(new FileInputStream(fileSettings)));
        String temp = null;
        String[] allSettings = new String[6];
        for (int i = 0; i < allSettings.length; i++) {
            temp = brSet.readLine();
            allSettings[i] = subStringForSettings(temp);
        }
        brSet.close();
        TASKNAME = allSettings[0];
        JOBFOLDER = allSettings[1];
        CLASSPATH = allSettings[2];
        TOTALNODES = Integer.parseInt(allSettings[3]);
        PROCPERNODE = Integer.parseInt(allSettings[4]);
        WALLTIME = Integer.parseInt(allSettings[5]);
        
    }

    public static String subStringForNode(String str) {
        int end = 0;
        for (int i = 0; i < str.length(); i++) {
            if (str.charAt(i) == ':') {
                end = i;
                break;
            }
        }
        str = str.substring(0, end);
        return str;
    }

    public static String subStringForSettings(String str) {
        int end = 0;
        for (int i = 0; i < str.length(); i++) {
            if (str.indexOf("=") == i) {
                end = ++i;
                break;
            }
        }
        str = str.substring(end);
        return str;
    }
    
    public static String[] createList(File listFile) {
        File exportFiles[] = listFile.listFiles();
        String[] names = new String[exportFiles.length];
        for (int i = 0; i < names.length; i++) {
            names[i] = exportFiles[i].getName();
        }
        return names;
    }

    
    /**
     * метод, генерирующий и исполняющий команду mpirun
     * @param np - общее количество ядер
     * @param npernode - количество ядер на одном узле
     * @param jobDir - путь к рабочей директории
     * @param hosts - путь к файлу с именами узлов
     * @param classPath - путь к классу параллельной программы
     * @return
     * @throws Exception 
     */
    public static Process makeMpirunCommand(String np, String npernode,String jobDir, String hosts, String classPath) throws Exception{
        String mpirun = "/gpfs/NETHOME/tamgu1/10p_soft/nodeBuild2/openmpi1.10/bin/mpirun";
//        String mpirun ="/gpfs/NETHOME/tamgu1/sirosh/openmpi/bin/mpirun";
//        String mpilib ="/gpfs/NETHOME/tamgu1/sirosh/openmpi/lib";
        int memory = 64000/Integer.parseInt(npernode);
        String java = "/gpfs/NETHOME/tamgu1/10p_soft/jdk-13.0.2/bin/java";
        ArrayList<String> mpirunStart = new ArrayList<String>();
        mpirunStart.add(mpirun);
        mpirunStart.add("-np");
        mpirunStart.add(np);
        mpirunStart.add("-npernode");
        mpirunStart.add(npernode);
        mpirunStart.add("--hostfile");
        mpirunStart.add(hosts);
        mpirunStart.add("-mca");
        mpirunStart.add("btl");
        mpirunStart.add("self,sm,openib");
        mpirunStart.add(java);
        mpirunStart.add("-Xmx"+memory+"m");
//        mpirunStart.add("-Djava.library.path="+mpilib);
        mpirunStart.add("-cp");
        int k = 0;
        for(int i =0; i < classPath.length();i++){
            if(classPath.charAt(i) == ' '){
                mpirunStart.add(classPath.substring(k, i));
                k = i+1;
            }
            if(i == classPath.length()-1){
               mpirunStart.add(classPath.substring(k, i+1)); 
            }
        }
        ProcessBuilder pBuilder = new ProcessBuilder(mpirunStart);
        File dir=new File(jobDir);
        pBuilder.directory(dir);
        Process process=pBuilder.start();
        return process;
    }
    /**
     * метод, запускающий параллельную задачу на кластере mvs10p
     * @param filePath - путь к файлу с настроками
     * @throws FileNotFoundException
     * @throws IOException
     * @throws Exception 
     */
    public static void launchTask(String filePath) throws FileNotFoundException, IOException, Exception{
       LaunchTaskOnCluster l = new LaunchTaskOnCluster(filePath);
        TASKNAME = MakeGetNodes.makeGetNodesCommand(String.valueOf(TOTALNODES), String.valueOf(WALLTIME), JOBFOLDER, TASKNAME);
        System.out.println("Path for Task Name = " +JOBFOLDER + "/" + TASKNAME);
        File outFromTask = new File(JOBFOLDER + "/" + TASKNAME + "/output");
        File hosts = new File(JOBFOLDER + "/" + TASKNAME + "/hosts");
        long waiting = 0;
        while (waiting == 0 && outFromTask.canRead()) {
            try {
                Thread.sleep(1000);
            } catch (InterruptedException ex) {
                Logger.getLogger(LaunchTaskOnCluster.class.getName()).log(Level.SEVERE, null, ex);
            }
            waiting = outFromTask.length();
        }
        BufferedReader br = new BufferedReader(new InputStreamReader(new FileInputStream(outFromTask)));
        BufferedWriter bw = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(hosts)));
        String temp2 = null;
        while ((temp2 = br.readLine()) != null) {

            if (temp2.charAt(0) == 'n') {
                bw.write(subStringForNode(temp2));
                bw.newLine();
            }

        }
        bw.close();
        br.close();
        int np = TOTALNODES*PROCPERNODE;
        Process child = makeMpirunCommand(String.valueOf(np), String.valueOf(PROCPERNODE), JOBFOLDER, hosts.getAbsolutePath(), CLASSPATH);
        InputStream stdout = child.getInputStream();
        InputStream stderr = child.getErrorStream();
        String line = null;
        BufferedReader outReader = new BufferedReader(new InputStreamReader(stdout));
        BufferedWriter out = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(JOBFOLDER + "/" + TASKNAME + "/out")));
        while ((line = outReader.readLine()) != null) {
            out.write(line);
            out.newLine();
        }
        outReader.close();
        out.close();
        BufferedReader errReader = new BufferedReader(new InputStreamReader(stderr));
        BufferedWriter err = new BufferedWriter(new OutputStreamWriter(new FileOutputStream(JOBFOLDER + "/" + TASKNAME + "/err")));
        while ((line = errReader.readLine()) != null) {
            err.write(line);
            err.newLine();
        }
        errReader.close();
        err.close();
     
    }

    public static void main(String[] args) throws IOException, Exception {
        launchTask(args[0]);
    }

}
