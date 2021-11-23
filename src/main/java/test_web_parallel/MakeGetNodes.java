/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package test_web_parallel;

import java.io.File;
import java.util.ArrayList;
import java.util.Scanner;

/**
 *
 * @author r1d1
 */
public class MakeGetNodes {

    
    
    /**
     * Checking for string is true job name in format 
     * "[job].[index]"
     * @param s -- string for checking
     * @return true\false
     */
    public static boolean stringIsCorrectJobName(String s, String taskName){
        //true task name can't be short 5 symbols
        if (s.length()<5){
            return false;
        }
        //first and last symbols is "
        if (s.charAt(0)!='"' || s.charAt(s.length()-1)!='"'){
            return false;
        }        
        int dotPos=-1;
        for (int i=1; i<s.length()-1; i++){
            //if there are more than one dots
            if (s.charAt(i)=='.' && dotPos!=-1){
                return false;
            }
            if (s.charAt(i)=='.'){
                dotPos=i;
            }
        }
        //If there no dots or after dot are not digits
        if (dotPos==-1 || dotPos==s.length()-2){
            return false;
        }
        //only digits after dot
        for (int i=dotPos+1; i<s.length()-1; i++){
            if (s.charAt(i)<'0' || s.charAt(i)>'9'){
                return false;
            }
        }
        if (dotPos-1 != taskName.length()){
            return false;
        }
        for (int i=1; i<dotPos; i++){
            if (s.charAt(i)!=taskName.charAt(i-1)){
                return false;
            }
        }
        return true;
    }
    
    /**
     * 
     * @param np -- total nodes
     * @param walltime -- job time
     * @param jobDir -- job directory, that will contains jobFolder with SUPPZ files
     * @return task name
     */
    public static String makeGetNodesCommand(String np, String walltime, String jobDir, String taskName) throws Exception{
        ProcessBuilder pBuilder=new ProcessBuilder("getnodes","-np",np,"-maxtime",walltime,taskName);
        File dir=new File(jobDir);
        pBuilder.directory(dir);
        Process process=pBuilder.start();
        process.waitFor();
        Scanner sc=new Scanner(process.getInputStream());        
        sc.useDelimiter(" ");
        while (sc.hasNext()){
            String curToken=sc.next();            
            if (stringIsCorrectJobName(curToken,taskName)){                
                return curToken.substring(1, curToken.length()-1);
            }
        }
        throw new Exception("Can't execute getnodes command!");        
    }
    
    public static void main(String[] args) throws Exception{     
        System.out.println(makeGetNodesCommand("1", "1", "/gpfs/NETHOME/tamgu1/gennadi", "test"));                
    }
    
}
