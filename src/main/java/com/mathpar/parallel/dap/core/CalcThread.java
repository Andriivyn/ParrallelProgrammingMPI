/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templatesamin.outputData
 * and open the template in the editor.
 */
package com.mathpar.parallel.dap.core;

import com.mathpar.log.MpiLogger;
import com.mathpar.number.Array;
import com.mathpar.number.Element;
import com.mathpar.number.Ring;

import java.util.*;

import liquibase.pro.packaged.L;
import mpi.MPI;
import mpi.MPIException;

/**
 * @author alla
 * @param
 */
public class CalcThread implements Runnable {

    private final static MpiLogger LOGGER = MpiLogger.getLogger(CalcThread.class);

    public Thread thread;
    Ring ring;

    private ArrayList<Amin> pine;
    ArrayList<Drop>[] vokzal;
    //список готових результатов для отправки родительским процессорам(отправляет диспетчер)
    volatile ArrayList<Drop> aerodromeResults;
    ArrayList<Drop> ownTrack;

    Element[] result;
    Drop currentDrop;
    static int myRank;

    volatile boolean finish;
    volatile boolean flToExit;
    volatile boolean IamFree;
    boolean changeTrack;
    long currentMemory;

    static long counterCycle = 0;

    long calcWorkTime = 0;
    long calcWaitTime = 0;

    public CalcThread(ArrayList<Amin> p,
                      ArrayList<Drop> ownTr, Ring ring) throws MPIException {
        thread = new Thread(this, "CalcThread");
        thread.setPriority(1);
        flToExit = false;
        //isEmptyVokzal = true;
        finish = false;
        pine = p;
        ownTrack = ownTr;
        this.ring = ring;
        aerodromeResults = new ArrayList<>();
        vokzal = new ArrayList[21];
        myRank = MPI.COMM_WORLD.getRank();
        IamFree = false;
        changeTrack = false;
        //LeavesStack = new Stack<>();

        for (int i = 0; i < vokzal.length; i++) {
            vokzal[i] = new ArrayList<Drop>();
        }

        thread.start();
    }

    public void DoneThread() {
        flToExit = true;
    }

    public void clear() {
        pine.clear();
        Arrays.stream(vokzal).forEach(ArrayList::clear);
        aerodromeResults.clear();
        ownTrack.clear();
    }


    public void putDropInTrack(Drop drop) {
        synchronized (ownTrack) {
            ownTrack.add(drop);
            if (drop.recNum < DispThread.trackLevel)
                DispThread.trackLevel = drop.recNum;
        }

    }
    public void putDropInVokzal(Drop drop) {
        //LOGGER.trace("put drop in vokzal rec " + drop.recNum);
      //  LOGGER.info("put drop in vokzal num = " + drop.number + ", id = " + drop.dropId + ", amin = " + drop.aminId);
        vokzal[drop.recNum].add(drop);
        if (drop.numberOfDaughterProc == -2) {
            drop.numberOfDaughterProc = -1;
        }

        if (drop.recNum > DispThread.myLevelH || DispThread.myLevelH == 20) {
            DispThread.myLevelH = drop.recNum;
        }

        if (DispThread.myLevel == 20 || DispThread.myLevel > DispThread.myLevelH) {
            DispThread.myLevel = DispThread.myLevelH;
        }
    }

    public void writeResultsToAmin(Drop drop) {
        //LOGGER.info("writeResultsToAmin + ");
        int aminId = drop.aminId;
        int dropId = drop.dropId;


        Amin amin = pine.get(aminId);
        //LOGGER.info("pine size ="+ pine.size());
        Drop aminDrop = Drop.getDropObject(amin.type, amin.config);
        aminDrop.key = amin.key;

        for (int i = 0; i < aminDrop.arcs[dropId + 1].length; i += 3) {

            int number = aminDrop.arcs[dropId + 1][i];
            int from = aminDrop.arcs[dropId + 1][i + 1];
            int to = aminDrop.arcs[dropId + 1][i + 2];

            if (aminDrop.arcs[number].length != 0) {
                Drop dependantDrop = amin.branch.get(number - 1);
                synchronized (dependantDrop) {

                    dependantDrop.inData[to] = drop.outData[from];
                    //LOGGER.info("dependantDrop = " + dependantDrop.type+ " to " +to+ " from "+from);
                    if (dependantDrop.hasFullInputData()){
                       // LOGGER.info("putDropInVokzal");
                        putDropInVokzal(dependantDrop);}
                }
            } else {
                //LOGGER.info("resultForOutFunction");
                amin.resultForOutFunction[to] = drop.outData[from];
               // LOGGER.info("amin = " + amin);
                if (amin.hasFullOutput()) {
                    //LOGGER.info("putResultsToAminOutput");
                    putResultsToAminOutput(amin);
                }
            }
        }
    }


    private void writeResultsAfterInpFunc(Drop drop, Amin curAmin, Element[] resInputFunc) {

        for (int i = 0; i < drop.arcs[0].length; i += 3) {
            int numOfDependantDrop = drop.arcs[0][i];
            int from = drop.arcs[0][i + 1];
            int to = drop.arcs[0][i + 2];
            Drop dependantDrop = curAmin.branch.get(numOfDependantDrop - 1);

           // LOGGER.info("dependantDrop = " + dependantDrop.type+ " to " +to+ " from "+from);
          //  LOGGER.info("dependantDrop.inData =  " + dependantDrop.inData);
            dependantDrop.inData[to] = resInputFunc[from];
           // LOGGER.info("dependantDrop.inData 2!!!="    + dependantDrop.inData);
           // LOGGER.info("dependantDrop.length = " + dependantDrop.inputDataLength );
            if (dependantDrop.hasFullInputData()) {
             //   LOGGER.info("putDropInVokzal");

                putDropInVokzal(dependantDrop);
                }
            }
    }

    private void addToAerodromeResults(Drop dropRes) {
       // LOGGER.info("addToAerodromeResults");
        synchronized (aerodromeResults) {
            // LOGGER.warn("put amin num = " + amin.aminIdInPine);
            aerodromeResults.add(dropRes);
        }
    }

    private void putResultsToAminOutput(Amin amin) {
      // LOGGER.info("putResultsToAminOutput");
        Drop drop = (Drop.getDropObject(amin.type, amin.config));
        drop.key = amin.key;
        drop.inData = amin.inputData;
        //if(amin.parentAmin!=-1)
        //    drop.outData = pine.get(amin.parentAmin).branch.get(amin.parentDrop).outData;
        amin.outputData = drop.outputFunction(amin.resultForOutFunction, ring);

        if (amin.parentAmin == -1 && myRank == 0 && Array.isEmptyArray(vokzal)) {
            //LOGGER.info("go finish");
            Drop resDrop = Drop.getDropObject(amin.type, amin.config);
            resDrop.outData = amin.outputData;
            finishWholeTask(resDrop);
        } else if (amin.parentProc != myRank) {
            Drop resultAmine = Drop.doNewDrop(amin.type, amin.key, amin.config, amin.aminIdInPine, amin.parentDrop,
                    amin.parentProc, amin.recNumb, amin.inputData);
            resultAmine.outData = amin.outputData;
            resultAmine.setNumbOfMyAmine(amin.parentAmin);
            addToAerodromeResults(resultAmine);
        } else {
            Drop dr = pine.get(amin.parentAmin).branch.get(amin.parentDrop);
            dr.outData = amin.outputData;
            putDropInTrack(dr);

            int aminIndex = pine.indexOf(amin);
            if (aminIndex != -1)
                pine.set(aminIndex, null);

            currentMemory = Runtime.getRuntime().totalMemory() - Runtime.getRuntime().freeMemory();
            if (currentMemory > DispThread.usedMemory)
                DispThread.usedMemory = currentMemory;
        }

    }

    private void finishWholeTask(Drop resDrop) {
        result = resDrop.recentCalc(ring);
        finish = true;
    }

    synchronized public Drop getTask(int thread) {
      // LOGGER.info("GET Task = "+thread);
        boolean empty = false;
     //   LOGGER.trace("mylevel " + DispThread.myLevel + " myLevelH = " + DispThread.myLevelH);

       // if(vokzal[DispThread.myLevelH].size() == 0) changeMyLevelH();
       // if(vokzal[DispThread.myLevel].size() == 0)changeMyLevel();
    //    LOGGER.trace("after mylevel " + DispThread.myLevel + " myLevelH = " + DispThread.myLevelH);
        if (DispThread.isEmptyVokzal()) {
      //     LOGGER.trace("Vokzal is empty");
            empty = true;
        }

        Drop drop = null;
        ArrayList<Drop> list;
        boolean getFromTrack = false;

        if (thread == 0) {
            if (ownTrack.size() != 0) {
                // for (int i = 0; i < ownTrack.size(); i++) {
                drop = ownTrack.get(0);
                deleteFromTrack(drop);
                // getFromTrack = true;
                //    break;
                // }
            } else if (/*!getFromTrack &&*/ !empty) {

             //   LOGGER.trace("after change DispThread.myLevelH = " + DispThread.myLevelH+ "empty vokzal " + Tools.isEmptyArray(vokzal));
                list = vokzal[DispThread.myLevelH];
                drop = list.get(0);
                list.remove(0);
                if (vokzal[DispThread.myLevelH].size() == 0) {
                    changeMyLevelH();
                }
                if (drop.numberOfDaughterProc == -1) {
                    drop.numberOfDaughterProc = myRank;
                }
            }
        } else if (thread == 1 && !empty) {
            list = vokzal[DispThread.myLevel];
            LOGGER.trace(String.format("vokzal level = " + DispThread.myLevel + DispThread.myLevelH
                    +" empty = " +Array.isEmptyArray(vokzal)) + "list size = "+ list.size());
            if(list.size() == 0)  { changeMyLevel();  list = vokzal[DispThread.myLevel];}
            if(DispThread.isEmptyVokzal()) return null;
            drop = list.get(0);
            list.remove(0);
            if (vokzal[DispThread.myLevel].size() == 0) {
                changeMyLevel();
            }
        }

        return drop;
    }


    private void changeMyLevelH() {
        while (DispThread.myLevelH != 0 && DispThread.myLevelH >= DispThread.myLevel && vokzal[DispThread.myLevelH].size() == 0) {
            DispThread.myLevelH--;
        }
        if (vokzal[DispThread.myLevelH].size() == 0) {
            DispThread.myLevelH = 20;
            DispThread.myLevel = 20;
        }
       // if(changeTrack) changeTrackLevel();


    }

    private void changeMyLevel() {
        DispThread.myLevel++;
        while (DispThread.myLevelH >= DispThread.myLevel && vokzal[DispThread.myLevel].size() == 0) {
            DispThread.myLevel++;
        }
        if (DispThread.myLevel > DispThread.myLevelH) {
            DispThread.myLevel = 20;
            DispThread.myLevelH = 20;
        }
        //if(changeTrack) changeTrackLevel();

    }

    private void changeTrackLevel() {
        //LOGGER.trace("ownTrack.size" + ownTrack.size());
        synchronized (ownTrack) {
            DispThread.trackLevel = ownTrack.size() == 0 ? 20 : ownTrack.stream().min(Comparator.comparing(Drop::getRecNum)).get().recNum;
            changeTrack = false;
        }
    }

    private void deleteFromTrack(Drop drop) {
        ownTrack.remove(drop);
        if(drop.recNum == DispThread.trackLevel)
            changeTrackLevel();
    }


    public void inputDataToAmin() throws MPIException {
        //LOGGER.info("inputDataToAmin = "+ currentDrop.type);
        Amin curAmin = null;
        Element[] resInputFunc;

        curAmin = new Amin(currentDrop, pine.size(), myRank);
        pine.add(curAmin);
        currentDrop.setNumbOfMyAmine(curAmin.aminIdInPine);
        //System.arraycopy(currentDrop.inData, 0, curAmin.inputData, 0, curAmin.inputData.size());
        curAmin.inputData = currentDrop.inData;
        resInputFunc = currentDrop.inputFunction(curAmin.inputData, curAmin , ring);

        writeResultsAfterInpFunc(currentDrop, curAmin, resInputFunc);

        currentDrop.independentCalc(ring, curAmin);
        //LOGGER.info("end = "+ currentDrop.type);

    }


    @Override
    public void run() {
        long calcWaitTimeStart, calcWaitTimeEnd, calcWorkTimeStart, calcWorkTimeEnd;
        calcWaitTimeStart = calcWaitTimeEnd = calcWorkTimeStart = calcWorkTimeEnd = System.currentTimeMillis();

        while (!flToExit) {
           // ++counterCycle;
          // if(myRank==0) LOGGER.info("in calc cycle owntrack size = "+ownTrack.size()+ " vokzalempty = " +  DispThread.isEmptyVokzal());
            if (ownTrack.size() == 0 && /*Tools.isEmptyArray(vokzal)*/DispThread.isEmptyVokzal()) {
                if (!IamFree) {
                    IamFree = true;
                }
                continue;
            } else {
                IamFree = false;
                calcWaitTimeEnd = System.currentTimeMillis();
                calcWaitTime+=calcWaitTimeEnd-calcWaitTimeStart;
            }
            try {
                calcWorkTimeStart = System.currentTimeMillis();
                ProcFunc();
                calcWorkTimeEnd = System.currentTimeMillis();
                calcWorkTime += calcWorkTimeEnd - calcWorkTimeStart;
                calcWaitTimeStart = System.currentTimeMillis();
            } catch (MPIException e) {
                e.printStackTrace();
            }
        }
    }
    private void ProcFunc() throws MPIException {
       // LOGGER.info("go to get task");
        currentDrop = getTask(0);
        if (currentDrop != null) {
          //  LOGGER.info("get drop number = " + currentDrop.number + " type = " + currentDrop.type + " proc = "+currentDrop.procId);
            //LOGGER.info("currentdrop out data = " + Array.toString(currentDrop.outData));
            if (!Array.isEmpty(currentDrop.outData)) {
             //   LOGGER.trace("Drop result");

                //LOGGER.info("amin = " + currentDrop.aminId);
                writeResultsToAmin(currentDrop);

            } else {
                if (currentDrop.isItLeaf()) {
                   // LOGGER.info("Drop is leaf " + currentDrop.aminId + " id = "+ currentDrop.dropId);
                    currentDrop.sequentialCalc(ring);

                   /* for (int i = 0; i <currentDrop.inputDataLength ; i++) {
                        LOGGER.trace("after currentDrop.INdata = " + currentDrop.inData[i]);

                    }*/

                   /* for (int i = 0; i <currentDrop.outData.length ; i++) {

                        LOGGER.trace("after currentDrop.outdata = " + currentDrop.outData[i]);
                    }*/

                    if (currentDrop.aminId == -1 && myRank == 0) {
                       // LOGGER.info("go to finish whole task");
                        finishWholeTask(currentDrop);
                    } else if (currentDrop.procId == myRank) {

                        writeResultsToAmin(currentDrop);
                       //LOGGER.info("after writeResultsToAmin");

                    } else {
                        //LOGGER.trace(" bef add aerodrome results");
                        addToAerodromeResults(currentDrop);
                    }
                  //  LOGGER.trace("after drop is leaf vokzal empty = " + Tools.isEmptyArray(vokzal));
                } else {
                  //LOGGER.info("drop is not a leaf");
                    inputDataToAmin();
                }
            }
        }

        return;
    }
}
