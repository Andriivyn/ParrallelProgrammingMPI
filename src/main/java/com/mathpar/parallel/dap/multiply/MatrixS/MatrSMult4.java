package com.mathpar.parallel.dap.multiply.MatrixS;

import com.mathpar.log.MpiLogger;
import com.mathpar.matrix.MatrixS;
import com.mathpar.number.Array;
import com.mathpar.number.Element;
import com.mathpar.number.Ring;
import com.mathpar.parallel.dap.core.Amin;
import com.mathpar.parallel.dap.core.Drop;
import com.mathpar.parallel.dap.ldumw.LdumwDto;

import java.util.ArrayList;

public class MatrSMult4 extends Drop {
    private final static MpiLogger LOGGER = MpiLogger.getLogger(MatrSMult4.class);
    protected static int leafSize = 4;


    private static int[][] _arcs = new int[][]{
            {1, 0, 0, 1, 4, 1, 1, 1, 2, 1, 6, 3, 2, 0, 0, 2, 5, 1, 2, 1, 2, 2, 7, 3,
                    3, 2, 0, 3, 4, 1, 3, 3, 2, 3, 6, 3, 4, 2, 0, 4, 5, 1, 4, 3, 2, 4, 7, 3},
            {5, 0, 0},
            {5, 0, 1},
            {5, 0, 2},
            {5, 0, 3},
            {}};

    public MatrSMult4() {

        //Дроп має тип 5
        type = 5;

        //унікальний номер дропа
        number = cnum++;
        arcs = _arcs;
    }

    //Розгортання аміну з дропами, відповідно до графу, для обрахунку поточного дропа.
    @Override
    public ArrayList<Drop> doAmin() {
        ArrayList<Drop> amin = new ArrayList<Drop>();

        amin.add(new MatrSMultiplyScalar());
        amin.add(new MatrSMultiplyScalar());
        amin.add(new MatrSMultiplyScalar());
        amin.add(new MatrSMultiplyScalar());

        return amin;
    }

    @Override
    public void setVars() {
        switch (key) {
            case (0):
            case (1):
            case (102):
            case (111):
            case (122):
            case (114):
            case (110):
            case (119):
            case (124):
            case (125):
            case (126):
            case (123): {
                inputDataLength = 2;
                outputDataLength = 1;
                resultForOutFunctionLength = 4;
                break;
            }

            case (2): {
                inputDataLength = 2;
                outputDataLength = 2;
                resultForOutFunctionLength = 5;
                break;
            }

            //drop 12
            case (3): {
                inputDataLength = 6;
                outputDataLength = 3;
                resultForOutFunctionLength = 5;
                break;
            }
            //drop 5
            case (4): {
                inputDataLength = 6;
                outputDataLength = 1;
                resultForOutFunctionLength = 5;
                break;
            }
            // LDUMW drop 3 and 4
            case (103): {
                inputDataLength = 3;
                outputDataLength = 2;
                resultForOutFunctionLength = 4;
                break;
            }
            case (104): {
                inputDataLength = 3;
                outputDataLength = 2;
                resultForOutFunctionLength = 5;
                break;
            }
            case (105): {
                inputDataLength = 2;
                outputDataLength = 1;
                resultForOutFunctionLength = 5;
                break;
            }
            case (107):
            case (115):
            case (120):
            case (116):
            case (121): {
                inputDataLength = 3;
                outputDataLength = 1;
                resultForOutFunctionLength = 4;
                break;
            }
            case (109): {
                inputDataLength = 5;
                outputDataLength = 2;
                resultForOutFunctionLength = 5;
                break;
            }
            case (112): {
                inputDataLength = 5;
                outputDataLength = 4;
                resultForOutFunctionLength = 8;
                break;
            }
            case (113): {
                inputDataLength = 4;
                outputDataLength = 1;
                resultForOutFunctionLength = 5;
                break;
            }
            case (118): {
                inputDataLength = 5;
                outputDataLength = 1;
                resultForOutFunctionLength = 5;
                break;
            }
        }

        inData = new Element[inputDataLength];
        outData = new Element[outputDataLength];

        //LOGGER.info(inData);
    }

    //Послідовний обрахунок листових вершин
    @Override
    public void sequentialCalc(Ring ring) {
        // LOGGER.info("in sequentialCalc indata = " + inData[0] + ",  "+inData[1]);

        MatrixS A;
        MatrixS B;
        switch (key) {
            case (0):
                outData[0] = inData[0].multiply(inData[1], ring);
                break;
            case (1):
                outData[0] = inData[0].multiply(inData[1], ring).negate(ring);
                break;
            case (2): {

                MatrixS b = ((MatrixS) inData[0]).transpose();
                outData[1] = b;
                MatrixS bbT = b.multiply((MatrixS) inData[0], ring);
                outData[0] = ((MatrixS) inData[1]).subtract(bbT, ring);

                break;
            }
            case (3): {
                A = (MatrixS) inData[0];
                B = (MatrixS) inData[1];
                MatrixS M22_2 = A.multiply(B, ring).divideByNumber(inData[4].
                        multiply(inData[4], ring), ring);

                Element M22_3 = inData[5].multiply(M22_2, ring);//temporary(use multiplyLeftI)
                Element ds = inData[2].multiply(inData[3], ring).divide(inData[4], ring);

                outData[0] = M22_2;
                outData[1] = ds;
                outData[2] = M22_3;

                break;
            }
            case (4): {
                A = (MatrixS) inData[0];
                B = (MatrixS) inData[1];
                MatrixS ET = ((MatrixS) inData[4]).transpose();
                Element d0 = inData[5];
                MatrixS Md = ((MatrixS) inData[2]).multiplyByNumber(inData[3], ring);
                MatrixS M22_1 = A.multiply(ET.multiply(B, ring), ring).negate(ring).
                        divideByNumber(d0, ring).add(Md, ring);
                outData[0] = M22_1;
                break;
            }
            case (102): {
                LdumwDto F11 = ((LdumwDto) inData[0]);
                MatrixS A12 = (MatrixS) inData[1];
                LOGGER.info("-------------------------102-------------------------------");

                MatrixS X_U2 = F11.J().multiply(F11.M(), ring).multiply(A12, ring).divideByNumber(F11.A_n(), ring);

                outData[0] = X_U2;
                break;
            }
            case (103): {
                LdumwDto F11 = ((LdumwDto) inData[0]);
                MatrixS A12 = (MatrixS) inData[1];
                Element a = inData[2];

                LOGGER.info("-------------------------103-------------------------------");

                MatrixS A12_0 = F11.M().multiply(A12, ring);
                MatrixS A12_2 = F11.Dbar().multiply(A12_0, ring).divideByNumber(a, ring);
                LOGGER.info(" F11.M()=" + F11.M());
                LOGGER.info(" A12=" + A12);
                LOGGER.info("-------------------------/103-------------------------------");
                outData[0] = A12_0;
                outData[1] = A12_2;
                break;
            }
            case (104): {
                LdumwDto F11 = ((LdumwDto) inData[0]);
                MatrixS A21 = (MatrixS) inData[1];
                Element a = inData[2];

                LOGGER.info("-------------------------104-------------------------------");

                MatrixS A21_0 = A21.multiply(F11.W(), ring);
                MatrixS A21_2 = A21_0.multiply(F11.Dbar(), ring).divideByNumber(a, ring);
                LOGGER.info(" S A21_0=" + A21_0);
                LOGGER.info(" S A21_2=" + A21_2);
                LOGGER.info("-------------------------/104-------------------------------");
                outData[0] = A21_0;
                outData[1] = A21_2;
                break;
            }
            case (105): {
                LdumwDto F11 = ((LdumwDto) inData[0]);
                MatrixS A21 = (MatrixS) inData[1];
                LOGGER.info("-------------------------105-------------------------------");
                MatrixS X_L3 = A21.multiply(F11.W(), ring).multiply(F11.I(), ring).divideByNumber(F11.A_n(), ring);
                LOGGER.info("-------------------------/105-------------------------------");
                outData[0] = X_L3;
                break;
            }
            case (107): {
                LdumwDto F11 = ((LdumwDto) inData[0]);
                MatrixS A21_0 = (MatrixS) inData[1];
                MatrixS A12_0 = (MatrixS) inData[2];
                LOGGER.info("-------------------------107-------------------------------");
                MatrixS A21_1 = A21_0.multiplyByNumber(F11.A_n(), ring).multiply(F11.Dhat(), ring);
                MatrixS A12_1 = F11.Dhat().multiplyByNumber(F11.A_n(), ring).multiply(A12_0, ring);

                MatrixS D11PLUS = F11.D().transpose();
                MatrixS A22_0 = A21_1.multiply(D11PLUS, ring).multiply(A12_1, ring);
                LOGGER.info("-------------------------/107-------------------------------");
                outData[0] = A22_0;
                break;
            }
            case (109): {
                LdumwDto F11 = ((LdumwDto) inData[0]);
                MatrixS A22 = (MatrixS) inData[1];
                MatrixS A22_0 = (MatrixS) inData[2];
                LdumwDto F21 = ((LdumwDto) inData[3]);
                Element a = inData[4];

                LOGGER.info("-------------------------109-------------------------------");
                Element ak2 = F11.A_n().multiply(F11.A_n(), ring);
                Element aak2 = a.multiply(ak2, ring);
                MatrixS A22_1left = A22.multiplyByNumber(aak2, ring).subtract(A22_0, ring);
                MatrixS A22_1 = A22_1left.divideByNumber(a.multiply(F11.A_n(), ring), ring);

                MatrixS X_A22_2 = F21.Dbar().multiply(F21.M(), ring).multiply(A22_1, ring);
                LOGGER.info("-------------------------/109-------------------------------");
                outData[0] = A22_1;
                outData[1] = X_A22_2;
                break;
            }
            case (110): {
                LdumwDto F11 = ((LdumwDto) inData[0]);
                LdumwDto F21 = ((LdumwDto) inData[1]);

                MatrixS UU = F21.U().multiply(F11.U(), ring);

                outData[0] = UU;
                break;
            }
            case (111): {
                LdumwDto F11 = ((LdumwDto) inData[0]);
                LdumwDto F21 = ((LdumwDto) inData[1]);

                MatrixS U1_m1 = F11.W()
                        .multiply(F11.Dhat(), ring)
                        .multiply(F21.W(), ring)
                        .multiply(F21.Dhat(), ring);

                outData[0] = U1_m1;
                break;
            }
            case (112): {
                MatrixS X_A22_2 = (MatrixS) inData[0];
                LdumwDto F12 = ((LdumwDto) inData[1]);
                LdumwDto F21 = ((LdumwDto) inData[2]);
                LdumwDto F11 = ((LdumwDto) inData[3]);
                Element a = inData[4];

                LOGGER.info("-------------------------112-------------------------------");


                MatrixS A22_2 = X_A22_2.multiply(F12.W(), ring).multiply(F12.Dbar(), ring);
                Element lambda = F21.A_n().divide(F11.A_n(), ring);
                Element as = lambda.multiply(F12.A_n(), ring);

                LOGGER.info("al=" + F21.A_n());
                LOGGER.info("ak=" + F11.A_n());
                LOGGER.info("am=" + F12.A_n());

                MatrixS I12lambda=(F12.I().multiplyByNumber(lambda, ring)).add(F12.Ibar(), ring);
                MatrixS invD12hat=  I12lambda.inverse(ring).multiply(F12.Dhat(), ring);

                Element ak2 = F11.A_n().multiply(F11.A_n(), ring);
                MatrixS A22_3 = A22_2.divideByNumber(ak2, ring).divideByNumber(a, ring);
                LOGGER.info("-------------------------/112-------------------------------");
                outData[0] = lambda;
                outData[1] = as;
                outData[2] = A22_3;
                outData[3] = invD12hat;
                break;
            }
            case (113): {
                LdumwDto F21 = ((LdumwDto) inData[0]);
                MatrixS A22_1 = (MatrixS) inData[1];
                MatrixS X_U2 = (MatrixS) inData[2];
                Element a = inData[3];

                LOGGER.info("-------------------------113-------------------------------");
                MatrixS U2 = F21.J().multiply(F21.M(), ring)
                        .multiply(A22_1, ring)
                        .divideByNumber(F21.A_n().multiply(a, ring), ring)
                        .add(X_U2, ring);

                LOGGER.info("-------------------------/113-------------------------------");
                outData[0] = U2;
                break;
            }
            case (114): {
                LdumwDto F21 = ((LdumwDto) inData[0]);
                MatrixS A22_1 = (MatrixS) inData[1];
                LOGGER.info("-------------------------114-------------------------------");
                MatrixS Y_L3 = F21.Dbar().multiply(F21.M(), ring).multiply(A22_1, ring);
                LOGGER.info("-------------------------/114-------------------------------");
                outData[0] = Y_L3;
                break;
            }
            case (115): {

                MatrixS invD12hat = (MatrixS) inData[0];
                LdumwDto F12 = ((LdumwDto) inData[1]);
                LdumwDto F11 = ((LdumwDto) inData[2]);

                LOGGER.info("-------------------------115-------------------------------");
                MatrixS L1_m1 = invD12hat.multiply(F12.M(), ring)
                        .multiply(F11.Dhat(), ring)
                        .multiply(F11.M(), ring);
                LOGGER.info("-------------------------/115-------------------------------");
                outData[0] = L1_m1;
                break;
            }
            case (116): {

                LdumwDto F11 = ((LdumwDto) inData[0]);
                LdumwDto F12 = ((LdumwDto) inData[1]);
                Element lambda = inData[2];


                LOGGER.info("-------------------------116-------------------------------");
                MatrixS I12lambda=(F12.I().multiplyByNumber(lambda, ring)).add(F12.Ibar(), ring);
                MatrixS L12tilde = F12.L().multiply(I12lambda, ring);
                MatrixS X_L = F11.L().multiply(L12tilde, ring);
                LOGGER.info("-------------------------/116-------------------------------");
                outData[0] = X_L;
                break;
            }
            case (118): {
                MatrixS Y_L3 = (MatrixS) inData[0];
                LdumwDto F12 = ((LdumwDto) inData[1]);
                LdumwDto F11 = ((LdumwDto) inData[2]);
                MatrixS X_L3 = (MatrixS) inData[3];
                Element a = inData[4];
                LOGGER.info("-------------------------118-------------------------------");
                MatrixS L_3 = Y_L3.multiply(F12.W(), ring)
                        .multiply(F12.I(), ring)
                        .divideByNumber(F12.A_n().multiply(F11.A_n(), ring).multiply(a, ring), ring)
                        .add(X_L3, ring);

                LOGGER.info("-------------------------/118-------------------------------");
                outData[0] = L_3;

                break;
            }
            case (119): {
                MatrixS U1_m1 = (MatrixS) inData[0];
                MatrixS U2 = (MatrixS) inData[1];
                LOGGER.info("-------------------------119-------------------------------");

                MatrixS X_W2 = U1_m1.multiply(U2, ring).negate(ring);
                outData[0] = X_W2;
                LOGGER.info("-------------------------/119-------------------------------");
                break;
            }
            case (120): {
                MatrixS invD12hat = (MatrixS) inData[0];
                LdumwDto F12 = ((LdumwDto) inData[1]);
                LdumwDto F22 = ((LdumwDto) inData[2]);
                LOGGER.info("-------------------------120-------------------------------");
                MatrixS U4_m1 = F12.W().multiply(invD12hat, ring).multiply(F22.W(), ring).multiply(F22.Dhat(), ring);
                LOGGER.info("-------------------------/120-------------------------------");
                outData[0] = U4_m1;

                break;
            }
            case (121): {
                LdumwDto F22 = ((LdumwDto) inData[0]);
                LdumwDto F12 = ((LdumwDto) inData[1]);
                Element lambda = inData[2];
                LOGGER.info("-------------------------121-------------------------------");
                MatrixS J12lambda = F12.J().multiplyByNumber(lambda, ring).add(F12.Jbar(), ring);
                MatrixS U12tilde = J12lambda.multiply(F12.U(), ring);
                MatrixS X_U = F22.U().multiply(U12tilde, ring);
                LOGGER.info("-------------------------/121-------------------------------");
                outData[0] = X_U;
                break;
            }
            case (122): {
                LdumwDto F22 = ((LdumwDto) inData[0]);
                LdumwDto F21 = ((LdumwDto) inData[1]);
                LOGGER.info("-------------------------122-------------------------------");
                MatrixS L4_m1 = F22.Dhat()
                        .multiply(F22.M(), ring)
                        .multiply(F21.Dhat(), ring)
                        .multiply(F21.M(), ring);
                LOGGER.info("-------------------------/122-------------------------------");
                outData[0] = L4_m1;
                break;
            }
            case (123): {
                LdumwDto F21 = ((LdumwDto) inData[0]);
                LdumwDto F22 = ((LdumwDto) inData[1]);
                LOGGER.info("-------------------------123-------------------------------");
                MatrixS LL = F21.L().multiply(F22.L(), ring);
                LOGGER.info("-------------------------/123-------------------------------");
                outData[0] = LL;
                break;
            }
            case (124): {
                MatrixS L4_m1 = (MatrixS) inData[0];
                MatrixS L3 = (MatrixS) inData[1];
                LOGGER.info("-------------------------124-------------------------------");
                MatrixS X_m12 = L4_m1.multiply(L3, ring).negate(ring);
                outData[0] = X_m12;
                LOGGER.info("-------------------------/124-------------------------------");
                break;
            }
            case (125): {
                MatrixS X_W21 = (MatrixS) inData[0];
                MatrixS U4_m1 = (MatrixS) inData[1];

                LOGGER.info("-------------------------125-------------------------------");
                outData[0] = X_W21.multiply(U4_m1, ring);

                LOGGER.info("-------------------------/125-------------------------------");
                break;
            }
            case (126): {
                MatrixS X_m12 = (MatrixS) inData[0];
                MatrixS L1_m1 = (MatrixS) inData[1];
                LOGGER.info("-------------------------126-------------------------------");
                MatrixS X_m2 = X_m12.multiply(L1_m1, ring);
                outData[0] = X_m2;
                LOGGER.info("-------------------------/126-------------------------------");
                break;
            }
        }
    }

    @Override
    //Вхідна функція дропа, розбиває вхідні дані на блоки.
    public MatrixS[] inputFunction(Element[] input, Amin amin, Ring ring) {

        MatrixS[] res = new MatrixS[8];

        MatrixS ms = null;
        MatrixS ms1 = null;

        if (input[0] instanceof MatrixS && input[1] instanceof MatrixS) {
            ms = (MatrixS) input[0];
            ms1 = (MatrixS) input[1];
        }

        switch (key) {
            case (2): {
                ms = ((MatrixS) input[0]).transpose();
                ms1 = (MatrixS) input[0];
                amin.resultForOutFunction[4] = ms;
                break;
            }
            case (4): {
                ms = (MatrixS) input[0];
                ms1 = (MatrixS) input[1];
                MatrixS E11T = ((MatrixS) input[4]).transpose();
                ms1 = E11T.multiply(ms1, ring);
                break;
            }
            case (102):
            case (103): {
                LdumwDto F11 = (LdumwDto) input[0];
                LOGGER.info("F11.M(): " + F11.M());
                ms = F11.M(); // M
                ms1 = (MatrixS) input[1]; // A12
                LOGGER.info("-------------------------102 103 input-------------------------------");
                break;
            }
            case (104):
            case (105): {
                LdumwDto F11 = (LdumwDto) input[0];
                ms = (MatrixS) input[1]; // A21
                ms1 = F11.W(); // W
                LOGGER.info("-------------------------104 105 input-------------------------------");
                break;
            }
            case (107): {
                LdumwDto F11 = ((LdumwDto) inData[0]);
                MatrixS A21_0 = (MatrixS) inData[1];
                MatrixS A12_0 = (MatrixS) inData[2];

                LOGGER.info("-------------------------107 input-------------------------------");

                MatrixS A21_1 = A21_0.multiplyByNumber(F11.A_n(), ring).multiply(F11.Dhat(), ring);
                MatrixS A12_1 = F11.Dhat().multiplyByNumber(F11.A_n(), ring).multiply(A12_0, ring);

                MatrixS D11PLUS = F11.D().transpose();

                ms = A21_1.multiply(D11PLUS, ring);
                ms1 = A12_1;
                break;
            }
            case (109): {
                LdumwDto F11 = ((LdumwDto) inData[0]);
                MatrixS A22 = (MatrixS) inData[1];
                MatrixS A22_0 = (MatrixS) inData[2];
                LdumwDto F21 = ((LdumwDto) inData[3]);
                Element a = inData[4];

                LOGGER.info("-------------------------109 input-------------------------------");
                Element ak2 = F11.A_n().multiply(F11.A_n(), ring);
                Element aak2 = a.multiply(ak2, ring);
                MatrixS A22_1left = A22.multiplyByNumber(aak2, ring).subtract(A22_0, ring);
                MatrixS A22_1 = A22_1left.divideByNumber(a.multiply(F11.A_n(), ring), ring);

                LOGGER.info("-------------------------/109 input-------------------------------");

                ms = F21.Dbar().multiply(F21.M(), ring);
                ms1 = A22_1;

                amin.resultForOutFunction[4] = A22_1;
                break;
            }
            case (110): {
                LdumwDto F11 = ((LdumwDto) inData[0]);
                LdumwDto F21 = ((LdumwDto) inData[1]);
                LOGGER.info("-------------------------110 input-------------------------------");
                ms = F21.U();
                ms1 = F11.U();
                break;
            }
            case (111): {
                LdumwDto F11 = ((LdumwDto) inData[0]);
                LdumwDto F21 = ((LdumwDto) inData[1]);
                LOGGER.info("-------------------------111 input-------------------------------");
                ms = F11.W().multiply(F11.Dhat(), ring);
                ms1 = F21.W().multiply(F21.Dhat(), ring);

                break;
            }
            case (112): {
                LdumwDto F12 = (LdumwDto) input[1];
                LOGGER.info("-------------------------112 input-------------------------------");
                ms = (MatrixS) input[0];
                ms1 = F12.W();
                break;
            }
            case (113): {
                LdumwDto F21 = ((LdumwDto) inData[0]);
                MatrixS A22_1 = (MatrixS) inData[1];
                LOGGER.info("-------------------------113 input-------------------------------");
                ms = F21.J().multiply(F21.M(), ring);
                ms1 = A22_1;
                break;
            }
            case (114): {
                LdumwDto F21 = ((LdumwDto) inData[0]);
                MatrixS A22_1 = (MatrixS) inData[1];
                LOGGER.info("-------------------------114 input-------------------------------");
                ms = F21.Dbar().multiply(F21.M(), ring);
                ms1 = A22_1;
                break;
            }
            case (115): {
                MatrixS invD12hat = (MatrixS) inData[0];
                LdumwDto F12 = ((LdumwDto) inData[1]);
                LdumwDto F11 = ((LdumwDto) inData[2]);
                LOGGER.info("-------------------------115 input-------------------------------");
                ms = invD12hat.multiply(F12.M(), ring);
                ms1 = F11.Dhat().multiply(F11.M(), ring);
                break;
            }
            case (116): {
                LdumwDto F11 = ((LdumwDto) inData[0]);
                LdumwDto F12 = ((LdumwDto) inData[1]);
                Element lambda = inData[2];

                LOGGER.info("-------------------------116 input-------------------------------");

                MatrixS I12lambda = F12.I().multiplyByNumber(lambda, ring).add(F12.Ibar(), ring);
                MatrixS L12tilde = F12.L().multiply(I12lambda, ring);
                ms = F11.L();
                ms1 = L12tilde;
                break;
            }
            case (118): {
                LdumwDto F12 = ((LdumwDto) inData[1]);
                LOGGER.info("-------------------------118 input-------------------------------");
                ms = (MatrixS) inData[0]; //Y_L3;
                ms1 = F12.W();
                break;
            }
            case (119): {
                ms = (MatrixS) inData[0];
                ms1 = (MatrixS) inData[1];
                LOGGER.info("-------------------------119 input-------------------------------");
                break;
            }
            case (120): {
                MatrixS invD12hat = (MatrixS) inData[0];
                LdumwDto F12 = ((LdumwDto) inData[1]);
                LdumwDto F22 = ((LdumwDto) inData[2]);
                LOGGER.info("-------------------------120 input-------------------------------");
                ms = F12.W().multiply(invD12hat, ring);
                ms1 = F22.W().multiply(F22.Dhat(), ring);
                break;
            }
            case (121): {
                LdumwDto F22 = ((LdumwDto) inData[0]);
                LdumwDto F12 = ((LdumwDto) inData[1]);
                Element lambda = inData[2];

                LOGGER.info("-------------------------121 input-------------------------------");

                MatrixS J12lambda = F12.J().multiplyByNumber(lambda, ring).add(F12.Jbar(), ring);
                MatrixS U12tilde = J12lambda.multiply(F12.U(), ring);

                ms = F22.U();
                ms1 = U12tilde;
                break;
            }
            case (122): {
                LdumwDto F22 = ((LdumwDto) inData[0]);
                LdumwDto F21 = ((LdumwDto) inData[1]);
                LOGGER.info("-------------------------122 input-------------------------------");
                ms = F22.Dhat().multiply(F22.M(), ring);
                ms1 = F21.Dhat().multiply(F21.M(), ring);
                break;
            }
            case (123): {
                LdumwDto F21 = ((LdumwDto) inData[0]);
                LdumwDto F22 = ((LdumwDto) inData[1]);
                LOGGER.info("-------------------------123 input-------------------------------");
                ms = F21.L();
                ms1 = F22.L();

                break;
            }
            case (124): {
                ms = (MatrixS) inData[0];
                ms1 = (MatrixS) inData[1];
                LOGGER.info("-------------------------124 input-------------------------------");
                break;
            }
            case (125): {
                ms = (MatrixS) inData[0];
                ms1 = (MatrixS) inData[1];
                LOGGER.info("-------------------------125 input-------------------------------");
                break;
            }
            case (126): {
                ms = (MatrixS) inData[0];
                ms1 = (MatrixS) inData[1];
                LOGGER.info("-------------------------126 input-------------------------------");
                break;
            }
        }
        Array.concatTwoArrays(ms.split(), ms1.split(), res);
        return res;

    }

    @Override
    public void independentCalc(Ring ring, Amin amin) {
        switch (key) {
            case (0):
            case (1):
            case (2):
            case (102):
            case (103):
            case (107):
            case (109):
            case (110):
            case (111):
            case (114):
            case (115):
            case (116):
            case (120):
            case (121):
            case (122):
            case (123):
            case (124):
            case (125):
            case (126):
            case (119):
                break;
            case (3): {
                Element ds = inData[2].multiply(inData[3], ring).divide(inData[4], ring);
                amin.resultForOutFunction[4] = ds;
                break;
            }
            case (4): {
                Element md = inData[2].multiply(inData[3], ring);
                amin.resultForOutFunction[4] = md;
                break;
            }
            case (104): {
                LdumwDto F11 = (LdumwDto) inData[0];
                Element a = inData[2];

                LOGGER.info("-------------------------104 indepen-------------------------------");


                amin.resultForOutFunction[4] = F11.Dbar().divideByNumber(a, ring);
                break;
            }
            case (105): {
                LdumwDto F11 = ((LdumwDto) inData[0]);

                LOGGER.info("-------------------------105 indepen-------------------------------");
                amin.resultForOutFunction[4] = F11.I().divideByNumber(F11.A_n(), ring); // F11.I / F11.ak
                break;
            }
            case (112): {
                LdumwDto F12 = ((LdumwDto) inData[1]);
                LdumwDto F21 = ((LdumwDto) inData[2]);
                LdumwDto F11 = ((LdumwDto) inData[3]);

                LOGGER.info("-------------------------112 indepen-------------------------------");

                LOGGER.info("al=" + F21.A_n());
                LOGGER.info("ak=" + F11.A_n());
                LOGGER.info("am=" + F12.A_n());


                Element lambda = F21.A_n().divide(F11.A_n(), ring);
                Element as = lambda.multiply(F12.A_n(), ring);

                MatrixS I12lambdaM2=(F12.I().divideByNumbertoFraction(lambda, ring)).add(F12.Ibar(), ring);
                MatrixS invD12hat=  I12lambdaM2.multiply(F12.Dhat(), ring);

                Element ak2 = F11.A_n().multiply(F11.A_n(), ring);

                amin.resultForOutFunction[4] = lambda;
                amin.resultForOutFunction[5] = as;
                amin.resultForOutFunction[6] = invD12hat;
                amin.resultForOutFunction[7] = ak2;
                break;
            }
            case (113): {
                LdumwDto F21 = ((LdumwDto) inData[0]);
                Element a = inData[3];


                LOGGER.info("-------------------------113 indepen-------------------------------");


                amin.resultForOutFunction[4] = F21.A_n().multiply(a, ring);

                break;
            }
            case (118): {
                LdumwDto F12 = ((LdumwDto) inData[1]);
                LdumwDto F11 = ((LdumwDto) inData[2]);

                Element a = inData[4];

                LOGGER.info("-------------------------118 indepen-------------------------------");

                MatrixS remainder = F12.I()
                        .divideByNumber(F12.A_n().multiply(F11.A_n(), ring).multiply(a, ring), ring);

                amin.resultForOutFunction[4] = remainder;
                break;
            }

        }
        return;
    }

    //Вихідна функція дропа, яка збирає блоки в результат
    @Override
    public Element[] outputFunction(Element[] input, com.mathpar.number.Ring ring) {
        MatrixS[] resmat = new MatrixS[4];
        for (int i = 0; i < 4; i++) {
            resmat[i] = (MatrixS) input[i];
        }

        Element[] res = new Element[outputDataLength];
        switch (key) {
            case (0):
            case (107):
            case (110):
            case (111):
            case (114):
            case (115):
            case (116):
            case (120):
            case (121):
            case (122):
            case (123):
            case (125):
            case (126):
                LOGGER.info("-------------------------0 output-------------------------------");
                res = new MatrixS[]{MatrixS.join(resmat)};
                break;
            case (1):
            case (119):
            case (124):
                LOGGER.info("------------------------- 1 output-------------------------------");
                res = new MatrixS[]{MatrixS.join(resmat).negate(ring)};
                break;
            case (2): {
                res = new Element[]{inData[1].add(MatrixS.join(resmat).negate(ring), ring), input[4]};
                break;
            }
            case (3): {

                MatrixS M22_2 = MatrixS.join(resmat).divideByNumber(inData[4].
                        multiply(inData[4], ring), ring);

                Element M22_3 = inData[5].multiply(M22_2, ring);//temporary(use multiplyLeftI)
                res = new Element[]{M22_2, input[4], M22_3};

                break;
            }
            case (4): {

                res = new MatrixS[]{MatrixS.join(resmat).negate(ring).
                        divideByNumber(inData[5], ring).add((MatrixS) input[4], ring)};
                break;
            }
            case (102): {
                LOGGER.info("-------------------------102 output-------------------------------");
                LdumwDto F11 = (LdumwDto) inData[0];

                res = new Element[]{
                        F11.J().multiply(MatrixS.join(resmat), ring).divideByNumber(F11.A_n(), ring)
                };
                break;
            }
            case (103): {
                LOGGER.info("-------------------------103 output-------------------------------");
                LdumwDto F11 = (LdumwDto) inData[0];
                Element a = inData[2];
                MatrixS A12_0 = MatrixS.join(resmat);
                MatrixS A12_2 = F11.Dbar().multiply(A12_0, ring).divideByNumber(a, ring);
                res = new Element[]{
                        A12_0, A12_2
                };
                break;
            }
            case (104): {
                LOGGER.info("-------------------------104 output-------------------------------");
                Element A21_0 = MatrixS.join(resmat);
                Element A21_2 = A21_0.multiply(input[4], ring);
                LOGGER.info(" O A21_0=" + A21_0);
                LOGGER.info(" O A21_2=" + A21_2);

                res = new Element[]{
                        A21_0, A21_2
                };
                break;
            }
            case (105): {
                LOGGER.info("-------------------------105 output-------------------------------");
                MatrixS A21xW21 = MatrixS.join(resmat);
                Element X_L3 = A21xW21.multiply(input[4], ring);

                res = new Element[]{
                        X_L3
                };
                break;
            }
            case (109): {
                LOGGER.info("-------------------------109 output-------------------------------");
                res = new Element[]{
                        input[4], MatrixS.join(resmat)
                };

                break;
            }

            case (112): {
                LOGGER.info("-------------------------112 output-------------------------------");
                LdumwDto F12 = (LdumwDto) inData[1];
                Element a = inData[4];
                MatrixS A22_2 = MatrixS.join(resmat).multiply(F12.Dbar(), ring);
                MatrixS A22_3 = A22_2.divideByNumber(input[7], ring).divideByNumber(a, ring);

                res = new Element[]{
                        input[4], input[5], A22_3, input[6]
                };
                break;
            }
            case (113): {
                LOGGER.info("-------------------------113 output-------------------------------");
                MatrixS matrixS = MatrixS.join(resmat);
                MatrixS X_U2 = (MatrixS) inData[2];
                Element ala = input[4];
                res = new Element[]{
                        matrixS.divideByNumber(ala, ring).add(X_U2, ring)
                };
                break;
            }
            case (118): {
                LOGGER.info("-------------------------118 output-------------------------------");
                Element matrixS = MatrixS.join(resmat);
                res = new Element[]{
                        matrixS.multiply(input[4], ring).add(inData[3], ring) // L3
                };

                break;
            }
        }

        return res;
    }

    //Перевіряє чи є дроп листовим
    @Override
    public boolean isItLeaf() {
        if (inData[0] instanceof LdumwDto) {
            LdumwDto ldumwDto = (LdumwDto) inData[0];
            return (ldumwDto.L().size <= leafSize);
        }
        MatrixS ms = (MatrixS) inData[0];
        return (ms.size <= leafSize);
    }

    @Override
    public void setLeafSize(int dataSize) {
        leafSize = dataSize;
    }
}
