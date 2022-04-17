package com.mathpar.parallel.dap.ldumw;

import com.mathpar.log.MpiLogger;
import com.mathpar.matrix.LDUMW;
import com.mathpar.matrix.MatrixS;
import com.mathpar.number.Element;
import com.mathpar.number.Fraction;
import com.mathpar.number.Ring;
import com.mathpar.parallel.dap.core.Amin;
import com.mathpar.parallel.dap.core.Drop;
import com.mathpar.parallel.dap.multiply.MatrixS.MatrSMult4;

import java.util.ArrayList;

import static com.mathpar.matrix.LDUMW.DtoUnit;

public class LdumwFact extends Drop {
    private static int leafSize = 4;
    private final static MpiLogger LOGGER = MpiLogger.getLogger(LdumwFact.class);
    private static int[][] arcs_ = new int[][]{ // TODO change the order of input params to use similar drops in different cases
            {1, 0, 0,    1, 4, 1,    2, 1, 1,    3, 1, 1,    3, 4, 2,    4, 2, 1,    4, 4, 2,    5, 2, 1,    9, 3, 1,    9, 4, 4,    12, 4, 4,    13, 4, 3,    18, 4, 4},   //27, 4, 8},//0. inputFunction
            {2, 0, 0,    3, 0, 0,    4, 0, 0,    5, 0, 0,    6, 1, 1,    7, 0, 0,    8, 1, 1,    9, 0, 0,    10, 0, 0,    11, 0, 0,    12, 0, 3,    15, 0, 2,    16, 0, 0,    18, 0, 2,    27, 0, 0},//1. F11 = LDU /* 1 крок - це все в одній змінній? */
            {13, 0, 2},//2. X_U2
            {6, 1, 0,    7, 0, 2},//3. A12_0  and  A12_2
            {7, 0, 1,    8, 1, 0},//4. A21_0  and  A21_2
            {18, 0, 3},//5. X_L3
            {12, 0, 1,    15, 0, 1,    16, 0, 1,    18, 0, 1,    20, 0, 1,    21, 0, 1,    27, 0, 2}, //6. F12(am)
            {9, 0, 2}, //7. A22_0
            {9, 0, 3,    10, 0, 1,    11, 0, 1,    12, 0, 2,    13, 0, 0,    14, 0, 0,    22, 0, 1,    23, 0, 0,    27, 0, 1}, //8. F21(al)
            {12, 1, 0,    13, 0, 1,    14, 0, 1}, //9. A22_1  and  X_A22_2
            {27, 0, 5}, //10. UU
            {19, 0, 0,     27, 0, 8}, //11. U1_m1
                {15, 3, 0,     16, 0, 2,     17, 1, 1,    17, 2, 0,    20, 3, 0,    21, 0, 2,    27, 0, 4,    27, 3, 17,    27, 1, 18}, //12. lambda  and  as  and  A22_3  and  invD12hat
            {19, 0, 1,    27, 0, 6}, //13. U2
            {18, 0, 0}, //14. Y_L3 (L3H2)
            {26, 0, 1,    27, 0, 13}, //15. L1_m1 ?? graphic has wrong number on input
            {27, 0, 16}, //16. X_L
            {20, 0, 2,    21, 0, 0,    22, 0, 0,    23, 0, 1,    27, 0, 3}, //17. F22
            {24, 0, 1,    27, 0, 15}, //18. L3
            {25, 0, 0}, //19. X_W21
            {25, 0, 1,    27, 0, 10}, //20. U4_m1
            {27, 0, 7}, //21. X_U
            {24, 0, 0,    27, 0, 11}, //22. L4_m1
            {27, 0, 14}, //23. LL
            {26, 0, 0}, //24. X_m12
            {27, 0, 9}, //25. X_W2
            {27, 0, 12}, //26. X_m2
            {}};//27. OutputFunction

    public LdumwFact() {
        arcs = arcs_;
        type = 23;
        number = cnum++;
        inputDataLength = 2;
        outputDataLength = 2;
        inData = new Element[inputDataLength];
        outData = new Element[outputDataLength];
        resultForOutFunctionLength = 19;
    }

    @Override
    public ArrayList<Drop> doAmin() {
        ArrayList<Drop> amin = new ArrayList<>();
        // step 1
        amin.add(new LdumwFact());
        // step 2
        amin.add(new MatrSMult4());
        amin.get(1).key = 102;
        // step 3
        amin.add(new MatrSMult4());
        amin.get(2).key = 103;
        // step 4
        amin.add(new MatrSMult4());
        amin.get(3).key = 104;
        // step 5
        amin.add(new MatrSMult4());
        amin.get(4).key = 105;
        // step 6
        amin.add(new LdumwFact());
        // step 7
        amin.add(new MatrSMult4());
        amin.get(6).key = 107;
        // step 8
        amin.add(new LdumwFact());
        // step 9
        amin.add(new MatrSMult4());
        amin.get(8).key = 109;
        // step 10
        amin.add(new MatrSMult4());
        amin.get(9).key = 110;
        // step 11
        amin.add(new MatrSMult4());
        amin.get(10).key = 111;
        // step 12
        amin.add(new MatrSMult4());
        amin.get(11).key = 112;
        // step 13
        amin.add(new MatrSMult4());
        amin.get(12).key = 113;
        // step 14
        amin.add(new MatrSMult4());
        amin.get(13).key = 114;
        // step 15
        amin.add(new MatrSMult4());
        amin.get(14).key = 115;
        // step 16
        amin.add(new MatrSMult4());
        amin.get(15).key = 116;
        // step 17
        amin.add(new LdumwFact());
        // step 18
        amin.add(new MatrSMult4());
        amin.get(17).key = 118;
        // step 19
        amin.add(new MatrSMult4());
        amin.get(18).key = 119;
        // step 20
        amin.add(new MatrSMult4());
        amin.get(19).key = 120;
        // step 21
        amin.add(new MatrSMult4());
        amin.get(20).key = 121;
        // step 22
        amin.add(new MatrSMult4());
        amin.get(21).key = 122;
        // step 23
        amin.add(new MatrSMult4());
        amin.get(22).key = 123;
        // step 24
        amin.add(new MatrSMult4());
        amin.get(23).key = 124;
        // step 25
        amin.add(new MatrSMult4());
        amin.get(24).key = 125;
        // step 26
        amin.add(new MatrSMult4());
        amin.get(25).key = 126;


        return amin;
    }


    @Override
    public void sequentialCalc(Ring ring) {
        MatrixS A = (MatrixS) inData[0];
        LOGGER.info("-------------------------sequentialCalc-------------------------------");
        LOGGER.info("DropID: " + dropId);
        LOGGER.info("A size: " + A.size);
        Element a = inData[1];

        LdumwDto FF = LDUMW.LDUMW(A, a, ring);
        LOGGER.info("FF is calculated");
        LOGGER.info("-------------------------/sequentialCalc-------------------------------");
        outData[0] = FF;
        outData[1] = FF.A_n();
    }

    @Override
    public Element[] inputFunction(Element[] input, Amin amin, Ring ring) {
        MatrixS A = (MatrixS) input[0];
        Element a = input[1];
        MatrixS[] split = A.split();
        Element[] res = new Element[5];
        System.arraycopy(split, 0, res, 0, 4);

//        LOGGER.info("inputFunction.a = " + A.size);
        res[4] = a;

        return res;
    }

    @Override
    public Element[] outputFunction(Element[] input, Ring ring) {
        LOGGER.info("Start outputFunction");

        LdumwDto F11 = (LdumwDto) input[0];
        LdumwDto F21 = (LdumwDto) input[1];
        LdumwDto F12 = (LdumwDto) input[2];
        LdumwDto F22 = (LdumwDto) input[3];
        Element lambda = input[4];

        MatrixS UU = (MatrixS) input[5];
        MatrixS U2 = (MatrixS) input[6];
        MatrixS X_U = (MatrixS) input[7];

        Element a = inData[1];

        MatrixS U1_m1 = (MatrixS) input[8];
        MatrixS X_W2 = (MatrixS) input[9];
        MatrixS U4_m1 = (MatrixS) input[10];

        MatrixS L4_m1 = (MatrixS) input[11];
        MatrixS X_m2 = (MatrixS) input[12];
        MatrixS L1_m1 = (MatrixS) input[13];

        MatrixS LL = (MatrixS) input[14];
        MatrixS L3 = (MatrixS) input[15];
        MatrixS X_L = (MatrixS) input[16];
        MatrixS invD12hat = (MatrixS) input[17];

        LOGGER.info("LL=" + LL);
        LOGGER.info("L3=" + L3);

        MatrixS L = MatrixS.join(new MatrixS[]{X_L, MatrixS.zeroMatrix(), L3, LL});
        MatrixS D = MatrixS.join(new MatrixS[]{
                F11.D(), F12.D().multiplyByNumber(lambda.multiply(lambda, ring), ring),
                F21.D(), F22.D()});
        MatrixS U = MatrixS.join(new MatrixS[]{UU, U2, MatrixS.zeroMatrix(), X_U});


        LdumwDto ldumw = new LdumwDto(L, D, U, F22.A_n());
        ldumw.IJMap(F22.A_n(), ring);

        // TODO algo by scheme
//        Element ar_m1 = F22.A_n().inverse(ring);
////        LOGGER.info("F22.A_n().inverse(ring)=" + F22.A_n().inverse(ring));
////        LOGGER.info("F22.A_n().negate(ring)=" + F22.A_n().negate(ring));
//        MatrixS DhatLeft = D.multiplyByNumber(a.multiply(ar_m1, ring), ring);
//        MatrixS DhatRight = ldumw.Dbar().multiplyByNumber(ar_m1, ring);
//        MatrixS Dhat = DhatLeft.add(DhatRight, ring);
//        MatrixS Dinv = Dhat.inverse(ring);
//        MatrixS M = Dinv.multiply(MatrixS.join(new MatrixS[]{L1_m1, MatrixS.zeroMatrix(), X_m2, L4_m1}), ring);
//
//
//        MatrixS W = MatrixS.join(new MatrixS[]{U1_m1, X_W2, MatrixS.zeroMatrix(), U4_m1}).multiply(Dinv, ring);
////
//
//        //TODO algo by standard
        Element am = F12.A_n();
        Element al = F21.A_n();
        Element ak = F11.A_n();
        Element as = input[18];
        Element a_n = ldumw.A_n();

        MatrixS L3prim=L3.negate(ring).multiply(invD12hat, ring).multiply(F12.M(), ring).multiply(F11.Dhat().multiply(F11.M(), ring), ring);
        MatrixS DhUnit= DtoUnit(D,ring.numberONE, ring).add(ldumw.Dbar(), ring).transpose();
        MatrixS[] Eprim=DhUnit.split();

        MatrixS U2prim=F11.W().multiply(F11.Dhat().multiply(F21.W().multiply(F21.Dhat().multiply(U2.negate(ring), ring), ring), ring), ring);
        // Du=
        MatrixS D11prim= DtoUnit(F11.D(),ring.numberONE,ring).add(F11. Dbar(), ring) ;
        MatrixS D12prim=DtoUnit(F21.D(), ak, ring).add(F21.Dbar().multiplyByNumber(a, ring), ring);
        MatrixS D21prim=DtoUnit(F12.D(), al ,ring).add(F12.Dbar().multiplyByNumber(a, ring), ring);
        MatrixS D22prim= DtoUnit(F22.D(), as,ring).add(F22.Dbar().multiplyByNumber(a, ring), ring);

        MatrixS V11A = F21.W().multiply(F21.Dbar(), ring).multiply(Eprim[0], ring);
        MatrixS V11B =F11.W().multiply(D11prim, ring).multiply(V11A, ring);
        MatrixS  V11= V11B.multiplyByNumber(new Fraction(a_n, ak.multiply(al, ring)), ring);
        MatrixS V12A= F21.W().multiply(D12prim, ring).multiply(Eprim[1], ring);
        MatrixS V12B= F11.W().multiply(F11.Dbar(), ring).multiply( V12A, ring);
        MatrixS  V12= V12B.multiplyByNumber(new Fraction(a_n, ak.multiply(al, ring).multiply(a , ring)), ring);
        MatrixS V21A =F12.W().multiply(D21prim, ring).multiply(F22.W(), ring).multiply(F22.Dbar(), ring);
        MatrixS V21B =V21A.multiply(Eprim[2], ring);
        MatrixS  V21= V21B.multiplyByNumber(new Fraction(ring.numberONE, am.multiply(a, ring)), ring);
        MatrixS V22A=F12.W().multiply(F12.Dbar(), ring).multiply(F22.W(), ring).multiply(D22prim, ring);
        MatrixS V22B=V22A.multiply(Eprim[3], ring);
        MatrixS  V22= V22B.multiplyByNumber(new Fraction(ring.numberONE, a.multiply(am, ring)), ring);


        MatrixS W=MatrixS.join(new MatrixS[]{V11.add(U2prim.multiply(V21, ring), ring),
                V12.add(U2prim.multiply(V22, ring), ring), V21, V22});
////        if(n>2){}
////     System.out.println( "V11="+V11+"; V12="+V12+"; V21="+V21+"; V22="+V22+"; U2prim="+U2prim+"; L3prim="+L3prim+";");
//
//
//
//
        MatrixS N11A =Eprim[0].multiply(F12.Dbar(), ring).multiply(F12.M(), ring);
        MatrixS  N11B =N11A.multiply(D11prim, ring).multiply(F11.M(), ring);
        MatrixS  N11= N11B.multiplyByNumber(new Fraction(a_n, ak.multiply(am, ring)), ring);
        MatrixS N21A=Eprim[2].multiply(D21prim, ring).multiply(F12.M(), ring);
        MatrixS N21B=N21A.multiply(F11.Dbar(), ring).multiply(F11.M(), ring);
        MatrixS    N21= N21B.multiplyByNumber(new Fraction(a_n, ak.multiply(am, ring).multiply(a , ring)), ring);
        MatrixS N12A =Eprim[1].multiply(F22.Dbar(), ring).multiply(F22.M(), ring);
        MatrixS N12B =N12A.multiply(D12prim, ring).multiply(F21.M(), ring);
        MatrixS  N12= N12B.multiplyByNumber(new Fraction(ring.numberONE, al.multiply(a, ring)), ring);
        MatrixS N22A=Eprim[3].multiply(D22prim, ring).multiply(F22.M(), ring);
        MatrixS N22B=N22A.multiply(F21.Dbar(), ring).multiply(F21.M(), ring);
        MatrixS  N22= N22B.multiplyByNumber(new Fraction(ring.numberONE, a.multiply(al, ring)), ring);

        MatrixS M =MatrixS.join(new MatrixS[]{N11.add(N12.multiply(L3prim, ring), ring),
                N12, N21.add(N22.multiply(L3prim, ring), ring), N22});


        ldumw.setM(M);
        ldumw.setW(W);
        LOGGER.info("End outputFunction");
        return new Element[]{ldumw, ldumw.A_n()};
    }

    @Override
    public boolean isItLeaf() {
        MatrixS ms = (MatrixS) inData[0];
        return (ms.size <= leafSize);
    }

    @Override
    public void setLeafSize(int dataSize) {
        leafSize = dataSize;
    }
}
