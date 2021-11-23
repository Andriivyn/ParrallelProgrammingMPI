package com.mathpar.matrix;
import com.mathpar.number.*;

import java.util.Random;

public class SVD {
    public static void main(String[] args) {
        Ring ring = new Ring("R64[]");
        double[][] m = generateMatrix(4);
        MatrixD A = new MatrixD(m, ring);
        MatrixD[] result = SVD(A, ring);
        MatrixD U = result[0];
        MatrixD D = result[1];
        MatrixD W = result[2];
        MatrixD B = U.multCU(D, ring).multCU(W, ring);
        MatrixD Check = B.subtract(A, ring);
        System.out.println("Check final SVD: max absolute value = " + Check.max(ring).value);
    }

    public static double[][] getPerfectMatrix() {
        return new double[][] {{1,1,1,1}, {1,1,2,1}, {1,1,1,2}, {1,2,1,1}};
    }

    public static double[][] generateMatrix(int N) {
        Random rand = new Random();
        double[][] res = new double[N][N];

        for(int i =0; i < N ; i++) {
            for (int j = 0; j < N; j++) {
                res[i][j] = rand.nextDouble() * 10.0;
            }
        }

        return res;
    }

    public static void print_matrix(MatrixD m) {
        for (int i = 0; i < m.M.length; ++i) {
            for (int j = 0; j < m.M.length; ++j)
                System.out.printf("%3.4f ", m.M[i][j].doubleValue());
            System.out.println("");
        }
        System.out.println("");
    }

    public static MatrixD[] SVD(MatrixD A, Ring ring) {
        if (A.rowNum() != A.colNum())
            return null;

        MatrixD[] result = biDiagonalize(A, ring);
        MatrixD[] diagonal = diagonalize(result[1], ring);
        result[0] = result[0].multCU(diagonal[0], ring);
        result[1] = diagonal[1];
        result[2] = diagonal[2].multCU(result[2], ring);

        return result;
    }

    public static MatrixD[] biDiagonalize(MatrixD A, Ring ring) {
        MatrixD I = new MatrixD(MatrixS.scalarMatrix(A.M.length, ring.numberONE, ring));
        int operationsAmount = A.M.length - 1;
        int i = 0;
        MatrixD Ai = A;
        MatrixD U = null;
        MatrixD W = null;

        while (i < operationsAmount) {
            VectorS x = Ai.takeColumn(i + 1);
            for (int k = 0; k < i; k++) {
                x.V[k] = ring.numberZERO;
            }

            if (!x.norm2(i + 1, ring).isZero(ring)) {
                Element norm_x = x.norm(ring);
                VectorS u = new VectorS(x.V.clone());

                if (u.V[i].value < 0)
                    u.V[i] = u.V[i].subtract(norm_x, ring);
                else
                    u.V[i] = u.V[i].add(norm_x, ring);

                Element e = x.V[i];
                if (x.V[i].value < 0)
                    e = x.V[i].multiply(new NumberR64(-1), ring);

                Element norm_Squ = norm_x.multiply(norm_x.add(e, ring), ring).multiply(new NumberR64(2), ring);
                VectorS ut = (VectorS) u.transpose(ring);
                NumberR64 norm_div = NumberR64.valueOf(2 / norm_Squ.value);

                if (U == null) {
                    MatrixD uut = (MatrixD) ut.multiply(u, ring);
                    MatrixD nd_uut = (MatrixD) norm_div.multiply(uut, ring);
                    U = I.subtract(nd_uut, ring);
                } else {
                    VectorS Put = (VectorS) U.multiply(ut, ring);
                    VectorS NdPut = (VectorS) norm_div.multiply(Put, ring);
                    MatrixD NdPutu = (MatrixD) NdPut.transpose(ring).multiply(u, ring);
                    U = U.subtract(NdPutu, ring);
                }

                Element dT = u.multiply(Ai, ring);
                Element nd_u = norm_div.multiply(ut, ring);
                Element nd_u_dT = nd_u.multiply(dT, ring);
                Ai = (MatrixD) Ai.subtract(nd_u_dT, ring);
            }

            VectorS x_r = new VectorS(Ai.takeRow(i + 1).V.clone());

            for (int k = 0; k < i + 1; k++) {
                x_r.V[k] = ring.numberZERO;
            }

            if (!x_r.norm2(i + 2, ring).isZero(ring)) {
                Element norm_x_r = x_r.norm(ring);
                VectorS u_r = new VectorS(x_r.V.clone());

                if (u_r.V[i + 1].value < 0)
                    u_r.V[i + 1] = u_r.V[i + 1].subtract(norm_x_r, ring);
                else
                    u_r.V[i + 1] = u_r.V[i + 1].add(norm_x_r, ring);

                Element e_r = x_r.V[i + 1];
                if (x_r.V[i + 1].value < 0)
                    e_r = x_r.V[i + 1].multiply(new NumberR64(-1), ring);
                Element norm_Squ_r = norm_x_r.multiply(norm_x_r.add(e_r, ring), ring).multiply(new NumberR64(2), ring);
                VectorS ut_r = (VectorS) u_r.transpose(ring);
                Element norm_div_r = NumberR64.valueOf(2 / norm_Squ_r.value);

                if (W == null) {
                    MatrixD utu_r = (MatrixD) ut_r.multiply(u_r, ring);
                    MatrixD nd_utu_r = (MatrixD) norm_div_r.multiply(utu_r, ring);
                    W = I.subtract(nd_utu_r, ring);
                } else {
                    VectorS Ndw = (VectorS) norm_div_r.multiply(u_r, ring);
                    VectorS wtQ = (VectorS) u_r.multiply(W, ring);
                    MatrixD NdwwtQ = (MatrixD) Ndw.transpose(ring).multiply(wtQ, ring);
                    W = W.subtract(NdwwtQ, ring);
                }

                VectorS d_r = (VectorS) Ai.multiply(ut_r, ring);
                VectorS nd_d_r = (VectorS) norm_div_r.multiply(d_r, ring);
                MatrixD nd_d_u_r = (MatrixD) nd_d_r.transpose(ring).multiply(u_r, ring);
                Ai = Ai.subtract(nd_d_u_r, ring);
            }

            i++;
        }

        return new MatrixD[]{U, Ai, W};
    }

// ********************* NIKOLAY *************************** 
    public static MatrixD[] diagonalize(MatrixD A, Ring r) {
        int n = A.rowNum();
        MatrixD left;
        MatrixD right;
        MatrixD tmp = A.copy();
        MatrixD L = MatrixD.ONE(n, r);
        MatrixD R = MatrixD.ONE(n, r);

        boolean side = true;
        int iteration = 1;

        while (!com.mathpar.students.OLD.savchenko.Utils.checkSecondDiagonalValues(tmp, n, r)) {
            if (side) {
                for (int i = 0; i < n - 1; ++i) {
                    if (!tmp.getElement(i, i + 1).isZero(r)) {
                        right = getGivensRotationMatrix(n, i, tmp.getElement(i, i), tmp.getElement(i, i + 1), r);
                        R = multLinesRight(R, right.transpose(r), i, r);
                        tmp = multSquareRight(tmp, right, i, r);
                    }
                }
            } else {
                for (int j = 0; j < n - 1; ++j) {
                    if (!tmp.getElement(j + 1, j).isZero(r)) {
                        left = getGivensRotationMatrix(n, j, tmp.getElement(j, j), tmp.getElement(j + 1, j), r);
                        L = multLinesLeft(left, L, j, r);
                        tmp = multSquareLeft(left.transpose(r), tmp, j, r);
                    }
                }
            }
            side = !side;
            iteration++;
        }
        return  new MatrixD[] {L, tmp, R};
    }

      public static MatrixD getGivensRotationMatrix(int n, int i, Element a, Element b, Ring ring) {
        MatrixD G = MatrixD.ONE(n, ring);
        if (b.isZero(ring)) {
            G.M[i][i] = ring.numberONE;
            G.M[i][i + 1] = ring.numberZERO;
            G.M[i + 1][i] = ring.numberZERO;
            G.M[i + 1][i + 1] = ring.numberONE;
        } else {
            Element r = a.pow(2, ring).add(b.pow(2, ring), ring).sqrt(ring);     // Math.sqrt(Math.pow(a, 2) + Math.pow(b, 2));
            Element c = a.divide(r, ring);                                            // a/r;
            Element s = b.divide(r, ring);                               // (-b)/r;
            if (c.isInfinite() || s.isInfinite()) {
                G.M[i][i] = ring.numberONE;
                G.M[i][i + 1] = ring.numberZERO;
                G.M[i + 1][i] = ring.numberZERO;
                G.M[i + 1][i + 1] = ring.numberONE;
            } else {
                G.M[i][i] = c;
                G.M[i][i + 1] = s.negate(ring);
                G.M[i + 1][i] = s;
                G.M[i + 1][i + 1] = c;
            }
        }
        return G;
    }

    public static boolean checkSecondDiagonalValues(MatrixD temp, int n, Ring ring) {
        for (int i = 0; i < (n-1); i++) {
            if (!temp.getElement(i, i+1).isZero(ring) || !temp.getElement(i+1, i).isZero(ring))
                return false;
        }
        return true;
    }


    public static MatrixD getTriangleMatrixNumber64(int n, int mod, Ring ring) {
        MatrixD L = new MatrixD(n, n, mod, ring);
        for (int i = 0; i < n; i++)
            for (int j = 0; j < n; j++)
                if (!(j<=i))
                    L.M[i][j] = ring.numberZERO;
        return L;
    }

    public static MatrixD[] getTwoGivensRotationMatrices(double a, double b, double d, int n, int i, int j, Ring ring) {
        MatrixD left = MatrixD.ONE(n, ring);
        MatrixD right = MatrixD.ONE(n, ring);
        double c = 0d;
        double s = 0d;
        double C = 0d;
        double S = 0d;
        double t = 0d;
        double T = 0d;

        if (a != 0d && d != 0d) {
//            System.out.println("a ≠ 0, d ≠ 0 \n");
            t = ((-1d*(b*b+d*d-a*a)) + Math.sqrt(Math.pow((b*b+d*d-a*a), 2) + 4d*a*a*b*b)) / (2d*a*b);
            T = (-1d/d)*(a*t + b);
        } else if (a != 0d && d == 0d) {
//            System.out.println("a ≠ 0, d = 0 \n");
            T = 0d;
            t = (-b)/a;
        } else if (a == 0d && d != 0d) {
//            System.out.println("a = 0, d ≠ 0 \n");
            t = 0d;
            T = (-b)/d;
        } else {                                    // a = 0 & d = 0
//            System.out.println("a = 0, d = 0 \n");
            c = 1d;
            s = 0d;
            C = 0d;
            S = 1d;
        }

        if (!(a == 0d && d == 0d)) {
            c = Math.sqrt(1d/(1d+t*t));                                 // c = Math.cos(Math.atan(t));
            s = t*c;                                                    // s = Math.sin(Math.atan(t));
            C = Math.sqrt(1d/(1d+T*T));                                 // C = Math.cos(Math.atan(T));
            S = T*C;                                                    // S = Math.sin(Math.atan(T));
        }

        left.M[i][i] = new NumberR64(c);
        left.M[i][j] = new NumberR64(-s);
        left.M[j][i] = new NumberR64(s);
        left.M[j][j] = new NumberR64(c);

        right.M[i][i] = new NumberR64(C);
        right.M[i][j] = new NumberR64(-S);
        right.M[j][i] = new NumberR64(S);
        right.M[j][j] = new NumberR64(C);

        return new MatrixD[]{left, right};
//            long [][] arr = {{1, 0}, {5, 2}};
//            MatrixD test = new MatrixD(arr, ring);
//            MatrixD[] lr = getTwoGivensRotationMatrices(test.getElement(0,0).doubleValue(), test.getElement(1, 0).doubleValue(),
//                    test.getElement(1,1).doubleValue(), 2, 0, 1, ring);
//            System.out.println("L * A * R = \n");
//            System.out.println(lr[0].multiplyMatr(test, ring).multiplyMatr(lr[1], ring).toString());
    }

    public static void removeNonDiagonalValues(MatrixD d, Ring ring) {
        for (int i = 0; i < d.M.length; i++) {
            for (int j = 0; j < d.M[0].length; j++) {
                if (i != j) {
                    d.M[i][j] = ring.numberZERO;
                }
            }
        }
    }

    public static boolean isPowerOfTwo(int number) {
        return number > 0 && ((number & (number - 1)) == 0);
    }

    public static MatrixD getSubMatrix(MatrixD matrix, int start_i, int end_i, int start_j, int end_j) {
        matrix = matrix.copy();
        int rowNum = end_i - start_i + 1;
        int colNum = end_j - start_j + 1;

        Element[][] e = new Element[rowNum][colNum];
        for (int i = start_i; i <= end_i; i++) {
            for (int j = start_j; j <= end_j; j++) {
                e[i-start_i][j-start_j] = matrix.getElement(i, j);
            }
        }

        return new MatrixD(e, 0);
    }

    public static MatrixD insertMatrixToMatrix(MatrixD matrix, MatrixD block, int i_start, int j_start){
        block = block.copy();
        MatrixD result = matrix.copy();

        for (int i = 0; i < block.rowNum(); i++) {
            for (int j = 0; j < block.colNum(); j++) {
                result.M[i+i_start][j+j_start] = block.getElement(i, j);
            }
        }

        return result;
    }

    public static void readBlock(MatrixD matrix, int iOffset, int jOffset, Element[][] elements, int h) {
        for (int i = iOffset; i < (h+iOffset); i++) {
            for (int j = jOffset; j < (h+jOffset); j++) {
                elements[i-iOffset][j-jOffset] = matrix.getElement(i, j);
            }
        }
    }

    public static MatrixD getBlock(MatrixD input, int block) {
        int rowNum = input.rowNum();
        int colNum = input.colNum();

        int i_start = 0, i_end = 0, j_start = 0, j_end = 0;

        if (block == 1) {
            i_start = rowNum / 4;
            i_end = rowNum - (rowNum / 4) - 1;
            j_start = 0;
            j_end = (colNum / 2) - 1;
        } else if (block == 2) {
            i_start = 0;
            i_end = (rowNum / 2) - 1;
            j_start = 0;
            j_end = (colNum / 2) - 1;
        } else if (block == 3) {
            i_start = rowNum / 2;
            i_end = rowNum - 1;
            j_start = colNum / 2;
            j_end = colNum - 1;
        } else if (block == 4) {
            i_start = rowNum / 4;
            i_end = rowNum - (rowNum / 4) - 1;
            j_start = colNum / 2;
            j_end = colNum - 1;
        }

        return getSubMatrix(input, i_start, i_end, j_start, j_end);
    }

    public static MatrixD block4_(MatrixD input, char b)  {
        if ((input.rowNum() != input.colNum()) || 
                !isPowerOfTwo(input.rowNum())) 
        {  return null;}
 

        MatrixD matrix = input.copy();
        int n = matrix.rowNum();
        int h = n/2;

        if (n == 1) {
            return matrix;
        } else {
            switch (b) {
                case 'A': return getSubMatrix(matrix, 0, h-1, 0, h-1);
                case 'B': return getSubMatrix(matrix, 0, h-1, h, n-1);
                case 'C': return getSubMatrix(matrix, h, n-1, 0, h-1);
                case 'D': return getSubMatrix(matrix, h, n-1, h, n-1);
                default: return matrix;
            }
        }
    }

    // В результате у матрицы B поменяется только ряд x и ряд y
    public static MatrixD multSquareRight(MatrixD a, MatrixD b, int pos, Ring ring) {

        Element c11 = a.M[pos][pos].multiply(b.M[pos][pos], ring);
        c11 = c11.add(a.M[pos][pos + 1].multiply(b.M[pos + 1][pos], ring), ring);

        Element c12 = a.M[pos][pos].multiply(b.M[pos][pos + 1], ring);
        c12 = c12.add(a.M[pos][pos + 1].multiply(b.M[pos + 1][pos + 1], ring), ring);

        Element c21 = a.M[pos + 1][pos].multiply(b.M[pos][pos], ring);
        c21 = c21.add(a.M[pos + 1][pos + 1].multiply(b.M[pos + 1][pos], ring), ring);

        Element c22 = a.M[pos + 1][pos].multiply(b.M[pos][pos + 1], ring);
        c22 = c22.add(a.M[pos + 1][pos + 1].multiply(b.M[pos + 1][pos + 1], ring), ring);

        a.M[pos][pos] =  c11;
        a.M[pos][pos + 1] =  c12;
        a.M[pos + 1][pos] =  c21;
        a.M[pos + 1][pos + 1] =  c22;

        return a;
    }

    public static MatrixD multSquareLeft(MatrixD a, MatrixD b, int pos, Ring ring) {

        Element c11 = a.M[pos][pos].multiply(b.M[pos][pos], ring);
        c11 = c11.add(a.M[pos][pos + 1].multiply(b.M[pos + 1][pos], ring), ring);

        Element c12 = a.M[pos][pos].multiply(b.M[pos][pos + 1], ring);
        c12 = c12.add(a.M[pos][pos + 1].multiply(b.M[pos + 1][pos + 1], ring), ring);

        Element c21 = a.M[pos + 1][pos].multiply(b.M[pos][pos], ring);
        c21 = c21.add(a.M[pos + 1][pos + 1].multiply(b.M[pos + 1][pos], ring), ring);

        Element c22 = a.M[pos + 1][pos].multiply(b.M[pos][pos + 1], ring);
        c22 = c22.add(a.M[pos + 1][pos + 1].multiply(b.M[pos + 1][pos + 1], ring), ring);

        b.M[pos][pos] =  c11;
        b.M[pos][pos + 1] =  c12;
        b.M[pos + 1][pos] =  c21;
        b.M[pos + 1][pos + 1] =  c22;

        return b;
    }

    public static MatrixD multLinesRight(MatrixD a, MatrixD b, int line, Ring ring) {
        int size = a.rowNum();
        Element[][] r = new Element[2][size];

        for (int i = 0; i < 2; ++i) {
            if (line + i < size) {

                for (int k = 0; k < size; ++k) {

                    Element v = a.M[k][0].multiply(b.M[0][line + i], ring);

                    for (int j = 1; j < size; ++j) {
                        v = v.add(a.M[k][j].multiply(b.M[j][line + i], ring), ring);
                    }

                    r[i][k] = v;
                }
            }
        }

        for (int i = 0; i < 2; ++i) {
            if (line + i < size) {
                for (int j = 0; j < size; ++j) {
                    a.M[j][i + line] = r[i][j];
                }
            }
        }
        return a;
    }

    public static MatrixD multLinesLeft(MatrixD a, MatrixD b, int line, Ring ring) {
        int size = a.rowNum();
        Element[][] r = new Element[2][size];

        for (int i = 0; i < 2; ++i) {
            if (line + i < size) {

                for (int k = 0; k < size; ++k) {

                    Element v = a.M[line + i][0].multiply(b.M[0][k], ring);

                    for (int j = 1; j < size; ++j) {
                        v = v.add(a.M[line + i][j].multiply(b.M[j][k], ring), ring);
                    }

                    r[i][k] = v;
                }
            }
        }

        for (int i = 0; i < 2; ++i) {
            if (line + i < size) {
                for (int j = 0; j < size; ++j) {
                    b.M[i + line][j] = r[i][j];
                }
            }
        }
        return b;
    }
        
}
