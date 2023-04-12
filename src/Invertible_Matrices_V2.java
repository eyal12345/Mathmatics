
public class Invertible_Matrices_V2 {

    /////////////////////////////////////////////// Print Methods /////////////////////////////////////////////////
    // display current status of the matrices M and InvM each time of iteration on an element
    public static void Print_Status_Matrices(float[][] M, float[][] invM) {
        for (int i = 0; i < M.length; i++) {
            for (int j = 0; j < M[0].length; j++) {
                if (j == M[0].length - 1) {
                    if ((Math.round(M[i][j] * 1000000.0) / 1000000.0) % 1 == 0) {
                        System.out.print((int) (Math.round(M[i][j] * 1000000.0) / 1000000.0) + " ");
                    } else {
                        System.out.print(Math.round(M[i][j] * 1000000.0) / 1000000.0 + " ");
                    }
                } else {
                    if ((Math.round(M[i][j] * 1000000.0) / 1000000.0) % 1 == 0) {
                        System.out.print((int) (Math.round(M[i][j] * 1000000.0) / 1000000.0) + " ,");
                    } else {
                        System.out.print(Math.round(M[i][j] * 1000000.0) / 1000000.0 + " ,");
                    }
                }
            }
            System.out.print("| ");
            for (int j = 0; j < invM[0].length; j++) {
                if (j == invM[0].length - 1) {
                    if ((Math.round(invM[i][j] * 1000000.0) / 1000000.0) % 1 == 0) {
                        System.out.print((int) (Math.round(invM[i][j] * 1000000.0) / 1000000.0) + " ");
                    } else {
                        System.out.print(Math.round(invM[i][j] * 1000000.0) / 1000000.0 + " ");
                    }
                } else {
                    if ((Math.round(invM[i][j] * 1000000.0) / 1000000.0) % 1 == 0) {
                        System.out.print((int) (Math.round(invM[i][j] * 1000000.0) / 1000000.0) + " ,");
                    } else {
                        System.out.print(Math.round(invM[i][j] * 1000000.0) / 1000000.0 + " ,");
                    }
                }
            }
            System.out.println();
        }
        System.out.println();
    }

    // display a matrix each current status
    public static void Print_Matrix(float[][] M) {
        for (int i = 0; i < M.length; i++) {
            for (int j = 0; j < M[0].length; j++) {
                if (j == M[0].length - 1) {
                    if ((Math.round(M[i][j] * 1000000.0) / 1000000.0) % 1 == 0) {
                        System.out.print((int) (Math.round(M[i][j] * 1000000.0) / 1000000.0) + " ");
                    } else {
                        System.out.print(Math.round(M[i][j] * 1000000.0) / 1000000.0 + " ");
                    }
                } else {
                    if ((Math.round(M[i][j] * 1000000.0) / 1000000.0) % 1 == 0) {
                        System.out.print((int) (Math.round(M[i][j] * 1000000.0) / 1000000.0) + " ,");
                    } else {
                        System.out.print(Math.round(M[i][j] * 1000000.0) / 1000000.0 + " ,");
                    }
                }
            }
            System.out.println();
        }
        System.out.println();
    }

    // display user interface by selection format for solution
    public static void User_Menu_Solution() {
        System.out.println("choose the character of format to representation of solution:");
        System.out.println("d. decimal");
        System.out.println("r. rational");
    }

    /////////////////////////////////////////// Auxiliary Operations /////////////////////////////////////////////
    // calculate multiplication between two matrices
    public static float[][] Mul_Mats(float[][] M1, float[][] M2) throws Exception {
        if (M1[0].length == M2.length) {
            float[][] M = new float[M1.length][M2[0].length];
            for (int i = 0; i < M1.length; i++) {
                for (int j = 0; j < M2[0].length; j++) {
                    for (int k = 0; k < M2.length; k++) {
                        M[i][j] += M1[i][k] * M2[k][j];
                    }
                }
            }
            return M;
        } else {
            throw new Exception("has no match between matrices dimensions");
        }
    }

    // replace between two rows in a matrix
    public static void Retreat_Rows_Matrix(float[][] M, int r1, int r2) {
        for (int j = 0; j < M[0].length; j++) {
            float t = M[r1][j];
            M[r1][j] = M[r2][j];
            M[r2][j] = t;
        }
    }

    // check if the specific row in the matrix is a unit vector
    public static boolean Is_Unit_Vector(float[][] M, int k) {
        int n = M.length;
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n && M[k][i] != 0; j++) {
                if (i != j && M[k][j] != 0) {
                    return false;
                }
            }
        }
        return true;
    }

    // check if a matrix is a unit matrix
    public static boolean Is_Unit_Matrix(float[][] M) {
        int n = M.length;
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                if (M[i][i] != 1 || (i != j && M[i][j] != 0)) {
                    return false;
                }
            }
        }
        return true;
    }

    // create a unit matrix with "n*n" size
    public static float[][] Unit_Matrix(int n) {
        float[][] I = new float[n][n];
        for (int i = 0; i < n; i++) {
            I[i][i] = 1;
        }
        return I;
    }

    // get the index starting from the specific column in the matrix which are him value not equal to 0
    public static int Get_Index_UnZero_Value(float[][] M, int k) {
        int n = M.length;
        for (int i = k + 1; i < n; i++) {
            if (M[i][k] != 0) {
                return i % n;
            }
        }
        return -1;
    }

    // duplicate the matrix values into a new matrix
    public static float[][] Copy_Matrix(float[][] M) {
        float[][] copyA = new float[M.length][M.length];
        for (int i = 0; i < M.length; i++) {
            for (int j = 0; j < M[0].length; j++) {
                copyA[i][j] = M[i][j];
            }
        }
        return copyA;
    }

    // convert a value to a format of rational number
    public static String convertDecimalToFraction(float x){
        if (x < 0) {
            return "-" + convertDecimalToFraction(-x);
        }
        float tolerance = (float) 1.0E-5, h1 = 1, h2 = 0, k1 = 0, k2 = 1, b = x;
        do {
            float a = (float) Math.floor(b);
            float aux = h1;
            h1 = a*h1 + h2;
            h2 = aux;
            aux = k1;
            k1 = a*k1 + k2;
            k2 = aux;
            b = 1 / (b - a);
        } while (Math.abs(x - h1/k1) > x*tolerance);
        if (k1 != 1) {
            return String.valueOf((int) h1 + "/" + (int) k1);
        } else {
            return String.valueOf((int) h1);
        }
    }

    //////////////////////////////////////////// Elementary Actions //////////////////////////////////////////////
    // show elementary actions for sum between rows in the matrices
    public static void Sum_Elementary_Action(float k, int j, int i) {
        if (k != 0) {
            int r = j + 1, c = i + 1;
            if (k > 0) {
                if (k % 1 == 0) {
                    if (k == 1) {
                        System.out.println("R" + r + " --> R" + r + " - R" + c + "\n");
                    } else {
                        System.out.println("R" + r + " --> R" + r + " - " + (int) k + "*R" + c + "\n");
                    }
                } else {
                    System.out.println("R" + r + " --> R" + r + " - " + k + "*R" + c + "\n");
                }
            } else {
                if (k % 1 == 0) {
                    if (k == -1) {
                        System.out.println("R" + r + " --> R" + r + " + R" + c + "\n");
                    } else {
                        System.out.println("R" + r + " --> R" + r + " + " + (int) (-k) + "*R" + c + "\n");
                    }
                } else {
                    System.out.println("R" + r + " --> R" + r + " + " + (-k) + "*R" + c + "\n");
                }
            }
        }
    }

    // show elementary actions for multiplication of a row in the matrices
    public static void Mul_Elementary_Action(float k, int j) {
        if (k != 1) {
            int r = j + 1;
            if (k % 1 == 0) {
                if (k == -1) {
                    System.out.println("R" + r + " --> - R" + r + "\n");
                } else {
                    System.out.println("R" + r + " --> " + (int) k + "*R" + r + "\n");
                }
            } else {
                System.out.println("R" + r + " --> " + k + "*R" + r + "\n");
            }
        }
    }

    /////////////////////////////////////////// Methods to Solution /////////////////////////////////////////////
    // invert the M matrix in parallel ranking by a multiplication of M in elementary matrix each iteration (iterative method)
    public static float[][] Invertible(float[][] M) throws Exception {
        int n = M.length;
        float[][] InvM = Unit_Matrix(n);
        int i = 0 ,j = 0;
        Print_Status_Matrices(M,InvM);
        while (!Is_Unit_Matrix(M)) {
            float[][] E = Unit_Matrix(n);
            if (M[i][i] == 0) {
                int k = Get_Index_UnZero_Value(M,i);
                System.out.println("R" + (i + 1) + " <--> R" + (k + 1) + "\n");
                Retreat_Rows_Matrix(E,i,k);
                M = Mul_Mats(E,M);
                InvM = Mul_Mats(E,InvM);
                Print_Status_Matrices(M,InvM);
            } else {
                if (i != j) {
                    E[j][i] -= (M[j][i] / M[i][i]);
                    Sum_Elementary_Action(-E[j][i],j,i);
                    M = Mul_Mats(E,M);
                    InvM = Mul_Mats(E,InvM);
                    M[j][i] = 0;
                    if (E[j][i] != 0) {
                        Print_Status_Matrices(M,InvM);
                    }
                } else if (Is_Unit_Vector(M,j)) {
                    E[j][j] = 1 / M[j][j];
                    Mul_Elementary_Action(E[j][j],j);
                    M = Mul_Mats(E,M);
                    InvM = Mul_Mats(E,InvM);
                    M[j][j] = 1;
                    if (E[j][j] != 1) {
                        Print_Status_Matrices(M,InvM);
                    }
                }
                if (j == n - 1) {
                    i = (i + 1) % n;
                }
                j = (j + 1) % n;
            }
        }
        return InvM;
    }

    // invert the M matrix in parallel ranking by a multiplication of M in elementary matrix each iteration (recursive method)
    public static float[][] Invertible_Rec(float[][] M, float[][] InvM, int i, int j) throws Exception {
        int n = M.length;
        if (Is_Unit_Matrix(M)) {
            return InvM;
        } else {
            float[][] E = Unit_Matrix(n);
            if (M[i][i] == 0) {
                int k = Get_Index_UnZero_Value(M,i);
                System.out.println("R" + (i + 1) + " <--> R" + (k + 1) + "\n");
                Retreat_Rows_Matrix(E,i,k);
                M = Mul_Mats(E,M);
                InvM = Mul_Mats(E,InvM);
                Print_Status_Matrices(M,InvM);
            } else {
                if (i != j) {
                    E[j][i] -= (M[j][i] / M[i][i]);
                    Sum_Elementary_Action(-E[j][i],j,i);
                    M = Mul_Mats(E,M);
                    InvM = Mul_Mats(E,InvM);
                    M[j][i] = 0;
                    if (E[j][i] != 0)
                        Print_Status_Matrices(M,InvM);
                } else if (Is_Unit_Vector(M,j)) {
                    E[j][j] = 1 / M[j][j];
                    Mul_Elementary_Action(E[j][j],j);
                    M = Mul_Mats(E,M);
                    InvM = Mul_Mats(E,InvM);
                    M[j][j] = 1;
                    if (E[j][j] != 1)
                        Print_Status_Matrices(M,InvM);
                }
                if (j == n - 1)
                    i = (i + 1) % n;
                j = (j + 1) % n;
            }
            return Invertible_Rec(M,InvM,i,j);
        }
    }

    ////////////////////////////////////////////// Run Progress ////////////////////////////////////////////////
    public static void main(String[] args) {
        float[][] A1 = {{2,1,-1},{-3,-1,2},{-2,1,2}};
        float[][] A2 = {{1,-1,-2},{2,-3,-5},{-1,3,5}};
        float[][] A3 = {{1,1,5},{1,2,7},{2,-1,4}};
        float[][] A4 = {{1,1,5},{-2,-2,-10},{6,6,30}};
        float[][] B = {{-2,2,2,-1,1},{-4,4,4,4,3},{2,3,2,3,2},{-3,-1,1,2,2},{5,5,3,5,5}};
        float[][] Inv;
        try {
            float[][] M = Copy_Matrix(A1);
//            Inv = Invertible(M);
            Print_Status_Matrices(M,Unit_Matrix(M.length));
            Inv = Invertible_Rec(M,Unit_Matrix(M.length) ,0 ,0);
            Print_Matrix(Inv);
        } catch (Exception e) {
            e.printStackTrace();
        }
    }
}
