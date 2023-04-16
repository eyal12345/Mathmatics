import java.util.Scanner;

public class Matrices_Decomposition {

    /////////////////////////////////////////////// Print Methods /////////////////////////////////////////////////
    // display the matrix M in the matrices format
    public static void Display_Exercise(float[][] M) {
        int n = M.length;
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                if ((Math.round(M[i][j] * 1000.0) / 1000.0) % 1 == 0) {
                    System.out.print((int) (Math.round(M[i][j] * 1000.0) / 1000.0));
                } if (j != n - 1) {
                    System.out.print(" ,");
                }
            }
            System.out.println();
        }
        System.out.println();
    }

    // display a matrix each current status
    public static void Print_Matrix(float[][] M, String fn) {
        int n = M.length;
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                if ((Math.round(M[i][j] * 1000.0) / 1000.0) % 1 == 0)
                    System.out.print((int)(Math.round(M[i][j] * 1000.0) / 1000.0));
                else if (fn.equals("d")) {
                    System.out.print(Math.round(M[i][j] * 1000.0) / 1000.0);
                } else if (fn.equals("r")) {
                    System.out.print(convertDecimalToFraction(M[i][j]));
                } if (j != n - 1) {
                    System.out.print(" ,");
                }
            }
            System.out.println();
        }
        System.out.println();
    }

    // display user interface by selection method for receive matrices
    public static void User_Menu_System_receive() {
        System.out.println("choose number method to solution:");
        System.out.println("1. LU decomposition by L and U multiplication (first method)");
        System.out.println("2. LL' decomposition by L and L' multiplication (first method)");
        System.out.println("3. LDL' decomposition by L, D and L' multiplication (first method)");
        System.out.println("4. LU decomposition by L and U multiplication (second method)");
        System.out.println("5. LL' decomposition by L and L' multiplication (second method)");
        System.out.println("6. LDL' decomposition by L, D and L' multiplication (second method)");
    }

    // display user interface by selection method for decompose matrices
    public static void User_Menu_System_decompose() {
        System.out.println("choose number method to solution:");
        System.out.println("1. get L and U matrices by decomposition of M (first method)");
        System.out.println("2. get L and L' matrices by decomposition of M (first method)");
        System.out.println("3. get L, D and L' matrices by decomposition of M (first method)");
        System.out.println("4. get L and U matrices by decomposition of M (second method)");
        System.out.println("5. get L and L' matrices by decomposition of M (second method)");
        System.out.println("6. get L, D and L' matrices by decomposition of M (second method)");
    }

    // display user interface by selection format for solution
    public static void User_Menu_Solution() {
        System.out.println("choose a character of format to representation of solution:");
        System.out.println("d. decimal");
        System.out.println("r. rational");
    }

    // display user interface by direction of the solution
    public static void User_Menu_Direction() {
        System.out.println("choose a character for direction of the solution:");
        System.out.println("d. get the components by M decomposition");
        System.out.println("r. get M by multiplication of the components");
    }

    /////////////////////////////////////////// Auxiliary Operations /////////////////////////////////////////////
    // check if the matrix is a square matrix
    public static boolean Is_Square_Matrix(float[][] M) {
        return M.length == M[0].length;
    }

    // check if the matrix is an upper triangular
    public static boolean Is_Upper_Triangular(float[][] M) {
        int n = M.length;
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < i; j++) {
                if (M[i][j] != 0) {
                    return false;
                }
            }
        }
        return true;
    }

    // check if the matrix is a lower triangular
    public static boolean Is_Lower_Triangular(float[][] M) {
        int n = M.length;
        for (int i = 0; i < n - 1; i++) {
            for (int j = i + 1; j < n; j++) {
                if (M[i][j] != 0) {
                    return false;
                }
            }
        }
        return true;
    }

    // check if the matrix is a symmetrical matrix
    public static boolean Is_Symmetrical_Matrix(float[][] M) {
        int n = M.length;
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                if (M[i][j] != M[j][i]) {
                    return false;
                }
            }
        }
        return true;
    }

    // check if all the values in the main diagonal which are positive
    public static boolean Is_Values_Positives(float[][] M) {
        int n = M.length;
        for (int i = 0; i < n; i++) {
            if (M[i][i] < 0) {
                return false;
            }
        }
        return true;
    }

    // check if all the values in the main diagonal which are equals to 1
    public static boolean Is_One_Slant(float[][] M) {
        int n = M.length;
        for (int i = 0; i < n; i++) {
            if (M[i][i] != 1) {
                return false;
            }
        }
        return true;
    }

    // calculate multiplication between two matrices
    public static float[][] Mul_Mats(float[][] M1, float[][] M2) {
        int n = M1.length;
        float[][] M = new float[n][n];
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                for (int k = 0; k < n; k++) {
                    M[i][j] += M1[i][k] * M2[k][j];
                }
            }
        }
        return M;
    }

    // replace between two rows in a matrix
    public static void Retreat_Rows_Matrix(float[][] M, int r1, int r2) {
        int n = M.length;
        for (int j = 0; j < n; j++) {
            float t = M[r1][j];
            M[r1][j] = M[r2][j];
            M[r2][j] = t;
        }
    }

    // calculate transpose matrix of a matrix
    public static float[][] Transpose(float[][] M) {
        int n = M.length;
        float[][] MT = new float[n][n];
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                MT[j][i] = M[i][j];
            }
        }
        return MT;
    }

    // check if a matrix which is transpose of another
    public static boolean isTranspose(float[][] M1, float[][] M2) {
        int n = M1.length;
        boolean flag = false;
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                if (M2[j][i] == M1[i][j]) {
                    flag = true;
                } else {
                    return false;
                }
            }
        }
        return flag;
    }

    // convert a value to a format of rational number
    public static String convertDecimalToFraction(float x) {
        if (x < 0) {
            return "-" + convertDecimalToFraction(-x);
        } else {
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
    }

    //////////////////////////////////////////// Elementary Actions //////////////////////////////////////////////
    // show elementary actions for sum between rows in the matrices
    public static void Sum_Elementary_Action(float k, int j, int i, String fn) {
        if (k != 0) {
            int r = j + 1, c = i + 1;
            if (k > 0) {
                if (k % 1 == 0) {
                    if (k == 1) {
                        System.out.println("R" + r + " --> R" + r + " - R" + c + "\n");
                    } else {
                        System.out.println("R" + r + " --> R" + r + " - " + (int) k + "*R" + c + "\n");
                    }
                } else if (fn.equals("d")) {
                    System.out.println("R" + r + " --> R" + r + " - " + k + "*R" + c + "\n");
                } else if (fn.equals("r")) {
                    System.out.println("R" + r + " --> R" + r + " - " + convertDecimalToFraction(k) + "*R" + c + "\n");
                }
            } else {
                if (k % 1 == 0) {
                    if (k == -1) {
                        System.out.println("R" + r + " --> R" + r + " + R" + c + "\n");
                    } else {
                        System.out.println("R" + r + " --> R" + r + " + " + (int) ((-1) * k) + "*R" + c + "\n");
                    }
                } else if (fn.equals("d")) {
                    System.out.println("R" + r + " --> R" + r + " + " + (-k) + "*R" + c + "\n");
                } else if (fn.equals("r")) {
                    System.out.println("R" + r + " --> R" + r + " + " + convertDecimalToFraction(-k) + "*R" + c + "\n");
                }
            }
        }
    }

    // show elementary actions for multiplication of a row in the matrices
    public static void Mul_Elementary_Action(float k, int j, String fn) {
        if (k != 1) {
            int r = j + 1;
            if (k % 1 == 0) {
                if (k == -1) {
                    System.out.println("R" + r + " --> - R" + r + "\n");
                } else {
                    System.out.println("R" + r + " --> " + (int) k + "*R" + r + "\n");
                }
            } else if (fn.equals("d")) {
                System.out.println("R" + r + " --> " + k + "*R" + r + "\n");
            } else if (fn.equals("r")) {
                System.out.println("R" + r + " --> " + convertDecimalToFraction(k) + "*R" + r + "\n");
            }
        }
    }

    /////////////////////////////////////// Methods to Solution (Receive M) /////////////////////////////////////////
    // get the LU decomposition by multiplication of L and U (first algorithm)
    public static float[][] From_LU_To_M_V1(float[][] L, float[][] U, String fn) throws Exception {
        System.out.println("U = ");
        Print_Matrix(U,fn);
        if (Is_Square_Matrix(L) && Is_Square_Matrix(U) && Is_Lower_Triangular(L) && Is_Upper_Triangular(U) && Is_One_Slant(L)) {
            return Mul_Mats(L,U);
        } else {
            throw new Exception("Not all conditions are held");
        }
    }

    // get the LL' decomposition by multiplication of L and L' (first algorithm)
    public static float[][] From_LLT_To_M_V1(float[][] L, float[][] M2, String fn) throws Exception {
        System.out.println("L' = ");
        Print_Matrix(M2,fn);
        if (Is_Square_Matrix(L) && Is_Lower_Triangular(L) && isTranspose(L,M2)) {
            return Mul_Mats(L,Transpose(L));
        } else {
            throw new Exception("Not all conditions are held");
        }
    }

    // get the LDL' decomposition by multiplication of L, D and L' (first algorithm)
    public static float[][] From_LDLT_To_M_V1(float[][] L, float[][] D, String fn) throws Exception {
        System.out.println("D = ");
        Print_Matrix(D,fn);
        System.out.println("L' = ");
        Print_Matrix(Transpose(L),fn);
        if (Is_Square_Matrix(L) && Is_Square_Matrix(D) && Is_Lower_Triangular(L) && Is_One_Slant(L) && Is_Lower_Triangular(D) && Is_Upper_Triangular(D)) {
            return Mul_Mats(Mul_Mats(L,D),Transpose(L));
        } else {
            throw new Exception("Not all conditions are held");
        }
    }

    // get the LU decomposition by multiplication of L and U (second algorithm)
    public static float[][] From_LU_To_M_V2(float[][] L, float[][] U, String fn) throws Exception {
        System.out.println("U = ");
        Print_Matrix(U,fn);
        if (Is_Square_Matrix(L) && Is_Square_Matrix(U) && Is_Lower_Triangular(L) && Is_Upper_Triangular(U) && Is_One_Slant(L)) {
            int n = L.length;
            float[][] M = new float[n][n];
            for (int i = 0; i < n; i++) {
                for (int j = 0; j < n; j++) {
                    int m = Math.min(i,j);
                    for (int k = 0; k <= m; k++) {
                        M[i][j] += L[i][k] * U[k][j];
                    }
                }
            }
            return M;
        } else {
            throw new Exception("Not all conditions are held");
        }
    }

    // get the LL' decomposition by multiplication of L and L' (second algorithm)
    public static float[][] From_LLT_To_M_V2(float[][] L, float[][] M2, String fn) throws Exception {
        System.out.println("L' = ");
        Print_Matrix(M2,fn);
        if (Is_Square_Matrix(L) && Is_Lower_Triangular(L) && isTranspose(L,M2)) {
            int n = L.length;
            float[][] M = new float[n][n];
            for (int i = 0; i < n; i++) {
                for (int j = 0; j < n; j++) {
                    int m = Math.min(i,j);
                    for (int k = 0; k <= m; k++) {
                        M[i][j] += L[i][k] * L[j][k];
                    }
                }
            }
            return M;
        } else {
            throw new Exception("Not all conditions are held");
        }
    }

    // get the LDL' decomposition by multiplication of L, D and L' (second algorithm)
    public static float[][] From_LDLT_To_M_V2(float[][] L, float[][] D, String fn) throws Exception {
        System.out.println("D = ");
        Print_Matrix(D,fn);
        System.out.println("L' = ");
        Print_Matrix(Transpose(L),fn);
        if (Is_Square_Matrix(L) && Is_Lower_Triangular(L) && Is_One_Slant(L) && Is_Square_Matrix(D) && Is_Lower_Triangular(D) && Is_Upper_Triangular(D)) {
            int n = L.length;
            float[][] M = new float[n][n];
            for (int i = 0; i < n; i++) {
                for (int j = 0; j < n; j++) {
                    int m = Math.min(i,j);
                    for (int k = 0; k <= m; k++) {
                        M[i][j] += L[i][k] * L[j][k] * D[k][k];
                    }
                }
            }
            return M;
        } else {
            throw new Exception("Not all conditions are held");
        }
    }

    ////////////////////////////////////// Methods to Solution (Decompose M) ////////////////////////////////////////
    // get the LU decomposition of M (first algorithm)
    public static float[][] From_M_To_LU_V1(float[][] M, String fn) throws Exception {
        if (Is_Square_Matrix(M)) {
            int n = M.length;
            float[][] L = new float[n][n];
            for (int i = 0; i < n; i++) {
                L[i][i] = 1;
                for (int j = i + 1; j < n; j++) {
                    if (M[i][i] == 0) {
                        System.out.println("R" + (i + 1) + " <--> R" + (j + 1) + "\n");
                        Retreat_Rows_Matrix(M,i,j);
                        Print_Matrix(M,fn);
                    }
                    L[j][i] = M[j][i] / M[i][i];
                    Sum_Elementary_Action(L[j][i],j,i,fn);
                    for (int k = 0; k < n; k++) {
                        M[j][k] -= M[i][k] * L[j][i];
                    }
                    M[j][i] = 0;
                    if (L[j][i] != 0) {
                        Print_Matrix(M,fn);
                    }
                }
            }
            System.out.println("L = ");
            Print_Matrix(L,fn);
            return M;
        } else {
            throw new Exception("The matrix is not square");
        }
    }

    // get the LL' decomposition of M (first algorithm)
    public static float[][] From_M_To_LLT_V1(float[][] M, String fn) throws Exception {
        if (Is_Square_Matrix(M) && Is_Symmetrical_Matrix(M) && Is_Values_Positives(M)) {
            int n = M.length;
            float[][] L = new float[n][n];
            for (int i = 0; i < n; i++) {
                for (int j = i + 1; j < n; j++) {
                    float c = M[j][i] / M[i][i];
                    Sum_Elementary_Action(c,j,i,fn);
                    for (int k = 0; k < n; k++) {
                        M[j][k] -= M[i][k] * c;
                    }
                    M[j][i] = 0;
                    if (c != 0) {
                        Print_Matrix(M,fn);
                    }
                }
                float c = (float) Math.sqrt(M[i][i]);
                Mul_Elementary_Action(1 / c,i,fn);
                for (int k = i; k < n; k++) {
                    M[i][k] /= c;
                    L[k][i] = M[i][k];
                }
                if (c != 1) {
                    Print_Matrix(M,fn);
                }
            }
            System.out.println("L = ");
            Print_Matrix(L,fn);
            return M;
        } else {
            throw new Exception("Not all conditions are held");
        }
    }

    // get the LDL' decomposition of M (first algorithm)
    public static float[][] From_M_To_LDLT_V1(float[][] M, String fn) throws Exception {
        if (Is_Square_Matrix(M) && Is_Symmetrical_Matrix(M)) {
            int n = M.length;
            float[][] L = new float[n][n];
            float[][] D = new float[n][n];
            for (int i = 0; i < n; i++) {
                for (int j = i + 1; j < n; j++) {
                    float c = M[j][i] / M[i][i];
                    Sum_Elementary_Action(c,j,i,fn);
                    for (int k = 0; k < n; k++) {
                        M[j][k] -= M[i][k] * c;
                    }
                    M[j][i] = 0;
                    if (c != 0) {
                        Print_Matrix(M,fn);
                    }
                }
                D[i][i] = M[i][i];
                Mul_Elementary_Action(1 / D[i][i],i,fn);
                for (int k = i; k < n; k++) {
                    M[i][k] /= D[i][i];
                    L[k][i] = M[i][k];
                }
                if (D[i][i] != 1) {
                    Print_Matrix(M,fn);
                }
            }
            System.out.println("L = ");
            Print_Matrix(L,fn);
            System.out.println("D = ");
            Print_Matrix(D,fn);
            return M;
        } else {
            throw new Exception("Not all conditions are held");
        }
    }

    // get the LU decomposition of M (second algorithm)
    public static float[][] From_M_To_LU_V2(float[][] M, String fn) throws Exception {
        if (Is_Square_Matrix(M)) {
            int n = M.length;
            float[][] L = new float[n][n];
            float[][] U = new float[n][n];
            for (int i = 0; i < n; i++) {
                L[i][i] = 1;
                for (int j = 0; j < n; j++) {
                    int m = Math.min(i,j);
                    for (int k = 0 ;k < m ;k++) {
                        M[i][j] -= L[i][k] * U[k][j];
                    }
                    if (i <= j) {
                        U[i][j] = M[i][j] / L[i][i];
                    } else {
                        L[i][j] = M[i][j] / U[j][j];
                    }
                }
            }
            System.out.println("L = ");
            Print_Matrix(L,fn);
            return M;
        } else {
            throw new Exception("The matrix is not square");
        }
    }

    // get the LL' decomposition of M (second algorithm)
    public static float[][] From_M_To_LLT_V2(float[][] M, String fn) throws Exception {
        if (Is_Square_Matrix(M) && Is_Symmetrical_Matrix(M) && Is_Values_Positives(M)) {
            int n = M.length;
            float[][] L = new float[n][n];
            for (int j = 0; j < n; j++) {
                for (int i = j; i < n; i++) {
                    for (int k = 0; k < j; k++) {
                        M[i][j] -= L[i][k] * L[j][k];
                    }
                    L[i][j] = M[i][j] / (float) Math.sqrt(M[j][j]);
                }
            }
            System.out.println("L = ");
            Print_Matrix(L,fn);
            return M;
        } else {
            throw new Exception("Not all conditions are held");
        }
    }

    // get the LDL' decomposition of M (second algorithm)
    public static float[][] From_M_To_LDLT_V2(float[][] M, String fn) throws Exception {
        if (Is_Square_Matrix(M) && Is_Symmetrical_Matrix(M)) {
            int n = M.length;
            float[][] L = new float[n][n];
            float[][] D = new float[n][n];
            for (int j = 0; j < n; j++) {
                L[j][j] = 1;
                for (int i = j; i < n; i++) {
                    for (int k = 0; k < j; k++) {
                        M[i][j] -= L[i][k] * L[j][k] * D[k][k];
                    }
                    L[i][j] = M[i][j] / M[j][j];
                }
                D[j][j] = M[j][j];
            }
            System.out.println("L = ");
            Print_Matrix(L,fn);
            System.out.println("D = ");
            Print_Matrix(D,fn);
            return M;
        } else {
            throw new Exception("Not all conditions are held");
        }
    }

    ///////////////////////////////////////////// User Interface ///////////////////////////////////////////////
    // get the matrix by a multiplication components
    public static void Receive_Matrix(float[][] L, float[][] M2) throws Exception {
        int m1 = L.length, n1 = L[0].length, m2 = M2.length, n2 = M2[0].length;
        if (m1 == m2 && n1 == n2) {
            System.out.println("receive the M matrix from the two next matrices (" + n1 + "*" + n1 + " size):");
            Display_Exercise(L);
            Display_Exercise(M2);
            Scanner sc = new Scanner(System.in);
            User_Menu_Solution();
            String fn = sc.next();
            if (fn.equals("d") || fn.equals("r")) {
                User_Menu_System_receive();
                int op = sc.nextInt();
                System.out.println("L = ");
                Print_Matrix(L,fn);
                float[][] M;
                switch (op) {
                    case 1:
                        M = From_LU_To_M_V1(L,M2,fn);
                        System.out.println("M = ");
                        Print_Matrix(M,fn);
                        break;
                    case 2:
                        M = From_LLT_To_M_V1(L,M2,fn);
                        System.out.println("M = ");
                        Print_Matrix(M,fn);
                        break;
                    case 3:
                        M = From_LDLT_To_M_V1(L,M2,fn);
                        System.out.println("M = ");
                        Print_Matrix(M,fn);
                        break;
                    case 4:
                        M = From_LU_To_M_V2(L,M2,fn);
                        System.out.println("M = ");
                        Print_Matrix(M,fn);
                        break;
                    case 5:
                        M = From_LLT_To_M_V2(L,M2,fn);
                        System.out.println("M = ");
                        Print_Matrix(M,fn);
                        break;
                    case 6:
                        M = From_LDLT_To_M_V2(L,M2,fn);
                        System.out.println("M = ");
                        Print_Matrix(M,fn);
                        break;
                    default:
                        throw new Exception("you entered an invalid number");
                }
            } else {
                throw new Exception("you entered invalid value for a representation elementary actions and solution");
            }
        } else {
            throw new Exception("your input does not meet the conditions for invertible matrices");
        }
    }

    // get the matrix components
    public static void Decompose_Matrix(float[][] M) throws Exception {
        int m = M.length, n = M[0].length;
        if (m == n) {
            System.out.println("decompose the next matrix (" + n + "*" + n + " size):");
            Display_Exercise(M);
            Scanner sc = new Scanner(System.in);
            User_Menu_Solution();
            String fn = sc.next();
            if (fn.equals("d") || fn.equals("r")) {
                User_Menu_System_decompose();
                int op = sc.nextInt();
                System.out.println("M = ");
                Print_Matrix(M,fn);
                float[][] M2;
                switch (op) {
                    case 1:
                        M2 = From_M_To_LU_V1(M,fn);
                        System.out.println("U = ");
                        Print_Matrix(M2,fn);
                        break;
                    case 2:
                        M2 = From_M_To_LLT_V1(M,fn);
                        System.out.println("L' = ");
                        Print_Matrix(M2,fn);
                        break;
                    case 3:
                        M2 = From_M_To_LDLT_V1(M,fn);
                        System.out.println("L' = ");
                        Print_Matrix(M2,fn);
                        break;
                    case 4:
                        M2 = From_M_To_LU_V2(M,fn);
                        System.out.println("U = ");
                        Print_Matrix(M2,fn);
                        break;
                    case 5:
                        M2 = From_M_To_LLT_V2(M,fn);
                        System.out.println("L' = ");
                        Print_Matrix(M2,fn);
                        break;
                    case 6:
                        M2 = From_M_To_LDLT_V2(M,fn);
                        System.out.println("L' = ");
                        Print_Matrix(M2,fn);
                        break;
                    default:
                        throw new Exception("you entered an invalid number");
                }
            } else {
                throw new Exception("you entered invalid value for a representation elementary actions and solution");
            }
        } else {
            throw new Exception("your input does not meet the conditions for invertible matrices");
        }
    }

    public static void main(String[] args) {
        float[][] LU = {{2,1,-1},{-3,-1,2},{-2,1,2}};
        float[][] LLT = {{4,12,-16},{12,37,-43},{-16,-43,98}};
        float[][] LLT2 = {{4,-2,2},{-2,17,11},{2,11,35}};
        float[][] LDLT = {{-1,-3,4},{-3,-5,32},{4,32,75}};
        float[][] L = {{1,0,0},{-2,1,0},{-1,4,1}};
        float[][] LTT = {{1,-2,-1},{0,1,4},{0,0,1}};
        float[][] U = {{2,1,-1},{0,3,1},{0,0,-1}};
        float[][] LT = {{2,0,0},{6,1,0},{-8,5,3}};
        float[][] D = {{-4,0,0},{0,-1,0},{0,0,9}};
        float[][] E = {{5}};
        try {
            Scanner sc = new Scanner(System.in);
            User_Menu_Direction();
            String de = sc.next();
            if (de.equals("d")) {
                Decompose_Matrix(LU);
            } else if (de.equals("r")) {
                Receive_Matrix(L,LTT);
            } else {
                throw new Exception("you entered invalid value for a direction of the solution");
            }
        } catch (Exception e) {
            e.printStackTrace();
        }
    }
}
