import java.util.Scanner;

public class Mat_Decomposition {

    public static void Print_Matrix(float[][] M) {
        for (int i = 0 ;i < M.length ;i++) {
            for (int j = 0 ;j < M[0].length ;j++) {
                if (j == M[0].length - 1){
                    if ((Math.round(M[i][j] * 1000000.0) / 1000000.0) % 1 == 0)
                        System.out.print((int)(Math.round(M[i][j] * 1000000.0) / 1000000.0) + " ");
                    else
                        System.out.print(Math.round(M[i][j] * 1000000.0) / 1000000.0 + " ");
                } else {
                    if ((Math.round(M[i][j] * 1000000.0) / 1000000.0) % 1 == 0)
                        System.out.print((int)(Math.round(M[i][j] * 1000000.0) / 1000000.0) + " ,");
                    else
                        System.out.print(Math.round(M[i][j] * 1000000.0) / 1000000.0 + " ,");
                }
            }
            System.out.println();
        }
        System.out.println();
    }

    private static boolean Is_Square_Matrix(float[][] M) {
        if (M.length != M[0].length)
            return false;
        return true;
    }

    private static boolean Is_Upper_Triangular(float[][] M) {
        int n = M.length;
        for (int i = 0 ;i < n ;i++)
            for (int j = 0 ;j < i ;j++)
                if (M[i][j] != 0)
                    return false;
        return true;
    }

    private static boolean Is_Lower_Triangular(float[][] M) {
        int n = M.length;
        for (int i = 0 ;i < n - 1 ;i++)
            for (int j = i + 1 ;j < n ;j++)
                if (M[i][j] != 0)
                    return false;
        return true;
    }

    private static boolean Is_Symmetrical_Matrix(float[][] M) {
        for (int i = 0 ;i < M.length ;i++)
            for (int j = 0 ;j < M[0].length ;j++)
                if (M[i][j] != M[j][i])
                    return false;
        return true;
    }

    private static boolean Is_Values_Positives(float[][] M) {
        for (int i = 0 ;i < M.length ;i++)
            if (M[i][i] < 0)
                return false;
        return true;
    }

    private static boolean Is_One_Slant(float[][] M) {
        int n = M.length;
        for (int i = 0 ;i < n ;i++)
            if (M[i][i] != 1)
                return false;
        return true;
    }

    public static float[][] Mult_Mats(float[][] M1 ,float[][] M2) throws Exception {
        if (M1[0].length == M2.length) {
            float[][] M = new float[M1.length][M2[0].length];
            for (int i = 0 ;i < M1.length ;i++)
                for (int j = 0 ;j < M2[0].length ;j++)
                    for (int k = 0 ;k < M2.length ;k++)
                        M[i][j] += M1[i][k] * M2[k][j];
            return M;
        } else
            throw new Exception ("The sizes Required is not Conditions");
    }

    private static void Retreat_Rows_Matrix(float[][] M ,int r1 ,int r2) {
        for (int j = 0 ;j < M[0].length ;j++) {
            float t = M[r1][j];
            M[r1][j] = M[r2][j];
            M[r2][j] = t;
        }
    }

    public static float[][] Transpose(float[][] M) {
        float[][] MT = new float[M[0].length][M.length];
        for (int i = 0 ;i < M.length ;i++)
            for (int j = 0 ;j < M[0].length ;j++)
                MT[j][i] = M[i][j];
        return MT;
    }

    public static void Sum_Elementaric_Action(float k ,int j ,int i) {
        if (k != 0) {
            int r = j + 1;
            int c = i + 1;
            if (k > 0) {
                if (k % 1 == 0)
                    if (k == 1)
                        System.out.println("R" + r + " --> R" + r + " - R" + c + "\n");
                    else
                        System.out.println("R" + r + " --> R" + r + " - " + (int)k + "*R" + c + "\n");
                else
                    System.out.println("R" + r + " --> R" + r + " - " + k + "*R" + c + "\n");
            } else {
                if (k % 1 == 0)
                    if (k == -1)
                        System.out.println("R" + r + " --> R" + r + " + R" + c + "\n");
                    else
                        System.out.println("R" + r + " --> R" + r + " + " + (int)((-1) * k) + "*R" + c + "\n");
                else
                    System.out.println("R" + r + " --> R" + r + " + " + (-1) * k + "*R" + c + "\n");
            }
        }
    }

    public static void Mult_Elementaric_Action(float k ,int j) {
        if (k != 1) {
            int r = j + 1;
            if (k % 1 == 0)
                if (k == -1)
                    System.out.println("R" + r + " --> - R" + r + "\n");
                else
                    System.out.println("R" + r + " --> " + (int)k + "*R" + r + "\n");
            else
                System.out.println("R" + r + " --> " + k + "*R" + r + "\n");
        }
    }

    public static void From_LU_To_M_V1(float[][] L ,float[][] U) throws Exception {
        System.out.println("L = ");
        Print_Matrix(L);
        System.out.println("U = ");
        Print_Matrix(U);
        if (Is_Square_Matrix(L) && Is_Square_Matrix(U) && Is_Lower_Triangular(L) && Is_Upper_Triangular(U) && Is_One_Slant(L)) {
            float[][] M = Mult_Mats(L,U);
            System.out.println("M = ");
            Print_Matrix(M);
        } else
            throw new Exception ("Not all conditions are held");
    }

    public static void From_M_To_LU_V1(float[][] M) throws Exception {
        System.out.println("M = ");
        Print_Matrix(M);
        if (Is_Square_Matrix(M)) {
            int n = M.length;
            float[][] L = new float[n][n];
            for (int i = 0 ;i < n ;i++) {
                L[i][i] = 1;
                for (int j = i + 1 ;j < n ;j++) {
                    if (M[i][i] == 0) {
                        System.out.println("R" + (i + 1) + " <--> R" + (j + 1) + "\n");
                        Retreat_Rows_Matrix(M,i,j);
                        Print_Matrix(M);
                    }
                    L[j][i] = M[j][i] / M[i][i];
                    Sum_Elementaric_Action(L[j][i],j,i);
                    for (int k = 0 ;k < n ;k++)
                        M[j][k] -= M[i][k] * L[j][i];
                    M[j][i] = 0;
                    if (L[j][i] != 0)
                        Print_Matrix(M);
                }
            }
            System.out.println("L = ");
            Print_Matrix(L);
            System.out.println("U = ");
            Print_Matrix(M);
        } else
            throw new Exception ("The matrix is not square");
    }

    public static void From_LLT_To_M_V1(float[][] L) throws Exception {
        System.out.println("L = ");
        Print_Matrix(L);
        if (Is_Square_Matrix(L) && Is_Lower_Triangular(L)) {
            float[][] M = Mult_Mats(L,Transpose(L));
            System.out.println("L* = ");
            Print_Matrix(Transpose(L));
            System.out.println("M = ");
            Print_Matrix(M);
        } else
            throw new Exception ("Not all conditions are held");
    }

    public static void From_M_To_LLT_V1(float[][] M) throws Exception {
        System.out.println("M = ");
        Print_Matrix(M);
        if (Is_Square_Matrix(M) && Is_Symmetrical_Matrix(M) && Is_Values_Positives(M)) {
            int n = M.length;
            float[][] L = new float[n][n];
            float c;
            for (int i = 0 ;i < n ;i++) {
                for (int j = i + 1 ;j < n ;j++) {
                    c = M[j][i] / M[i][i];
                    Sum_Elementaric_Action(c,j,i);
                    for (int k = 0 ;k < n ;k++)
                        M[j][k] -= M[i][k] * c;
                    M[j][i] = 0;
                    if (c != 0)
                        Print_Matrix(M);
                }
                c = (float) Math.sqrt(M[i][i]);
                Mult_Elementaric_Action(1 / c,i);
                for (int k = i ;k < n ;k++) {
                    M[i][k] /= c;
                    L[k][i] = M[i][k];
                }
                if (c != 1)
                    Print_Matrix(M);
            }
            System.out.println("L = ");
            Print_Matrix(L);
            System.out.println("L* = ");
            Print_Matrix(M);
        } else
            throw new Exception ("Not all conditions are held");
    }

    public static void From_LDLT_To_M_V1(float[][] L ,float[][] D) throws Exception {
        System.out.println("L = ");
        Print_Matrix(L);
        System.out.println("D = ");
        Print_Matrix(D);
        if (Is_Square_Matrix(L) && Is_Lower_Triangular(L) && Is_One_Slant(L) && Is_Square_Matrix(D) && Is_Lower_Triangular(D) && Is_Upper_Triangular(D)) {
            float[][] M = Mult_Mats(Mult_Mats(L,D),Transpose(L));
            System.out.println("L* = ");
            Print_Matrix(Transpose(L));
            System.out.println("M = ");
            Print_Matrix(M);
        } else
            throw new Exception ("Not all conditions are held");
    }

    public static void From_M_To_LDLT_V1(float[][] M) throws Exception {
        System.out.println("M = ");
        Print_Matrix(M);
        if (Is_Square_Matrix(M) && Is_Symmetrical_Matrix(M)) {
            int n = M.length;
            float[][] L = new float[n][n];
            float[][] D = new float[n][n];
            for (int i = 0 ;i < n ;i++) {
                for (int j = i + 1 ;j < n ;j++) {
                    float c = M[j][i] / M[i][i];
                    Sum_Elementaric_Action(c,j,i);
                    for (int k = 0 ;k < n ;k++)
                        M[j][k] -= M[i][k] * c;
                    M[j][i] = 0;
                    if (c != 0)
                        Print_Matrix(M);
                }
                D[i][i] = M[i][i];
                Mult_Elementaric_Action(1 / D[i][i],i);
                for (int k = i ;k < n ;k++) {
                    M[i][k] /= D[i][i];
                    L[k][i] = M[i][k];
                }
                if (D[i][i] != 1)
                    Print_Matrix(M);
            }
            System.out.println("L = ");
            Print_Matrix(L);
            System.out.println("D = ");
            Print_Matrix(D);
            System.out.println("L* = ");
            Print_Matrix(M);
        } else
            throw new Exception ("Not all conditions are held");
    }

    public static void From_LU_To_M_V2(float[][] L ,float[][] U) throws Exception {
        System.out.println("L = ");
        Print_Matrix(L);
        System.out.println("U = ");
        Print_Matrix(U);
        if (Is_Square_Matrix(L) && Is_Square_Matrix(U) && Is_Lower_Triangular(L) && Is_Upper_Triangular(U) && Is_One_Slant(L)) {
            int n = L.length;
            float[][] M = new float[n][n];
            for (int i = 0 ;i < n ;i++) {
                for (int j = 0; j < n; j++) {
                    int m = Math.min(i, j);
                    for (int k = 0; k <= m; k++)
                        M[i][j] += L[i][k] * U[k][j];
                }
            }
            System.out.println("M = ");
            Print_Matrix(M);
        }  else
            throw new Exception ("Not all conditions are held");
    }

    public static void From_M_To_LU_V2(float[][] M) throws Exception {
        System.out.println("M = ");
        Print_Matrix(M);
        if (Is_Square_Matrix(M)) {
            int n = M.length;
            float[][] L = new float[n][n];
            float[][] U = new float[n][n];
            for (int i = 0 ;i < n ;i++) {
                L[i][i] = 1;
                for (int j = 0 ;j < n ;j++) {
                    int m = Math.min(i,j);
                    for (int k = 0 ;k < m ;k++)
                        M[i][j] -= L[i][k] * U[k][j];
                    if (i <= j)
                        U[i][j] = M[i][j] / L[i][i];
                    else
                        L[i][j] = M[i][j] / U[j][j];
                }
            }
            System.out.println("L = ");
            Print_Matrix(L);
            System.out.println("U = ");
            Print_Matrix(U);
        } else
            throw new Exception ("The matrix is not square");
    }

    public static void From_LLT_To_M_V2(float[][] L) throws Exception {
        System.out.println("L = ");
        Print_Matrix(L);
        if (Is_Square_Matrix(L) && Is_Lower_Triangular(L)) {
            int n = L.length;
            float[][] M = new float[n][n];
            for (int i = 0 ;i < n ;i++) {
                for (int j = 0; j < n; j++) {
                    int m = Math.min(i,j);
                    for (int k = 0; k <= m; k++)
                        M[i][j] += L[i][k] * L[j][k];
                }
            }
            System.out.println("L* = ");
            Print_Matrix(Transpose(L));
            System.out.println("M = ");
            Print_Matrix(M);
        } else
            throw new Exception ("Not all conditions are held");
    }

    public static void From_M_To_LLT_V2(float[][] M) throws Exception {
        System.out.println("M = ");
        Print_Matrix(M);
        if (Is_Square_Matrix(M) && Is_Symmetrical_Matrix(M) && Is_Values_Positives(M)) {
            int n = M.length;
            float[][] L = new float[n][n];
            for (int j = 0 ;j < n ;j++) {
                for (int i = j ;i < n ;i++) {
                    for (int k = 0 ;k < j ;k++)
                        M[i][j] -= L[i][k] * L[j][k];
                    L[i][j] = M[i][j] / (float) Math.sqrt(M[j][j]);
                }
            }
            System.out.println("L = ");
            Print_Matrix(L);
            System.out.println("L* = ");
            Print_Matrix(Transpose(L));
        } else
            throw new Exception ("Not all conditions are held");
    }

    public static void From_LDLT_To_M_V2(float[][] L ,float[][] D) throws Exception {
        System.out.println("L = ");
        Print_Matrix(L);
        System.out.println("D = ");
        Print_Matrix(D);
        if (Is_Square_Matrix(L) && Is_Lower_Triangular(L) && Is_One_Slant(L) && Is_Square_Matrix(D) && Is_Lower_Triangular(D) && Is_Upper_Triangular(D)) {
            int n = L.length;
            float[][] M = new float[n][n];
            for (int i = 0 ;i < n ;i++) {
                for (int j = 0; j < n; j++) {
                    int m = Math.min(i,j);
                    for (int k = 0; k <= m; k++)
                        M[i][j] += L[i][k] * L[j][k] * D[k][k];
                }
            }
            System.out.println("L* = ");
            Print_Matrix(Transpose(L));
            System.out.println("M = ");
            Print_Matrix(M);
        } else
            throw new Exception ("Not all conditions are held");
    }

    public static void From_M_To_LDLT_V2(float[][] M) throws Exception {
        System.out.println("M = ");
        Print_Matrix(M);
        if (Is_Square_Matrix(M) && Is_Symmetrical_Matrix(M)) {
            int n = M.length;
            float[][] L = new float[n][n];
            float[][] D = new float[n][n];
            for (int j = 0 ;j < n ;j++) {
                L[j][j] = 1;
                for (int i = j ;i < n ;i++) {
                    for (int k = 0 ;k < j ;k++)
                        M[i][j] -= L[i][k] * L[j][k] * D[k][k];
                    L[i][j] = M[i][j] / M[j][j];
                }
                D[j][j] = M[j][j];
            }
            System.out.println("L = ");
            Print_Matrix(L);
            System.out.println("D = ");
            Print_Matrix(D);
            System.out.println("L* = ");
            Print_Matrix(Transpose(L));
        } else
            throw new Exception ("Not all conditions are held");
    }

    public static void main(String[] args) {
        float[][] LU = {{2, 1, -1}, {-3, -1, 2}, {-2, 1, 2}};
        float[][] LLT = {{4,12,-16},{12,37,-43},{-16,-43,98}};
        float[][] LLT2 = {{4,-2,2},{-2,17,11},{2,11,35}};
        float[][] LDLT = {{-1,-3,4},{-3,-5,32},{4,32,75}};
        float[][] L = {{1, 0, 0}, {-2, 1, 0}, {-1, 4, 1}};
        float[][] U = {{2, 1, -1}, {0, 3, 1}, {0, 0, -1}};
        float[][] LT = {{2, 0, 0}, {6, 1, 0}, {-8, 5, 3}};
        float[][] D = {{-4, 0, 0}, {0, -1, 0}, {0, 0, 9}};
        try {
            Scanner sc = new Scanner(System.in);
            int op = sc.nextInt();
            switch (op) {
                case 1:
                    From_LU_To_M_V1(L,U);
                    break;
                case 2:
                    From_M_To_LU_V1(LU);
                    break;
                case 3:
                    From_LLT_To_M_V1(LT);
                    break;
                case 4:
                    From_M_To_LLT_V1(LLT);
                    break;
                case 5:
                    From_LDLT_To_M_V1(L,D);
                    break;
                case 6:
                    From_M_To_LDLT_V1(LDLT);
                    break;
                case 7:
                    From_LU_To_M_V2(L,U);
                    break;
                case 8:
                    From_M_To_LU_V2(LU);
                    break;
                case 9:
                    From_LLT_To_M_V2(LT);
                    break;
                case 10:
                    From_M_To_LLT_V2(LLT);
                    break;
                case 11:
                    From_LDLT_To_M_V2(L,D);
                    break;
                case 12:
                    From_M_To_LDLT_V2(LDLT);
                    break;
            }
        } catch (Exception e) {
            e.printStackTrace();
        }
    }
}
