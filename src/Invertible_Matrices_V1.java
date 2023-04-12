import java.util.Scanner;

public class Invertible_Matrices_V1 {

    private static void Print_Status_Matrices(float[][] M, float[][] InvM) {
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
            for (int j = 0; j < InvM[0].length; j++) {
                if (j == InvM[0].length - 1) {
                    if ((Math.round(InvM[i][j] * 1000000.0) / 1000000.0) % 1 == 0) {
                        System.out.print((int) (Math.round(InvM[i][j] * 1000000.0) / 1000000.0) + " ");
                    } else {
                        System.out.print(Math.round(InvM[i][j] * 1000000.0) / 1000000.0 + " ");
                    }
                } else {
                    if ((Math.round(InvM[i][j] * 1000000.0) / 1000000.0) % 1 == 0) {
                        System.out.print((int) (Math.round(InvM[i][j] * 1000000.0) / 1000000.0) + " ,");
                    } else {
                        System.out.print(Math.round(InvM[i][j] * 1000000.0) / 1000000.0 + " ,");
                    }
                }
            }
            System.out.println();
        }
        System.out.println();
    }

    private static void Retreat_Rows_Matrices(float[][] M, float[][] invM, int r1, int r2) {
        for (int j = 0; j < M[0].length; j++) {
            float t = M[r1][j];
            M[r1][j] = M[r2][j];
            M[r2][j] = t;
            float inv_t = invM[r1][j];
            invM[r1][j] = invM[r2][j];
            invM[r2][j] = inv_t;
        }
    }

    public static float[][] Unit_Matrix(int n) {
        float[][] I = new float[n][n];
        for (int i = 0; i < n; i++) {
            I[i][i] = 1;
        }
        return I;
    }

    public static float Determinant(float[][] M) {
        if (M.length == 1) {
            return M[0][0];
        }
        float sum = 0;
        for (int i = 0; i < M.length; i++) {
            sum += Math.pow(-1,i) * M[0][i] * Determinant(Sub_Matrix(M,0,i));
        }
        return sum;
    }

    public static float[][] Sub_Matrix(float[][] M, int x, int y) {
        int n = M.length, p = 0, q = 0;
        float[][] subM = new float[n - 1][n - 1];
        for (int i = 0; i < M.length; i++) {
            for (int j = 0; j < M[0].length; j++) {
                if (i != x && j != y) {
                    subM[p][q] = M[i][j];
                    q++;
                    if (q == subM[0].length) {
                        p++;
                        q = 0;
                    }
                }
            }
        }
        return subM;
    }

    public static float[][] Adjoint(float[][] M) {
        int n = M.length;
        float[][] Adj = new float[n][n];
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                Adj[j][i] = (float) Math.pow(-1,i + j) * Determinant(Sub_Matrix(M,i,j));
            }
        }
        return Adj;
    }

    private static int Get_Index_UnZero_Value(float[][] M, int k) {
        int n = M.length;
        for (int i = 0; i < n; i++) {
            if (M[k][i] != 0) {
                return i % n;
            }
        }
        return -1;
    }

    private static boolean Is_Unit_Vector(float[][] M, int k) {
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

    public static float[][] Mul_Const_Matrix(float k, float[][] M) {
        for (int i = 0; i < M.length; i++) {
            for (int j = 0; j < M[0].length; j++) {
                M[i][j] = k * M[i][j];
            }
        }
        return M;
    }

    private static boolean Is_Upper_Triangular(float[][] M) {
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

    private static boolean Is_Lower_Triangular(float[][] M) {
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

    private static void Sum_Elementary_Action(float k, int j, int i) {
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
                        System.out.println("R" + r + " --> R" + r + " + " + (int) ((-1) * k) + "*R" + c + "\n");
                    }
                } else {
                    System.out.println("R" + r + " --> R" + r + " + " + (-1) * k + "*R" + c + "\n");
                }
            }
        }
    }

    private static void Mul_Elementary_Action(float k, int j) {
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

    private static float[][] Upper_Ranking_Method(float[][] M, float[][] InvM) throws Exception {
        float det = Determinant(M);
        if (det != 0) {
            if (Is_Upper_Triangular(M) && Is_Lower_Triangular(M)) {
                System.out.println("M is already parallel triangular so now will be change directly to I:");
            } else if (Is_Upper_Triangular(M) && !Is_Lower_Triangular(M)) {
                System.out.println("M is already upper triangular so now we'll go directly to the lower ranking:");
                return Lower_Ranking_Method(M,InvM);
            } else if (!Is_Upper_Triangular(M) && Is_Lower_Triangular(M)) {
                System.out.println("transform L matrix to I by an elementary actions:");
            } else {
                System.out.println("transform M matrix to U by an elementary actions:");
            }
            int n = M.length;
            Print_Status_Matrices(M,InvM);
            for (int i = 0; i < n; i++) {
                if (Is_Lower_Triangular(M)) {
                    float c = 1 / M[i][i];
                    Mul_Elementary_Action(c,i);
                    for (int j = 0; j < n; j++) {
                        InvM[i][j] /= M[i][i];
                    }
                    M[i][i] = 1;
                    if (c != 1) {
                        Print_Status_Matrices(M,InvM);
                    }
                }
                for (int j = i + 1; j < n; j++) {
                    if (M[i][i] == 0) {
                        System.out.println("R" + (i + 1) + " <--> R" + (j + 1) + "\n");
                        Retreat_Rows_Matrices(M,InvM,i,j);
                        Print_Status_Matrices(M,InvM);
                    }
                    float c = M[j][i] / M[i][i];
                    Sum_Elementary_Action(c,j,i);
                    for (int k = 0; k < n; k++) {
                        M[j][k] -= M[i][k] * c;
                        InvM[j][k] -= InvM[i][k] * c;
                    }
                    M[j][i] = 0;
                    if (c != 0) {
                        Print_Status_Matrices(M,InvM);
                    }
                }
            }
            if (!Is_Upper_Triangular(M) || !Is_Lower_Triangular(M)) {
                System.out.print("and Then ");
                return Lower_Ranking_Method(M,InvM);
            }
            return InvM;
        } else {
            throw new Exception("this is a singular matrix");
        }
    }

    private static float[][] Lower_Ranking_Method(float[][] M, float[][] InvM) throws Exception {
        float det = Determinant(M);
        if (det != 0) {
            if (Is_Upper_Triangular(M) && Is_Lower_Triangular(M)) {
                System.out.println("M is already parallel triangular so now will be change directly to I:");
            } else if (!Is_Upper_Triangular(M) && Is_Lower_Triangular(M)) {
                System.out.println("M is already lower triangular so now we'll go directly to the upper ranking:");
                return Upper_Ranking_Method(M,InvM);
            } else if (Is_Upper_Triangular(M) && !Is_Lower_Triangular(M)) {
                System.out.println("transform U matrix to I by an elementary actions:");
            } else {
                System.out.println("transform M matrix to L by an elementary actions:");
            }
            int n = M.length;
            Print_Status_Matrices(M,InvM);
            for (int i = n - 1; i >= 0; i--) {
                if (Is_Upper_Triangular(M)) {
                    float c = 1 / M[i][i];
                    Mul_Elementary_Action(c,i);
                    for (int j = 0; j < n; j++) {
                        InvM[i][j] /= M[i][i];
                    }
                    M[i][i] = 1;
                    if (c != 1) {
                        Print_Status_Matrices(M,InvM);
                    }
                }
                for (int j = i - 1; j >= 0; j--) {
                    if (M[i][i] == 0) {
                        System.out.println("R" + (i + 1) + " <--> R" + (j + 1) + "\n");
                        Retreat_Rows_Matrices(M,InvM,i,j);
                        Print_Status_Matrices(M,InvM);
                    }
                    float c = M[j][i] / M[i][i];
                    Sum_Elementary_Action(c,j,i);
                    for (int k = n - 1; k >= 0; k--) {
                        M[j][k] -= M[i][k] * c;
                        InvM[j][k] -= InvM[i][k] * c;
                    }
                    M[j][i] = 0;
                    if (c != 0) {
                        Print_Status_Matrices(M,InvM);
                    }
                }
            }
            if (!Is_Upper_Triangular(M) || !Is_Lower_Triangular(M)) {
                System.out.print("and Then ");
                return Upper_Ranking_Method(M,InvM);
            }
            return InvM;
        } else {
            throw new Exception("this is a singular matrix");
        }
    }

    public static float[][] Invertible(float[][] M) throws Exception {
        float det = Determinant(M);
        if (det != 0) {
            System.out.println("matrix ranking:");
            int n = M.length;
            float[][] InvM = Unit_Matrix(n);
            Print_Status_Matrices(M,InvM);
            for (int i = 0; i < n; i++) {
                if (M[i][i] == 0) {
                    int r = Get_Index_UnZero_Value(M,i);
                    System.out.println("R" + (i + 1) + " <--> R" + (r + 1) + "\n");
                    Retreat_Rows_Matrices(M,InvM,i,r);
                    Print_Status_Matrices(M,InvM);
                }
                for (int j = 0; j < n; j++) {
                    if (i != j) {
                        float c = M[j][i] / M[i][i];
                        Sum_Elementary_Action(c,j,i);
                        for (int k = 0; k < n; k++) {
                            M[j][k] -= M[i][k] * c;
                            InvM[j][k] -= InvM[i][k] * c;
                        }
                        M[j][i] = 0;
                        if (c != 0) {
                            Print_Status_Matrices(M,InvM);
                        }
                    }
                    if (Is_Unit_Vector(M,j)) {
                        float c = 1 / M[j][j];
                        Mul_Elementary_Action(c,j);
                        for (int k = 0; k < n; k++) {
                            InvM[j][k] /= M[j][j];
                        }
                        M[j][j] = 1;
                        if (c != 1) {
                            Print_Status_Matrices(M,InvM);
                        }
                    }
                }
            }
            return InvM;
        } else {
            throw new Exception("this is a singular matrix");
        }
    }

    public static float[][] Invertible_Direct(float[][] M) throws Exception {
        float det = Determinant(M);
        if (det != 0) {
            int n = M.length;
            float[][] invM = new float[n][n];
            float[][] adj = Adjoint(M);
            for (int i = 0; i < n; i++) {
                for (int j = 0; j < n; j++) {
                    invM[i][j] = (1 / det) * adj[i][j];
                }
            }
            return invM;
        } else {
            throw new Exception("this is a singular matrix");
        }
    }

    public static float[][] Invertible_Direct_Plus(float[][] M) throws Exception {
        float det = Determinant(M);
        if (det != 0) {
            return Mul_Const_Matrix(1 / det,Adjoint(M));
        } else {
            throw new Exception("this is a singular matrix");
        }
    }

    public static void Correct_Checking_Matrix(float[][] M) throws Exception {
        System.out.println("find the next matrix's invertible:");
        Print_Matrix(M);
        Scanner sc = new Scanner(System.in);
        int op = sc.nextInt();
        switch (op) {
            case 1:
                M = Lower_Ranking_Method(M ,Unit_Matrix(M.length));
                Print_Matrix(M);
                break;
            case 2:
                M = Upper_Ranking_Method(M ,Unit_Matrix(M.length));
                Print_Matrix(M);
                break;
            case 3:
                M = Lower_Ranking_Method(Unit_Matrix(M.length), Upper_Ranking_Method(M,Unit_Matrix(M.length)));
                Print_Matrix(M);
                break;
            case 4:
                M = Upper_Ranking_Method(Unit_Matrix(M.length), Lower_Ranking_Method(M,Unit_Matrix(M.length)));
                Print_Matrix(M);
                break;
            case 5:
                M = Lower_Ranking_Method(Unit_Matrix(M.length), Lower_Ranking_Method(M,Unit_Matrix(M.length)));
                Print_Matrix(M);
                break;
            case 6:
                M = Upper_Ranking_Method(Unit_Matrix(M.length), Upper_Ranking_Method(M,Unit_Matrix(M.length)));
                Print_Matrix(M);
                break;
            case 7:
                M = Invertible(M);
                Print_Matrix(M);
                break;
            case 8:
                M = Invertible(Invertible(M));
                Print_Matrix(M);
                break;
            case 9:
                M = Upper_Ranking_Method(Invertible(M),Unit_Matrix(M.length));
                Print_Matrix(M);
                break;
            case 10:
                M = Lower_Ranking_Method(Invertible(M),Unit_Matrix(M.length));
                Print_Matrix(M);
                break;
            case 11:
                M = Invertible(Upper_Ranking_Method(M,Unit_Matrix(M.length)));
                Print_Matrix(M);
                break;
            case 12:
                M = Invertible(Lower_Ranking_Method(M,Unit_Matrix(M.length)));
                Print_Matrix(M);
                break;
            case 13:
                M = Invertible_Direct(M);
                Print_Matrix(M);
                break;
            case 14:
                M = Invertible_Direct_Plus(M);
                Print_Matrix(M);
                break;
        }
    }

    public static void main(String[] args) {
        float[][] A1 = {{2,1,-1},{-3,-1,2},{-2,1,2}};
        float[][] A2 = {{1,-1,-2},{2,-3,-5},{-1,3,5}};
        float[][] A3 = {{1,1,5},{1,2,7},{2,-1,4}};
        float[][] A4 = {{1,2,-1},{2,4,-3},{-1,-2,0}};
        float[][] B = {{-2,2,2,-1,1},{-4,4,4,4,3},{2,3,2,3,2},{-3,-1,1,2,2},{5,5,3,5,5}};
        float[][] C = {{-2,3,3,-2},{-1,4,2,-2},{1,3,1,3},{-3,-2,4,-5}};
        float[][] D = {{2,10,7,-5,0,8},{-12,-8,-9,-1,-11,-5},{-1,8,0,-4,-9,11},{12,6,-10,5,-11,-12},{12,11,-10,-2,-4,-12},{5,-5,3,6,8,-10}};
        float[][] E = {{1,1,-2,0,-1,-1},{0,-2,-1,2,-2,2},{0,0,-1,-1,-1,1},{0,0,-2,2,-2,0},{-2,-2,0,0,-2,2},{1,0,-1,0,0,2}};
        float[][] F = {{-2,0,-4,-2,5,-6},{6,-1,-2,0,1,0},{6,0,0,-2,5,2},{-4,0,6,4,4,6},{-4,-2,3,-6,4,-4},{5,1,-5,-1,5,-2}};
        float[][] G = {{-1,0,4,-1,-2,5,-1},{6,-7,7,1,-7,-6,5},{0,2,6,1,-7,6,1},{2,2,0,-2,-6,6,-5},{3,3,-3,1,7,-1,-2},{1,-5,1,-4,3,1,-4},{7,-5,5,0,-4,-4,1}};
        try {
            Correct_Checking_Matrix(A1);
        } catch (Exception e) {
            e.printStackTrace();
        }
    }
}
