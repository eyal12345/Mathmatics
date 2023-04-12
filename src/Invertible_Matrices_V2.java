
public class Invertible_Matrices_V2 {

    private static void Print_State(float[][] M ,float[][] invM) {
        for (int i = 0 ;i < M.length ;i++) {
            for (int j = 0 ;j < M[0].length ;j++) {
                if (j == M[0].length - 1) {
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
            System.out.print("| ");
            for (int j = 0 ;j < invM[0].length ;j++) {
                if (j == invM[0].length - 1) {
                    if ((Math.round(invM[i][j] * 1000000.0) / 1000000.0) % 1 == 0)
                        System.out.print((int)(Math.round(invM[i][j] * 1000000.0) / 1000000.0) + " ");
                    else
                        System.out.print(Math.round(invM[i][j] * 1000000.0) / 1000000.0 + " ");
                } else {
                    if ((Math.round(invM[i][j] * 1000000.0) / 1000000.0) % 1 == 0)
                        System.out.print((int)(Math.round(invM[i][j] * 1000000.0) / 1000000.0) + " ,");
                    else
                        System.out.print(Math.round(invM[i][j] * 1000000.0) / 1000000.0 + " ,");
                }
            }
            System.out.println();
        }
        System.out.println();
    }

    public static void Print_Matrix(float[][] A) {
        for (int i = 0 ;i < A.length ;i++) {
            for (int j = 0 ;j < A[0].length ;j++) {
                if (j == A[0].length - 1){
                    if ((Math.round(A[i][j] * 1000000.0) / 1000000.0) % 1 == 0)
                        System.out.print((int)(Math.round(A[i][j] * 1000000.0) / 1000000.0) + " ");
                    else
                        System.out.print(Math.round(A[i][j] * 1000000.0) / 1000000.0 + " ");
                } else {
                    if ((Math.round(A[i][j] * 1000000.0) / 1000000.0) % 1 == 0)
                        System.out.print((int)(Math.round(A[i][j] * 1000000.0) / 1000000.0) + " ,");
                    else
                        System.out.print(Math.round(A[i][j] * 1000000.0) / 1000000.0 + " ,");
                }
            }
            System.out.println();
        }
        System.out.println();
    }

    private static float[][] Mult_Mats(float[][] M1 ,float[][] M2) throws Exception {
        if (M1[0].length == M2.length) {
            float[][] M = new float[M1.length][M2[0].length];
            for (int i = 0 ;i < M1.length ;i++)
                for (int j = 0 ;j < M2[0].length ;j++)
                    for (int k = 0 ;k < M2.length ;k++)
                        M[i][j] += M1[i][k] * M2[k][j];
            return M;
        } else
            throw new Exception ("Error: The sizes Required is not Conditions");
    }

    private static void Retreat_Rows_Matrix(float[][] M ,int r1 ,int r2) {
        for (int j = 0 ;j < M[0].length ;j++) {
            float t = M[r1][j];
            M[r1][j] = M[r2][j];
            M[r2][j] = t;
        }
    }

    private static boolean Is_Unit_Vector(float[][] A ,int k) {
        int n = A.length;
        for (int i = 0 ;i < n ;i++)
            for (int j = 0 ;j < n && A[k][i] != 0 ;j++)
                if (i != j && A[k][j] != 0)
                    return false;
        return true;
    }

    private static boolean Is_Unit_Matrix(float[][] M) {
        int n = M.length;
        for (int i = 0 ;i < n ;i++)
            for (int j = 0 ;j < n ;j++)
                if (M[i][i] != 1 || (i != j && M[i][j] != 0))
                    return false;
        return true;
    }

    public static float[][] Unit_Matrix(int n) {
        float[][] I = new float[n][n];
        for (int i = 0 ;i < n ;i++)
            I[i][i] = 1;
        return I;
    }

    private static int Get_Index_UnZero(float[][] A ,int k) {
        int n = A.length;
        for (int i = k + 1 ;i < n ;i++)
            if (A[i][k] != 0)
                return i % n;
        return -1;
    }

    private static float[][] Copy_Matrix(float[][] A) {
        float[][] copyA = new float[A.length][A.length];
        for (int i = 0 ;i < A.length ;i++)
            for (int j = 0 ;j < A[0].length ;j++)
                copyA[i][j] = A[i][j];
        return copyA;
    }

    private static void Sum_Elementaric_Action(float k ,int j ,int i) {
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
                        System.out.println("R" + r + " --> R" + r + " + " + (int)(-k) + "*R" + c + "\n");
                else
                    System.out.println("R" + r + " --> R" + r + " + " + (-k) + "*R" + c + "\n");
            }
        }
    }

    private static void Mult_Elementaric_Action(float k ,int j) {
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

    public static float[][] Inv_Iter(float[][] M) throws Exception {
        int n = M.length;
        float[][] InvM = Unit_Matrix(n);
        int i = 0 ,j = 0;
        Print_State(M,InvM);
        while (!Is_Unit_Matrix(M)) {
            float[][] E = Unit_Matrix(n);
            if (M[i][i] == 0) {
                int k = Get_Index_UnZero(M,i);
                System.out.println("R" + (i + 1) + " <--> R" + (k + 1) + "\n");
                Retreat_Rows_Matrix(E,i,k);
                M = Mult_Mats(E,M);
                InvM = Mult_Mats(E,InvM);
                Print_State(M,InvM);
            } else {
                if (i != j) {
                    E[j][i] -= (M[j][i] / M[i][i]);
                    Sum_Elementaric_Action(-E[j][i],j,i);
                    M = Mult_Mats(E,M);
                    InvM = Mult_Mats(E,InvM);
                    M[j][i] = 0;
                    if (E[j][i] != 0)
                        Print_State(M,InvM);
                } else {
                    if (Is_Unit_Vector(M,j)) {
                        E[j][j] = 1 / M[j][j];
                        Mult_Elementaric_Action(E[j][j],j);
                        M = Mult_Mats(E,M);
                        InvM = Mult_Mats(E,InvM);
                        M[j][j] = 1;
                        if (E[j][j] != 1)
                            Print_State(M,InvM);
                    }
                }
                if (j == n - 1)
                    i = (i + 1) % n;
                j = (j + 1) % n;
            }
        }
        return InvM;
    }

    public static float[][] Inv_Rec(float[][] M ,float[][] InvM ,int i ,int j) throws Exception {
        int n = M.length;
        if (Is_Unit_Matrix(M))
            return InvM;
        else {
            float[][] E = Unit_Matrix(n);
            if (M[i][i] == 0) {
                int k = Get_Index_UnZero(M,i);
                System.out.println("R" + (i + 1) + " <--> R" + (k + 1) + "\n");
                Retreat_Rows_Matrix(E,i,k);
                M = Mult_Mats(E,M);
                InvM = Mult_Mats(E,InvM);
                Print_State(M,InvM);
            } else {
                if (i != j) {
                    E[j][i] -= (M[j][i] / M[i][i]);
                    Sum_Elementaric_Action(-E[j][i],j,i);
                    M = Mult_Mats(E,M);
                    InvM = Mult_Mats(E,InvM);
                    M[j][i] = 0;
                    if (E[j][i] != 0)
                        Print_State(M,InvM);
                } else {
                    if (Is_Unit_Vector(M,j)) {
                        E[j][j] = 1 / M[j][j];
                        Mult_Elementaric_Action(E[j][j],j);
                        M = Mult_Mats(E,M);
                        InvM = Mult_Mats(E,InvM);
                        M[j][j] = 1;
                        if (E[j][j] != 1)
                            Print_State(M,InvM);
                    }
                }
                if (j == n - 1)
                    i = (i + 1) % n;
                j = (j + 1) % n;
            }
            return Inv_Rec(M,InvM,i,j);
        }
    }

    public static void main(String[] args) {
        float[][] A1 = {{2,1,-1},{-3,-1,2},{-2,1,2}};
        float[][] A2 = {{1,-1,-2},{2,-3,-5},{-1,3,5}};
        float[][] A3 = {{1,1,5},{1,2,7},{2,-1,4}};
        float[][] A4 = {{1,1,5},{-2,-2,-10},{6,6,30}};
        float[][] B = {{-2,2,2,-1,1},{-4,4,4,4,3},{2,3,2,3,2},{-3,-1,1,2,2},{5,5,3,5,5}};
        float[][] Inv;
        try {
            float[][] A = Copy_Matrix(A1);
            Inv = Inv_Iter(A);
            //Print_State(A ,Unit_Matrix(A.length));
            //Inv = Inv_Rec(A ,Unit_Matrix(A.length) ,0 ,0);
            Print_Matrix(Inv);
        } catch (Exception e) {
            e.printStackTrace();
        }
    }
}
