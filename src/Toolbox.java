
public class Toolbox {

    private static void Print_Matrix(float[][] A) {
        int n = A.length;
        for (int i = 0 ;i < n ;i++) {
            for (int j = 0 ;j < n ;j++)
                if (j == n - 1) {
                    if ((Math.round(A[i][j] * 10000.0) / 10000.0) % 1 == 0)
                        System.out.print((int)(Math.round(A[i][j] * 10000.0) / 10000.0) + " \n");
                    else
                        System.out.print(Math.round(A[i][j] * 10000.0) / 10000.0 + " \n");
                } else {
                    if ((Math.round(A[i][j] * 10000.0) / 10000.0) % 1 == 0)
                        System.out.print((int)(Math.round(A[i][j] * 10000.0) / 10000.0) + ",");
                    else
                        System.out.print(Math.round(A[i][j] * 10000.0) / 10000.0 + ",");
                }
        }
        System.out.println();
    }

    private static void Print_State(float[][] A ,float[] b) {
        int n = A.length;
        for (int i = 0 ;i < n ;i++) {
            for (int j = 0 ;j < n ;j++) {
                if (j == n - 1) {
                    if ((Math.round(A[i][j] * 10000.0) / 10000.0) % 1 == 0)
                        System.out.print((int)(Math.round(A[i][j] * 10000.0) / 10000.0) + " ");
                    else
                        System.out.print(Math.round(A[i][j] * 10000.0) / 10000.0 + " ");
                } else {
                    if ((Math.round(A[i][j] * 10000.0) / 10000.0) % 1 == 0)
                        System.out.print((int)(Math.round(A[i][j] * 10000.0) / 10000.0) + " ,");
                    else
                        System.out.print(Math.round(A[i][j] * 10000.0) / 10000.0 + " ,");
                }
            }
            if ((Math.round(b[i] * 10000.0) / 10000.0) % 1 == 0)
                System.out.print("| " + (int)(Math.round(b[i] * 10000.0) / 10000.0) + " ,");
            else
                System.out.print("| " + Math.round(b[i] * 10000.0) / 10000.0 + " ,");
        }
        System.out.println();
    }

    // calculate multipication between two matrixs provided that M1's length column is equal to M2's length row
    private static float[][] Mult_Mats(float[][] M1 ,float[][] M2) throws Exception {
        if (M1[0].length == M2.length) {
            float[][] A = new float[M1.length][M2[0].length];
            for (int i = 0 ;i < M1.length ;i++)
                for (int j = 0; j < M2[0].length; j++)
                    for (int k = 0; k < M2.length; k++)
                        A[i][j] += M1[i][k] * M2[k][j];
            return A;
        } else
            throw new Exception ("Error: The sizes Required is not Conditions");
    }

    // calculate ranking of M matrix
    public static float[][] Ranking_Matrix(float[][] M) {
        int n = M.length;
        float[][] rankM = Copy_Matrix(M);
        for (int i = 0 ;i < n - 1 ;i++) {
            for (int j = i + 1 ;j < n ;j++) {
                if (rankM[i][i] == 0)
                    Retreat_Rows_Matrix(rankM, i, (i + 1) % n);
                float c = rankM[j][i] / rankM[i][i];
                for (int k = 0 ;k < n ;k++)
                    rankM[j][k] = rankM[j][k] - rankM[i][k] * c;
            }
        }
        return rankM;
    }

    // retreat between rows "r1" and "r2" in matrix M
    private static void Retreat_Rows_Matrix(float[][] M ,int r1 ,int r2) {
        int n = M.length;
        for (int j = 0 ;j < n ;j++) {
            float t = M[r1][j];
            M[r1][j] = M[r2][j];
            M[r2][j] = t;
        }
    }

    // create unit matrix in "n" size
    public static float[][] Create_Unit_Matrix(int n) {
        float[][] I = new float[n][n];
        for (int i = 0 ;i < n ;i++)
            I[i][i] = 1;
        return I;
    }

    // create zero vector in "n" size
    private static float[] Create_Zero_Vector(int n) {
        float[] x = new float[n];
        for (int k = 0 ;k < n ;k++)
            x[k] = 0;
        return x;
    }

    // duplicate matrix value of M to new matrix
    private static float[][] Copy_Matrix(float[][] M) {
        int n = M.length;
        float[][] copyM = new float[n][n];
        for (int i = 0 ;i < n ;i++)
            for (int j = 0 ;j < n ;j++)
                copyM[i][j] = M[i][j];
        return copyM;
    }

    // get the index with number row is not zero value
    private static int Get_Index_UnZero_Value(float[][] M ,int k) {
        int n = M.length;
        for (int i = k + 1 ;i < n ;i++)
            if (M[i][k] != 0)
                return i % n;
        return -1;
    }

    // check if in the M matrix is unit matrix
    private static boolean Is_Unit_Matrix(float[][] M) {
        int n = M.length;
        for (int i = 0 ;i < n ;i++)
            for (int j = 0 ;j < n ;j++)
                if (M[i][i] != 1 || (i != j && M[i][j] != 0))
                    return false;
        return true;
    }

    // check if in the M matrix exist row that unit vector
    private static boolean Is_Unit_Vector(float[][] M ,int k) {
        int n = M.length;
        for (int i = 0 ;i < n ;i++)
            for (int j = 0 ;j < n && M[k][i] != 0 ;j++)
                if (i != j && M[k][j] != 0)
                    return false;
        return true;
    }

    // check if M is zero matrix
    private static boolean Is_Zero_Matrix(float[][] M) {
        int n = M.length;
        for (int i = 0 ;i < n ;i++)
            for (int j = 0; j < n; j++)
                if (M[i][j] != 0)
                    return false;
        return true;
    }

    // check if b is zero vector
    private static boolean Is_Zero_Vector(float[] b) {
        int n = b.length;
        for (int i = 0; i < n; i++)
            if (b[i] != 0)
                return false;
        return true;
    }

    // check if vector values is equal to each other
    private static boolean Is_Equals_Values(float[] r) {
        int n = r.length;
        for (int i = 0 ;i < n - 1 ;i++)
            if (r[i] != r[i + 1])
                return false;
        return true;
    }

    // check if two vectors in the M matrix is linear dependents
    private static boolean Is_Linear_Dependent(float[][] M) {
        int n = M.length;
        float[] R = new float[n];
        for (int i1 = 0 ;i1 < n - 1 ;i1++) {
            for (int i2 = i1 + 1 ;i2 < n ;i2++) {
                for (int j = 0 ;j < n ;j++)
                    R[j] = M[i1][j] / M[i2][j];
                if (Is_Equals_Values(R))
                    return true;
            }
        }
        return false;
    }

    // check if matrix M is square matrix
    private static boolean Is_Square_Matrix(float[][] M) {
        if (M.length != M[0].length)
            return false;
        return true;
    }

    // check if M matrix is upper-triangular
    private static boolean Is_Upper_Triangular(float[][] M) {
        int n = M.length;
        for (int i = 0 ;i < n ;i++)
            for (int j = 0; j < i; j++)
                if (M[i][j] != 0)
                    return false;
        return true;
    }

    // check if M matrix is lower-triangular
    private static boolean Is_Lower_Triangular(float[][] M) {
        int n = M.length;
        for (int i = 0 ;i < n - 1 ;i++)
            for (int j = i + 1; j < n; j++)
                if (M[i][j] != 0)
                    return false;
        return true;
    }

    // check if eagen values in matrix M is positives
    private static boolean Is_Values_Positives(float[][] M) {
        int n = M.length;
        for (int i = 0 ;i < n ;i++)
            if (M[i][i] < 0)
                return false;
        return true;
    }

    // check if in matrix M main diagonal values is 1
    private static boolean Is_One_Slant(float[][] M) {
        int n = M.length;
        for (int i = 0 ;i < n ;i++)
            if (M[i][i] != 1)
                return false;
        return true;
    }

    // multipicate const "k" in matrix M
    public static float[][] Mult_Const_Matrix(float k ,float[][] M) {
        for (int i = 0 ;i < M.length ;i++)
            for (int j = 0 ;j < M[0].length ;j++)
                M[i][j] = k * M[i][j];
        return M;
    }

    // calculate determinant of matrix M
    public static float Determinant(float[][] M) {
        int n = M.length;
        if (n == 1)
            return M[0][0];
        else {
            float sum = 0;
            for (int i = 0; i < n; i++)
                sum += Math.pow(-1, i) * M[0][i] * Determinant(Sub_Matrix(M, 0, i));
            return sum;
        }
    }

    // calculate sub-matrix of M by cutting row "x" and column "y"
    public static float[][] Sub_Matrix(float[][] M ,int x ,int y) {
        int n = M.length ,p = 0 ,q = 0;
        float[][] subM = new float[n - 1][n - 1];
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                if (i != x && j != y) {
                    subM[p][q] = M[i][j];
                    q++;
                    if (q == n - 1) {
                        p++;
                        q = 0;
                    }
                }
            }
        }
        return subM;
    }

    // calculate sub-vector of v by without index "x"
    public static float[] Sub_Vector(float[] v ,int x) {
        int n = v.length;
        float[] subv = new float[n - 1];
        for (int i = 0 ;i < n ;i++)
            if (i != x)
                subv[i] = v[i];
        return subv;
    }

    // calculate transpose matrix of M
    public static float[][] Transpose(float[][] M) {
        int n = M.length;
        float[][] MT = new float[n][n];
        for (int i = 0 ;i < n ;i++)
            for (int j = 0 ;j < n ;j++)
                MT[j][i] = M[i][j];
        return MT;
    }

    // calculate adjoint matrix of M
    public static float[][] Adjoint(float[][] M) {
        int n = M.length;
        float[][] AdjM = new float[n][n];
        for (int i = 0 ;i < n ;i++)
            for (int j = 0; j < n; j++)
                AdjM[j][i] = (float) (Math.pow(-1, i + j) * Determinant(Sub_Matrix(M, i, j)));
        return AdjM;
    }

    // calculate invertible matrix of M
    public static float[][] Invertible(float[][] M) {
        int n = M.length;
        float[][] invM = new float[n][n];
        float[][] adjM = Adjoint(M);
        float detM = Determinant(M);
        for (int i = 0 ;i < n ;i++)
            for (int j = 0; j < n; j++)
                invM[i][j] = (1 / detM) * adjM[i][j];
        return invM;
    }

    public static void Random_Systems(int n) {
        try {
            float[][] A = new float[n][n];
            float[] b = new float[n];
            int low = 10*n + 1;
            for (int k = 1 ;k <= 20 ;k++) {
                System.out.println("Equation Number " + k + ":");
                for (int i = 0 ;i < n ;i++) {
                    for (int j = 0 ;j < n ;j++) {
                        A[i][j] = (int)((Math.random()*22) - 11);
                    }
                    b[i] = (int)((Math.random()*2*low) - low);
                }
                Print_State(A,b);
                System.out.println();
            }
        } catch (Exception ex) {
            ex.printStackTrace();
        }
    }

    private static String IntegerToRomanNumeral(int input) {
        if (input < 1 || input > 99)
            return "Invalid Roman Number Value";
        String s = "";
        while (input >= 10) {
            s += "X";
            input -= 10;
        }
        while (input >= 9) {
            s += "IX";
            input -= 9;
        }
        while (input >= 5) {
            s += "V";
            input -= 5;
        }
        while (input >= 4) {
            s += "IV";
            input -= 4;
        }
        while (input >= 1) {
            s += "I";
            input -= 1;
        }
        return s;
    }

    private static void Random_Determinants(int n ,int r) {
        float[][] A = new float[n][n];
        float det = -1;
        int count = 1;
        while (det != 0) {
            System.out.println("System Number " + count);
            for (int i = 0 ;i < n ;i++)
                for (int j = 0 ;j < n ;j++)
                    A[i][j] = (int)((Math.random()*2*(r + 1)) - (r + 1));
            det = Determinant(A);
            Print_Matrix(A);
            System.out.println("det = " + (int)det);
            count++;
        }
    }

    private static float[][] resizeArray (float[][] oldArray, int addrows) {
        int oldSize = java.lang.reflect.Array.getLength(oldArray);
        Class elementType = oldArray.getClass().getComponentType();
        int newSize = oldSize + addrows;
        float[][] newArray = (float[][]) java.lang.reflect.Array.newInstance(elementType, newSize);
        int preserveLength = Math.min(oldSize, newSize);
        if (preserveLength > 0)
            System.arraycopy(oldArray, 0, newArray, 0, preserveLength);
        return newArray;
    }

    // solve system 3x3 in O(1) complexity
    public static void Solve_System_Direct(float[][] A ,float[] b) {
        float det = A[0][0]*A[1][1]*A[2][2] - A[0][0]*A[1][2]*A[2][1] - A[0][1]*A[1][0]*A[2][2] + A[0][2]*A[1][0]*A[2][1] + A[0][1]*A[1][2]*A[2][0] - A[0][2]*A[1][1]*A[2][0];
        if (det != 0) {
            float det1 = b[0]*A[1][1]*A[2][2] - b[0]*A[1][2]*A[2][1] - A[0][1]*b[1]*A[2][2] + A[0][2]*b[1]*A[2][1] + A[0][1]*A[1][2]*b[2] - A[0][2]*A[1][1]*b[2];
            System.out.println(det1 / det);
            float det2 = A[0][0]*b[1]*A[2][2] - A[0][0]*A[1][2]*b[2] - b[0]*A[1][0]*A[2][2] + A[0][2]*A[1][0]*b[2] + b[0]*A[1][2]*A[2][0] - A[0][2]*b[1]*A[2][0];
            System.out.println(det2 / det);
            float det3 = A[0][0]*A[1][1]*b[2] - A[0][0]*b[1]*A[2][1] - A[0][1]*A[1][0]*b[2] + b[0]*A[1][0]*A[2][1] + A[0][1]*b[1]*A[2][0] - b[0]*A[1][1]*A[2][0];
            System.out.println(det3 / det);
        }
    }

    public static void main(String[] args) {
        Random_Determinants(8,4);
    }
}
