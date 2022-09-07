import java.io.IOException;
import java.util.Scanner;
import java.util.Vector;

public class Matrix_Aritmetics {

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

    public static float[][] Create_Matrix(int n){
        Scanner sc = new Scanner(System.in);
        float[][] M = new float[n][n];
        System.out.println("insert values to matrix:");
        for(int i = 0 ;i < n ;i++)
            for(int j = 0 ;j < n ;j++){
                System.out.print("M["+i+"]["+j+"]->");
                M[i][j] = sc.nextFloat();
            }
        return M;
    }

    public static float[][] Unit_Matrix(int n) {
        float[][] I = new float[n][n];
        for (int i = 0 ;i < I.length ;i++)
            for (int j = 0 ;j < I[0].length ;j++) {
                if (i == j)
                    I[i][j] = 1;
                else
                    I[i][j] = 0;
            }
        return I;
    }

    public static float[][] Transpose(float[][] M) {
        float[][] MT = new float[M[0].length][M.length];
        for (int i = 0 ;i < M.length ;i++)
            for (int j = 0 ;j < M[0].length ;j++)
                MT[j][i] = M[i][j];
        return MT;
    }

    public static float Determinant(float[][] M) {
        int n = M.length;
        if (n == 1)
            return M[0][0];
        float det = 0;
        for (int i = 0 ;i < n ;i++)
            det += Math.pow(-1,i) * M[0][i] * Determinant(Sub_Matrix(M,0,i));
        return det;
    }

    public static float[][] Sub_Matrix(float[][] M, int x, int y) {
        int n = M.length ,p = 0 ,q = 0;
        float[][] subA = new float[n - 1][n - 1];
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                if (i != x && j != y) {
                    subA[p][q] = M[i][j];
                    q++;
                    if (q == n - 1) {
                        p++;
                        q = 0;
                    }
                }
            }
        }
        return subA;
    }

    public static float[][] Invertible(float[][] M) {
        float[][] invM = new float[M.length][M.length];
        float[][] adj = Adjoint(M);
        for (int i = 0 ;i < M.length ;i++)
            for (int j = 0 ;j < M.length ;j++)
                invM[i][j] = (1 / Determinant(M)) * adj[i][j];
        return invM;
    }

    public static float[][] Adjoint(float[][] M) {
        int n = M.length;
        float[][] Adj = new float[n][n];
        for (int i = 0 ;i < n ;i++)
            for (int j = 0 ;j < n ;j++)
                Adj[j][i] = (float)(Math.pow(-1,i + j) * Determinant(Sub_Matrix(M,i,j)));
        return Adj;
    }

    public static float Trace_Matrix(float[][] M) {
        float sum = 0;
        for (int i = 0 ;i < M.length ;i++)
            sum += M[i][i];
        return sum;
    }

    public static float Norm_Standard(float[][] M) {
        float sum = 0;
        float norm;
        for (int i = 0 ;i < M.length ;i++)
            for (int j = 0 ;j < M[0].length ;j++)
                sum += + Math.pow(M[i][j], 2);
        norm = (float) Math.sqrt(sum);
        return norm;
    }

    public static float Norm_Infinite(float[][] M) {
        float max = Integer.MIN_VALUE;
        for (int i = 0 ;i < M.length ;i++) {
            float sum = 0;
            for (int j = 0 ;j < M[0].length ;j++)
                sum += Math.abs(M[i][j]);
            if (max < sum)
                max = sum;
        }
        return max;
    }

    public static float Norm_One(float[][] M) {
        float max = Integer.MIN_VALUE;
        for (int j = 0 ;j < M[0].length ;j++) {
            float sum = 0;
            for (int i = 0 ;i < M.length ;i++)
                sum += Math.abs(M[i][j]);
            if (max < sum)
                max = sum;
        }
        return max;
    }

    public static float Condition_Number(float[][] M) {
        return Norm_Standard(M) * Norm_Standard(Invertible(M));
    }

    public static float Condition_Number_Infinite(float[][] M) { return Norm_Infinite(M) * Norm_Infinite(Invertible(M)); }

    public static float Condition_Number_One(float[][] M) {
        return Norm_One(M) * Norm_One(Invertible(M));
    }

    public static double[][] Add_Or_Sub_Mats(double[][] M1, double[][] M2, char oper) throws Exception {
        if (M1.length == M2.length && M1[0].length == M2[0].length) {
            int r = M1.length;
            int c = M1[0].length;
            double[][] M = new double[r][c];
            for (int i = 0 ;i < r ;i++)
                for (int j = 0 ;j < c ;j++){
                    if(oper == '+')
                        M[i][j] = M1[i][j] + M2[i][j];
                    if(oper == '-')
                        M[i][j] = M1[i][j] - M2[i][j];
                }
            return M;
        } else {
            throw new Exception ("Error: The sizes of matrix unlike");
        }
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

    private static boolean Is_Square_Matrix(float[][] M) {
        return M.length == M[0].length;
    }

    public static Vector<Float> Eigen_Values(float[][] M) {
        Vector<Float> v = null;
        for (int i = 0 ;i < M.length ;i++)
            v.add(M[i][i]);
        return v;
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

    private static boolean Is_Transpose_Matrix(float[][] M1, float[][] M2) {
        for (int i = 0 ;i < M1.length ;i++)
            for (int j = 0 ;j < M1[0].length ;j++)
                if (M1[i][j] != M2[j][i])
                    return false;
        return true;
    }

    private static boolean Is_One_Slant(float[][] M) {
        boolean flag = true;
        for (int i = 0 ;i < M.length ;i++) {
            if (M[i][i] == 1)
                flag = true;
            else {
                flag = false;
                break;
            }
        }
        return flag;
    }

    public static float[][] Mult_Const_Matrix(float k ,float[][] M) {
        for (int i = 0 ;i < M.length ;i++)
            for (int j = 0 ;j < M[0].length ;j++)
                M[i][j] = k * M[i][j];
        return M;
    }

    public static float[][] Random_Matrix(int low ,int high ,int n) {
        float[][] M = new float[n][n];
        for (int i = 0 ;i < M.length ;i++)
            for (int j = 0 ;j < M[0].length ;j++)
                M[i][j] = (int)(Math.random() * (high - low + 1)) + low;
        return M;
    }

    public static double Average_Values(float[][] M) {
        float average ,count = 0;
        int N = M.length * M[0].length;
        for (int i = 0 ;i < M.length ;i++)
            for (int j = 0 ;j < M[0].length ;j++)
                count += M[i][j];
        average = count / N;
        return average;
    }

    public static float[][] Sort_Rows(float[][] M) {
        float[] helper = new float[M[0].length];
        for (int i = 0 ;i < M.length ;i++) {
            for (int j = 0 ;j < M[0].length ;j++)
                helper[j] = M[i][j];
            helper = Selection_Sort(helper);
            for (int j = 0 ;j < M[0].length ;j++)
                M[i][j] = helper[j];
        }
        return M;
    }

    public static float[][] Sort_Cols(float[][] M) {
        float[] helper = new float[M.length];
        for (int j = 0 ;j < M[0].length ;j++) {
            for (int i = 0 ;i < M.length ;i++)
                helper[i] = M[i][j];
            helper = Selection_Sort(helper);
            for (int i = 0 ;i < M.length ;i++)
                M[i][j] = helper[i];
        }
        return M;
    }

    public static float[][] Sort_Rows_Cols(float[][] M) {
        M = Sort_Rows(M);
        M = Sort_Cols(M);
        return M;
    }

    public static float[][] Sort_Cols_Rows(float[][] M) {
        M = Sort_Cols(M);
        M = Sort_Rows(M);
        return M;
    }

    private static float[] Selection_Sort(float[] V) {
        for (int i = 0 ;i < V.length - 1 ;i++){
            int index = i;
            for (int j = i + 1 ;j < V.length ;j++)
                if (V[j] < V[index])
                    index = j;

            float smallerNumber = V[index];
            V[index] = V[i];
            V[i] = smallerNumber;
        }
        return V;
    }

    public static void Calc_Sub_Matrixs(float[][] M) throws Exception {
        if(Is_Square_Matrix(M)){
            int n = M.length;
            Print_Matrix(M);
            System.out.println();
            for(int i = 0 ;i < n ;i++) {
                for(int j = 0 ;j < n ;j++) {
                    float[][] sub = Sub_Matrix(M,i,j);
                    Print_Matrix(sub);
                    System.out.println();
                }
            }
        } else {
            throw new Exception("The Matrix is not Square!!!");
        }
    }

    //////////////////////////////////

    private static boolean Is_Orthogonal(float[][] M) throws Exception {
        int n = M.length;
        if (Mult_Mats(M ,Transpose(M)) == Unit_Matrix(n))
            return true;
        return false;
    }

    public static double normVector(Vector<Double> v) {
        double sum = 0;
        for(int i = 0 ;i < v.size() ;i++)
            sum += Math.pow(v.elementAt(i), 2);
        return Math.sqrt(sum);
    }

    /*public static Vector<Vector<Double>> Gram_Schmidt(Vector<Vector<Float>> vSet){
        int n = vSet.size();
        float sum = 0;
        float norm;
        Vector<Vector<Float>> uSet = new Vector<Vector<Float>>(n);
        Vector<Vector<Float>> e = new Vector<Vector<Float>>(n);
        for(int i = 0 ;i < vSet.size() ;i++){
            norm = normVector(vSet.elementAt(i));
            Vector<Float> helper = vSet.elementAt(i);
            if(i == 0){
                for(int j = 0 ;j < helper.size() ;j++){
                    helper.elementAt(j) = helper.elementAt(j)/norm;
                }
            }else{

            }
            uSet.add(vSet.elementAt(i));
        }
        return uSet;
    }*/

    public static float Determinant_Up_To_3x3(float[][] M) {
        int n = M.length;
        float det = (n > 1) ? 0 : M[0][0];
        for (int i = 0 ;i < n ;i++) {
            float m_d = 1 ,s_d = -1;
            for (int j = 0 ;j < n ;j++) {
                int dl = j - i ,dr = j + i;
                if (dl < 0 && dr >= n) {
                    m_d *= M[dr - n][n + dl];
                    s_d *= M[2*n - 1 - dr][j];
                } else if (dl < 0) {
                    m_d *= M[dr][n + dl];
                    s_d *= M[n - 1 - dr][j];
                } else if (dr >= n) {
                    m_d *= M[dr - n][dl];
                    s_d *= M[2*n - 1 - dr][j];
                } else {
                    m_d *= M[dr][dl];
                    s_d *= M[n - 1 - dr][j];
                }
            }
            det += m_d + s_d;
        }
        return det;
    }

    private static boolean Is_Linear_Dependent(float[][] M) {
        boolean flag = true;
        int n = M[0].length;
        float[] rats = new float[n];
        for (int i1 = 0 ;i1 < n - 1 ;i1++) {
            float[] v1 = Get_Row_From_Matrix(M,i1);
            for (int i2 = i1 + 1 ;i2 < n ;i2++) {
                float[] v2 = Get_Row_From_Matrix(M,i2);
                for (int j = 0 ;j < n ;j++)
                    rats[j] = v1[j] / v2[j];
                if (Is_Equals_Values(rats)) {
                    flag = true;
                    break;
                }
                else
                    flag = false;
            }
            if (flag)
                break;
        }
        return flag;
    }

    private static float[] Get_Row_From_Matrix(float[][] M ,int i) {
        int n = M[0].length;
        float[] v = new float[n];
        for (int j = 0 ;j < n ;j++)
            v[j] = M[i][j];
        return v;
    }

    private static boolean Is_Equals_Values(float[] v) {
        for(int i = 0 ;i < v.length - 1 ;i++)
            if(v[i] != v[i + 1])
                return false;
        return true;
    }

    public static void main(String[] args) {
        Scanner sc = new Scanner(System.in);
        System.out.println("Press Number Function to Perform:");
        System.out.println("1.Invertible Matrix");
        System.out.println("2.Adjoint Matrix");
        System.out.println("3.Transpose Matrix");
        System.out.println("4.Determinant of Matrix");
        System.out.println("5.Matrixs Multiplication");
        System.out.println("6.Norm of Matrix (Normal ,Infinite Norm ,One Norm)");
        System.out.println("7.Condition Number of Matrix (Normal ,Infinite Norm ,One Norm)");
        System.out.println("8.Sub Matrixs for Matrix");
        System.out.println("9.Special Method for Determinant up to 3x3");
        System.out.println("10.Linear Dependent Vectors in Matrix");
        int func = sc.nextInt();
        int n;
        if (func == 1) {
            System.out.print("Enter Number Variables:");
            n = sc.nextInt();
        } else {
            System.out.print("Enter Size Matrixs:");
            n = sc.nextInt();
        }
        float[][] M = Create_Matrix(n);
        switch (func){
            case 1:
                Print_Matrix(Invertible(M));
                break;
            case 2:
                Print_Matrix(Adjoint(M));
                break;
            case 3:
                Print_Matrix(Transpose(M));
                break;
            case 4:
                System.out.println(Determinant(M));
                break;
            case 5:
                try {
                    float[][] M1 = new float[n][n];
                    float[][] M2 = new float[n][n];
                    System.out.println("enter values to first matrix:");
                    for(int i = 0 ;i < M1.length ;i++)
                        for(int j = 0 ;j < M1[0].length ;j++){
                            System.out.print("M1["+i+"]["+j+"]->");
                            M1[i][j] = sc.nextFloat();
                        }
                    System.out.println("enter values to second matrix:");
                    for(int i = 0 ;i < M2.length ;i++)
                        for(int j = 0 ;j < M2[0].length ;j++){
                            System.out.print("M2["+i+"]["+j+"]->");
                            M2[i][j] = sc.nextFloat();
                        }
                    Print_Matrix(Mult_Mats(M1,M2));
                    break;
                } catch (Exception e1) {
                    e1.printStackTrace();
                }
            case 6:
                float norm = Norm_Standard(M);
                float normInf = Norm_Infinite(M);
                float normOne = Norm_One(M);
                System.out.println("Normal: " + norm);
                System.out.println("Infinite Norm : " + normInf);
                System.out.println("One Norm : " + normOne);
                break;
            case 7:
                float conNum = Condition_Number(M);
                float conInf = Condition_Number_Infinite(M);
                float conOne = Condition_Number_One(M);
                System.out.println("Normal: " + conNum);
                System.out.println("Infinite Norm : " + conInf);
                System.out.println("One Norm : " + conOne);
                break;
            case 8:
                try {
                    Calc_Sub_Matrixs(M);
                } catch (Exception e) {
                    e.printStackTrace();
                }
                break;
            case 9:
                System.out.println(Determinant_Up_To_3x3(M));
                break;
            case 10:
                System.out.println(Is_Linear_Dependent(M));
                break;
            default:
                try {
                    throw new IOException("error: indFunc out of the range");
                } catch (IOException e) {
                    e.printStackTrace();
                }
        }
        // |{{1,-3,1},{2,4,-5},{5,5,1}}| = 100
    }
}
