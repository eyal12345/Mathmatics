import java.util.Vector;

public class Ranking {

    public static void Print_Equation(float[][] A ,float[] b) {
        if(A[0].length == 1){
            System.out.println(A[0][0] + "*x = " + b[0]);
        }else{
            for (int i = 0 ;i < A.length ;i++) {
                for (int j = 0 ;j < A[0].length ;j++) {
                    System.out.print(A[i][j] + "*x" + (j + 1));
                    if(j == A[0].length - 1)
                        System.out.print(" = ");
                    else
                        System.out.print(" + ");
                }
                System.out.println(b[i]);
            }
        }
    }

    private static void Print_System(float[][] mat ,float[] vec) {
        for (int i = 0 ;i < mat.length ;i++) {
            for (int j = 0 ;j < mat[0].length ;j++) {
                System.out.print(Math.round(mat[i][j] * 10000.0) / 10000.0 + " ");
            }
            System.out.println("| " + Math.round(vec[i] * 10000.0) / 10000.0);
        }
        System.out.println();
    }

    private static boolean Is_Square_Matrix(float[][] mat) {
        if (mat.length != mat[0].length)
            return false;
        return true;
    }

    private static boolean Is_Length_Equals(float[][] mat ,float[] arr){
        return mat.length == arr.length && mat[0].length == arr.length;
    }

    private static boolean Is_Upper_Triangular(float[][] mat){
        for (int i = 0 ;i < mat.length ;i++)
            for (int j = 0 ;j < mat[0].length ;j++)
                if (i > j)
                    if (mat[i][j] != 0)
                        return false;
        return true;
    }

    private static boolean Is_Lower_Triangular(float[][] mat){
        for (int i = 0 ;i < mat.length ;i++)
            for (int j = 0 ;j < mat[0].length ;j++)
                if (i < j)
                    if (mat[i][j] != 0)
                        return false;
        return true;
    }

    public static float Determine(float[][] mat) {
        if (mat.length == 1)
            return mat[0][0];
        float sum = 0;
        int k = 1;
        for (int i = 0 ;i < mat.length ;i++) {
            sum += mat[0][i] * Determine(Sub_Matrix(mat, 0, i)) * k;
            k = k * (-1);
        }
        return sum;
    }

    public static float[][] Sub_Matrix(float[][] mat, int x, int y) {
        int n = mat.length ,p = 0 ,q = 0;
        float[][] subMat = new float[n - 1][n - 1];
        for (int i = 0 ;i < mat.length ;i++)
            for (int j = 0 ;j < mat[0].length ;j++)
                if (i != x && j != y) {
                    subMat[p][q] = mat[i][j];
                    q++;
                    if (q == subMat[0].length) {
                        p++;
                        q = 0;
                    }
                }
        return subMat;
    }

    /////////////////////////////////////////////////////////////////

    private static float[] Upper_Ranking(float[][] mat ,float[] vec){
        System.out.println("Upper Ranking:");
        System.out.println();
        int n = mat.length;
        boolean flag = false;
        Print_System(mat,vec);
        for (int i = 0 ;i < n - 1 ;i++) {
            if(Is_Lower_Triangular(mat)){
                System.out.println("R" + (i + 1) + " --> " + (1/mat[i][i]) + "*R" + (i + 1));
                System.out.println();
                vec[i] /= mat[i][i];
                mat[i][i] = 1;
                Print_System(mat,vec);
            }
            flag = true;
            for (int j = i + 1 ;j < n ;j++) {
                System.out.println("row = " + (j + 1) + " ,col = " + (i + 1));
                if(mat[i][i] == 0){
                    System.out.println("R" + (i + 1) + " <--> R" + (j + 1));
                    Retreat_Rows(mat,i,j);
                    Swap(vec,i,j);
                }
                float c = mat[j][i] / mat[i][i];
                if(c > 0)
                    System.out.println("R" + (j + 1) + " --> R" + (j + 1) + " - " + c + "*R" + (i + 1));
                if(c < 0)
                    System.out.println("R" + (j + 1) + " --> R" + (j + 1) + " + " + (-1)*c + "*R" + (i + 1));
                System.out.println();
                for (int k = 0 ;k < n ;k++) {
                    mat[j][k] = mat[j][k] - mat[i][k] * c;
                    if(flag){
                        vec[j] = vec[j] - vec[i] * c;
                        flag = false;
                    }
                }
                Print_System(mat,vec);
                flag = true;
            }
        }
        if(Is_Upper_Triangular(mat) && Is_Lower_Triangular(mat)){
            System.out.println("R" + (n) + " --> " + (1/mat[n - 1][n - 1]) + "*R" + (n));
            System.out.println();
            vec[n - 1] /= mat[n - 1][n - 1];
            mat[n - 1][n - 1] = 1;
            Print_System(mat,vec);
        }else{
            vec = Lower_Ranking(mat, vec);
        }
        return vec;
    }

    private static float[] Lower_Ranking(float[][] mat ,float[] vec){
        System.out.println("Lower Ranking:");
        System.out.println();
        int n = mat.length;
        boolean flag = false;
        Print_System(mat,vec);
        for (int i = n - 1 ;i > 0 ;i--) {
            if(Is_Upper_Triangular(mat)){
                System.out.println("R" + (i + 1) + " --> " + (1/mat[i][i]) + "*R" + (i + 1));
                System.out.println();
                vec[i] /= mat[i][i];
                mat[i][i] = 1;
                Print_System(mat,vec);
            }
            flag = true;
            for (int j = i - 1 ;j >= 0 ;j--) {
                System.out.println("row = " + (j + 1) + " ,col = " + (i + 1));
                if(mat[i][i] == 0){
                    System.out.println("R" + (i + 1) + " <--> R" + (j + 1));
                    Retreat_Rows(mat,i,j);
                    Swap(vec,i,j);
                }
                float c = mat[j][i] / mat[i][i];
                if(c > 0)
                    System.out.println("R" + (j + 1) + " --> R" + (j + 1) + " - " + c + "*R" + (i + 1));
                if(c < 0)
                    System.out.println("R" + (j + 1) + " --> R" + (j + 1) + " + " + (-1)*c + "*R" + (i + 1));
                System.out.println();
                for (int k = n - 1 ;k >= 0 ;k--) {
                    mat[j][k] = mat[j][k] - mat[i][k] * c;
                    if(flag){
                        vec[j] = vec[j] - vec[i] * c;
                        flag = false;
                    }
                }
                Print_System(mat,vec);
                flag = true;
            }
        }
        if(Is_Upper_Triangular(mat) && Is_Lower_Triangular(mat)){
            System.out.println("R" + 1 + " --> " + (1/mat[0][0]) + "*R" + 1);
            System.out.println();
            vec[0] /= mat[0][0];
            mat[0][0] = 1;
            Print_System(mat,vec);
        }else{
            vec = Upper_Ranking(mat, vec);
        }
        return vec;
    }

    /////////////////////////////////////////////////////////////////

    private static void Swap(float[]vec ,int i ,int j){
        float t = vec[i];
        vec[i] = vec[j];
        vec[j] = t;
    }

    private static void Retreat_Rows(float[][] mat ,int r1 ,int r2){
        for(int j = 0 ;j < mat[0].length ;j++){
            float t = mat[r1][j];
            mat[r1][j] = mat[r2][j];
            mat[r2][j] = t;
        }
    }

    private static float[][] Copy_Matrix(float[][] mat){
        float[][] copy = new float[mat.length][mat.length];
        for (int i = 0 ;i < mat.length ;i++)
            for (int j = 0 ;j < mat[0].length ;j++)
                copy[i][j] = mat[i][j];
        return copy;
    }

    public static float[][] Ranking_Matrix(float[][] mat){
        int n = mat.length;
        float[][] matRank = Copy_Matrix(mat);
        for (int i = 0 ;i < n - 1 ;i++) {
            for (int j = i + 1 ;j < n ;j++) {
                float c = matRank[j][i] / matRank[i][i];
                for (int k = 0 ;k < n ;k++) {
                    matRank[j][k] = matRank[j][k] - matRank[i][k] * c;
                }
            }
        }
        return matRank;
    }

    public static float[] Sub_Vector(float[] vec, int x){
        int n = vec.length;
        float[] subVec = new float[n - 1];
        for (int i = 0 ;i < vec.length ;i++)
            if (i != x){
                subVec[i] = vec[i];
            }
        return subVec;
    }

    private static boolean Is_Zero_Vector(float[] vec){
        int n = vec.length;
        boolean flag = true;
        for(int i = 0; i < n; i++){
            if(vec[i] == 0) flag = true;
            else{
                flag = false;
                break;
            }
        }
        return flag;
    }

    private static float[] Ranking_Method(float[][] A ,float[] b){
        System.out.println("Transform A Matrix To I (Unit Matrix) By Ranking:");
        System.out.println();
        int n = A.length;
        boolean flag = false;
        Print_System(A,b);
        for (int i = 0 ;i < n ;i++) {
            for (int j = 0 ;j < n ;j++) {
                if(A[i][i] == 0){
                    System.out.println("R" + (i + 1) + " <--> R" + (j + i + 2));
                    Retreat_Rows(A,i,j + i + 1);
                    Swap(b,i,j + i + 1);
                    System.out.println();
                    Print_System(A,b);
                }
                if(i != j){
                    flag = true;
                    System.out.println("row = " + (j + 1) + " ,col = " + (i + 1));
                    float c = A[j][i] / A[i][i];
                    if(c > 0)
                        System.out.println("R" + (j + 1) + " --> R" + (j + 1) + " - " + c + "*R" + (i + 1));
                    if(c < 0)
                        System.out.println("R" + (j + 1) + " --> R" + (j + 1) + " + " + (-1)*c + "*R" + (i + 1));
                    System.out.println();
                    for (int k = 0 ;k < n ;k++) {
                        A[j][k] = A[j][k] - A[i][k] * c;
                        if(i == k && A[j][k] != 0)
                            A[j][k] = 0;
                        if(flag){
                            b[j] = b[j] - b[i] * c;
                            flag = false;
                        }
                    }
                    Print_System(A,b);
                }else {
                    continue;
                }
            }
        }
        if(Is_Upper_Triangular(A) && Is_Lower_Triangular(A)){
            for(int i = 0 ;i < n ;i++){
                System.out.println("R" + (i + 1) + " --> " + (1/A[i][i]) + "*R" + (i + 1));
                System.out.println();
                b[i] /= A[i][i];
                A[i][i] = 1;
                Print_System(A,b);
            }
        }
        return b;
    }

    public static void Check_System(float[][]A ,float[]b) throws Exception{
        if(Is_Square_Matrix(A) && Is_Length_Equals(A,b)){
            int n = b.length;
            Vector<Float> v = new Vector<Float>();
            if(n > 1){
                if(Determine(A) == 0) {
                    if(Is_Zero_Vector(b)){
                        /* solve equation Ax = 0 (x is vector variables ,
						b is zero vector and A is matrix coefficients)*/
                        System.out.println("det(A) = 0 ,but system homogeneus");
                        A = Ranking_Matrix(A);
                        for(int i = 0 ;i < n - 1 ;i++)
                            b[i] = (-1) * A[i][n - 1];
                        System.out.println("solve the next system: (suppose that: x" + (n) + " = " + 1 + ")");
                        System.out.println();
                        b = Ranking_Method(Sub_Matrix(A,n - 1,n - 1),Sub_Vector(b,n - 1));
                        for(int i = 0 ;i < n - 1 ;i++)
                            v.add((float) (Math.round(b[i] *10000.0)/10000.0));
                        v.add((float) 1.0);
                        System.out.print("Solution is: x = k * " + v + " For any k is Real");
                    }else {
                        throw new Exception("not exists solutions: det(A) = 0");
                    }
                }else{
                    b = Ranking_Method(A,b);
                    for(int i = 0 ;i < n ;i++)
                        v.add((float) (Math.round(b[i] *10000.0)/10000.0));
                    System.out.print("Solution is: x = " + v);
                }
            }else if(n == 1 && A[0][0] != 0){
                // solve equation ax = b (x is variable and a,b fixed)
                System.out.print("Solution is: x = " + b[0]/A[0][0]);
            }else{
                // checking if 0x = c ? (c fixed)
                if(b[0] == 0 && A[0][0] == 0)
                    throw new Exception("exists infinite number of solutions to the equation");
                else
                    throw new Exception("not exists any solutions to the equation");
            }
        }else{
            throw new Exception("Error: The input is not exist conditions to linear equations");
        }
    }

    public static void main(String[] args) {
        float[][] A = {{2,1,-1},{-3,-1,2},{-2,1,2}};
        float[] a = {8,-11,-3};
        float[][] B = {{-2,2,2,-1,1},{-4,4,4,4,3},{2,3,2,3,2},{-3,-1,1,2,2},{5,5,3,5,5}};
        float[] b = {1,31,15,8,52};
        float[][] C = {{-2,3,3,-2},{-1,4,2,-2},{1,3,1,3},{-3,-2,4,-5}};
        float[] c = {8,5,19,-19};
        float[][] D = {{2,10,7,-5,0,8},{-12,-8,-9,-1,-11,-5},{-1,8,0,-4,-9,11},{12,6,-10,5,-11,-12},{12,11,-10,-2,-4,-12},{5,-5,3,6,8,-10}};
        float[] d = {-45,-13,-23,36,-43,28};
        float[][] F = {{-2,0,-4,-2,5,-6},{6,-1,-2,0,1,0},{6,0,0,-2,5,2},{-4,0,6,4,4,6},{-4,-2,3,-6,4,-4},{5,1,-5,-1,5,-2}};
        float[] f = {0,0,0,0,0,0};
        System.out.println("Solve The Next Exercise:");
        Print_Equation(F,f);
        System.out.println();
        try {
            Check_System(F,f);
        } catch (Exception e) {
            e.printStackTrace();
        }
    }
}
