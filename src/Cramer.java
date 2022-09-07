import java.util.Vector;

public class Cramer {

    public static void Print_Matrix(float[][] mat) {
        for (int i = 0 ;i < mat.length ;i++) {
            for (int j = 0 ;j < mat[0].length ;j++) {
                System.out.print(Math.round(mat[i][j] * 10000.0) / 10000.0 + " ");
            }
            System.out.println();
        }
    }

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

    private static boolean Is_Square_Matrix(float[][] mat) {
        if (mat.length != mat[0].length)
            return false;
        return true;
    }

    private static boolean Is_Length_Equals(float[][] mat ,float[] arr){
        return mat.length == arr.length && mat[0].length == arr.length;
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

    public static Vector<Float> Cramer_Method(float[][] A, float[] b){
        /* solve equation Ax = b (x is vector variables ,
		b is vector leading coefficients and A is matrix coefficients)*/
        int n = A.length;
        Vector<Float> v = new Vector<Float>();
        float det = Determine(A);
        Print_Matrix(A);
        System.out.println();
        System.out.println("det(A) = " + det);
        System.out.println();
        float[] dets = new float[n];
        float[] helper = new float[n];
        for(int j = 0 ;j < n ;j++){
            for(int i = 0 ;i < n ;i++){
                helper[i] = A[i][j];
                A[i][j] = b[i];
            }
            dets[j] = Determine(A);
            Print_Matrix(A);
            System.out.println();
            System.out.println("det(A" + (j + 1) + ") = " + dets[j]);
            for(int i = 0 ;i < n ;i++)
                A[i][j] = helper[i];
            v.add(dets[j]/det);
            System.out.println("x" + (j + 1) + " = " + v.elementAt(j));
            System.out.println();
        }
        return v;
    }

    public static void Check_System(float[][] A, float[] b) throws Exception{
        if(Is_Square_Matrix(A) && Is_Length_Equals(A,b)){
            int n = A.length;
            if(n > 1){
                if(Determine(A) == 0) {
                    if(Is_Zero_Vector(b)) {
                        /* solve equation Ax = 0 (x is vector variables ,
						b is zero vector and A is matrix coefficients)*/
                        System.out.println("det(A) = 0 ,but system homogeneus");
                        A = Ranking_Matrix(A);
                        for(int i = 0 ;i < n - 1 ;i++)
                            b[i] = (-1) * A[i][n - 1];
                        A = Sub_Matrix(A,n - 1,n - 1);
                        b = Sub_Vector(b,n - 1);
                        Vector<Float> v = Cramer_Method(A,b);
                        System.out.println("x" + (n) + " = " + 1.0 + " (Posting Guess)");
                        v.add((float) 1.0);
                        System.out.println();
                        System.out.print("Solution is: x = k * " + v + " For any k is Real");
                    }else{
                        throw new Exception("not exists solutions: det(A) = 0");
                    }
                }else{
					/* solve equation Ax = b (x is vector variables ,
						b is vector leading coefficients and A is matrix coefficients)*/
                    Vector<Float> v = Cramer_Method(A,b);
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
        float[][] D = {{6,1,3,3,-11,1},{11,-6,5,11,0,-4},{-2,-2,-4,2,-3,0},{2,12,-1,-7,3,0},{5,-11,-11,8,-8,-2},{3,2,2,-1,-1,1}};
        float[] d = {-1,41,-30,-8,-27,18};
        float[][] F = {{-2,0,-4,-2,5,-6},{6,-1,-2,0,1,0},{6,0,0,-2,5,2},{-4,0,6,4,4,6},{-4,-2,3,-6,4,-4},{5,1,-5,-1,5,-2}};
        float[] f = {0,0,0,0,0,0};
        System.out.println("Solve The Next Exercise:");
        Print_Equation(B,b);
        System.out.println();
        try {
            Check_System(B,b);
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

}
