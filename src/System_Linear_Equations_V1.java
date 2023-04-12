import java.util.Scanner;
import java.util.Arrays;

public class System_Linear_Equations_V1 {

    /////////////////////////////////////////////// Print Methods /////////////////////////////////////////////////
    // display the system Ax = b in the linear equations format
    public static void Display_Exercise(float[][] A, float[] b) {
        int n = A.length;
        if (n == 1) {
            System.out.println("solve the next equation (" + 1 + " variable):");
            if (A[0][0] % 1 == 0) {
                if (A[0][0] == 1) {
                    if (b[0] % 1 == 0) {
                        System.out.println("x = " + (int) b[0]);
                    } else {
                        System.out.println("x = " + b[0]);
                    }
                } else if (A[0][0] == -1) {
                    if (b[0] % 1 == 0) {
                        System.out.println("-x = " + (int) b[0]);
                    } else {
                        System.out.println("-x = " + b[0]);
                    }
                } else {
                    if (b[0] % 1 == 0) {
                        System.out.println((int) A[0][0] + "x = " + (int) b[0]);
                    } else {
                        System.out.println((int) A[0][0] + "x = " + b[0]);
                    }
                }
            } else {
                if (b[0] % 1 == 0) {
                    System.out.println(A[0][0] + "x = " + (int) b[0]);
                } else {
                    System.out.println(A[0][0] + "x = " + b[0]);
                }
            }
        } else {
            System.out.println("solve the next system (" + n + " variables):");
            for (int i = 0; i < n; i++) {
                System.out.print("eq" + (i + 1) + ": ");
                for (int j = 0; j < n; j++) {
                    if (A[i][j] > 0) {
                        System.out.print("+ ");
                        if (Math.abs(A[i][j]) == 1) {
                            System.out.print("  x" + (j + 1) + " ");
                        } else if (A[i][j] % 1 == 0) {
                            System.out.print((int) A[i][j] + "*x" + (j + 1) + " ");
                        } else {
                            System.out.print(A[i][j] + "*x" + (j + 1) + " ");
                        }
                    } else if (A[i][j] < 0) {
                        System.out.print("- ");
                        if (Math.abs(A[i][j]) == 1) {
                            System.out.print("  x" + (j + 1) + " ");
                        } else if (A[i][j] % 1 == 0) {
                            System.out.print((int) Math.abs(A[i][j]) + "*x" + (j + 1) + " ");
                        } else {
                            System.out.print(Math.abs(A[i][j]) + "*x" + (j + 1) + " ");
                        }
                    } else {
                        System.out.print(" ".repeat(7));
                    }
                }
                if (b[i] % 1 == 0) {
                    System.out.println(" = " + (int) (Math.round(b[i] * 1000.0) / 1000.0));
                } else {
                    System.out.println(" = " + (Math.round(b[i] * 1000.0) / 1000.0));
                }
            }
        }
        System.out.println();
    }

    // display current status of the system Ax = b each time of iteration on an element
    public static void Print_Status_System(float[][] A, float[] b, String fn) {
        int n = A.length;
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                if ((Math.round(A[i][j] * 1000.0) / 1000.0) % 1 == 0) {
                    System.out.print((int) (Math.round(A[i][j] * 1000.0) / 1000.0) + " ");
                } else if (fn.equals("d")) {
                    System.out.print(Math.round(A[i][j] * 1000.0) / 1000.0 + " ");
                } else if (fn.equals("r")) {
                    System.out.print(convertDecimalToFraction(A[i][j]) + " ");
                } if (j != n - 1) {
                    System.out.print(",");
                }
            }
            System.out.print("| ");
            if ((Math.round(b[i] * 1000.0) / 1000.0) % 1 == 0) {
                System.out.print((int) (Math.round(b[i] * 1000.0) / 1000.0));
            } else if (fn.equals("d")) {
                System.out.print(Math.round(b[i] * 1000.0) / 1000.0);
            } else if (fn.equals("r")) {
                System.out.print(convertDecimalToFraction(b[i]));
            }
            System.out.println();
        }
        System.out.println();
    }

    // display user interface by selection method for solution
    public static void User_Menu_System() {
        System.out.println("choose number method to solution:");
        System.out.println("1. numeric method");
        System.out.println("2. cramer method");
        System.out.println("3. forward backward method");
        System.out.println("4. upper -> lower ranking method");
        System.out.println("5. lower -> upper ranking method");
        System.out.println("6. parallel ranking method (iterative method)");
        System.out.println("7. parallel ranking method (recursive method)");
    }

    // display user interface by selection format for solution
    public static void User_Menu_Solution() {
        System.out.println("choose the character of format to representation of solution:");
        System.out.println("d. decimal");
        System.out.println("r. rational");
    }

    // display a matrix each current status
    public static void Print_Matrix(float[][] A, String fn) {
        int n = A.length;
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                if ((Math.round(A[i][j] * 1000.0) / 1000.0) % 1 == 0) {
                    System.out.print((int) (Math.round(A[i][j] * 1000.0) / 1000.0) + " ");
                } else if (fn.equals("d")) {
                    System.out.print(Math.round(A[i][j] * 1000.0) / 1000.0 + " ");
                } else if (fn.equals("r")) {
                    System.out.print(convertDecimalToFraction(A[i][j]) + " ");
                } if (j != n - 1) {
                    System.out.print(",");
                }
            }
            System.out.println();
        }
        System.out.println();
    }

    // display a vector each current status
    public static void Print_Vector(float[] b, String fn) {
        int n = b.length;
        for (int i = 0; i < n; i++) {
            if ((Math.round(b[i] * 1000.0) / 1000.0) % 1 == 0) {
                System.out.println((int) (Math.round(b[i] * 1000.0) / 1000.0));
            } else if (fn.equals("d")) {
                System.out.println(Math.round(b[i] * 1000.0) / 1000.0);
            } else if (fn.equals("r")) {
                System.out.println(convertDecimalToFraction(b[i]));
            }
        }
        System.out.println();
    }

    // show the resulting solution as a vector representation
    public static void Print_Solution(float[] x, String fn) {
        int n = x.length;
        if (n == 1) {
            if ((Math.round(x[0] * 1000.0) / 1000.0) % 1 == 0) {
                System.out.println("x = " + (int) (Math.round(x[0] * 1000.0) / 1000.0));
            } else if (fn.equals("d")) {
                System.out.println("x = " + Math.round(x[0] * 1000.0) / 1000.0);
            } else if (fn.equals("r")) {
                System.out.println("x = " + convertDecimalToFraction(x[0]));
            }
        } else {
            String s = "( ";
            for (int i = 0; i < n; i++) {
                if (i == n - 1) {
                    if ((Math.round(x[i] * 1000.0) / 1000.0) % 1 == 0) {
                        s += (int) (Math.round(x[i] * 1000.0) / 1000.0) + " )";
                    } else if (fn.equals("d")) {
                        s += Math.round(x[i] * 1000.0) / 1000.0 + " )";
                    } else if (fn.equals("r")) {
                        s += convertDecimalToFraction(x[i]) + " )";
                    }
                } else {
                    if ((Math.round(x[i] * 1000.0) / 1000.0) % 1 == 0) {
                        s += (int) (Math.round(x[i] * 1000.0) / 1000.0) + " , ";
                    } else if (fn.equals("d")) {
                        s += Math.round(x[i] * 1000.0) / 1000.0 + " , ";
                    } else if (fn.equals("r")) {
                        s += convertDecimalToFraction(x[i]) + " , ";
                    }
                }
            }
            System.out.println(s);
        }
    }

    /////////////////////////////////////////// Auxiliary Operations /////////////////////////////////////////////
    // check if a matrix is zero matrix
    public static boolean Is_Zero_Matrix(float[][] A) {
        int n = A.length;
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                if (A[i][j] != 0) {
                    return false;
                }
            }
        }
        return true;
    }

    // check if a vector is zero vector
    public static boolean Is_Zero_Vector(float[] b) {
        int n = b.length;
        for (int i = 0; i < n; i++) {
            if (b[i] != 0) {
                return false;
            }
        }
        return true;
    }

    // check if a matrix is a unit matrix
    public static boolean Is_Unit_Matrix(float[][] A) {
        int n = A.length;
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                if (A[i][i] != 1 || (i != j && A[i][j] != 0)) {
                    return false;
                }
            }
        }
        return true;
    }

    // check if the specific row in the matrix is a unit vector
    public static boolean Is_Unit_Vector(float[][] A, int k) {
        int n = A.length;
        boolean flag = true;
        for (int i = 0; i < n && flag; i++) {
            for (int j = 0; j < n && A[k][i] != 0; j++) {
                if (j != i && A[k][j] != 0) {
                    flag = false;
                    break;
                }
            }
        }
        return flag;
    }

    // check if all vector values is equal to each other
    public static boolean Is_Equals_Values(float[] r) {
        for (int i = 0; i < r.length - 1; i++) {
            if (r[i] != r[i + 1]) {
                return false;
            }
        }
        return true;
    }

    // check if the matrix is an upper triangular
    public static boolean Is_Upper_Triangular(float[][] A) {
        int n = A.length;
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < i; j++) {
                if (A[i][j] != 0) {
                    return false;
                }
            }
        }
        return true;
    }

    // check if the matrix is a lower triangular
    public static boolean Is_Lower_Triangular(float[][] A) {
        int n = A.length;
        for (int i = 0; i < n - 1; i++) {
            for (int j = i + 1; j < n; j++) {
                if (A[i][j] != 0) {
                    return false;
                }
            }
        }
        return true;
    }

    // check if exists two vectors in the matrix which are linearly dependent
    public static boolean Is_Linear_Dependent_Rows(float[][] A) {
        int n = A.length;
        float[] R = new float[n];
        for (int i1 = 0; i1 < n - 1; i1++) {
            for (int i2 = i1 + 1; i2 < n; i2++) {
                for (int j = 0; j < n; j++) {
                    R[j] = A[i1][j] / A[i2][j];
                }
                if (Is_Equals_Values(R)) {
                    return true;
                }
            }
        }
        return false;
    }

    // calculate determinant of a matrix
    public static float Determinant(float[][] A) {
        int n = A.length;
        if (n == 1) {
            return A[0][0];
        }
        float sum = 0;
        for (int i = 0; i < n; i++) {
            sum += A[0][i] * Determinant(Sub_Matrix(A,0,i)) * Math.pow(-1,i);
        }
        return sum;
    }

    // calculate of sub-matrix from a matrix by cutting row "x" and column "y"
    public static float[][] Sub_Matrix(float[][] A, int x, int y) {
        int n = A.length ,p = 0 ,q = 0;
        float[][] subA = new float[n - 1][n - 1];
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                if (i != x && j != y) {
                    subA[p][q] = A[i][j];
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

    // calculate of sub-vector from a vector without index "x"
    public static float[] Sub_Vector(float[] b, int x) {
        int n = b.length;
        float[] subb = new float[n - 1];
        for (int i = 0; i < n; i++) {
            if (i != x) {
                subb[i] = b[i];
            }
        }
        return subb;
    }

    // calculate adjoint matrix of a matrix
    public static float[][] Adjoint(float[][] A) {
        int n = A.length;
        float[][] Adj = new float[n][n];
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                Adj[j][i] = (float) (Math.pow(-1,i + j) * Determinant(Sub_Matrix(A,i,j)));
            }
        }
        return Adj;
    }

    // calculate invertible matrix of a matrix
    public static float[][] Invertible(float[][] A) {
        int n = A.length;
        float[][] invA = new float[n][n];
        float[][] adj = Adjoint(A);
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                invA[i][j] = (1 / Determinant(A)) * adj[i][j];
            }
        }
        return invA;
    }

    // calculate multiplication between two matrices provided that M1's length column is equal to M2's length row
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

    // calculate ranking of a matrix by upper triangular
    public static float[][] Upper_Ranking_Matrix(float[][] A, String fn) {
        int n = A.length;
        Print_Matrix(A,fn);
        float[][] rA = Copy_Matrix(A);
        for (int i = 0; i < n - 1; i++) {
            for (int j = i + 1; j < n; j++) {
                if (rA[i][i] == 0) {
                    System.out.println("R" + (i + 1) + " <--> R" + (i + 1) % n + 1 + "\n");
                    Retreat_Rows_Matrix(rA,i,(i + 1) % n);
                }
                if (rA[j][i] != 0) {
                    float c = rA[j][i] / rA[i][i];
                    Sum_Elementary_Action(c,j,i,fn);
                    for (int k = 0; k < n; k++) {
                        rA[j][k] = rA[j][k] - rA[i][k] * c;
                    }
                    Print_Matrix(rA,fn);
                }
            }
        }
        return rA;
    }

    // replace between two rows in a system Ax = b
    public static void Retreat_Rows_System(float[][] A, float[] b, int r1, int r2) {
        int n = A.length;
        for (int j = 0; j < n; j++) {
            float t = A[r1][j];
            A[r1][j] = A[r2][j];
            A[r2][j] = t;
        }
        float t = b[r1];
        b[r1] = b[r2];
        b[r2] = t;
    }

    // replace between two rows in a matrix
    public static void Retreat_Rows_Matrix(float[][] A, int r1, int r2) {
        int n = A.length;
        for (int j = 0; j < n; j++) {
            float t = A[r1][j];
            A[r1][j] = A[r2][j];
            A[r2][j] = t;
        }
    }

    // duplicate the matrix values into a new matrix
    public static float[][] Copy_Matrix(float[][] A) {
        int n = A.length;
        float[][] copyA = new float[n][n];
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                copyA[i][j] = A[i][j];
            }
        }
        return copyA;
    }

    // get the index starting from the specific column in the matrix which are him value not equal to 0
    public static int Get_Index_UnZero_Value(float[][] A, int k) {
        int n = A.length;
        for (int i = k + 1; i < n + k; i++) {
            if (A[i % n][k] != 0) {
                return i % n;
            }
        }
        return -1;
    }

    // create zero vector with "n" cells
    public static float[] Zero_Vector(int n) {
        float[] x = new float[n];
        for (int k = 0; k < n; k++) {
            x[k] = 0;
        }
        return x;
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
    // show elementary actions for sum between rows in the system
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
                        System.out.println("R" + r + " --> R" + r + " + " + (int) (-k) + "*R" + c + "\n");
                    }
                } else if (fn.equals("d")) {
                    System.out.println("R" + r + " --> R" + r + " + " + (-k) + "*R" + c + "\n");
                } else if (fn.equals("r")) {
                    System.out.println("R" + r + " --> R" + r + " + " + convertDecimalToFraction(-k) + "*R" + c + "\n");
                }
            }
        }
    }

    // show elementary actions for multiplication of a row in the system
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

    // show elementary actions for replace between rows in the system
    public static void Retreat_Elementary_Action(int i, int j) {
        int r1 = i + 1, r2 = j + 1;
        if (r1 <= r2) {
            System.out.println("R" + r1 + " <--> R" + r2 + "\n");
        } else {
            System.out.println("R" + r2 + " <--> R" + r1 + "\n");
        }
    }

    /////////////////////////////////////////// Methods to Solution /////////////////////////////////////////////
    // solve system of linear equations Ax = b by invertible multiplication method: x = inv(A)*b
    public static float[] Invertible_Method(float[][] A, float[] b, String fn) {
        int n = A.length;
        float[] x = new float[n];
        float[][] invA = Invertible(A);
        System.out.println("the solution is: inv(A)*b:");
        Print_Matrix(invA,fn);
        Print_Vector(b,fn);
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                x[i] += b[j] * invA[i][j];
            }
            x[i] = (float)(Math.round(x[i] * 10000.0) / 10000.0);
        }
        return x;
    }

    // solve system of linear equations Ax = b by cramer method: x[i] = det(A[i])/det(A) for 1 <= i <= n
    public static float[] Cramer_Method(float[][] A, float[] b, String fn) {
        int n = A.length;
        float[] x = new float[n];
        float det = Determinant(A);
        Print_Matrix(A,fn);
        System.out.println("|A| = " + det);
        float[] h = new float[n];
        for (int j = 0; j < n; j++) {
            for (int i = 0; i < n; i++) {
                h[i] = A[i][j];
                A[i][j] = b[i];
            }
            float detj = Determinant(A);
            Print_Matrix(A,fn);
            System.out.println("|A" + (j + 1) + "| = " + detj);
            for (int i = 0; i < n; i++) {
                A[i][j] = h[i];
            }
            x[j] = detj / det;
            System.out.println("x" + (j + 1) + " = " + x[j] + "\n");
        }
        return x;
    }

    // solve system of linear equations Ax = b by forward backward method: Ly = b and then Ux = y
    public static float[] Forward_Backward_Method(float[][] A, float[] b, String fn) throws Exception {
        int n = b.length;
        System.out.println("first, we will upper ranking of A:");
        float[][] U = Upper_Ranking_Matrix(A,fn);
        Print_Matrix(U,fn);
        System.out.println("second, we will calculate lower ranking of A:");
        float[][] L = Mul_Mats(A, Invertible(U));
        Print_Matrix(L,fn);
        for (int i = 0; i < n - 1; i++) {
            for (int j = i + 1; j < n; j++) {
                if (L[i][i] == 0) {
                    Retreat_Rows_System(L,b,i,j);
                    break;
                }
            }
        }
        float[] x = new float[n];
        float[] y = new float[n];
        System.out.println("third, we will solve forward system Ly = b:");
        Print_Status_System(L,b,fn);
        for (int i = 0; i < n; i++) {
            y[i] = b[i];
            for (int j = 0; j < i; j++) {
                y[i] -= L[i][j] * y[j];
            }
            y[i] /= L[i][i];
        }
        System.out.println("finally, we will solve backward system Ux = y:");
        Print_Status_System(U,y,fn);
        for (int i = n - 1; i >= 0; i--) {
            x[i] = y[i];
            for (int j = i + 1; j < n; j++) {
                x[i] -= U[i][j] * x[j];
            }
            x[i] /= U[i][i];
            x[i] = (float)(Math.round(x[i] * 10000.0) / 10000.0);
        }
        return x;
    }

    // solve system of linear equations Ax = b by an upper ranking and then a lower ranking
    public static float[] Upper_Ranking_Method(float[][] A, float[] b, String fn) throws Exception {
        float det = Determinant(A);
        if (det != 0) {
            if (Is_Upper_Triangular(A) && Is_Lower_Triangular(A)) {
                System.out.println("A is already parallel triangular so now will be change directly to I:");
            } else if (Is_Upper_Triangular(A) && !Is_Lower_Triangular(A)) {
                System.out.println("A is already upper triangular so now we'll go directly to the lower ranking:");
                return Lower_Ranking_Method(A,b,fn);
            } else if (!Is_Upper_Triangular(A) && Is_Lower_Triangular(A)) {
                System.out.println("transform L matrix to I by a elementary actions:");
            } else {
                System.out.println("transform A matrix to U by a elementary actions:");
            }
            int n = A.length;
            for (int i = 0; i < n; i++) {
                if (A[i][i] == 0) {
                    int r = Get_Index_UnZero_Value(A,i);
                    Retreat_Elementary_Action(i,r);
                    Retreat_Rows_System(A,b,i,r);
                    Print_Status_System(A,b,fn);
                }
                if (Is_Unit_Vector(A,i)) {
                    float c = 1 / A[i][i];
                    Mul_Elementary_Action(c,i,fn);
                    b[i] /= A[i][i];
                    A[i][i] = 1;
                    if (c != 1) {
                        Print_Status_System(A,b,fn);
                    }
                }
                for (int j = i + 1; j < n; j++) {
                    float c = A[j][i] / A[i][i];
                    Sum_Elementary_Action(c,j,i,fn);
                    for (int k = 0; k < n; k++) {
                        A[j][k] -= A[i][k] * c;
                    }
                    A[j][i] = 0;
                    b[j] -= b[i] * c;
                    if (c != 0) {
                        Print_Status_System(A,b,fn);
                    }
                }
            }
            if (!Is_Upper_Triangular(A) || !Is_Lower_Triangular(A)) {
                System.out.print("and than ");
                return Lower_Ranking_Method(A,b,fn);
            }
            return b;
        } else {
            throw new Exception("this is a singular matrix");
        }
    }

    // solve system of linear equations Ax = b by a lower ranking and then an upper ranking
    public static float[] Lower_Ranking_Method(float[][] A, float[] b, String fn) throws Exception {
        float det = Determinant(A);
        if (det != 0) {
            if (Is_Upper_Triangular(A) && Is_Lower_Triangular(A)) {
                System.out.println("A is already parallel triangular so now will be change directly to I:");
            } else if (!Is_Upper_Triangular(A) && Is_Lower_Triangular(A)) {
                System.out.println("A is already lower triangular so now we'll go directly to the upper ranking:");
                return Upper_Ranking_Method(A,b,fn);
            } else if (Is_Upper_Triangular(A) && !Is_Lower_Triangular(A)) {
                System.out.println("transform U matrix to I by a elementary actions:");
            } else {
                System.out.println("transform A matrix to L by a elementary actions:");
            }
            int n = A.length;
            for (int i = n - 1; i >= 0; i--) {
                if (A[i][i] == 0) {
                    int r = n - 1 - Get_Index_UnZero_Value(A,i);
                    Retreat_Elementary_Action(i,r);
                    Retreat_Rows_System(A,b,i,r);
                    Print_Status_System(A,b,fn);
                }
                if (Is_Unit_Vector(A,i)) {
                    float c = 1 / A[i][i];
                    Mul_Elementary_Action(c,i,fn);
                    b[i] /= A[i][i];
                    A[i][i] = 1;
                    if (c != 1) {
                        Print_Status_System(A,b,fn);
                    }
                }
                for (int j = i - 1; j >= 0; j--) {
                    float c = A[j][i] / A[i][i];
                    Sum_Elementary_Action(c,j,i,fn);
                    for (int k = n - 1; k >= 0; k--) {
                        A[j][k] -= A[i][k] * c;
                    }
                    A[j][i] = 0;
                    b[j] -= b[i] * c;
                    if (c != 0) {
                        Print_Status_System(A,b,fn);
                    }
                }
            }
            if (!Is_Upper_Triangular(A) || !Is_Lower_Triangular(A)) {
                System.out.print("and than ");
                return Upper_Ranking_Method(A,b,fn);
            }
            return b;
        } else {
            throw new Exception("this is a singular matrix");
        }
    }

    // solve system of linear equations Ax = b by parallel ranking (iterative method)
    public static float[] Parallel_Ranking_Method(float[][] A, float[] b, String fn) throws Exception {
        float det = Determinant(A);
        if (det != 0) {
            System.out.println("transform A matrix to I by a elementary actions:");
            int n = A.length;
            float c;
            for (int i = 0; i < n; i++) {
                if (A[i][i] == 0) {
                    int r = Get_Index_UnZero_Value(A,i);
                    Retreat_Elementary_Action(i,r);
                    Retreat_Rows_System(A,b,i,r);
                    Print_Status_System(A,b,fn);
                }
                for (int j = 0; j < n; j++) {
                    if (i != j) {
                        c = A[j][i] / A[i][i];
                        Sum_Elementary_Action(c,j,i,fn);
                        for (int k = 0; k < n; k++) {
                            A[j][k] -= A[i][k] * c;
                        }
                        A[j][i] = 0;
                        b[j] -= b[i] * c;
                        if (c != 0) {
                            Print_Status_System(A,b,fn);
                        }
                    }
                    if (Is_Unit_Vector(A,j)) {
                        c = 1 / A[j][j];
                        Mul_Elementary_Action(c,j,fn);
                        b[j] /= A[j][j];
                        A[j][j] = 1;
                        if (c != 1) {
                            Print_Status_System(A,b,fn);
                        }
                    }
                }
            }
            return b;
        } else {
            throw new Exception("this is a singular matrix");
        }
    }

    // solve system of linear equations Ax = b by a parallel ranking (recursive method)
    public static float[] Parallel_Ranking_Method_Rec(float[][] A, float[] b, int i, int j, String fn) throws Exception {
        float det = Determinant(A);
        if (det != 0) {
            int n = A.length;
            if (Is_Unit_Matrix(A)) {
                return b;
            } else {
                if (A[i][i] == 0) {
                    int r = Get_Index_UnZero_Value(A,i);
                    Retreat_Elementary_Action(i,r);
                    Retreat_Rows_System(A,b,i,r);
                    Print_Status_System(A,b,fn);
                }
                if (i != j) {
                    float c = A[j][i] / A[i][i];
                    Sum_Elementary_Action(c,j,i,fn);
                    for (int k = 0; k < n; k++) {
                        A[j][k] -= A[i][k] * c;
                    }
                    A[j][i] = 0;
                    b[j] -= b[i] * c;
                    if (c != 0) {
                        Print_Status_System(A,b,fn);
                    }
                }
                if (Is_Unit_Vector(A,j)) {
                    float c = 1 / A[j][j];
                    Mul_Elementary_Action(c,j,fn);
                    b[j] /= A[j][j];
                    A[j][j] = 1;
                    if (c != 1) {
                        Print_Status_System(A,b,fn);
                    }
                }
                if (j == n - 1) {
                    i = (i + 1) % n;
                }
                j = (j + 1) % n;
                return Parallel_Ranking_Method_Rec(A,b,i,j,fn);
            }
        } else {
            throw new Exception("this is a singular matrix");
        }
    }

    ///////////////////////////////////////////// User Interface ///////////////////////////////////////////////
    // choose action in order to solve a system Ax = b
    public static float[] Solve_System(float[][] A, float[] b, int op, String fn) throws Exception {
        Print_Status_System(A,b,fn);
        switch (op) {
            case 1:
                return Invertible_Method(A,b,fn);
            case 2:
                return Cramer_Method(A,b,fn);
            case 3:
                return Forward_Backward_Method(A,b,fn);
            case 4:
                return Upper_Ranking_Method(A,b,fn);
            case 5:
                return Lower_Ranking_Method(A,b,fn);
            case 6:
                return Parallel_Ranking_Method(A,b,fn);
            case 7:
                return Parallel_Ranking_Method_Rec(A,b,0,0,fn);
            default:
                throw new Exception("you entered an invalid number");
        }
    }

    //////////////////////////////////////////// Check Structure //////////////////////////////////////////////
    // check the system of linear equations Ax = b and solve in the method requested
    public static void Check_System(float[][] A, float[] b) throws Exception {
        if (A.length == A[0].length && A.length == b.length) {
            Scanner sc = new Scanner(System.in);
            Display_Exercise(A,b);
            int n = A.length;
            if (n > 1) { // R2 Space or higher
                User_Menu_System();
                int op = sc.nextInt();
                User_Menu_Solution();
                String fn = sc.next();
                float det = Determinant(A);
                if (det == 0) {
                    if (Is_Zero_Vector(b)) {
                        System.out.println("|b| = 0 therefore system is homogeneous");
                        if (Is_Zero_Matrix(A)) {
                            System.out.print("exists an infinite number of solutions to the system in the space R" + (n));
                        } else if (Is_Linear_Dependent_Rows(A)) {
                            System.out.print("exist a single solution for the system which is: x = ");
                            Print_Solution(Zero_Vector(n),fn);
                        } else {
                            System.out.println("ranking A matrix by upper ranking:");
                            A = Upper_Ranking_Matrix(A,fn);
                            for (int i = 0; i < n - 1; i++) {
                                b[i] -= A[i][n - 1];
                            }
                            System.out.println("suppose that x" + (n) + " = " + 1 + " and we will move to right wing");
                            A = Sub_Matrix(A,n - 1,n - 1);
                            b = Sub_Vector(b,n - 1);
                            Display_Exercise(A,b);
                            float[] x = Solve_System(A,b,op,fn);
                            System.out.println("now adding the x" + (n) + " to vector solution");
                            x = Arrays.copyOf(x,n);
                            x[n - 1] = (float) 1.0;
                            System.out.print("the solution is an infinite set of linearly dependent vectors in the base ");
                            Print_Solution(x,fn);
                        }
                    } else {
                        System.out.print("does not an exist solutions because |A| is 0 and |b| is not 0");
                    }
                } else {
                    float[] x = Solve_System(A,b,op,fn);
                    System.out.print("exist a single solution for the system which is: x = ");
                    Print_Solution(x,fn);
                }
            } else { // R1 Space
                if (A[0][0] == 0) {
                    if (b[0] == 0) {
                        System.out.print("exists an infinite number of solutions to the equation in the space R1");
                    } else {
                        System.out.print("does not an exists any solution to the equation");
                    }
                } else {
                    float c = b[0] / A[0][0];
                    if (c % 1 == 0) {
                        System.out.println("exist a single solution for the equation which is: x = " + (int) c);
                    } else {
                        System.out.println("exist a single solution for the equation which is: x = " + c);
                    }
                }
            }
        } else {
            throw new Exception("your input does not meet the conditions for system of linear equations");
        }
    }

    ////////////////////////////////////////////// Run Progress ////////////////////////////////////////////////
    public static void main(String[] args) {
        float[][] A1 = {{(float) -0.5}};
        float[] a1 = {4};
        float[][] B1 = {{0}};
        float[] b1 = {4};
        float[][] A2 = {{1,-4},{-2,8}};
        float[] a2 = {5,-10};
        float[][] B2 = {{0,1},{1,1}};
        float[] b2 = {1,1};
        float[][] A3 = {{2,1,-1},{-3,-1,2},{-2,1,2}};
        float[] a3 = {8,-11,-3};
        float[][] B3 = {{1,1,2},{1,2,3},{2,2,4}};
        float[] b3 = {0,0,0};
        float[][] C3 = {{1,1,5},{1,2,7},{2,-1,4}};
        float[] c3 = {0,0,0};
        float[][] D3 = {{2,1,1},{-3,-2,-1},{4,-1,5}};
        float[] d3 = {1,-1,5};
        float[][] E3 = {{2,1,3},{10,-1,1},{4,-1,-1}};
        float[] e3 = {10,26,8};
        float[][] F3 = {{0,0,0},{0,0,0},{0,0,0}};
        float[] f3 = {0,0,0};
        float[][] G3 = {{2,1,-1},{-3,-1,2},{-2,1,2}};
        float[] g3 = {0,0,0};
        float[][] A4 = {{-2,3,3,-2},{-1,4,2,-2},{1,3,1,3},{-3,-2,4,-5}};
        float[] a4 = {8,5,19,-19};
        float[][] A5 = {{-2,2,2,-1,1},{-4,4,4,4,3},{2,3,2,3,2},{-3,-1,1,2,2},{5,5,3,5,5}};
        float[] a5 = {1,31,15,8,52};
        float[][] A6 = {{6,1,3,3,-11,1},{11,-6,5,11,0,-4},{-2,-2,-4,2,-3,0},{2,12,-1,-7,3,0},{5,-11,-11,8,-8,-2},{3,2,2,-1,-1,1}};
        float[] a6 = {-1,41,-30,-8,-27,18};
        float[][] B6 = {{-2,0,-4,-2,5,-6},{6,-1,-2,0,1,0},{6,0,0,-2,5,2},{-2,0,3,2,2,3},{-4,-2,3,-6,4,-4},{5,1,-5,-1,5,-2}};
        float[] b6 = {0,0,0,0,0,0};
        float[][] A7 = {{-1,0,4,-1,-2,5,-1},{6,-7,7,1,-7,-6,5},{0,2,6,1,-7,6,1},{2,2,0,-2,-6,6,-5},{3,3,-3,1,7,-1,-2},{1,-5,1,-4,3,1,-4},{7,-5,5,0,-4,-4,1}};
        float[] a7 = {6,-12,8,20,-3,-4,-6};
        float[][] B7 = {{2,3,1,-4,0,-3,0},{-3,1,1,1,0,-4,-1},{0,1,0,-2,1,-1,1},{-4,1,-3,1,0,-2,1},{1,-3,0,-2,-4,1,0},{1,-2,3,0,-4,-2,-4},{0,4,-4,-2,-3,-2,3}};
        float[] b7 = {0,0,0,0,0,0,0};
        try {
            Check_System(A4,a4);
        } catch (Exception ex) {
            ex.printStackTrace();
        }
    }
}
