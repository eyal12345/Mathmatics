import java.util.Scanner;

public class System_Linear_Equations_V2 {

    /////////////////////////////////////////////// Print Methods /////////////////////////////////////////////////
    // display the system Ax = b in the linear equations format
    public static void Display_Exercise(float[][] A, float[] b) {
        int n = A.length;
        if (n == 1) {
            System.out.println("Solve The Next Exercise (" + 1 + " Variable):");
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
                    System.out.println(" = " + (int) (Math.round(b[i] * 10000.0) / 10000.0));
                } else {
                    System.out.println(" = " + (Math.round(b[i] * 10000.0) / 10000.0));
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
        System.out.println("1. elementaric actions in iterative method");
        System.out.println("2. elementaric actions in recursive method");
        System.out.println("3. cramer in iterative method");
        System.out.println("4. cramer in recursive method");
        System.out.println("5. get single value from system by cramer method");
        System.out.println("6. investigation system linear equations");
    }

    // display user interface by selection format for solution
    public static void User_Menu_Solution() {
        System.out.println("choose the character of format to representation of solution:");
        System.out.println("d. decimal");
        System.out.println("r. rational");
    }

    // display a vector each current status
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
            String s = "x = ( ";
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
    // calculate multiplication between two matrices
    public static float[][] Mul_Mats(float[][] M1, float[][] M2) throws Exception {
        if (M1[0].length == M2.length) {
            float[][] A = new float[M1.length][M2[0].length];
            for (int i = 0; i < M1.length; i++) {
                for (int j = 0; j < M2[0].length; j++) {
                    for (int k = 0; k < M2.length; k++) {
                        A[i][j] += M1[i][k] * M2[k][j];
                    }
                }
            }
            return A;
        } else {
            throw new Exception("has no match between matrices dimensions");
        }
    }

    // replace between two rows in a system Ax = b
    public static void Retreat_Rows_System(float[][] A, float[] b, int r1, int r2) {
        for (int j = 0; j < A[0].length; j++) {
            float t = A[r1][j];
            A[r1][j] = A[r2][j];
            A[r2][j] = t;
        }
        float t = b[r1];
        b[r1] = b[r2];
        b[r2] = t;
    }

    // check if the specific row in the matrix is a unit vector
    public static boolean Is_Unit_Vector(float[][] A, int k) {
        int n = A.length;
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n && A[k][i] != 0; j++) {
                if (i != j && A[k][j] != 0) {
                    return false;
                }
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

    // create a unit matrix with "n*n" size
    public static float[][] Unit_Matrix(int n) {
        float[][] I = new float[n][n];
        for (int i = 0; i < n; i++) {
            I[i][i] = 1;
        }
        return I;
    }

    // get the index starting from the specific column in the matrix which are him value not equal to 0
    public static int Get_Index_UnZero_Value(float[][] A, int k) {
        int n = A.length;
        for (int i = k + 1; i < n; i++) {
            if (A[i][k] != 0) {
                return i % n;
            }
        }
        return -1;
    }

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

    // check if all vector values is equal to each other
    public static boolean Is_Equals_Values(float[] r) {
        int n = r.length;
        for (int i = 0; i < n - 1; i++) {
            if (r[i] != r[i + 1]) {
                return false;
            }
        }
        return true;
    }

    // check if exists two vectors in the matrix which are linearly dependent
    public static boolean Is_Linear_Dependent_Rows(float[][] A) {
        int n = A.length;
        float[] R = new float[n];
        for (int i1 = 0; i1 < n - 1; i1++) {
            float[] v1 = Get_Row_From_Matrix(A,i1);
            for (int i2 = i1 + 1; i2 < n; i2++) {
                float[] v2 = Get_Row_From_Matrix(A,i2);
                for (int j = 0; j < n; j++) {
                    R[j] = v1[j] / v2[j];
                }
                if (Is_Equals_Values(R)) {
                    return true;
                }
            }
        }
        return false;
    }

    // create zero vector with "n" cells
    public static float[] Zero_Vector(int n) {
        float[] x = new float[n];
        for (int k = 0; k < n; k++) {
            x[k] = 0;
        }
        return x;
    }

    // import the specific row from the matrix
    public static float[] Get_Row_From_Matrix(float[][] A, int i) {
        int n = A.length;
        float[] v = new float[n];
        for (int j = 0; j < n; j++) {
            v[j] = A[i][j];
        }
        return v;
    }

    // calculate determinant of a matrix
    public static float Determinant(float[][] A) {
        int n = A.length;
        if (n == 1) {
            return A[0][0];
        }
        float det = 0;
        for (int i = 0; i < n; i++) {
            det += Math.pow(-1,i) * A[0][i] * Determinant(Sub_Matrix(A,0,i));
        }
        return det;
    }

    // calculate of sub-matrix from a matrix by cutting row "x" and column "y"
    public static float[][] Sub_Matrix(float[][] A, int x, int y) {
        int n = A.length, p = 0, q = 0;
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

    // show elementary actions for multiplication of a row in the system
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
    // solve system of linear equations Ax = b by multiplication elementary matrix each iteration (iterative method)
    public static float[] Elementary_Method(float[][] A, float[] b, String fn) throws Exception {
        int n = A.length, i = 0, j = 0;
        while (!Is_Unit_Matrix(A)) {
            float[][] E = Unit_Matrix(n);
            if (A[i][i] == 0) {
                if (n == 1) {
                    if (b[0] == 0) {
                        throw new Exception("exists an infinite number of solutions to the equation in the space R1");
                    } else {
                        throw new Exception("does not an exists any solution to the equation");
                    }
                } else {
                    int k = Get_Index_UnZero_Value(A,i);
                    System.out.println("R" + (i + 1) + " <--> R" + (k + 1) + "\n");
                    Retreat_Rows_System(E,b,i,k);
                    A = Mul_Mats(E,A);
                    Print_Status_System(A,b,fn);
                }
            } else {
                if (i != j) {
                    E[j][i] -= (A[j][i] / A[i][i]);
                    Sum_Elementary_Action(-E[j][i],j,i);
                    b[j] += b[i] * E[j][i];
                    A = Mul_Mats(E,A);
                    A[j][i] = 0;
                    if (E[j][i] != 0) {
                        Print_Status_System(A,b,fn);
                    }
                } else if (Is_Unit_Vector(A,j)) {
                    E[j][j] = 1 / A[j][j];
                    Mul_Elementary_Action(E[j][j],j);
                    b[j] *= E[j][j];
                    A = Mul_Mats(E,A);
                    A[j][j] = 1;
                    if (E[j][j] != 1) {
                        Print_Status_System(A,b,fn);
                    }
                }
                if (j == n - 1) {
                    i = (i + 1) % n;
                }
                j = (j + 1) % n;
            }
        }
        return b;
    }

    // solve system of linear equations Ax = b by multiplication elementary matrix each iteration (recursive method)
    public static float[] Elementary_Method_Rec(float[][] A, float[] b, int i, int j, String fn) throws Exception {
        int n = A.length;
        if (Is_Unit_Matrix(A)) {
            return b;
        } else {
            float[][] E = Unit_Matrix(n);
            if (A[i][i] == 0) {
                if (n == 1) {
                    if (b[0] == 0) {
                        throw new Exception("exists an infinite number of solutions to the equation in the space R1");
                    } else {
                        throw new Exception("does not an exists any solution to the equation");
                    }
                } else {
                    int k = Get_Index_UnZero_Value(A,i);
                    System.out.println("R" + (i + 1) + " <--> R" + (k + 1) + "\n");
                    Retreat_Rows_System(E,b,i,k);
                    A = Mul_Mats(E,A);
                    Print_Status_System(A,b,fn);
                }
            } else {
                if (i != j) {
                    E[j][i] -= (A[j][i] / A[i][i]);
                    Sum_Elementary_Action(-E[j][i],j,i);
                    b[j] += b[i] * E[j][i];
                    A = Mul_Mats(E,A);
                    A[j][i] = 0;
                    if (E[j][i] != 0) {
                        Print_Status_System(A,b,fn);
                    }
                } else if (Is_Unit_Vector(A,j)) {
                    E[j][j] = 1 / A[j][j];
                    Mul_Elementary_Action(E[j][j],j);
                    b[j] *= E[j][j];
                    A = Mul_Mats(E,A);
                    A[j][j] = 1;
                    if (E[j][j] != 1) {
                        Print_Status_System(A,b,fn);
                    }
                }
                if (j == n - 1) {
                    i = (i + 1) % n;
                }
                j = (j + 1) % n;
            }
            return Elementary_Method_Rec(A,b,i,j,fn);
        }
    }

    // solve system of linear equations Ax = b by cramer method (iterative method)
    public static float[] Cramer_Method(float[][] A, float[] b) throws Exception {
        int n = A.length;
        float det = Determinant(A);
        float[] x = Zero_Vector(n);
        if (det != 0) {
            for (int i = 0; i < n; i++) {
                float sum = 0;
                for (int j = 0; j < n; j++) {
                    sum += Math.pow(-1,i + j) * b[j] * Determinant(Sub_Matrix(A,j,i));
                }
                x[i] = sum / det;
            }
            return x;
        } else {
            if (Is_Zero_Vector(b)) {
                throw new Exception("infinity solutions");
            }
            throw new Exception("has no solution because det(A) = 0");
        }
    }

    // solve system of linear equations Ax = b by cramer method (recursive method)
    public static float[] Cramer_Method_Rec(float[][] A, float[] b, float[] x, int i) throws Exception {
        int n = A.length;
        if (i == n) {
            return x;
        } else {
            float det = Determinant(A);
            if (det != 0 && i < n) {
                float sum = 0;
                for (int j = 0; j < n; j++) {
                    sum += Math.pow(-1,i + j) * b[j] * Determinant(Sub_Matrix(A,j,i));
                }
                x[i] = sum / det;
                return Cramer_Method_Rec(A,b,x,i + 1);
            } else {
                if (Is_Zero_Vector(b)) {
                    throw new Exception("infinity solutions");
                }
                throw new Exception("has no solution because det(A) = 0");
            }
        }
    }

    // get one specific value from vector of solution a system Ax = b
    public static void Get_Single_Value(float[][] A, float[] b, int i, String fn) throws Exception {
        int n = A.length;
        float det = Determinant(A), val;
        if (det != 0) {
            if (n == 1) {
                val = b[0] / det;
                if (val % 1 == 0) {
                    System.out.println("x = " + (int) val);
                } else if (fn.equals("d")) {
                    System.out.println("x = " + val);
                } else if (fn.equals("r")) {
                    System.out.println("x = " + convertDecimalToFraction(val));
                }
            } else {
                float sum = 0;
                for (int j = 1; j <= n; j++) {
                    sum += Math.pow(-1,i + j) * b[j - 1] * Determinant(Sub_Matrix(A,j - 1,i - 1));
                }
                val = sum / det;
                if (val % 1 == 0) {
                    System.out.println("x" + i + " = " + (int) val);
                } else if (fn.equals("d")) {
                    System.out.println("x" + i + " = " + val);
                } else if (fn.equals("r")) {
                    System.out.println("x" + i + " = " + convertDecimalToFraction(val));
                }
            }
        } else {
            throw new Exception("has no solution because det(A) = 0");
        }
    }

    // investigate system of linear equations Ax = b without solve it
    public static void System_Investigation(float[][] A, float[] b, String fn) throws Exception {
        if (A.length == A[0].length && A.length == b.length) {
            int n = b.length;
            if (n > 1) { // R2 Space or higher
                System.out.println("system of linear equations from the shape Ax = b:");
                if (Determinant(A) == 0) {
                    if (Is_Zero_Vector(b)) {
                        if (Is_Zero_Matrix(A)) {
                            System.out.print("exists an infinite number of solutions to the system in a space R" + n);
                        } else if (Is_Linear_Dependent_Rows(A)) {
                            System.out.print("exist a single solution for the system which is: x = ");
                            Print_Solution(Zero_Vector(n),fn);
                        } else {
                            System.out.print("the solution is an infinite set of linearly dependent vectors in the base x");
                        }
                    } else {
                        System.out.print("does not an exists solutions because det(A) = 0");
                    }
                } else {
                    System.out.print("vector x is a single solution for system Ax = b");
                }
            } else { // R1 Space
                System.out.println("equation from the shape ax = b:");
                if (A[0][0] == 0) {
                    if (b[0] == 0) {
                        System.out.print("exists an infinite number of solutions to the equation in the space R1");
                    } else {
                        System.out.print("does not an exists any solution to the equation");
                    }
                } else {
                    System.out.print("exist a single solution for the equation which is: x = b/a");
                }
            }
        } else {
            throw new Exception("your input does not meet the conditions for system of linear equations");
        }
    }

    ///////////////////////////////////////////// User Interface ///////////////////////////////////////////////
    // choose action in order to solve a system Ax = b
    public static void Solve_System(float[][] A, float[] b, int op, String fn) throws Exception {
        Display_Exercise(A,b);
        Print_Status_System(A,b,fn);
        float[] x;
        if (op == 1) {
            x = Elementary_Method(A,b,fn);
            Print_Solution(x,fn);
        } else if (op == 2) {
            x = Elementary_Method_Rec(A,b,0,0,fn);
            Print_Solution(x,fn);
        } else if (op == 3) {
            x = Cramer_Method(A,b);
            Print_Solution(x,fn);
        } else if (op == 4) {
            x = Cramer_Method_Rec(A,b,Zero_Vector(A.length),0);
            Print_Solution(x,fn);
        } else if (op == 5) {
            int n = A.length;
            if (n == 1) { // R1 Space
                System.out.println("the single value is: ");
                Get_Single_Value(A,b,0,fn);
            } else { // R2 Space or higher
                System.out.println("enter number between 1 to " + n + ":");
                Scanner sc = new Scanner(System.in);
                int i = sc.nextInt();
                if (i >= 1 && i <= n) {
                    System.out.println("the single value in location " + i + " is: ");
                    Get_Single_Value(A,b,i,fn);
                } else {
                    throw new Exception("you entered an invalid index");
                }
            }
        } else if (op == 6) {
            System_Investigation(A,b,fn);
        } else {
            throw new Exception("you entered an invalid number");
        }
    }

    ////////////////////////////////////////////// Run Progress ////////////////////////////////////////////////
    public static void main(String[] args) {
        //////////////////////////////////////////////////////////////////////////////////////////////
        float[][] A1 = {{-1,0,4,-1,-2,5,-1},{6,-7,7,1,-7,-6,5},{0,2,6,1,-7,6,1},{2,2,0,-2,-6,6,-5},{3,3,-3,1,7,-1,-2},{1,-5,1,-4,3,1,-4},{7,-5,5,0,-4,-4,1}};
        float[] b1 = {6,-12,8,20,-3,-4,-6};
        // x = (2, 17, 5, -29, -1, -7, 8)
        float[][] A2 = {{7,4,4,-1,-3,-1,-3},{-5,-3,4,2,-4,-1,7},{-5,7,-7,-1,-4,0,0},{6,3,5,3,-5,0,0},{3,1,5,-1,7,-5,-1},{1,0,7,-5,-2,-2,-7},{-1,0,3,3,7,-2,-5}};
        float[] b2 = {0,1,5,6,-18,7,14};
        // x = (-22, 57, 53, -31, 41, 113, 27)
        float[][] A3 = {{6,-5,4,6,6,5,1},{4,6,2,1,-5,-6,4},{0,-2,4,-7,1,7,-7},{4,-5,-3,3,0,4,-4},{-7,-2,-5,-7,4,-4,-1},{7,-1,3,-4,-4,3,-6},{-7,-6,-3,-3,-4,6,-7}};
        float[] b3 = {-10,13,-6,-20,-9,-12,15};
        // x = (316.9, -186, -509.8, -695.4, -300.2, 792.7, 1207.9)
        float[][] A4 = {{-6,0,-3,-5,1,-3,6},{-2,1,1,-2,-7,6,-2},{-7,-6,-2,-5,-4,5,-1},{-6,-3,4,-4,-5,7,-1},{-7,0,-1,-2,4,-6,-7},{-6,-4,4,-5,-2,5,3},{3,5,-1,5,2,-5,-5}};
        float[] b4 = {-18,-7,13,16,-20,-17,-4};
        // x = (1588, -4485.5, 3378.5, -183.5, -7366, -7791, 453.5)
        float[][] A5 = {{7,6,4,2,-2,-7,5},{0,-3,7,-7,4,0,-7},{5,-3,4,-5,0,4,0},{4,3,-5,0,4,-2,-3},{6,-3,6,-3,-4,2,5},{-7,1,-7,6,-2,3,2},{6,-2,-5,-1,-3,4,-2}};
        float[] b5 = {-8,-8,10,14,1,8,2};
        // x = (13.6, -16.4, 13.8, 40, 18.6, 9.4, -7.4)
        float[][] A6 = {{-3,6,-3,3,-4,6,6},{4,0,2,-4,-6,-6,-3},{-2,-2,-7,5,3,3,2},{-5,0,2,0,3,2,4},{-6,0,0,7,-3,4,6},{5,-5,2,1,7,-2,-5},{-3,-2,6,5,0,-1,7}};
        float[] b6 = {14,-6,20,6,16,-10,-15};
        // x = (12.5, -31.5, -1, -19.5, -12.5, 28, 13)
        float[][] A7 = {{6,-3,7,-6,-5,-2,3},{5,4,7,6,6,2,-5},{-4,-2,-1,-5,-2,3,-5},{-7,1,2,-1,-5,0,6},{7,-2,0,7,-2,4,-2},{-1,-3,-4,-4,0,0,2},{7,-6,-6,-4,6,-1,2}};
        float[] b7 = {-2,-9,15,0,18,-13,-15};
        // x = (-12.4, -30, 9.7, 14.6, 1.4, -18.2, -9.1)
        float[][] A8 = {{0,-3,-1,5,7,0,0},{6,-5,1,6,0,5,-7},{-1,-5,7,-5,-7,4,0},{3,-5,4,7,6,0,3},{-4,-1,4,-6,-4,-6,-3},{2,6,0,4,-1,-1,-1},{-7,0,7,-5,0,-2,-3}};
        float[] b8 = {-5,7,18,-9,-12,11,-6};
        // x = (-6.7, -0.85, -1.65, 6.7, -6.1, 2.35, 1.05)
        float[][] A9 = {{0,0,3,4,-1,6,-3},{7,-4,0,-7,1,-7,-4},{-7,3,-3,5,1,3,4},{-6,5,3,0,6,5,-5},{2,-5,1,-6,-7,7,-6},{6,7,-3,2,-6,3,3},{2,2,-7,-3,3,-6,0}};
        float[] b9 = {3,5,2,-13,20,2,5};
        // x = (-20.65, 11.8, 2.2, 13.55, -33.25, -29.75, -29.15)
        float[][] A10 = {{0,-5,3,3,-4,7,-6},{-2,1,3,-7,6,7,-6},{0,1,0,-2,5,5,2},{-6,6,6,-3,-4,1,-5},{-3,1,1,-3,-5,6,0},{4,2,-4,0,-1,0,0},{1,0,-1,-1,4,-4,5}};
        float[] b10 = {-15,-2,17,1,-10,15,5};
        // x = (3.11, 4.87, 1.4, 2.58, 1.58, 1.45, 1.07)
        float[][] A11 = {{4,-5,1,1,-1,-6,7},{-5,0,1,-7,-2,-5,-5},{0,-1,5,0,7,6,0},{6,-3,1,0,-3,-3,6},{-5,-5,4,-4,4,-4,-5},{0,1,-2,6,4,-7,6},{-1,-6,2,-2,-4,1,7}};
        float[] b11 = {-13,19,-8,-19,16,11,-20};
        // x = (-2.1, 1.68, 1.46, 1.58, -0.3, -1.92, -1.58)
        float[][] A12 = {{-3,-6,0,-1,-5,-2,-4},{3,-2,-2,3,5,5,0},{5,7,4,4,7,-1,-2},{3,0,0,-7,3,7,-6},{-6,0,0,-2,-7,-4,-6},{0,1,7,-3,4,5,3},{0,7,2,0,2,-1,-3}};
        float[] b12 = {2,-3,-10,5,-7,7,-11};
        // x = (3.19, -0.21, 1.33, -0.21, -3.58, 1.64, 1.13)
        float[][] A13 = {{-6,-6,-4,4,7,-6,4},{-2,1,0,7,-7,-2,3},{-2,6,0,-7,7,0,-1},{5,7,7,0,-3,7,-4},{7,-2,4,-3,7,6,-6},{-1,-3,-4,-7,3,-2,0},{0,7,1,-6,0,2,-1}};
        float[] b13 = {10,-4,15,-13,16,5,-20};
        // x = (67.26, -2.56, 78.06, -26.1, 6.2, -87.76, 61.22)
        float[][] A14 = {{4,0,0,1,-3,4,-7},{5,4,4,4,0,7,0},{4,-2,0,-3,0,-6,3},{-3,-2,0,4,-3,3,4},{-6,4,5,-2,7,-6,-2},{0,4,-7,-6,-3,-4,-2},{3,6,1,0,0,3,0}};
        float[] b14 = {12,8,17,-17,-15,16,-12};
        // x = (15.12, 10.07, -32.94, 55.46, 17.66, -28.28, -8.88)
        float[][] A15 = {{0,-4,-2,7,7,1,0},{0,2,3,1,2,3,0},{0,-1,-1,-5,2,-3,-1},{2,1,3,0,-5,1,-6},{2,4,0,4,-3,2,4},{4,2,6,3,0,-4,-2},{6,2,5,-3,-3,-4,-2}};
        float[] b15 = {3,8,9,11,14,18,18};
        // x = (4.3, 9, -5.2, 0.4, 3.8, -0.8, -4.8)
        //////////////////////////////////////////////////////////////////////////////////////////////
        float[][] X2 = {{1,-1,-2},{2,-3,-5},{-1,3,5}};
        float[] x2 = {1,0,-3};
        float[][] X3 = {{1,1,5},{1,2,7},{2,-1,4}};
        float[] x3 = {0,0,0};
        float[][] X4 = {{2,-1},{3,2}};
        float[] x4 = {5,4};
        float[][] X5 = {{1,1,2},{1,2,3},{2,2,4}};
        float[] x5 = {0,0,0};
        float[][] S1 = {{(float) -0.5}};
        float[] s1 = {4};
        try {
            Scanner sc = new Scanner(System.in);
            User_Menu_System();
            int op = sc.nextInt();
            User_Menu_Solution();
            String fn = sc.next();
            Solve_System(A15,b15,op,fn);
        } catch (Exception e) {
            e.printStackTrace();
        }
    }
}
