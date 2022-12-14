import java.text.DecimalFormat;
import java.util.Scanner;
import java.util.Vector;

public class System_Linear_Equations_V3 {

    // display system Ax = b in linear equations format
    private static void Display_Exercise(float[][] A ,float[] b) {
        int m = A.length ,n = A[0].length;
        if (n == 1) {
            System.out.println("Solve The Next Equation Under Space R" + 1 + ":");
            if (A[0][0] % 1 == 0) {
                if (A[0][0] == 1) {
                    if (b[0] % 1 == 0)
                        System.out.println("x = " + (int) b[0]);
                    else
                        System.out.println("x = " + b[0]);
                } else if (A[0][0] == -1) {
                    if (b[0] % 1 == 0)
                        System.out.println("-x = " + (int) b[0]);
                    else
                        System.out.println("-x = " + b[0]);
                } else {
                    if (b[0] % 1 == 0)
                        System.out.println((int)A[0][0] + "x = " + (int) b[0]);
                    else
                        System.out.println((int)A[0][0] + "x = " + b[0]);
                }
            } else {
                if (b[0] % 1 == 0)
                    System.out.println(A[0][0] + "x = " + (int)b[0]);
                else
                    System.out.println(A[0][0] + "x = " + b[0]);
            }
        } else {
            if (m == 1)
                System.out.println("Solve The Next Equation Under Space R" + n + ":");
            else
                System.out.println("Solve The Next System Under Space R" + n + ":");
            for (int i = 0 ;i < m ;i++) {
                boolean flag = true;
                for (int j = 0 ;j < n ;j++) {
                    if (A[i][j] > 0) {
                        if (flag) {
                            System.out.print("  ");
                            flag = false;
                        } else
                            System.out.print(" + ");
                        if (Math.abs(A[i][j]) == 1)
                            System.out.print("x" + (j + 1));
                        else
                            System.out.print((int) A[i][j] + "*x" + (j + 1));
                    } else if (A[i][j] < 0) {
                        if (flag) {
                            System.out.print("- ");
                            flag = false;
                        } else
                            System.out.print(" - ");
                        if (Math.abs(A[i][j]) == 1)
                            System.out.print("x" + (j + 1));
                        else
                            System.out.print((int) Math.abs(A[i][j]) + "*x" + (j + 1));
                    }
                }
                System.out.println(" = " + (int) b[i]);
            }
        }
        System.out.println();
    }

    private static void Print_State(float[][] A ,float[][] b ,String ch) throws Exception {
        if (!ch.equals("d") && !ch.equals("r"))
            throw new Exception("enter error character for representation solution");
        else {
            int m = A.length ,n = A[0].length;
            for (int i = 0 ;i < m ;i++) {
                for (int j = 0 ;j < n ;j++) {
                    if ((Math.round(A[i][j] * 1000.0) / 1000.0) % 1 == 0)
                        System.out.print((int)(Math.round(A[i][j] * 1000.0) / 1000.0) + " ");
                    else if (ch.equals("d"))
                        System.out.print(Math.round(A[i][j] * 1000.0) / 1000.0 + " ");
                    else if (ch.equals("r"))
                        System.out.print(convertDecimalToFraction(A[i][j]) + " ");
                    if (j != n - 1)
                        System.out.print(",");
                }
                System.out.print("| ");
                for (int j = 0 ;j < n + 1 ;j++) {
                    if ((Math.round(b[i][j] * 1000.0) / 1000.0) % 1 == 0)
                        System.out.print((int)(Math.round(b[i][j] * 1000.0) / 1000.0) + " ");
                    else if (ch.equals("d"))
                        System.out.print(Math.round(b[i][j] * 1000.0) / 1000.0 + " ");
                    else if (ch.equals("r"))
                        System.out.print(convertDecimalToFraction(b[i][j]) + " ");
                    if (Is_Zero_Col(b,j + 1))
                        break;
                    if (j != n)
                        System.out.print(",");
                }
                System.out.println();
            }
            System.out.println();
        }
    }

    // display select interface for users
    private static void User_Menu_System() {
        System.out.println("choose number method to solution:");
        System.out.println("1. upper -> lower ranking method");
        System.out.println("2. lower -> upper ranking method");
        System.out.println("3. parallel ranking method (iterative method)");
        System.out.println("4. parallel ranking method (recursive method)");
    }

    private static void User_Menu_Solution() {
        System.out.println("choose representation format for solution:");
        System.out.println("d. decimal");
        System.out.println("r. rational");
    }

    private static void Sum_Elementaric_Action(float k ,int j ,int i ,String ch) throws Exception {
        if (!ch.equals("d") && !ch.equals("r"))
            throw new Exception("enter error character for representation solution");
        else {
            if (k != 0) {
                int r = j + 1 ,c = i + 1;
                if (k > 0) {
                    if (k % 1 == 0)
                        if (k == 1)
                            System.out.println("R" + r + " --> R" + r + " - R" + c + "\n");
                        else
                            System.out.println("R" + r + " --> R" + r + " - " + (int)k + "*R" + c + "\n");
                    else if (ch.equals("d")) {
                        k = (float) (Math.round(k * 10000.0) / 10000.0);
                        if (k % 1 == 0)
                            System.out.println("R" + r + " --> R" + r + " - " + (int)k + "*R" + c + "\n");
                        else
                            System.out.println("R" + r + " --> R" + r + " - " + k + "*R" + c + "\n");
                    } else if (ch.equals("r")) {
                        String v = convertDecimalToFraction(k);
                        if (!v.equals("1"))
                            System.out.println("R" + r + " --> R" + r + " - " + convertDecimalToFraction(k) + "*R" + c + "\n");
                        else
                            System.out.println("R" + r + " --> R" + r + " - R" + c + "\n");
                    }
                } else {
                    if (k % 1 == 0)
                        if (k == -1)
                            System.out.println("R" + r + " --> R" + r + " + R" + c + "\n");
                        else
                            System.out.println("R" + r + " --> R" + r + " + " + (int)(-k) + "*R" + c + "\n");
                    else if (ch.equals("d")) {
                        k = (float) (Math.round(k * 10000.0) / 10000.0);
                        if (k % 1 == 0)
                            System.out.println("R" + r + " --> R" + r + " + " + (int)(-k) + "*R" + c + "\n");
                        else
                            System.out.println("R" + r + " --> R" + r + " + " + (-k) + "*R" + c + "\n");
                    } else if (ch.equals("r")) {
                        String v = convertDecimalToFraction(-k);
                        if (!v.equals("-1"))
                            System.out.println("R" + r + " --> R" + r + " + " + convertDecimalToFraction(-k) + "*R" + c + "\n");
                        else
                            System.out.println("R" + r + " --> R" + r + " + R" + c + "\n");
                    }
                }
            }
        }
    }

    private static void Mult_Elementaric_Action(float k ,int j ,String ch) throws Exception {
        if (!ch.equals("d") && !ch.equals("r"))
            throw new Exception("enter error character for representation solution");
        else {
            if (k != 1) {
                int r = j + 1;
                if (k % 1 == 0)
                    if (k == -1)
                        System.out.println("R" + r + " --> - R" + r + "\n");
                    else
                        System.out.println("R" + r + " --> " + (int)k + "*R" + r + "\n");
                else if (ch.equals("d")) {
                    k = (float) (Math.round(k * 10000.0) / 10000.0);
                    if (k % 1 == 0)
                        System.out.println("R" + r + " --> " + (int)k + "*R" + r + "\n");
                    else
                        System.out.println("R" + r + " --> " + k + "*R" + r + "\n");
                } else if (ch.equals("r"))
                    System.out.println("R" + r + " --> " + convertDecimalToFraction(k) + "*R" + r + "\n");
            }
        }
    }

    private static void Retreat_Elementaric_Action(int i , int j) {
        int r1 = i + 1 ,r2 = j + 1;
        if (r1 <= r2)
            System.out.println("R" + r1 + " <--> R" + r2 + "\n");
        else
            System.out.println("R" + r2 + " <--> R" + r1 + "\n");
    }

    private static void Print_Solution(float[][] x ,String ch) throws Exception {
        if (!ch.equals("d") && !ch.equals("r"))
            throw new Exception("enter error character for representation solution");
        else {
            int m = x.length;
            String s = "";
            if (Is_Zero_Col(x,0) && !Is_Zero_Col(x,1)) { // ??u
                int n = x[0].length;
                System.out.print("solution is infinite set of vectors that linear dependent in line: x = ");
                for (int t = 1 ;t < n && !Is_Zero_Col(x,t) ;t++) {
                    s = (Is_Zero_Col(x,2)) ? s + "??*( " : s + "??" + t + "*( ";
                    for (int i = 0 ;i < m ;i++) {
                        if (i == m - 1) {
                            if ((Math.round(x[i][t] * 1000.0) / 1000.0) % 1 == 0)
                                s += (int) (Math.round(x[i][t] * 1000.0) / 1000.0) + " )";
                            else if (ch.equals("d"))
                                s += Math.round(x[i][t] * 1000.0) / 1000.0 + " )";
                            else if (ch.equals("r"))
                                s += convertDecimalToFraction(x[i][t]) + " )";
                        } else {
                            if ((Math.round(x[i][t] * 1000.0) / 1000.0) % 1 == 0)
                                s += (int) (Math.round(x[i][t] * 1000.0) / 1000.0) + " , ";
                            else if (ch.equals("d"))
                                s += Math.round(x[i][t] * 1000.0) / 1000.0 + " , ";
                            else if (ch.equals("r"))
                                s += convertDecimalToFraction(x[i][t]) + " , ";
                        }
                    }
                    if (m <= n && !Is_Zero_Col(x,t + 1))
                        s += " + ";
                }
            } else if (!Is_Zero_Col(x,0) && !Is_Zero_Col(x,1)) { // x0 + ??u
                System.out.print("solution is infinite set of vectors that linear dependent in line: x = ");
                s += "( ";
                for (int i = 0 ;i < m ;i++)
                    if (i == m - 1) {
                        if ((Math.round(x[i][0] * 1000.0) / 1000.0) % 1 == 0)
                            s += (int) (Math.round(x[i][0] * 1000.0) / 1000.0) + " )";
                        else if (ch.equals("d"))
                            s += Math.round(x[i][0] * 1000.0) / 1000.0 + " )";
                        else if (ch.equals("r"))
                            s += convertDecimalToFraction(x[i][0]) + " )";
                    } else {
                        if ((Math.round(x[i][0] * 1000.0) / 1000.0) % 1 == 0)
                            s += (int) (Math.round(x[i][0] * 1000.0) / 1000.0) + " , ";
                        else if (ch.equals("d"))
                            s += Math.round(x[i][0] * 1000.0) / 1000.0 + " , ";
                        else if (ch.equals("r"))
                            s += convertDecimalToFraction(x[i][0]) + " , ";
                    }
                for (int t = 1 ;t < m ;t++) {
                    if (!Is_Zero_Col(x,t)) {
                        s = (Is_Zero_Col(x,2)) ? s + " + ??*( " : s + " + ??" + t + "*( ";
                        for (int i = 0 ;i < m ;i++)
                            if (i == m - 1) {
                                if ((Math.round(x[i][t] * 1000.0) / 1000.0) % 1 == 0)
                                    s += (int) (Math.round(x[i][t] * 1000.0) / 1000.0) + " )";
                                else if (ch.equals("d"))
                                    s += Math.round(x[i][t] * 1000.0) / 1000.0 + " )";
                                else if (ch.equals("r"))
                                    s += convertDecimalToFraction(x[i][t]) + " )";
                            } else {
                                if ((Math.round(x[i][t] * 1000.0) / 1000.0) % 1 == 0)
                                    s += (int) (Math.round(x[i][t] * 1000.0) / 1000.0) + " , ";
                                else if (ch.equals("d"))
                                    s += Math.round(x[i][t] * 1000.0) / 1000.0 + " , ";
                                else if (ch.equals("r"))
                                    s += convertDecimalToFraction(x[i][t]) + " , ";
                            }
                    }
                }
            } else { // x0
                System.out.print("exist single solution for the system that is: x = ");
                s += "( ";
                for (int i = 0 ;i < m ;i++) {
                    x[i][0] = roundTwoDecimals(x[i][0]);
                    if (i == m - 1) {
                        if ((Math.round(x[i][0] * 1000.0) / 1000.0) % 1 == 0)
                            s += (int) (Math.round(x[i][0] * 1000.0) / 1000.0) + " )";
                        else if (ch.equals("d"))
                            s += Math.round(x[i][0] * 1000.0) / 1000.0 + " )";
                        else if (ch.equals("r"))
                            s += convertDecimalToFraction(x[i][0]) + " )";
                    } else {
                        if ((Math.round(x[i][0] * 1000.0) / 1000.0) % 1 == 0)
                            s += (int) (Math.round(x[i][0] * 1000.0) / 1000.0) + " , ";
                        else if (ch.equals("d"))
                            s += Math.round(x[i][0] * 1000.0) / 1000.0 + " , ";
                        else if (ch.equals("r"))
                            s += convertDecimalToFraction(x[i][0]) + " , ";
                    }
                }
            }
            System.out.println(s);
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

    private static float[][] Increase_Rows_in_Matrix(float[][] A ,int m) {
        int n = A[0].length;
        float[][] nA = new float[n][n];
        for (int i = 0 ;i < m ;i++)
            for (int j = 0; j < n; j++)
                nA[i][j] = A[i][j];
        return nA;
    }

    private static float[][] Increase_Rows_in_Vector(float[] b ,int m) {
        int n = b.length;
        float[][] nb = new float[m][1];
        for (int i = 0 ;i < n ;i++)
            nb[i][0] = (i < m) ? b[i] : 0;
        return nb;
    }

    private static float[][] Increase_Cols_in_Vector(float[][] b) {
        int m = b.length ,n = b[0].length;
        float[][] nb = new float[m][n + 1];
        for (int i = 0 ;i < m ;i++)
            for (int j = 0 ;j < n ;j++)
                nb[i][j] = b[i][j];
        return nb;
    }

    private static float[][] Decrease_Cols_in_Vector(float[][] b ,int k) {
        int m = b.length ,n = b[0].length;
        float[][] nb = new float[m][n - 1];
        for (int i = 0 ;i < m ;i++)
            for (int j = 0 ;j < n - 1 ;j++)
                nb[i][j] = b[i][j] - b[i][n - 1] * (b[k][j] / b[k][n - 1]);
        return nb;
    }

    // check if b is zero vector
    private static boolean Is_Zero_Vector(float[] b) {
        int n = b.length;
        for (int i = 0 ;i < n ;i++)
            if (b[i] != 0)
                return false;
        return true;
    }

    // check if exist zero row in A Matrix
    private static boolean Is_Zero_Row(float[][] A ,int r) {
        int n = A[0].length;
        for (int j = 0 ;j < n ;j++)
            if (A[r][j] != 0)
                return false;
        return true;
    }

    // check if exist zero col in b Vector
    private static boolean Is_Zero_Col(float[][] b ,int c) {
        int m = b.length ,n = b[0].length;
        if (c >= n)
            return true;
        for (int i = 0 ;i < m ;i++)
            if (b[i][c] != 0)
                return false;
        return true;
    }

    // check if A matrix is upper-triangular
    private static boolean Is_Upper_Triangular(float[][] A) {
        int n = A.length;
        for (int i = 0 ;i < n ;i++)
            for (int j = 0; j < i; j++)
                if (A[i][j] != 0)
                    return false;
        return true;
    }

    // check if A matrix is lower-triangular
    private static boolean Is_Lower_Triangular(float[][] A) {
        int n = A.length;
        for (int i = 0 ;i < n - 1 ;i++)
            for (int j = i + 1; j < n; j++)
                if (A[i][j] != 0)
                    return false;
        return true;
    }

    private static int Get_Linear_Dependent_Columns(float[][] A) {
        int m = A.length ,n = A[0].length;
        for (int c1 = 0 ;c1 < n - 1 ;c1++) {
            for (int c2 = c1 + 1 ;c2 < n ;c2++) {
                Vector<Float> C  = new Vector<Float>();
                for (int r = 0 ;r < m ;r++)
                    if (A[r][c1] != 0 || A[r][c2] != 0)
                        C.add(A[r][c1] / A[r][c2]);
                if (Is_Equals_Values(C) && C.size() > 1)
                    return c1;
            }
        }
        return -1;
    }

    // check if two vectors in the A matrix is linear dependents
    private static boolean Is_Linear_Dependent_Rows(float[][] A) {
        int m = A.length ,n = A[0].length;
        for (int r1 = 0 ;r1 < m - 1 ;r1++) {
            for (int r2 = r1 + 1 ;r2 < m ;r2++) {
                if (!Is_Unit_Vector(A,r1) && (!Is_Unit_Vector(A,r2))) {
                    Vector<Float> R  = new Vector<Float>();
                    for (int j = 0 ;j < n ;j++) {
                        if (A[r1][j] != 0 || A[r2][j] != 0) {
                            R.add(A[r1][j] / A[r2][j]);
                            if (r2 == r1 + 1 && Is_Equals_Values(R) && R.size() >= 2)
                                return true;
                        }
                    }
                }
            }
        }
        return false;
    }

    // check if two vectors in the A matrix is linear independents
    private static boolean Is_Linear_Independent_System(float[][] A ,float[] b) {
        int m = A.length ,n = A[0].length;
        for (int r = 0 ;r < m ;r++)
            if (Is_Zero_Row(A,r) && b[r] != 0)
                return true;
        for (int r1 = 0 ;r1 < m - 1 ;r1++) {
            for (int r2 = r1 + 1 ;r2 < m ;r2++) {
                Vector<Float> R  = new Vector<Float>();
                for (int j = 0 ;j < n ;j++)
                    if (A[r1][j] != 0 || A[r2][j] != 0)
                        R.add(A[r1][j] / A[r2][j]);
                if (Is_Equals_Values(R) && (b[r1] != 0 || b[r2] != 0) && (b[r1] / b[r2] != R.get(0)))
                    return true;
            }
        }
        return false;
    }

    // check if vector values is equal to each other
    private static boolean Is_Equals_Values(Vector<Float> R) {
        int n = R.size();
        for (int i = 0 ;i < n - 1 ;i++)
            if (!R.get(i).equals(R.get(i + 1)))
                return false;
        return true;
    }

    // get row from A
    private static float[] Get_Row_from_Matrix(float[][] A , int r) {
        int n = A[0].length;
        float[] v = new float[n];
        for (int j = 0 ;j < n ;j++)
            v[j] = A[r][j];
        return v;
    }

    private static int Get_Index_Row_from_Matrix(float[][] A ,int r) {
        int n = A[0].length;
        float[] v = Get_Row_from_Matrix(A,r);
        v[r] = (Is_Zero_Vector(v)) ? 1 : v[r];
        int c1 = Get_Index_for_Unit_Vector(v);
        for (int i = 0 ;i < n ;i++) {
            int c2 = Get_Index_for_Unit_Vector(Get_Row_from_Matrix(A,i));
            if (c1 == c2 && c1 != -1)
                return i;
        }
        return -1;
    }

    private static boolean Is_Exist_Vector(float[][] A ,int r) {
        int n = A[0].length;
        float[] v = Get_Row_from_Matrix(A,r);
        v[r] = (Is_Zero_Vector(v)) ? 1 : v[r];
        int c1 = Get_Index_for_Unit_Vector(v);
        for (int i = 0 ;i < n ;i++) {
            if (i != r && Is_Unit_Vector(A,r) && Is_Unit_Vector(A,i)) {
                int c2 = Get_Index_for_Unit_Vector(Get_Row_from_Matrix(A,i));
                if (c1 == c2 && c1 != -1)
                    return true;
            }
        }
        return false;
    }

    // check if A is zero matrix
    private static boolean Is_Unit_Matrix(float[][] A) {
        int n = A.length;
        for (int i = 0 ;i < n ;i++)
            for (int j = 0 ;j < n ;j++)
                if (A[i][i] != 1 || (i != j && A[i][j] != 0))
                    return false;
        return true;
    }

    // check if in the A matrix exist row that unit vector
    private static boolean Is_Unit_Vector(float[][] A , int i) {
        int n = A[0].length;
        int c = 0;
        for (int j = 0 ;j < n ;j++)
            if (A[i][j] != 0)
                c++;
        if (c <= 1)
            return true;
        return false;
    }

    // get index column for unit vector in A matrix
    private static int Get_Index_for_Unit_Vector(float[] v) {
        int n = v.length;
        for (int c = 0 ;c < n ;c++)
            if (v[c] != 0)
                return c;
        return -1;
    }

    private static void Retreat_Rows_System(float[][] A ,float[][] b ,int r1 ,int r2) {
        int n = A[0].length ,m = b[0].length;
        for (int j = 0 ;j < n ;j++) {
            float k = A[r1][j];
            A[r1][j] = A[r2][j];
            A[r2][j] = k;
        }
        for (int j = 0 ;j < m ;j++) {
            if (!Is_Zero_Col(b,j)) {
                float k = b[r1][j];
                b[r1][j] = b[r2][j];
                b[r2][j] = k;
            }
        }
    }

    // get the index with number row is not zero value
    private static int Get_Index_UnZero_Value(float[][] A ,int k ,boolean flag) {
        int n = A.length;
        if (flag) {
            for (int i = k + 1 ;i < n + k ;i++)
                if (A[i % n][k] != 0)
                    return i % n;
        } else {
            for (int i = n + k - 1 ;i > k - 1 ;i--)
                if (A[i % n][k] != 0)
                    return i % n;
        }
        return -1;
    }

    private static int Intersection_Zero_Row_Col(float[][] A ,int r) {
        int n = A[0].length;
        for (int c = 0 ;c < n ;c++)
            if (Is_Zero_Col(A,c) && Is_Zero_Row(A,r))
                return c;
        return -1;
    }

    // duplicate matrix value of M to new matrix
    private static float[][] Copy_Matrix(float[][] A) {
        int m = A.length ,n = A[0].length;
        float[][] copyM = new float[m][n];
        for (int i = 0 ;i < m ;i++)
            for (int j = 0 ;j < n ;j++)
                copyM[i][j] = A[i][j];
        return copyM;
    }

    private static String convertDecimalToFraction(float x){
        if (x < 0)
            return "-" + convertDecimalToFraction(-x);
        float tolerance = (float) 1.0E-5;
        float h1 = 1 ,h2 = 0;
        float k1 = 0 ,k2 = 1;
        float b = x;
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
        if (k1 != 1)
            return String.valueOf((int)h1 + "/" + (int)k1);
        else
            return String.valueOf((int)h1);
    }

    private static float roundTwoDecimals(float d) {
        DecimalFormat twoDForm = new DecimalFormat("#.###");
        return Float.valueOf(twoDForm.format(d));
    }

    private static float[][] Upper_Ranking_Method(float[][] A ,float[][] b ,int t ,String ch) throws Exception {
        int n = A.length;
        for (int i = 0 ;i < n ;i++) {
            A[i][i] = (A[i][i] >= -0.0001 && A[i][i] <= 0.0001) ? 0 : A[i][i];
            if (Is_Zero_Row(A,i) && !Is_Linear_Dependent_Rows(A) && !Is_Exist_Vector(A,i)) {
                System.out.println("define new column in b when x" + (i + 1) + " is free variable in R" + n + ":");
                b = Increase_Cols_in_Vector(b);
                A[i][i] = 1;
                b[i][++t] = 1;
                Print_State(A,b,ch);
            }
            if (Is_Zero_Row(A,i) && !Is_Linear_Dependent_Rows(A)) {
                int d1 = Intersection_Zero_Row_Col(A,i);
                int d2 = Get_Linear_Dependent_Columns(A);
                if (d1 != -1) {
                    System.out.println("define new column in b when x" + (d1 + 1) + " is free variable in R" + n + ":");
                    b = Increase_Cols_in_Vector(b);
                    A[i][d1] = 1;
                    b[i][++t] = 1;
                    Print_State(A,b,ch);
                } else if (d2 != -1) {
                    System.out.println("define new column in b when x" + (d2 + 1) + " is free variable in R" + n + ":");
                    b = Increase_Cols_in_Vector(b);
                    A[i][d2] = 1;
                    b[i][++t] = 1;
                    Print_State(A,b,ch);
                } else if (!Is_Exist_Vector(A,i)) {
                    System.out.println("define new column in b when x" + (i + 1) + " is free variable in R" + n + ":");
                    b = Increase_Cols_in_Vector(b);
                    A[i][i] = 1;
                    b[i][++t] = 1;
                    Print_State(A,b,ch);
                }
            }
            if (A[i][i] == 0) {
                int r = Get_Index_UnZero_Value(A,i,true);
                int l = Get_Index_Row_from_Matrix(A,i);
                if (Is_Exist_Vector(A,i) && l < i)
                    r = Get_Index_UnZero_Value(A,i,false);
                if (r >= 0 && r < n && r != i) {
                    A[r][i] = (A[r][i] >= -0.0001 && A[r][i] <= 0.0001) ? 0 : A[r][i];
                    Retreat_Elementaric_Action(i,r);
                    Retreat_Rows_System(A,b,i,r);
                    Print_State(A,b,ch);
                }
            }
            for (int j = i + 1 ;j < n ;j++) {
                if (A[i][i] != 0 && A[j][i] != 0) {
                    float c = A[j][i] / A[i][i];
                    Sum_Elementaric_Action(c,j,i,ch);
                    for (int k = 0 ;k < n ;k++) {
                        A[j][k] -= A[i][k] * c;
                        if (k <= t) {
                            b[j][k] -= b[i][k] * c;
                            b[j][k] = (b[j][k] >= -0.0001 && b[j][k] <= 0.0001) ? 0 : b[j][k];
                        }
                    }
                    A[j][i] = (A[j][i] >= -0.0001 && A[j][i] <= 0.0001) ? 0 : A[j][i];
                    Print_State(A,b,ch);
                }
                A[j][j] = (A[j][j] >= -0.0001 && A[j][j] <= 0.0001) ? 0 : A[j][j];
                if (Is_Unit_Vector(A,j)) {
                    int d = Get_Index_for_Unit_Vector(Get_Row_from_Matrix(A,j));
                    if (d != -1 && A[j][d] != 0 && A[j][d] != 1) {
                        float c = 1 / A[j][d];
                        Mult_Elementaric_Action(c,j,ch);
                        for (int k = 0 ;k <= t ;k++)
                            b[j][k] /= A[j][d];
                        b[j][0] = (float) (Math.round(b[j][0] * 10000.0) / 10000.0);
                        A[j][d] = 1;
                        Print_State(A,b,ch);
                    }
                }
                if (Is_Zero_Row(A,j) && !Is_Zero_Row(b,j)) {
                    System.out.println("we got certain legality for free variables ,so:\n" +
                            " * express last free variable created by other elements.\n" +
                            " * putting him the rest of the elements in each line.\n" +
                            " * lower a last column from b.");
                    b = Decrease_Cols_in_Vector(b, j);
                    Print_State(A,b,ch);
                    t--;
                }
            }
            if (Is_Upper_Triangular(A) && !Is_Lower_Triangular(A))
                return Lower_Ranking_Method(A,b,t,ch);
        }
        if (!Is_Lower_Triangular(A))
            return Lower_Ranking_Method(A,b,t,ch);
        return b;
    }

    private static float[][] Lower_Ranking_Method(float[][] A ,float[][] b ,int t ,String ch) throws Exception {
        int n = A.length;
        for (int i = n - 1 ;i >= 0 ;i--) {
            A[i][i] = (A[i][i] >= -0.0001 && A[i][i] <= 0.0001) ? 0 : A[i][i];
            if (Is_Zero_Row(A,i) && !Is_Linear_Dependent_Rows(A) && !Is_Exist_Vector(A,i)) {
                System.out.println("define new column in b when x" + (i + 1) + " is free variable in R" + n + ":");
                b = Increase_Cols_in_Vector(b);
                A[i][i] = 1;
                b[i][++t] = 1;
                Print_State(A,b,ch);
            }
            if (Is_Zero_Row(A,i) && !Is_Linear_Dependent_Rows(A)) {
                int d1 = Intersection_Zero_Row_Col(A,i);
                int d2 = Get_Linear_Dependent_Columns(A);
                if (d1 != -1) {
                    System.out.println("define new column in b when x" + (d1 + 1) + " is free variable in R" + n + ":");
                    b = Increase_Cols_in_Vector(b);
                    A[i][d1] = 1;
                    b[i][++t] = 1;
                    Print_State(A,b,ch);
                } else if (d2 != -1) {
                    System.out.println("define new column in b when x" + (d2 + 1) + " is free variable in R" + n + ":");
                    b = Increase_Cols_in_Vector(b);
                    A[i][d2] = 1;
                    b[i][++t] = 1;
                    Print_State(A,b,ch);
                } else if (!Is_Exist_Vector(A,i)) {
                    System.out.println("define new column in b when x" + (i + 1) + " is free variable in R" + n + ":");
                    b = Increase_Cols_in_Vector(b);
                    A[i][i] = 1;
                    b[i][++t] = 1;
                    Print_State(A,b,ch);
                }
            }
            if (A[i][i] == 0) {
                int r = Get_Index_UnZero_Value(A,i,false);
                int l = Get_Index_Row_from_Matrix(A,i);
                if (Is_Exist_Vector(A,i) && l < i)
                    r = Get_Index_UnZero_Value(A,i,true);
                if (r >= 0 && r < n && r != i) {
                    A[r][i] = (A[r][i] >= -0.0001 && A[r][i] <= 0.0001) ? 0 : A[r][i];
                    Retreat_Elementaric_Action(i,r);
                    Retreat_Rows_System(A,b,i,r);
                    Print_State(A,b,ch);
                }
            }
            for (int j = i - 1 ;j >= 0 ;j--) {
                if (A[i][i] != 0 && A[j][i] != 0) {
                    float c = A[j][i] / A[i][i];
                    Sum_Elementaric_Action(c,j,i,ch);
                    for (int k = n - 1 ;k >= 0 ;k--) {
                        A[j][k] -= A[i][k] * c;
                        if (k <= t) {
                            b[j][k] -= b[i][k] * c;
                            b[j][k] = (b[j][k] >= -0.0001 && b[j][k] <= 0.0001) ? 0 : b[j][k];
                        }
                    }
                    A[j][i] = (A[j][i] >= -0.0001 && A[j][i] <= 0.0001) ? 0 : A[j][i];
                    Print_State(A,b,ch);
                }
                A[j][j] = (A[j][j] >= -0.0001 && A[j][j] <= 0.0001) ? 0 : A[j][j];
                if (Is_Unit_Vector(A,j)) {
                    int d = Get_Index_for_Unit_Vector(Get_Row_from_Matrix(A,j));
                    if (d != -1 && A[j][d] != 0 && A[j][d] != 1) {
                        float c = 1 / A[j][d];
                        Mult_Elementaric_Action(c,j,ch);
                        for (int k = 0 ;k <= t ;k++)
                            b[j][k] /= A[j][d];
                        b[j][0] = (float) (Math.round(b[j][0] * 10000.0) / 10000.0);
                        A[j][d] = 1;
                        Print_State(A,b,ch);
                    }
                }
                if (Is_Zero_Row(A,j) && !Is_Zero_Row(b,j)) {
                    System.out.println("we got certain legality for free variables ,so:\n" +
                            " * express last free variable created by other elements.\n" +
                            " * putting him the rest of the elements in each line.\n" +
                            " * lower a last column from b.");
                    b = Decrease_Cols_in_Vector(b, j);
                    Print_State(A,b,ch);
                    t--;
                }
            }
            if (!Is_Upper_Triangular(A) && Is_Lower_Triangular(A))
                return Upper_Ranking_Method(A,b,t,ch);
        }
        if (!Is_Upper_Triangular(A))
            return Upper_Ranking_Method(A,b,t,ch);
        return b;
    }

    private static float[][] Parallel_Ranking_Method(float[][] A ,float[][] b ,String ch) throws Exception {
        System.out.println("transform A matrix to I by elementary actions:");
        int t = 0;
        while (!Is_Unit_Matrix(A)) {
            int n = A.length;
            for (int i = 0 ;i < n ;i++) {
                A[i][i] = (A[i][i] >= -0.0001 && A[i][i] <= 0.0001) ? 0 : A[i][i];
                if (Is_Zero_Row(A,i) && !Is_Linear_Dependent_Rows(A)) {
                    int d1 = Intersection_Zero_Row_Col(A,i);
                    int d2 = Get_Linear_Dependent_Columns(A);
                    if (d1 != -1) {
                        System.out.println("define new column in b when x" + (d1 + 1) + " is free variable in R" + n + ":");
                        b = Increase_Cols_in_Vector(b);
                        A[i][d1] = 1;
                        b[i][++t] = 1;
                        Print_State(A,b,ch);
                    } else if (d2 != -1) {
                        System.out.println("define new column in b when x" + (d2 + 1) + " is free variable in R" + n + ":");
                        b = Increase_Cols_in_Vector(b);
                        A[i][d2] = 1;
                        b[i][++t] = 1;
                        Print_State(A,b,ch);
                    } else if (!Is_Exist_Vector(A,i)) {
                        System.out.println("define new column in b when x" + (i + 1) + " is free variable in R" + n + ":");
                        b = Increase_Cols_in_Vector(b);
                        A[i][i] = 1;
                        b[i][++t] = 1;
                        Print_State(A,b,ch);
                    }
                }
                if (A[i][i] == 0) {
                    int r = Get_Index_UnZero_Value(A,i,true);
                    if (r >= 0 && r < n && r != i) {
                        A[r][i] = (A[r][i] >= -0.0001 && A[r][i] <= 0.0001) ? 0 : A[r][i];
                        Retreat_Elementaric_Action(i,r);
                        Retreat_Rows_System(A,b,i,r);
                        Print_State(A,b,ch);
                    }
                }
                for (int j = 0 ;j < n ;j++) {
                    if (i != j && A[i][i] != 0 && A[j][i] != 0) {
                        float c = A[j][i] / A[i][i];
                        Sum_Elementaric_Action(c,j,i,ch);
                        for (int k = 0 ;k < n ;k++) {
                            A[j][k] -= A[i][k] * c;
                            if (k <= t)
                                b[j][k] -= b[i][k] * c;
                        }
                        A[j][i] = (A[j][i] >= -0.0001 && A[j][i] <= 0.0001) ? 0 : A[j][i];
                        Print_State(A,b,ch);
                    }
                    A[j][j] = (A[j][j] >= -0.0001 && A[j][j] <= 0.0001) ? 0 : A[j][j];
                    if (Is_Unit_Vector(A,j)) {
                        int d = Get_Index_for_Unit_Vector(Get_Row_from_Matrix(A,j));
                        if (d != -1 && A[j][d] != 0 && A[j][d] != 1) {
                            float c = 1 / A[j][d];
                            Mult_Elementaric_Action(c,j,ch);
                            for (int k = 0 ;k <= t ;k++)
                                b[j][k] /= A[j][d];
                            b[j][0] = (float) (Math.round(b[j][0] * 10000.0) / 10000.0);
                            A[j][d] = 1;
                            Print_State(A,b,ch);
                        }
                    }
                }
            }
        }
        return b;
    }

    private static float[][] Parallel_Ranking_Method_Rec(float[][] A ,float[][] b ,int i ,int j ,int t ,String ch) throws Exception {
        if (Is_Unit_Matrix(A)) {
            return b;
        } else {
            int n = A.length;
            A[i][i] = (A[i][i] >= -0.0001 && A[i][i] <= 0.0001) ? 0 : A[i][i];
            if (Is_Zero_Row(A,i) && !Is_Linear_Dependent_Rows(A)) {
                int d1 = Intersection_Zero_Row_Col(A,i);
                int d2 = Get_Linear_Dependent_Columns(A);
                if (d1 != -1) {
                    System.out.println("define new column in b when x" + (d1 + 1) + " is free variable in R" + n + ":");
                    b = Increase_Cols_in_Vector(b);
                    A[i][d1] = 1;
                    b[i][++t] = 1;
                    Print_State(A,b,ch);
                } else if (d2 != -1) {
                    System.out.println("define new column in b when x" + (d2 + 1) + " is free variable in R" + n + ":");
                    b = Increase_Cols_in_Vector(b);
                    A[i][d2] = 1;
                    b[i][++t] = 1;
                    Print_State(A,b,ch);
                } else if (!Is_Exist_Vector(A,i)) {
                    System.out.println("define new column in b when x" + (i + 1) + " is free variable in R" + n + ":");
                    b = Increase_Cols_in_Vector(b);
                    A[i][i] = 1;
                    b[i][++t] = 1;
                    Print_State(A,b,ch);
                }
            }
            if (A[i][i] == 0) {
                int r = Get_Index_UnZero_Value(A,i,true);
                if (r >= 0 && r < n && r != i) {
                    A[r][i] = (A[r][i] >= -0.0001 && A[r][i] <= 0.0001) ? 0 : A[r][i];
                    Retreat_Elementaric_Action(i,r);
                    Retreat_Rows_System(A,b,i,r);
                    Print_State(A,b,ch);
                }
            }
            if (i != j && A[i][i] != 0 && A[j][i] != 0) {
                float c = A[j][i] / A[i][i];
                Sum_Elementaric_Action(c,j,i,ch);
                for (int k = 0 ;k < n ;k++) {
                    A[j][k] -= A[i][k] * c;
                    if (k <= t)
                        b[j][k] -= b[i][k] * c;
                }
                A[j][i] = (A[j][i] >= -0.0001 && A[j][i] <= 0.0001) ? 0 : A[j][i];
                Print_State(A,b,ch);
            }
            A[j][j] = (A[j][j] >= -0.0001 && A[j][j] <= 0.0001) ? 0 : A[j][j];
            if (Is_Unit_Vector(A,j)) {
                int d = Get_Index_for_Unit_Vector(Get_Row_from_Matrix(A,j));
                if (d != -1 && A[j][d] != 0 && A[j][d] != 1) {
                    float c = 1 / A[j][d];
                    Mult_Elementaric_Action(c,j,ch);
                    for (int k = 0 ;k <= t ;k++)
                        b[j][k] /= A[j][d];
                    b[j][0] = (float) (Math.round(b[j][0] * 10000.0) / 10000.0);
                    A[j][d] = 1;
                    Print_State(A,b,ch);
                }
            }
            if (j == n - 1)
                i = (i + 1) % n;
            j = (j + 1) % n;
            return Parallel_Ranking_Method_Rec(A,b,i,j,t,ch);
        }
    }

    // check linear equations system Ax = b and solve in method requested
    private static void Check_System(float[][] A ,float[] b) throws Exception {
        int m = A.length ,n = A[0].length ,k = b.length;
        if (m <= n && m == k) {
            Scanner sc = new Scanner(System.in);
            Display_Exercise(A,b);
            if (n > 1) {
                if (Is_Linear_Independent_System(A,b))
                    System.out.println("does not exists solutions for this system");
                else {
                    User_Menu_System();
                    int op = sc.nextInt();
                    A = (m < n) ? Increase_Rows_in_Matrix(A,m) : A;
                    float[][] bt = Increase_Rows_in_Vector(b,n);
                    User_Menu_Solution();
                    String ch = sc.next();
                    Print_State(A,bt,ch);
                    float[][] x;
                    if (op == 1) {
                        x = Upper_Ranking_Method(A,bt,0,ch);
                        Print_Solution(x,ch);
                    } else if (op == 2) {
                        x = Lower_Ranking_Method(A,bt,0,ch);
                        Print_Solution(x,ch);
                    } else if (op == 3) {
                        x = Parallel_Ranking_Method(A,bt,ch);
                        Print_Solution(x,ch);
                    } else if (op == 4) {
                        x = Parallel_Ranking_Method_Rec(A,bt,0,0,0,ch);
                        Print_Solution(x,ch);
                    } else
                        throw new Exception("was inserted error value for option number solution");
                }
            } else { // R1 Space
                User_Menu_Solution();
                String ch = sc.next();
                if (A[0][0] == 0) {
                    if (b[0] == 0)
                        System.out.println("exists infinite number of solutions for the equation that is: x = ?? (free variable) when ?? ??? R");
                    else
                        System.out.println("does not exists any solution to the equation");
                } else {
                    if ((b[0] / A[0][0]) % 1 == 0)
                        System.out.println("exist single solution for the equation that is: x = " + (int)(b[0] / A[0][0]));
                    else if (ch.equals("d"))
                        System.out.println("exist single solution for the equation that is: x = " + b[0] / A[0][0]);
                    else if (ch.equals("r"))
                        System.out.println("exist single solution for the equation that is: x = " + convertDecimalToFraction(b[0] / A[0][0]));
                }
            }
        } else
            throw new Exception("the input is not exist conditions for system linear equations");
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

    // test main
    public static void main(String[] args) {
        float[][] A71 = {{0,-4,-2,7,7,1,0},{0,2,3,1,2,3,0},{0,-1,-1,-5,2,-3,-1},{2,1,3,0,-5,1,-6},{2,4,0,4,-3,2,4},{4,2,6,3,0,-4,-2},{6,2,5,-3,-3,-4,-2}};
        float[] b71 = {3,8,9,11,14,18,18};
        // x = (4.3, 9, -5.2, 0.4, 3.8, -0.8, -4.8)
        float[][] A72 = {{2,3,1,-4,0,-3,0},{-3,1,1,1,0,-4,-1},{0,1,0,-2,1,-1,1},{-4,1,-3,1,0,-2,1},{1,-3,0,-2,-4,1,0},{1,-2,3,0,-4,-2,-4},{0,4,-4,-2,-3,-2,3}};
        float[] b72 = {0,0,0,0,0,0,0};
        // upper & parallel: x = ??*( 2.5 , -2.25 , -1.5 , 1.25 , 1 , -2.75 , 1 ) = ??*( 5/2 , -9/4 , -3/2 , 5/4 , 1 , -11/4 , 1 )
        // lower: x = ??*( 1 , -0.9 , -0.6 , 0.5 , 0.4 , -1.1 , 0.4 ) = ??*( 1 , -9/10 , -3/5 , 1/2 , 2/5 , -11/10 , 2/5 )
        float[][] A81 = {{4,4,-4,-1,-4,-1,4,-1},{0,4,-4,4,0,3,1,4},{-4,3,-3,0,0,1,-2,-2},{0,3,-2,2,4,0,0,0},{1,-1,2,0,-3,2,2,0},{3,-2,0,-1,0,4,-4,1},{-2,3,-1,3,-1,1,4,4},{-4,-1,0,-4,-4,-3,-1,-4}};
        float[] b81 = {0,0,0,0,0,0,0,0};
        // upper & parallel: x = ??*( 1.5 , 6 , 5.5 , 1.5 , -2.5 , -2.5 , -4.5 , 1 ) = ??*( 3/2 , 6 , 11/2 , 3/2 , -5/2 , -5/2 , -9/2 , 1 )
        // lower: x = ??*( 1 , 4 , 3.667 , 1 , -1.667 , -1.667 , -3 , 0.667 ) = ??*( 1 , 4 , 11/3 , 1 , -5/3 , -5/3 , -3 , 2/3 )
        float[][] A82 = {{2,0,2,0,-1,1,-3,2},{-1,-3,-1,-3,-4,0,-2,-1},{4,3,-1,2,4,0,0,0},{0,-3,1,0,0,2,-1,4},{2,-3,-2,4,2,-1,2,-4},{-3,-2,1,2,3,1,4,-2},{0,0,-4,-4,2,-3,-2,3},{0,-3,-1,3,0,-3,1,-3}};
        float[] b82 = {0,0,0,0,0,0,0,0};
        // upper & parallel: x = ??*( 9 , -3 , 7 , -8 , -1 , -5 , 10 , 1 ) = ??*( 9 , -3 , 7 , -8 , -1 , -5 , 10 , 1 )
        // lower: x = ??*( 1 , -0.333 , 0.778 , -0.889 , -0.111 , -0.556 , 1.111 , 0.111 ) = ??*( 1 , -1/3 , 7/9 , -8/9 , -1/9 , -5/9 , 10/9 , 1/9 )
        try {
            Check_System(A81,b81);
        } catch (Exception ex) {
            ex.printStackTrace();
        }
    }
}
