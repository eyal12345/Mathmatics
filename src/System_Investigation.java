public class System_Investigation {

    /////////////////////////////////////////////// Print Methods /////////////////////////////////////////////////
    // display the system Ax = b in the linear equations format
    public static void Display_Exercise(float[][] A, float[] b) {
        int n = A.length;
        if (n == 1) {
            System.out.println("investigate the next equation (" + 1 + " variable):");
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
            System.out.println("investigate the next system (" + n + " variables):");
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

    // check if all vector values is equal to each other
    public static boolean Is_Equals_Values(float[] r) {
        for (int i = 0; i < r.length - 1; i++) {
            if (r[i] != r[i + 1]) {
                return false;
            }
        }
        return true;
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

    //////////////////////////////////////////////// Investigation //////////////////////////////////////////////////
    // investigate system of linear equations Ax = b without solve it
    public static void Investigation(float[][] A, float[] b) throws Exception {
        if (A.length == A[0].length && A.length == b.length) {
            int n = b.length;
            Display_Exercise(A,b);
            if (n > 1) { // R2 Space or higher
                System.out.println("system of linear equations from the shape Ax = b");
                if (Determinant(A) == 0) {
                    if (Is_Zero_Vector(b)) {
                        if (Is_Zero_Matrix(A)) {
                            System.out.print("exists an infinite number of solutions to the system in a space R" + n);
                        } else if (Is_Linear_Dependent_Rows(A)) {
                            System.out.print("exist a single solution for the system which is: x = ");
                        } else {
                            System.out.print("the solution is an infinite set of linearly dependent vectors in the base x");
                        }
                    } else {
                        System.out.print("does not an exists solutions because det(A) = 0");
                    }
                } else {
                    System.out.print("exist a single solution for a system");
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
            Investigation(A4,a4);
        } catch (Exception ex) {
            ex.printStackTrace();
        }
    }
}
