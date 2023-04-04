public class Determinant_Calculate {
    public static float Determinant(float[][] A) {
        int n = A.length;
        if (n == 1) {
            return A[0][0];
        } else {
            float sum = 0;
            for (int i = 0; i < n; i++)
                sum += A[0][i] * Determinant(Sub_Matrix(A, 0, i)) * Math.pow(-1, i);
            return sum;
        }
    }
    public static float[][] Sub_Matrix(float[][] A ,int x ,int y) {
        int n = A.length ,p = 0 ,q = 0;
        float[][] subA = new float[n - 1][n - 1];
        for (int i = 0 ;i < n ;i++) {
            for (int j = 0 ;j < n ;j++) {
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
    public static float Determinant_Up_To_3x3(float[][] M) {
        int n = M.length;
        float det = (n > 1) ? 0 : M[0][0];
        for (int i = 0 ;i < n ;i++) {
            float m_d = 1 ,s_d = 1;
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
            det += (m_d - s_d);
        }
        return det;
    }
    public static float Determinant_Up_To_3x3_V2(float[][] M) {
        int n = M.length;
        float det = (n > 1) ? 0 : M[0][0];
        for (int i = 0 ;i < n ;i++) {
            float m_d = 1 ,s_d = 1;
            for (int j = 0 ;j < n ;j++) {
                int dl = n + j - i ,dr = n + j + i;
                m_d *= M[dr % n][dl % n];
                s_d *= M[(3*n - 1 - dr) % n][j];
            }
            det += (m_d - s_d);
        }
        return det;
    }
    public static void main(String[] args) {
        float[][] A1 = {{1}}; // det = 1
        float[][] A2 = {{1,-3},{2,4}}; // det = 10
        float[][] A3 = {{1,-3,1},{2,4,-5},{5,5,1}}; // det = 100
        float[][] A41 = {{2,-1,3,3},{4,-2,6,5},{6,-3,-1,2},{8,-4,-4,-5}}; // det = 0
        float[][] A42 = {{1,4,-3,-3},{3,0,-1,-1},{5,4,-4,-3},{2,1,-1,-1}}; // det = 4
        System.out.println(Determinant(A41));
        System.out.println(Determinant(A42));
    }
}
