namespace Com_Methods {

    class LU_Method {
        public Matrix LU;

        public void Initialize(Matrix A) {

            LU = new Matrix(A.M, A.N);
            LU.Copy(A);

            //Составление U (верхней треугольной матрицы)
            for (int i = 0; i < LU.M - 1; i++){

                for (int j = i + 1; j < LU.M; j++) {
                    double multiplyer;

                    multiplyer = LU.Elem[j][i] / LU.Elem[i][i];

                    LU.Elem[j][i] = 0;

                    for (int k = i + 1; k < A.N; k++) {
                        LU.Elem[j][k] -= multiplyer * LU.Elem[i][k];
                    }
                }

            }

            //Составление L (нижней треугольной матрицы с единичной диагональю)
            for (int i = 1; i < LU.M; i++) {
                
                for (int j = 0; j < i; j++) {

                    double sum = 0;
                    if (i > 1) {
                        for (int k = 0; k < j; k++) {
                            sum += LU.Elem[i][k]*LU.Elem[k][j];
                        }
                    }

                    LU.Elem[i][j] = (A.Elem[i][j] - sum) / LU.Elem[j][j];
                }

            }
        }


        public void Solver(Vector F, Vector RES) {

            Vector Y = new Vector(F.N);

            Matrix L = new Matrix(LU.M, LU.N);

            for (int i = 0; i < L.M; i++)
            {
                for (int j = 0; j < L.N; j++)
                {
                    if (i == j)
                        L.Elem[i][j] = 1;
                    else if (i < j)
                        L.Elem[i][j] = 0;
                    else
                        L.Elem[i][j] = LU.Elem[i][j];
                }
            }

            Substitution_Method.Direct_Row_Substitution(L, F, Y);            
            Substitution_Method.Back_Row_Substitution(LU, Y, RES);
        }

    }
}