using System;

namespace Com_Methods {
    class Gauss_Method {
        public static int Find_Main_Element(Matrix A, int j)
        {
            int Index = j;
            for (int i = j + 1; i < A.M; i++)
                if (Math.Abs(A.Elem[i][j]) > Math.Abs(A.Elem[Index][j])) Index = i;
            if (Math.Abs(A.Elem[Index][j]) < CONST.EPS) throw new Exception("Gauss_Method: degenerate matrix...");
            return Index;
        }


        public void Direct_Way(Matrix A, Vector F)
        {
            double help;
            
            for (int i = 0; i < A.M - 1; i++)
            {
                int I = Find_Main_Element(A, i);

                if (I != i)
                {
                    var Help = A.Elem[I];
                    A.Elem[I] = A.Elem[i];
                    A.Elem[i] = Help;

                    help = F.Elem[i];
                    F.Elem[i] = F.Elem[I];
                    F.Elem[I] = help;
                }

                for (int j = i + 1; j < A.M; j++)
                {
                    help = A.Elem[j][i] / A.Elem[i][i];
                    A.Elem[j][i] = 0;
                    for (int k = i + 1; k < A.N; k++)
                    {
                        A.Elem[j][k] -= help * A.Elem[i][k];
                    }
                    F.Elem[j] -= help * F.Elem[i];
                }
            }
        }


        public Vector Solver(Matrix A, Vector F) {

            Direct_Way(A, F);

            Vector RES = new Vector(F.N);

            Substitution_Method.Back_Row_Substitution(A, F, RES);

            return RES;
        }
    }
}