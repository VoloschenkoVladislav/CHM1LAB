using System;
using System.Collections.Generic;

namespace Com_Methods {

    class QR_Method {
        public Matrix Q, R;

        public void Gramm_Schmidt_Initialize(Matrix A)
        {
            Q = new Matrix(A.M, A.M);
            for (int i = 0; i < Q.M; i++)
            {
                for (int j = 0; j < Q.N; j++)
                {
                    if (i == j)
                        Q.Elem[i][j] = 1;
                    else
                        Q.Elem[i][j] = 0;
                }
            }
            R = new Matrix(A.M, A.N);

            Vector q = new Vector(A.M);

            for (int j = 0; j < A.N; j++)
            {
                for (int elem = 0; elem < A.N; elem++) q.Elem[elem] = A.Elem[elem][j];


                for (int i = 0; i < j; i++)
                {
                    for (int k = 0; k < q.N; k++) R.Elem[i][j] += q.Elem[k] * Q.Elem[k][i];

                    for (int k = 0; k < q.N; k++) q.Elem[k] -= R.Elem[i][j] * Q.Elem[k][i];
                }

                R.Elem[j][j] = q.Norma();

                if (R.Elem[j][j] < CONST.EPS) return;

                for (int i = 0; i < A.M; i++) Q.Elem[i][j] = q.Elem[i] / R.Elem[j][j];
            }
        }

        public void Givens_Initialize(Matrix A) {

            Q = new Matrix(A.M, A.N);
            for (int i = 0; i < Q.M; i++)
            {
                for (int j = 0; j < Q.N; j++)
                {
                    if (i == j)
                        Q.Elem[i][j] = 1;
                    else
                        Q.Elem[i][j] = 0;
                }
            }
            R = new Matrix(A.M, A.N);
            R.Copy(A);

            double sinO = 0, cosO = 0;
            double between_num1, between_num2;

            for (int j = 0; j < A.N - 1; j++)
            {
                for (int i = j + 1; i < A.M; i++)
                {
                    if (Math.Abs(R.Elem[i][j]) > CONST.EPS)
                    {
                        double divider;
                        divider = Math.Sqrt(Math.Pow(R.Elem[i][j], 2) + Math.Pow(R.Elem[j][j], 2));
                        cosO = R.Elem[j][j] / divider;
                        sinO = R.Elem[i][j] / divider;

                        for (int k = j; k < R.N; k++)
                        {
                            between_num1 = cosO * R.Elem[j][k] + sinO * R.Elem[i][k];
                            between_num2 = cosO * R.Elem[i][k] - sinO * R.Elem[j][k];
                            R.Elem[j][k] = between_num1;
                            R.Elem[i][k] = between_num2;
                        }

                        for (int k = 0; k < Q.M; k++)
                        {
                            between_num1 = cosO * Q.Elem[k][j] + sinO * Q.Elem[k][i];
                            between_num2 = cosO * Q.Elem[k][i] - sinO * Q.Elem[k][j];
                            Q.Elem[k][j] = between_num1;
                            Q.Elem[k][i] = between_num2;
                        }
                    }
                }
            }          
        }

        public void Householder_Initialize(Matrix A)
        {
            Q = new Matrix(A.M, A.N);
            for (int i = 0; i < Q.M; i++)
            {
                for (int j = 0; j < Q.N; j++)
                {
                    if (i == j)
                        Q.Elem[i][j] = 1;
                    else
                        Q.Elem[i][j] = 0;
                }
            }
            R = new Matrix(A.M, A.N);
            R.Copy(A);

            Vector p = new Vector(A.M);

            for (int i = 0; i < A.N - 1; i++)
            {
                double norma = 0;
                for (int i_ = i; i_ < A.M; i_++) norma += Math.Pow(R.Elem[i_][i], 2);

                if ( Math.Abs(norma - Math.Pow(R.Elem[i][i], 2) ) > CONST.EPS)
                {
                    double b;
                    b = Math.Sign(-R.Elem[i][i]) * Math.Sqrt(norma);

                    double mu;

                    mu = 1 / b / (b - R.Elem[i][i]);

                    for (int i_ = 0; i_ < R.M; i_++)
                    {
                        p.Elem[i_] = 0;
                        if (i_ >= i)
                            p.Elem[i_] = R.Elem[i_][i];
                    }

                    p.Elem[i] -= b;

                    for (int m = i; m < R.N; m++)
                    {
                        norma = 0;
                        for (int n = i; n < R.M; n++) norma += R.Elem[n][m] * p.Elem[n];
                        norma *= mu;
                        for (int n = i; n < R.M; n++) R.Elem[n][m] -= norma * p.Elem[n];
                    }

                    for (int m = 0; m < R.N; m++)
                    {
                        norma = 0;
                        for (int n = i; n < Q.M; n++) norma += Q.Elem[m][n] * p.Elem[n];
                        norma *= mu;
                        for (int n = i; n < Q.M; n++) Q.Elem[m][n] -= norma * p.Elem[n];
                    }
                }
            }

            Q.TransposeMatrix();
        }


        public void Solver(Vector F, Vector RES)
        {
            Vector Y;

            Matrix Qt = Q.TransposeMatrix();

            Y = Qt * F;

            Substitution_Method.Back_Row_Substitution(R, Y, RES);
        }

    }
}