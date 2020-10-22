using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace Com_Methods
{

    class CONST
    {
        public const double EPS = 10e-6;
    }
    class Program
    { 
        static void Main(string[] args)
        {
            List<int> testSet = new List<int> { 100, 200, 300, 400, 500 };

            //метод Гаусса
            foreach (int range in testSet)
            {
                Vector X = new Vector(range);
                for (int i = 0; i < range; i++)
                {
                    X.Elem[i] = 1;
                }

                Matrix A = new Matrix(range, range);
                A.GenerateRandomDate();

                Vector F = A * X;


                Stopwatch stopwatch = new Stopwatch();

                Gauss_Method slau = new Gauss_Method();
                stopwatch.Start();
                Vector RES = slau.Solver(A, F);
                stopwatch.Stop();
                TimeSpan ts = stopwatch.Elapsed;
                string time = String.Format("{0:00}.{1:000}", ts.Seconds, ts.Milliseconds);

                //Норма единичного вектора, коим является вектор X, равна единице
                //Потому делитель отсутствует
                double error = (RES - X).Norma();

                Console.WriteLine("Solved with Gauss method at dim(A) = " + range + ":");
                Console.WriteLine(String.Format("Time: {0}", time));
                Console.WriteLine(String.Format("Error: {0}\n", error));
            }

            //LU-разложение
            foreach (int range in testSet)
            {
                Vector X = new Vector(range);
                for (int i = 0; i < range; i++)
                {
                    X.Elem[i] = 1;
                }

                Matrix A = new Matrix(range, range);
                A.GenerateRandomDate();

                Vector F = A * X;


                Stopwatch stopwatch = new Stopwatch();

                LU_Method slau = new LU_Method();
                Vector RES = new Vector(F.N);
                stopwatch.Start();
                slau.Initialize(A);
                slau.Solver(F, RES);
                stopwatch.Stop();
                TimeSpan ts = stopwatch.Elapsed;
                string time = String.Format("{0:00}.{1:000}", ts.Seconds, ts.Milliseconds);

                //Норма единичного вектора, коим является вектор X, равна единице
                //Потому делитель отсутствует
                double error = (RES - X).Norma();

                Console.WriteLine("Solved with LU-method at dim(A) = " + range + ":");
                Console.WriteLine(String.Format("Time: {0}", time));
                Console.WriteLine(String.Format("Error: {0}\n", error));
            }


            //QR-разложение (Грамм-Шмидт)
            foreach (int range in testSet)
            {
                Vector X = new Vector(range);
                for (int i = 0; i < range; i++)
                {
                    X.Elem[i] = 1;
                }

                Matrix A = new Matrix(range, range);
                A.GenerateRandomDate();

                Vector F = A * X;


                Stopwatch stopwatch = new Stopwatch();

                QR_Method slau = new QR_Method();
                Vector RES = new Vector(F.N);
                stopwatch.Start();
                slau.Gramm_Schmidt_Initialize(A);
                slau.Solver(F, RES);
                stopwatch.Stop();
                TimeSpan ts = stopwatch.Elapsed;
                string time = String.Format("{0:00}.{1:000}", ts.Seconds, ts.Milliseconds);

                //Норма единичного вектора, коим является вектор X, равна единице
                //Потому делитель отсутствует
                double error = (RES - X).Norma();

                Console.WriteLine("Solved with QR-method (Gramm-Schmidt) at dim(A) = " + range + ":");
                Console.WriteLine(String.Format("Time: {0}", time));
                Console.WriteLine(String.Format("Error: {0}\n", error));
            }


            //QR-разложение (Гивенс)
            foreach (int range in testSet)
            {
                Vector X = new Vector(range);
                for (int i = 0; i < range; i++)
                {
                    X.Elem[i] = 1;
                }

                Matrix A = new Matrix(range, range);
                A.GenerateRandomDate();

                Vector F = A * X;


                Stopwatch stopwatch = new Stopwatch();

                QR_Method slau = new QR_Method();
                Vector RES = new Vector(F.N);
                stopwatch.Start();
                slau.Givens_Initialize(A);
                slau.Solver(F, RES);
                stopwatch.Stop();
                TimeSpan ts = stopwatch.Elapsed;
                string time = String.Format("{0:00}.{1:000}", ts.Seconds, ts.Milliseconds);

                //Норма единичного вектора, коим является вектор X, равна единице
                //Потому делитель отсутствует
                double error = (RES - X).Norma();

                Console.WriteLine("Solved with QR-method (Givens) at dim(A) = " + range + ":");
                Console.WriteLine(String.Format("Time: {0}", time));
                Console.WriteLine(String.Format("Error: {0}\n", error));
            }

            //QR-разложение (Хаусхолдер)
            foreach (int range in testSet)
            {
                Vector X = new Vector(range);
                for (int i = 0; i < range; i++)
                {
                    X.Elem[i] = 1;
                }

                Matrix A = new Matrix(range, range);
                A.GenerateRandomDate();

                Vector F = A * X;


                Stopwatch stopwatch = new Stopwatch();

                QR_Method slau = new QR_Method();
                Vector RES = new Vector(F.N);
                stopwatch.Start();
                slau.Householder_Initialize(A);
                slau.Solver(F, RES);
                stopwatch.Stop();
                TimeSpan ts = stopwatch.Elapsed;
                string time = String.Format("{0:00}.{1:000}", ts.Seconds, ts.Milliseconds);

                //Норма единичного вектора, коим является вектор X, равна единице
                //Потому делитель отсутствует
                double error = (RES - X).Norma();

                Console.WriteLine("Solved with QR-method (Householder) at dim(A) = " + range + ":");
                Console.WriteLine(String.Format("Time: {0}", time));
                Console.WriteLine(String.Format("Error: {0}\n", error));
            }

            Console.ReadLine();
        }
    }
}
