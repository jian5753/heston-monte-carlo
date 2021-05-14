using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

using DFinNR;

namespace ConsoleApplication1
{
    class Program
    {
        static Matrix corrToCov(Matrix corrMtrx, double[] variance)
        {

            Matrix covMtrx = new Matrix(corrMtrx);
            for (int rowIdx = 0; rowIdx < corrMtrx.rows(); rowIdx++)
            {
                for (int colIdx = 0; colIdx < corrMtrx.columns(); colIdx++)
                {
                    covMtrx[rowIdx, colIdx] *= Math.Sqrt(variance[rowIdx]);
                    covMtrx[rowIdx, colIdx] *= Math.Sqrt(variance[colIdx]);
                }
            }
            return covMtrx;
        }
        static void Main(string[] args)
        {
            Random rv = new Random(1234);
            #region test1
            /*
             normSample test = new normSample(2, 1000);
            Console.WriteLine(test.showPara());
            test.draw();

            Console.WriteLine("corr to cov");
            double[] data = new double[4] { 1, -0.3, -0.3, 1};
            Matrix testCorr = new Matrix(data, 2, 2);
            double[] var = new double[2] { 0.2, 0.5 };
            Matrix testCov = corrToCov(testCorr, var);
            Console.Write(testCov.ToString());
            Console.ReadLine();

            Console.WriteLine("inverse cholesky");
            test.inverseCholesky();
            Matrix adjustedCovMtrx = Matrix.sampleCovMtrx(test.data);
            Console.WriteLine($"{adjustedCovMtrx}");
            Console.ReadLine();

            Console.WriteLine("correlated sample");
            Matrix upTri = testCov.choleskyDecomp();
            Console.WriteLine(upTri.ToString());
            Matrix dotted = Matrix.transpose(upTri) * test.data;
            Console.WriteLine(Matrix.sampleVar(dotted).ToString());
            Console.ReadLine();

            Console.WriteLine("check stats");
            //Console.WriteLine($"{dotted.getRow(0)[9999]}");
            double testMean = Matrix.sampleMean(dotted.getRow(0));
            Console.WriteLine($"{testMean}");
            double testCov2 = Matrix.sampleCov(dotted.getRow(0), dotted.getRow(1));
            Console.WriteLine($"{testCov2}");
            Matrix testCovMtrx = Matrix.sampleCovMtrx(dotted);
            Console.WriteLine(testCovMtrx.ToString());
            Console.ReadLine();
             */
            #endregion
            
            //double[,] test = new double[2, 365];
            


            #region paras
            double s0 = 101.52;
            double strike = 100.0;
            double rf = 0.001521;

            double v0 = 0.00770547621786487;
            double kappa = 2.20366282736578;
            double theta = 0.0164951784035976;
            double sigma = 0.33220849746904;
            double rho = -0.277814270110106;

            double deltat = 1.0 / 365;
            #endregion

            Matrix test = new Matrix(2, 365);
            for (int t = 0; t < 365; t++)
            {
                for (int pathIdx = 0; pathIdx < 2; pathIdx++)
                {
                    test[pathIdx, t] = DStat.N_Inv(rv.NextDouble());
                }
            }
            

            double[] data = new double[4] { 1, rho, rho, 1 };
            Matrix testCorr = new Matrix(data, 2, 2);
            double[] var = new double[2] { 1, 1 };
            Matrix testCov = corrToCov(testCorr, var);

            Matrix upTri = testCov.choleskyDecomp();
            Console.WriteLine(upTri.ToString());
            Matrix dotted = Matrix.transpose(upTri) * test;

            // double v1 = Math.Max(s0 + kappa * (theta - s0) * deltat + sigma * Math.Sqrt(v0 * deltat) * test[0, 0], 0);
            //double s1 = s0 + rf * s0 * deltat + Math.Sqrt(v0 * deltat) * s0 * test[0, 0];
            double s1 = s0 * Math.Exp((rf - 0.5 * v0) * deltat + Math.Sqrt(v0 * deltat) * dotted[0, 0]);
            double v1 = Math.Max(v0 + kappa * (theta - v0) * deltat + sigma * Math.Sqrt(v0) * Math.Sqrt(deltat) * dotted[1, 0], 0);
            Console.ReadLine();
            /*
            vol = Math.Sqrt(Math.Abs(Var));
            //if (Var < 0.0) Var = 0.0;
            //vol = Math.Sqrt(Var);
            vol2 = Sigma * vol;
            mu = r - y - 0.5 * vol * vol;
            nu = Kappa * (Theta - vol * vol);

            S = S * Math.Exp(mu * dt + vol * N1 * Sqrtdt);
            Var = vol * vol + nu * dt + vol2 * Sqrtdt * N2
             */

        }
    }
}
