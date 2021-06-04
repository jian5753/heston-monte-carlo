using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace MT_hestonSim
{
    class MonteCarloSimulation_hestonModel
    {
        private VanillaOption option;
        private double s0;
        private double var0;
        private double T;
        private double rf;

        private double rho;
        private double kappa;
        private double theta;
        private double sigma;

        private int pathLen;

        #region constructor
        public MonteCarloSimulation_hestonModel()
        {
            s0 = 0; var0 = 0; T = 0; rf = 0;
            rho = 0; kappa = 0; theta = 0; sigma = 0;
            pathLen = 0;
        }

        public MonteCarloSimulation_hestonModel(
            VanillaOption option,
            double rho, double kappa, double theta, double sigma, int pathLen)
        {
            this.option = option;
            s0 = option.getS0(); var0 = option.getVar0(); T = option.getT(); rf = option.getRf();
            this.rho = rho; this.kappa = kappa; this.theta = theta; this.sigma = sigma;
            this.pathLen = pathLen;
        }

        #endregion

        public double[] drawSPath(Random rv)
        {
            double[,] corrData = { { 1, rho }, { rho, 1 } };
            Mtrx testCorr = new Mtrx(2, 2, ref corrData);

            Zmtrx test = new Zmtrx(2, pathLen, rv);
            Mtrx upTri = testCorr.choleskyDecomp();
            Mtrx dotted = upTri.T().dot(test);

            SVpath hestonPath = new SVpath(
                ref dotted,
                option.getS0(), option.getVar0(),
                kappa, theta, sigma,
                option.getRf(), option.getT() / pathLen);
            return hestonPath.getSpath();

        }

        public double[] drawVPath(Random rv)
        {
            double[,] corrData = { { 1, rho }, { rho, 1 } };
            Mtrx testCorr = new Mtrx(2, 2, ref corrData);

            Zmtrx test = new Zmtrx(2, pathLen, rv);
            Mtrx upTri = testCorr.choleskyDecomp();
            Mtrx dotted = upTri.T().dot(test);

            SVpath hestonPath = new SVpath(
                ref dotted,
                option.getS0(), option.getVar0(),
                kappa, theta, sigma,
                option.getRf(), option.getT() / pathLen);
            return hestonPath.getVPath();
        }

        public double[][] drawSandVPath(Random rv)
        {
            double[,] corrData = { { 1, rho }, { rho, 1 } };
            Mtrx testCorr = new Mtrx(2, 2, ref corrData);

            Zmtrx test = new Zmtrx(2, pathLen, rv);
            Mtrx upTri = testCorr.choleskyDecomp();
            Mtrx dotted = upTri.T().dot(test);

            SVpath hestonPath = new SVpath(
                ref dotted,
                option.getS0(), option.getVar0(),
                kappa, theta, sigma,
                option.getRf(), option.getT());
            return hestonPath.getSandVPath();
        }

        public double[] drawSt(int pathCnt, Random rv)
        {
            double[] StArr = new double[pathCnt];

            double[,] corrData = { { 1, rho }, { rho, 1 } };
            Mtrx testCorr = new Mtrx(2, 2, ref corrData);

            Parallel.For(0, pathCnt, i =>
            {

                Zmtrx test = new Zmtrx(2, pathLen, rv);
                Mtrx upTri = testCorr.choleskyDecomp();
                Mtrx dotted = upTri.T().dot(test);

                SVpath hestonPath = new SVpath
                (
                    ref dotted,
                    option.getS0(), option.getVar0(),
                    kappa, theta, sigma,
                    option.getRf(), option.getT() / pathLen
                );

                StArr[i] = hestonPath.getSt();
            });
            return StArr;
        }
        public double[] drawSt(int pathCnt)
        {
            double[] StArr = new double[pathCnt];

            double[,] corrData = { { 1, rho }, { rho, 1 } };
            Mtrx testCorr = new Mtrx(2, 2, ref corrData);

            Parallel.For(0, pathCnt, i =>
            {
                /* 
                  fix rv for test
                 */
                Random rv = new Random(1234);
                Zmtrx test = new Zmtrx(2, pathLen);
                Mtrx upTri = testCorr.choleskyDecomp();
                Mtrx dotted = upTri.T().dot(test);

                SVpath hestonPath = new SVpath
                (
                    ref dotted,
                    option.getS0(), option.getVar0(),
                    kappa, theta, sigma,
                    option.getRf(), option.getT() / pathLen
                );

                StArr[i] = hestonPath.getSt();

            });
            return StArr;
        }

        public double meanPrice(double[] stArr, int pathCnt)
        {
            double ans = 0;
            for (int i = 0; i < pathCnt; i++)
            {
                ans += option.payoff(stArr[i]);
            }
            ans /= pathCnt;
            return ans * Math.Exp(-rf * T);
        }

        public double meanPrice(Random rv, int pathCnt)
        {
            double[] stArr = drawSt(pathCnt, rv);
            return meanPrice(stArr, pathCnt);
        }
    }
}
