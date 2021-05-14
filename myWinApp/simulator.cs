using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using DFinNR;

namespace myWinApp
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

        private double pathLen;

        public MonteCarloSimulation_hestonModel(){
            s0 = 0; var0 = 0; T = 0; rf = 0;
            rho = 0; kappa = 0; theta = 0; sigma = 0;
            pathLen = 0;
        }

        public MonteCarloSimulation_hestonModel(
            VanillaOption option,
            double rho, double kappa, double theta, double sigma, double pathLen)
        {
            this.option = option;
            s0 = option.getS0(); var0 = option.getVar0(); T = option.getT(); rf = option.getRf();
            this.rho = rho; this.kappa = kappa; this.theta = theta; this.sigma = sigma;
            this.pathLen = pathLen;
        }

        public double[] drawSPath(Random rv)
        {
            double[,] corrData = { { 1, rho }, { rho, 1 } };
            Mtrx testCorr = new Mtrx(2, 2, ref corrData);

            Zmtrx test = new Zmtrx(2, 365, rv);
            Mtrx upTri = testCorr.choleskyDecomp();
            Mtrx dotted = upTri.T().dot(test);

            SVpath hestonPath = new SVpath(
                ref dotted,
                option.getS0(), option.getVar0(),
                kappa, theta, sigma,
                option.getRf(), option.getT() / pathLen);
            return hestonPath.getSpath();

        }

        public double[] drawSt(int pathCnt, Random rv)
        {
            double[] StArr = new double[pathCnt];

            double[,] corrData = { { 1, rho }, { rho, 1 } };
            Mtrx testCorr = new Mtrx(2, 2, ref corrData);

            for (int i = 0; i < pathCnt; i++)
            {
                Zmtrx test = new Zmtrx(2, 365, rv);
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
            }
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

    /*
    class Simulator 
    {
        private int pathCnt; private int pathLen;
        private double s0; private double v0;
        private double kappa; private double theta; private double sigma; private double rho;
        private double rf; private double T; private double dt;

        public Simulator()
        {
            pathCnt = 128 * 128;
            pathLen = 365;
            s0 = 101.52;
            //k = 100.0;
            rf = 0.001521;

            v0 = 0.00770547621786487;
            kappa = 2.20366282736578;
            theta = 0.0164951784035976;
            sigma = 0.33220849746904;
            rho = -0.277814270110106;

            dt = 1.0 / pathLen;
        }

        public void sim()
        {
            Random rv = new Random(1234);
            for (int pathIdx = 0; pathIdx < pathCnt; pathIdx++)
            {
                #region draw z
                double[,] z = new double[2, pathLen];

                for (int timeIdx = 0; timeIdx < pathLen; timeIdx++)
                {
                    z[0, timeIdx] = DStat.N_Inv(rv.NextDouble());
                    z[1, timeIdx] = DStat.N_Inv(rv.NextDouble());
                }
                #endregion
            }
        }
    }
    */

}
