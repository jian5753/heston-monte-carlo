using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace myWinApp
{
    
    class VanillaOption
    {
        protected double s0; protected double var0; protected double k;
        protected double T; protected double rf;

        public VanillaOption()
        {
            s0 = 0.0; var0 = 0.0; k = 0.0; T = 0.0; rf = 0.0;
        }

        public VanillaOption(double s0, double var0, double k, double T, double rf)
        {
            this.s0 = s0; this.var0 = var0; this.k = k; this.T = T; this.rf = rf;
        }

        public double getS0() { return s0; }
        public double getVar0() { return var0; }
        public double getRf() { return rf; }
        public double getT() { return T; }

        public virtual double payoff(double St)
        {
            return 0.0;
        }
    }

    class VanillaCall : VanillaOption
    {
        public VanillaCall(
            double s0, double var0, double k, double T, double rf
            ):base(s0, var0, k, T, rf){}

        public VanillaCall() :base() { }
        public override double payoff(double st)
        {
            return Math.Max(st - k, 0);
        }
    }

    class VanillaPut : VanillaOption
    {
        public VanillaPut(
            double s0, double v0, double k, double T, double rf
            ) : base(s0, v0, k, T, rf) { }

        public VanillaPut() : base() { }

        public override double payoff(double st)
        {
            return Math.Max(k - st, 0);
        }
    }

    class Call
    {
        private int trialCnt; private int pathCnt; private int pathLen;
        private double s0; private double v0; private double k;
        private double kappa; private double theta; private double sigma; private double rho;
        private double rf; private double T; private double deltat; 
        private Simulation sims;
        private double price;
        public Call(int trialCnt, int pathCnt, int pathLen,
            double s0, double v0, double k,
            double kappa, double theta, double sigma, double rho,
            double rf, double T)
        {
            this.trialCnt = trialCnt; this.pathCnt = pathCnt; this.pathLen = pathLen;
            this.s0 = s0; this.v0 = v0; this.k = k;
            this.kappa = kappa; this.theta = theta;
            this.sigma = sigma; this.rho = rho;
            this.rf = rf; this.T = T;
            this.deltat = T / 365.0;
            this.price = -1;
            
        }

        public double payoffFcn(double St)
        {
            return Math.Max(St - k, 0);
        }

        public void sim()
        {
            sims = new Simulation(
                trialCnt, pathCnt, pathLen, s0, v0, kappa, theta, sigma, rho, rf, deltat
            );

            double[] mean_eachTrial = new double[trialCnt];
            price = 0;
            for (int trialIdx = 0; trialIdx < trialCnt; trialIdx++)
            {
                double payoff = 0;
                for (int pathIdx = 0; pathIdx < pathCnt; pathIdx++)
                {
                    double st = sims[trialIdx, pathIdx].getSpath()[pathLen];
                    double tempPayoff = payoffFcn(st);
                    payoff += payoffFcn(st);
                }
                mean_eachTrial[trialIdx] = payoff / pathCnt;
                price += payoff / pathCnt;
            }
            price /= trialCnt;
            price *= Math.Exp(-rf * pathLen * deltat);
        }

        public Simulation getSim()
        {
            return sims;
        }

        public double getPrice()
        {
            return price;
        }
    }
}
