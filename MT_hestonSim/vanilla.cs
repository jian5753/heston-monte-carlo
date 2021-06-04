using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace MT_hestonSim
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
            ) : base(s0, var0, k, T, rf) { }

        public VanillaCall() : base() { }
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
}
