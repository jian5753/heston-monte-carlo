using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Data;
using System.Drawing;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows.Forms;
using System.Diagnostics;

namespace MT_hestonSim
{
    public partial class MT_hestonSim_form : Form
    {
        private double s0;
        private double k;
        private double var0;
        private double rf;
        private double T;

        private double rho;
        private double kappa;
        private double theta;
        private double sigma;

        private int seed;
        private int pathCnt;

        public MT_hestonSim_form()
        {
            InitializeComponent();
            #region default parameters;
            s0 = 101.52;
            k = 100.0;
            var0 = 0.00770547621786487;
            rf = 0.001521;
            T = 0.01;

            rho = -0.9;
            kappa = 1.5;
            theta = 0.04;
            sigma = 0.3;

            seed = 1234;
            pathCnt = 10000;
            #endregion
        }

        private void startTest_Click(object sender, EventArgs e)
        {
            VanillaCall testCall = new VanillaCall(s0, var0, k, T, rf);
            Random rv = new Random(1234);
            int pathLen = (int)(365 * T);
            MonteCarloSimulation_hestonModel simForCall =
                new MonteCarloSimulation_hestonModel(testCall, rho, kappa, theta, sigma, pathLen);

            Stopwatch SW = new Stopwatch();
            SW.Start();
            double[] StArr = simForCall.drawSt(2, rv);
            SW.Stop();
            double t0 = SW.ElapsedMilliseconds;
  
        }
    }
}
