using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Data;
using System.Drawing;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows.Forms;
using System.Windows.Forms.DataVisualization.Charting;
using System.Diagnostics;
using DFinNR;

namespace myWinApp
{
    public partial class Form1 : Form
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

        public Form1()
        {
            InitializeComponent();
        }

        private void button2_Click(object sender, EventArgs e)
        {
            #region default parameters;
            double s0 = 101.52;
            double k = 100.0;
            double var0 = 0.00770547621786487;
            double rf = 0.001521;
            double T = 1.0;

            double rho = -0.9;
            double kappa = 1.5;
            double theta = 0.04;
            double sigma = 0.3;

            int seed = 1234;
            int pathCnt = 10000;
            #endregion

            #region parse input
            try { s0 = double.Parse(textBox_s0.Text); } catch { };
            try { k = double.Parse(textBox_k.Text); } catch { };
            try { var0 = double.Parse(textBox_var0.Text); } catch { };
            try { T = double.Parse(textBox_T.Text); } catch { };
            try { rf = double.Parse(textBox_rf.Text); } catch { };

            try { rho = double.Parse(textBox_rho.Text); } catch { };
            try { kappa = double.Parse(textBox_kappa.Text); } catch { };
            try { theta = double.Parse(textBox_theta.Text); } catch { };
            try { sigma = double.Parse(textBox_sigma.Text); } catch { };

            try { pathCnt = int.Parse(textBox_pathCnt.Text); } catch { };
            try { seed = int.Parse(textBox_seed.Text); } catch { seed = 0; };
            #endregion

            VanillaCall testCall = new VanillaCall(s0, var0, k, T, rf);
            VanillaPut testPut = new VanillaPut(s0, var0, k, T, rf);
            MonteCarloSimulation_hestonModel simForCall = 
                new MonteCarloSimulation_hestonModel(testCall, rho, kappa, theta, sigma, 365 * T);
            MonteCarloSimulation_hestonModel simForPut =
                new MonteCarloSimulation_hestonModel(testPut, rho, kappa, theta, sigma, 365 * T);

            Stopwatch SW = new Stopwatch();
            SW.Start();
            Random rv;
            if (seed == 0) { rv = new Random(); }
            else { rv = new Random(seed); };

            double[] stArr = simForCall.drawSt(pathCnt, rv);
            double callPrice = simForCall.meanPrice(stArr, pathCnt);
            double putPrice = simForPut.meanPrice(stArr, pathCnt);
            textBox_callPrice.Text = callPrice.ToString("F4");
            textBox_putPrice.Text = putPrice.ToString("F4");


            SW.Stop();
            double timeConsumption = SW.ElapsedMilliseconds;
            msgBox.Text += $"done. time consumption {timeConsumption} ms.";
            msgBox.Text += $"call price: {textBox_callPrice.Text}, put price: {textBox_putPrice.Text}. \n";

        }

        private void clear_Click(object sender, EventArgs e)
        {
            msgBox.Text = "";
        }
    }
}
