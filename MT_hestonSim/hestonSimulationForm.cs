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

namespace MT_hestonSim
{
    public partial class HestonSimuoationForm : Form
    {
        #region attributes
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
        #endregion

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

        public HestonSimuoationForm()
        {
            InitializeComponent();
            #region default parameters;
            s0 = 101.52;
            k = 100.0;
            var0 = 0.00770547621786487;
            rf = 0.001521;
            T = 1.0;

            rho = -0.9;
            kappa = 1.5;
            theta = 0.04;
            sigma = 0.3;

            seed = 1234;
            pathCnt = 10000;
            #endregion

            #region chart

            #region spath chart
            chart_sPath.Series[0].Points.AddXY(0, 60);
            chart_sPath.Series[0].Points.AddXY(365, 60);
            chart_sPath.ChartAreas[0].AxisX.Maximum = 365.0;
            #endregion

            #region vpath chart
            chart_vPath.Series[0].Points.AddXY(0, theta);
            chart_vPath.Series[0].Points.AddXY(365, theta);
            chart_vPath.ChartAreas[0].AxisY.Maximum = 0.08;
            chart_vPath.ChartAreas[0].AxisX.Maximum = 365.0;
            #endregion

            #endregion
        }

        private void button2_Click(object sender, EventArgs e)
        {
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

            int pathLen = (int)(365 * T);
            #endregion

            VanillaCall testCall = new VanillaCall(s0, var0, k, T, rf);
            VanillaPut testPut = new VanillaPut(s0, var0, k, T, rf);
            MonteCarloSimulation_hestonModel simForCall = 
                new MonteCarloSimulation_hestonModel(testCall, rho, kappa, theta, sigma, pathLen);
            MonteCarloSimulation_hestonModel simForPut =
                new MonteCarloSimulation_hestonModel(testPut, rho, kappa, theta, sigma, pathLen);

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

        private void button_drawPath_Click(object sender, EventArgs e)
        {
            #region rv
            Random rv = new Random();
            //if (seed == 0) { rv = new Random(); }
            //else { rv = new Random(seed); };
            #endregion
            chart_sPath.Series[0].Points.Clear();
            chart_vPath.Series[0].Points.Clear();

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

            int pathLen = (int)(365 * T);
            #endregion

            #region draw path 
            VanillaOption forPath = new VanillaOption(s0, var0, k, T, rf);
            MonteCarloSimulation_hestonModel pathSim
                = new MonteCarloSimulation_hestonModel(forPath, rho, kappa, theta, sigma, pathLen);
            double[][] sAndvPath = pathSim.drawSandVPath(rv);
            #endregion

            for (int i = 0; i < pathLen; i++){
                chart_sPath.Series[0].Points.AddY(sAndvPath[0][i]);
                chart_vPath.Series[0].Points.AddY(sAndvPath[1][i]);
            }
            chart_sPath.ChartAreas[0].AxisY.Maximum = sAndvPath[0].Max() * 1.1;
            chart_sPath.ChartAreas[0].AxisY.Minimum = sAndvPath[0].Min() * 0.9;
            chart_vPath.ChartAreas[0].AxisY.Maximum = sAndvPath[1].Max() * 1.1;
            chart_vPath.ChartAreas[0].AxisY.Minimum = sAndvPath[1].Min() * 0.9;

            chart_sPath.ChartAreas[0].AxisX.Maximum = pathLen;
            chart_vPath.ChartAreas[0].AxisX.Maximum = pathLen;


        }

        private void simulationPara_Enter(object sender, EventArgs e)
        {

        }
    }
}
