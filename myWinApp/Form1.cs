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

        VanillaCall theCall;

        int trialCnt = 1; int pathCnt = 10000; int pathLen = 365;

        double s0 = 101.52;
        double k = 100.0;
        double var0 = 0.00770547621786487;
        double rf = 0.001521;
        double T = 1.0;

        double kappa = 2.20366282736578;
        double theta = 0.0164951784035976;
        double sigma = 0.33220849746904;
        double rho = -0.277814270110106;
        
        double dt = 1.0 / 365.0;
        
        private void button1_Click_1(object sender, EventArgs e)
        {
            #region old code
            /*
            normSample test = new normSample(2, 1000);
            test.draw();

            try { rho = double.Parse(textBox_rho.Text); }
            catch { rho = 0; }

            double[] corrData = new double[4] { 1, rho, rho, 1 };
            Matrix testCorr = new Matrix(corrData, 2, 2);

            double[] var = new double[2] { 1, 1 };
            Matrix testCov = Matrix.corrToCov(testCorr, var);

            test.inverseCholesky();

            Matrix upTri = testCov.choleskyDecomp();
            Matrix dotted = Matrix.transpose(upTri) * test.data;

            double testMean = Matrix.sampleMean(dotted.getRow(0));
            double testCov2 = Matrix.sampleCov(dotted.getRow(0), dotted.getRow(1));
            Matrix testCovMtrx = Matrix.sampleCovMtrx(dotted);
            richTextBox1.Text = testCovMtrx.ToString();

            

            hestonSVpaths t1 = new hestonSVpaths(ref dotted, s0, v0, kappa, theta, sigma, rf, dt);

            Simulation_singleTrial t2 = new Simulation_singleTrial(
                5, 1000, s0, v0, kappa, theta, sigma, rho, rf, dt
            );

            Simulation t3 = new Simulation(
                5, 5, 1000, s0, v0, kappa, theta, sigma, rho, rf, dt
                );
            t1.drawVpath(chart1.Series["volatility"]);
            //t1.drawSpath(chart1.Series["stock"]);
            //t2.drawSingleSpath(chart1.Series["stock"], 2);
            t3[2, 3].drawSpath(chart1.Series["stock"]);
            //chart1.ChartAreas[0].AxisY.Minimum = t1.getSpath().Min() * 0.9;
            chart1.ChartAreas[0].AxisY2.Minimum = t1.getVpath().Min() * 0.9;
            */
            #endregion
            
            Random rv = new Random(1234);
            

            try { rho = double.Parse(textBox_rho_old.Text); }
            catch { rho = 0; }
            rho = -0.277814270110106;

            double[,] corrData =  { { 1, rho } , { rho, 1 }  };
            Mtrx testCorr = new Mtrx(2, 2, ref corrData);

            Stopwatch SW = new Stopwatch();
            SW.Start();

            int pathCnt = 10000;
            double[] StArr = new double[pathCnt];
            for(int i = 0; i < pathCnt; i++)
            {
                Zmtrx test = new Zmtrx(2, 365, rv);
                Mtrx upTri = testCorr.choleskyDecomp();
                Mtrx dotted = upTri.T().dot(test);

                SVpath testPath = new SVpath(ref dotted, s0, var0, kappa, theta, sigma, rf, dt);
                StArr[i] = testPath.getSt();
            }
            SW.Stop();
            double timeConsumption = SW.ElapsedMilliseconds;
        }

        private void pricingButton_Click(object sender, EventArgs e)
        {
            /*
            Call testCall = new Call(trialCnt, pathCnt, pathLen, s0, v0, k,
                kappa, theta, sigma, rho, rf, 1);
            testCall.sim();
            priceResult.Text = testCall.getPrice().ToString("F4");*/
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

            //double[] sPath = testModel.drawSPath(new Random(1234));

            Stopwatch SW = new Stopwatch();
            SW.Start();

            //double[] stArr = simForCall.drawSt(50000, new Random(1234));

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
