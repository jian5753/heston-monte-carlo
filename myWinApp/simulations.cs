using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows.Forms.DataVisualization.Charting;
using DFinNR;

namespace myWinApp
{
    

    class Simulation_singleTrial
    {
        private hestonSVpaths[] pathi;

        public Simulation_singleTrial(
            int pathCnt, int pathLen, // trial parameters
            double s0, double v0,
            double kappa, double theta, double sigma, double rho, // process parameters
            double rf, double deltat
        )
        {
            this.pathi = new hestonSVpaths[pathCnt];
            for (int i =0; i < pathCnt; i++)
            {
                normSample zTemp = new normSample(2, pathLen);
                zTemp.draw();
                //zTemp.inverseCholesky();

                double[] corrData = new double[4] { 1, rho, rho, 1 };
                Matrix corrMtrx = new Matrix(corrData, 2, 2);
                double[] var = new double[2] { 1, 1 };
                Matrix covMtrx = Matrix.corrToCov(corrMtrx, var);

                
                Matrix upTri = covMtrx.choleskyDecomp();
                Matrix dotted = Matrix.transpose(upTri) * zTemp.data;

                this.pathi[i] = new hestonSVpaths(ref dotted, s0, v0, kappa, theta, sigma, rf, deltat);
            }
        }

        #region methods

        
        public hestonSVpaths this[int i]
        {
            get { return this.pathi[i]; }
            set { this.pathi[i] = value; }
        }

        public void drawSingleSpath(Series target, int pathId)
        {
            this[pathId].drawSpath(target);
        }
        public void drawSingleVpath(Series target, int pathId)
        {
            this[pathId].drawVpath(target);
        }
        #endregion
    }
    class Simulation
    {
        private int trialCnt;
        private Simulation_singleTrial[] trials;

        public Simulation(
            int trialCnt, int pathCnt, int pathLen,
            double s0, double v0, double kappa, double theta, double sigma, double rho,
            double rf, double deltat
            )
        {
            this.trials = new Simulation_singleTrial[trialCnt];
            for (int i = 0; i < trialCnt; i++)
            {
                trials[i] = new Simulation_singleTrial(
                    pathCnt, pathLen, s0, v0, kappa, theta, sigma, rho, rf, deltat
                );
            }
        }

        public hestonSVpaths this[int trialIdx, int pathIdx]
        {
            get { return this.trials[trialIdx][pathIdx]; }
            set { this.trials[trialIdx][pathIdx] = value; }
        }

    }


}
