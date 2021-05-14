using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows.Forms.DataVisualization.Charting;
using DFinNR;

namespace myWinApp
{
    public class hestonSVpaths
    {
        private Matrix source;
        private int length;
        private double[] vPath;
        private double[] sPath;
        private double deltat;

        public hestonSVpaths(ref Matrix source,
            double s0, double v0,
            double kappa, double theta, double sigma,
            double rf,
            double deltat)
        {
            this.source = source;
            this.length = source.columns();
            double sqrtdt = Math.Sqrt(deltat);
            vPath = new double[this.length+1];
            sPath = new double[this.length+1];
            this.deltat = deltat;

            vPath[0] = v0;
            sPath[0] = s0;

            for (int t = 0; t < this.length; t++)
            {
                double Zt_1 = this.source.row(0)[t];
                double Zt_2 = this.source.row(1)[t];
                /* double s1 = s0 * Math.Exp((rf - 0.5 * v0) * deltat + Math.Sqrt(v0 * deltat) * dotted[0, 0]);
                double v1 = Math.Max(v0 + kappa * (theta - v0) * deltat + sigma * Math.Sqrt(v0) * Math.Sqrt(deltat) * dotted[1, 0], 0);*/
                sPath[t + 1] = sPath[t] * Math.Exp((rf - 0.5 * vPath[t]) * deltat + Math.Sqrt(vPath[t]) * sqrtdt * Zt_1);
                vPath[t + 1] = Math.Max(vPath[t] + kappa * (theta - vPath[t]) * deltat + sigma * Math.Sqrt(vPath[t]) * sqrtdt * Zt_2, 0);
            }
        }
        public double[] getVpath()
        {
            return this.vPath;
        }

        public double[] getSpath()
        {
            return this.sPath;
        }

        public void drawVpath(Series target)
        {
            target.Points.Clear();
            for(int t = 0; t < this.length; t++)
            {
                target.Points.AddXY(t, this.vPath[t]);
            }
        }

        public void drawSpath(Series target)
        {
            target.Points.Clear();
            for(int t = 0; t < this.length; t++)
            {
                target.Points.AddXY(t, this.sPath[t]);
            }
        }

    }
}
