using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading;
using System.Threading.Tasks;
using DFinNR;

namespace MT_hestonSim
{
    class Mtrx
    {
        protected int rowCnt;
        protected int colCnt;
        public double[,] data;

        public Mtrx(int rowCnt, int colCnt)
        {
            this.rowCnt = rowCnt;
            this.colCnt = colCnt;
            data = new double[rowCnt, colCnt];
        }
        public Mtrx(int rowCnt, int colCnt, ref double[,] data)
        {
            this.rowCnt = rowCnt;
            this.colCnt = colCnt;
            this.data = data;
        }

        #region query
        public int getColCnt()
        {
            return colCnt;
        }
        public int getRowCnt()
        {
            return rowCnt;
        }

        public double[] getRow(int rowIdx)
        {
            double[] ans = new double[colCnt];
            Parallel.For(0, colCnt, colIdx =>
           {
               ans[colIdx] = data[rowIdx, colIdx];
               //Thread.Sleep(1);
           });

            return ans;
        }

        public double[] getRow_sThread(int rowIdx)
        {
            double[] ans = new double[colCnt];
            for (int colIdx = 0; colIdx < colCnt; colIdx++)
            {
                ans[colIdx] = data[rowIdx, colIdx];
                Thread.Sleep(1);
            }
            return ans;
        }

        public double[] getCol(int colIdx)
        {
            double[] ans = new double[rowCnt];
            Parallel.For(0, rowCnt, rowIdx => {
                ans[rowIdx] = data[rowIdx, colIdx];
            });
            return ans;
        }
        #endregion

        #region operation overload
        public double this[int i, int j]
        {
            get { return data[i, j]; }
            set { this.data[i, j] = value; }
        }
        #endregion

        #region algebra operation
        public Mtrx dot(Mtrx rightMtrx)
        {
            Mtrx ans = new Mtrx(this.rowCnt, rightMtrx.getColCnt());
            /*
            Parallel.For(0, rowCnt, rowIdx =>
            {
                Parallel.For(0, colCnt, colIdx =>
                {
                    double[] vec1 = this.getRow(rowIdx);
                    double[] vec2 = rightMtrx.getCol(colIdx);
                    ans[rowIdx, colIdx] = 0;

                    Parallel.For(0, colCnt, k =>
                    {
                        ans[rowIdx, colIdx] += vec1[k] * vec2[k];
                    });
                });
            });
            */
            

            for (int rowIdx = 0; rowIdx < rowCnt; rowIdx++)
            {
                for (int colIdx = 0; colIdx < rightMtrx.getColCnt(); colIdx++)
                {
                    double[] vec1 = this.getRow(rowIdx);
                    double[] vec2 = rightMtrx.getCol(colIdx);
                    ans[rowIdx, colIdx] = 0;
                    for (int k = 0; k < colCnt; k++)
                    {
                        ans[rowIdx, colIdx] += vec1[k] * vec2[k];
                    }
                }
            }
            
            return ans;
            
        }

        public Mtrx T()
        {
            Mtrx ans = new Mtrx(colCnt, rowCnt);
            for (int rowIdx = 0; rowIdx < rowCnt; rowIdx++)
            {
                for (int colIdx = 0; colIdx < colCnt; colIdx++)
                {
                    ans.data[colIdx, rowIdx] = this.data[rowIdx, colIdx];
                }
            }
            return ans;
        }

        public Mtrx choleskyDecomp()
        {
            Mtrx theAns = new Mtrx(rowCnt, this.colCnt);

            //the [1, 1] data
            theAns.data[0, 0] = Math.Sqrt(data[0, 0]);
            //the first row
            for (int colIdx = 1; colIdx < this.colCnt; colIdx++)
            {
                theAns[0, colIdx] = data[0, colIdx] / theAns[0, 0];
            }

            for (int rowIdx = 1; rowIdx < rowCnt; rowIdx++)
            {
                for (int colIdx = rowIdx; colIdx < this.colCnt; colIdx++)
                {
                    //Console.WriteLine($"{rowIdx}, {colIdx}");
                    theAns[rowIdx, colIdx] = this[rowIdx, colIdx];
                    if (rowIdx == colIdx)
                    {

                        for (int k = 0; k < colIdx; k++)
                        {
                            //Console.WriteLine(theAns[k, colIdx]);
                            theAns[rowIdx, colIdx] -= theAns[k, colIdx] * theAns[k, colIdx];
                        }
                        theAns[rowIdx, colIdx] = Math.Sqrt(theAns[rowIdx, colIdx]);
                    }
                    else
                    {
                        for (int k = 0; k < colIdx; k++)
                        {
                            theAns[rowIdx, colIdx] -= theAns[k, rowIdx] * theAns[k, colIdx];
                        }
                        theAns[rowIdx, colIdx] /= theAns[rowIdx, colIdx];
                    }
                }
            }

            return theAns;
        }
        #endregion

        #region statistics
        public Mtrx mean_throughRow()
        {
            Mtrx ans = new Mtrx(rowCnt, 1);
            Parallel.For(0, rowCnt, rowIdx =>
            {
                double rowMean = 0;
                for(int colIdx = 0; colIdx< colCnt; colIdx++)
                {
                    rowMean += data[rowIdx, colIdx];
                    //Thread.Sleep(1);
                }
                ans[rowIdx, 0] = rowMean / colCnt;
            });
            return ans;
        }

        public Mtrx mean_throughRow_sThread()
        {
            Mtrx ans = new Mtrx(rowCnt, 1);
            for (int rowIdx = 0; rowIdx < rowCnt; rowIdx++)
            {
                double rowMean = 0;
                for (int colIdx = 0; colIdx < colCnt; colIdx++)
                {
                    rowMean += data[rowIdx, colIdx];
                    Thread.Sleep(1);
                }
                ans[rowIdx, 0] = rowMean / colCnt;
            }
            return ans;
        }

        public Mtrx sampleCov_rowAsRV()
        {
            Mtrx ans = new Mtrx(rowCnt, rowCnt);
            Mtrx means = mean_throughRow();
            Parallel.For(0, rowCnt, rowIdx =>
            {
                Parallel.For(0, rowCnt, colIdx =>
                {
                    double xy = 0;
                    for (int k = 0; k < colCnt; k++)
                    {
                        xy += data[rowIdx, k] * data[colIdx, k];
                    }
                    xy /= colCnt;
                    ans[rowIdx, colIdx] = xy - means[rowIdx, 0] * means[colIdx, 0];
                });
            });
            return ans;
        }

        public Mtrx sampleCov_rowAsRV_sthread()
        {
            Mtrx ans = new Mtrx(rowCnt, rowCnt);
            Mtrx means = mean_throughRow();
            for (int rowIdx = 0; rowIdx < rowCnt; rowIdx++)
            {
                for (int colIdx = 0; colIdx < rowCnt; colIdx++)
                {
                    double xy = 0;
                    for (int k = 0; k < colCnt; k++)
                    {
                        xy += data[rowIdx, k] * data[colIdx, k];
                    }
                    xy /= colCnt;
                    ans[rowIdx, colIdx] = xy - means[rowIdx, 0] * means[colIdx, 0];
                }
            }
            return ans;
        }
        #endregion
    }
    class Zmtrx : Mtrx
    {
        private Random rv;

        #region constructor;
        public Zmtrx(int rowCnt, int colCnt, Random rv) : base(rowCnt, colCnt)
        {
            this.rv = rv;
            #region draw z
            Parallel.For(0, colCnt, colIdx =>
            {
                Parallel.For(0, rowCnt, rowIdx =>
                {
                    //&data;
                        double temp = DStat.N_Inv(rv.NextDouble());
                        data[rowIdx, colIdx] = temp;
                    //Thread.Sleep(1);
                });
            });

            /*
            for (int colIdx = 0; colIdx < colCnt; colIdx++)
            {
                for (int rowIdx = 0; rowIdx < rowCnt; rowIdx++)
                {
                    data[rowIdx, colIdx] = DStat.N_Inv(rv.NextDouble());
                }
            }
            */
            #endregion
        }

        public Zmtrx(int rowCnt, int colCnt) : base(rowCnt, colCnt)
        {
            this.rv = new Random();
            #region draw z
            Parallel.For(0, colCnt, colIdx =>
            {
                Parallel.For(0, rowCnt, rowIdx =>
                {
                    data[rowIdx, colIdx] = DStat.N_Inv(rv.NextDouble());
                    //Thread.Sleep(1);
                });
            });
            #endregion
        }

        public Zmtrx(Mtrx source) : base(source.getRowCnt(), source.getColCnt())
        {
            data = source.data;
        }
        #endregion

        #region algebra operation

        #endregion

        #region operation overload

        #endregion
    }

    class SVpath
    {
        private int length;
        private double[] vPath;
        private double[] sPath;
        private double T;
        public SVpath(ref Mtrx source,
            double s0, double v0,
            double kappa, double theta, double sigma,
            double rf, double T
        )
        {
            this.length = source.getColCnt();
            this.T = T;
            double deltaT = T / this.length;
            double sqrtdt = Math.Sqrt(deltaT);
            vPath = new double[this.length + 1];
            sPath = new double[this.length + 1];

            vPath[0] = v0;
            sPath[0] = s0;

            for (int t = 0; t < this.length; t++)
            {
                sPath[t + 1] = sPath[t] * Math.Exp((rf - 0.5 * vPath[t]) * deltaT + Math.Sqrt(vPath[t]) * sqrtdt * source[0, t]);
                vPath[t + 1] = Math.Max(vPath[t] + kappa * (theta - vPath[t]) * deltaT + sigma * Math.Sqrt(vPath[t]) * sqrtdt * source[1, t], 0);
            }
        }

        public double getSt()
        {
            return sPath[length];
        }

        public double[] getSpath()
        {
            return sPath;
        }

        public double[] getVPath()
        {
            return vPath;
        }

        public double[][] getSandVPath()
        {
            double[][] ans = new double[2][];
            ans[0] = sPath;
            ans[1] = vPath;
            return ans;
        }
    }
}
