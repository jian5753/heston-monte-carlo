using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace DFinNR
{
    public class normSample
    {
        private int pathCnt;
        private int length;
        private float mean;
        private float sigma;
        public Matrix data;

        public normSample(int pathCnt, int length, float mean, float sigma)
        {
            this.pathCnt = pathCnt;
            this.length = length;
            this.mean = mean;
            this.sigma = sigma;
        }

        public normSample(int _pathCnt, int _length)
        {
            this.pathCnt = _pathCnt;
            this.length = _length;
            this.mean = 0;
            this.sigma = 1;     
        }

        public string showPara()
        {
            return $"{this.pathCnt}, {this.length}, {this.mean}, {this.sigma}";
        }
        public void draw()
        {
            Random rv = new Random();
            this.data = new Matrix(this.pathCnt, this.length);
            for (int timeIdx = 0; timeIdx < this.length; timeIdx++)
            {
                for (int pathIdx = 0; pathIdx < this.pathCnt; pathIdx++)
                {
                    this.data[pathIdx, timeIdx] = DStat.N_Inv(rv.NextDouble());
                }
            }
        }
        
        public void inverseCholesky()
        {
            Matrix covMtrxCap = Matrix.sampleCovMtrx(this.data);
            Matrix dTri = Matrix.transpose(covMtrxCap.choleskyDecomp());
            Matrix invDtri = Matrix.inverse(dTri);
            this.data = invDtri * this.data;
        }
       
    }
}
