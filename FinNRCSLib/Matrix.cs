using System;
using System.Linq;
using System.Collections.Generic;

// Modified By Andy Dong

namespace DFinNR
{
    //! %Matrix used in linear algebra.
    /*! This class implements the concept of Matrix as used in linear
        algebra. As such, it is <b>not</b> meant to be used as a
        container.
    */
    public struct Matrix 
    {
        #region properties
        private int rows_, columns_;
        private double[] data_;
        public double this[int i] { get { return data_[i]; } set { data_[i] = value; } }

        public int rows() { return rows_; }
        public int columns() { return columns_; }
        public bool empty() { return rows_ == 0 || columns_ == 0; }

        public double this[int i, int j] 
        { get { return data_[i * columns_ + j]; } set { data_[i * columns_ + j] = value; } }
        
        public Vector row(int r) 
        {
            Vector result = new Vector(columns_);
            for (int i = 0; i < columns_; i++)
                result[i] = this[r, i];
            return result; 
        }
        
        public Vector column(int c) 
        {
            Vector result = new Vector(rows_);
            for (int i = 0; i < rows_; i++)
                result[i] = this[i, c];
            return result;
        }
        
        public Vector diagonal() 
        {
            int arraySize = Math.Min(rows(), columns());
            Vector tmp = new Vector(arraySize);
            for(int i = 0; i < arraySize; i++)
                tmp[i] = this[i, i];
            return tmp;
        }

        public double[] GetData()
        {
            return data_;
        }

        public Vector GetRange(int start, int length) 
        {
            return new Vector(data_.Skip(start).Take(length).ToList());
        }
        #endregion
        /* addd by K.L. Huang 20210505*/
        #region query
        public Matrix getRow(int rowIdx)
        {
            double[] result = new double[this.columns_];
            for(int i = 0; i < this.columns_; i++)
            {
                result[i] = this[rowIdx, i];
            }
            Matrix theRow = new Matrix(result, 1, this.columns_);
            return theRow;
        }

        public Matrix getCol(int colIdx)
        {
            double[] result = new double[this.rows_];
            for (int i = 0; i < this.rows_; i++)
            {
                result[i] = this[i, colIdx];
            }
            Matrix theCol = new Matrix(result, this.rows_, 1);
            return theCol;
        }
        #endregion


        #region Constructors
        //! creates a null matrix
        // public Matrix() : base(0) { rows_ = 0; columns_ = 0; }

        //! creates a matrix with the given dimensions
        public Matrix(int rows, int columns) 
        {
            data_ = new double[rows * columns];
            rows_ = rows;
            columns_ = columns;
        }

        //! creates the matrix and fills it with <tt>value</tt>
        public Matrix(int rows, int columns, double value) 
        {
            data_ = new double[rows * columns];
            for (int i = 0; i < data_.Length; i++)
                data_[i] = value;
            rows_ = rows;
            columns_ = columns;
        }

        public Matrix(Matrix from) 
        {
            data_ = !from.empty() ? (double[])from.data_.Clone() : null;
            rows_ = from.rows_;
            columns_ = from.columns_;
        }

        public Matrix(double[] elements, int row, int col)
        {
            data_ = new double[row * col];
            for (int i = 0; i < data_.Length; i++)
                data_[i] = elements[i];
            rows_ = row;
            columns_ = col;
        }
        
	    #endregion
    
        #region Algebraic operators
        
        /*! \pre all matrices involved in an algebraic expression must have the same size. */
        public static Matrix operator +(Matrix m1, Matrix m2) 
        { return operMatrix(ref m1, ref m2, (x, y) => x + y); }
        public static Matrix operator -(Matrix m1, Matrix m2) 
        { return operMatrix(ref m1, ref m2, (x, y) => x - y); }
        public static Matrix operator *(double value, Matrix m1) 
        { return operValue(ref m1, value, (x, y) => x * y); }
        public static Matrix operator /(double value, Matrix m1) 
        { return operValue(ref m1, value, (x, y) => x / y); }
        public static Matrix operator *(Matrix m1, double value) 
        { return operValue(ref m1, value, (x, y) => x * y); }
        public static Matrix operator /(Matrix m1, double value) 
        { return operValue(ref m1, value, (x, y) => x / y); }
        
        private static Matrix operMatrix(ref Matrix m1, ref Matrix m2, 
            Func<double, double, double> func) 
        {
            if (!(m1.rows_ == m2.rows_ && m1.columns_ == m2.columns_))
                throw new ApplicationException("operation on matrices with different sizes (" 
                    +m2.rows_+"x" + m2.columns_ + ", " + m1.rows_ + "x" + m1.columns_ + ")");

            Matrix result = new Matrix(m1.rows_, m1.columns_);
            for (int i = 0; i < m1.rows_; i++)
                for (int j = 0; j < m1.columns_; j++)
                    result[i, j] = func(m1[i, j], m2[i, j]);
            return result;
        }
        
        private static Matrix operValue(ref Matrix m1, double value, 
            Func<double, double, double> func) 
        {
            Matrix result = new Matrix(m1.rows_, m1.columns_);
            for (int i = 0; i < m1.rows_; i++)
                for (int j = 0; j < m1.columns_; j++)
                    result[i, j] = func(m1[i, j], value);
            return result;
        }

        public static Vector operator *(Vector v, Matrix m) 
        {
            if (!(v.Count == m.rows()))
                throw new ApplicationException("vectors and matrices with different sizes ("
                    +v.Count+", " + m.rows() + "x" + m.columns() + ") cannot be multiplied");
            Vector result = new Vector(m.columns());
            for (int i = 0; i < result.Count; i++)
                result[i] = v * m.column(i);
            return result;
        }
        
        /*! \relates Matrix */
        public static Vector operator *(Matrix m, Vector v) 
        {
            if (!(v.Count == m.columns()))
                throw new ApplicationException("vectors and matrices with different sizes ("
                    +v.Count+", " + m.rows() + "x" + m.columns() + ") cannot be multiplied");
            Vector result = new Vector(m.rows());
            for (int i = 0; i < result.Count; i++)
                result[i] = m.row(i) * v;
            return result;
        }
        
        /*! \relates Matrix */
        public static Matrix operator *(Matrix m1, Matrix m2) 
        {
            if (!(m1.columns() == m2.rows()))
                throw new ApplicationException("matrices with different sizes (" +
                       m1.rows() + "x" + m1.columns() + ", " +
                       m2.rows() + "x" + m2.columns() + ") cannot be multiplied");
            Matrix result = new Matrix(m1.rows(), m2.columns());
            for (int i = 0; i < result.rows(); i++)
                for (int j = 0; j < result.columns(); j++)
                    result[i, j] = m1.row(i) * m2.column(j);
            return result;
        }
        #endregion

        #region stats
        /* add by K.L. huang in 20210505*/
        public static Matrix sampleVar(Matrix m)
        {
            Matrix theAns = new Matrix(m.rows_, 1);
            for (int rowIdx = 0; rowIdx < m.rows_; rowIdx++)
            {
                double squaredMean = 0;
                double mean = 0;
                for(int colIdx = 0; colIdx < m.columns_; colIdx++)
                {
                    squaredMean += Math.Pow(m[rowIdx, colIdx], 2);
                    mean += m[rowIdx, colIdx];
                }

                squaredMean /= m.columns_;
                mean /= m.columns_;
                theAns[rowIdx, 0] = squaredMean - Math.Pow(mean, 2);
            }
            return theAns;
        }
        
        public static double sampleMean(Matrix m1)
        {
            double mean = 0;
            for(int colIdx = 0; colIdx < m1.columns_; colIdx++)
            {
                mean += m1[0, colIdx];
            }
            return mean / m1.columns_;
        }
        public static double sampleCov(Matrix m1, Matrix m2)
        {
            double xy = 0;
            for(int colIdx = 0; colIdx < m1.columns_; colIdx++)
            {
                xy += m1[0, colIdx] * m2[0, colIdx];
            }
            return (xy/m1.columns_) - sampleMean(m1) * sampleMean(m2);
        }
        public static Matrix sampleCovMtrx(Matrix m)
        {
            Matrix theAns = new Matrix(m.rows_, m.rows_);
            for (int rowIdx_cov = 0; rowIdx_cov < m.rows_; rowIdx_cov++)
            {
                for (int colIdx_cov = 0; colIdx_cov < m.rows_; colIdx_cov++)
                {
                    theAns[rowIdx_cov, colIdx_cov] = sampleCov(
                        m.getRow(rowIdx_cov),
                        m.getRow(colIdx_cov));
                }
            }
            return theAns;
        }
        #endregion

        public static Matrix transpose(Matrix m) 
        {
            Matrix result = new Matrix(m.columns(),m.rows());
            for (int i=0; i<m.rows(); i++)
                for (int j=0; j<m.columns();j++)
                    result[j,i] = m[i,j];
            return result;
        }

        public static Matrix outerProduct(List<double> v1begin, List<double> v2begin) 
        {
            int size1 = v1begin.Count;
            if (!(size1>0)) throw new ApplicationException("null first vector");

            int size2 = v2begin.Count;
            if(!(size2>0)) throw new ApplicationException("null second vector");

            Matrix result = new Matrix(size1, size2);

            for (int i=0; i<v1begin.Count; i++)
                for(int j=0; j<v2begin.Count; j++)
                    result[i,j] = v1begin[i] * v2begin[j];
            return result;
        }

        /*! \relates Matrix */
        public static Matrix inverse(Matrix m)
        {
            int numColumns = m.columns_;

            double[] elements = new double[m.rows_ * m.columns_];
            for (int s = 0; s < m.rows_; s++)
                for (int t = 0; t < m.columns_; t++)
                    elements[t + s * m.columns_] = m[s, t];

            int i, j, k, l, u, v;            
            double d = 0, p = 0;

            int[] pnRow = new int[numColumns];
            int[] pnCol = new int[numColumns];

            for (k = 0; k <= numColumns - 1; k++)
            {
                d = 0.0;
                for (i = k; i <= numColumns - 1; i++)
                {
                    for (j = k; j <= numColumns - 1; j++)
                    {
                        l = i * numColumns + j; p = Math.Abs(elements[l]);
                        if (p > d)
                        {
                            d = p;
                            pnRow[k] = i;
                            pnCol[k] = j;
                        }
                    }
                }
                
                if (d == 0.0)
                {
                    return m;
                }

                if (pnRow[k] != k)
                {
                    for (j = 0; j <= numColumns - 1; j++)
                    {
                        u = k * numColumns + j;
                        v = pnRow[k] * numColumns + j;
                        p = elements[u];
                        elements[u] = elements[v];
                        elements[v] = p;
                    }
                }

                if (pnCol[k] != k)
                {
                    for (i = 0; i <= numColumns - 1; i++)
                    {
                        u = i * numColumns + k;
                        v = i * numColumns + pnCol[k];
                        p = elements[u];
                        elements[u] = elements[v];
                        elements[v] = p;
                    }
                }

                l = k * numColumns + k;
                elements[l] = 1.0 / elements[l];
                for (j = 0; j <= numColumns - 1; j++)
                {
                    if (j != k)
                    {
                        u = k * numColumns + j;
                        elements[u] = elements[u] * elements[l];
                    }
                }

                for (i = 0; i <= numColumns - 1; i++)
                {
                    if (i != k)
                    {
                        for (j = 0; j <= numColumns - 1; j++)
                        {
                            if (j != k)
                            {
                                u = i * numColumns + j;
                                elements[u] = elements[u] - elements[i * numColumns + k] 
                                    * elements[k * numColumns + j];
                            }
                        }
                    }
                }

                for (i = 0; i <= numColumns - 1; i++)
                {
                    if (i != k)
                    {
                        u = i * numColumns + k;
                        elements[u] = -elements[u] * elements[l];
                    }
                }
            }
            
            for (k = numColumns - 1; k >= 0; k--)
            {
                if (pnCol[k] != k)
                {
                    for (j = 0; j <= numColumns - 1; j++)
                    {
                        u = k * numColumns + j;
                        v = pnCol[k] * numColumns + j;
                        p = elements[u];
                        elements[u] = elements[v];
                        elements[v] = p;
                    }
                }

                if (pnRow[k] != k)
                {
                    for (i = 0; i <= numColumns - 1; i++)
                    {
                        u = i * numColumns + k;
                        v = i * numColumns + pnRow[k];
                        p = elements[u];
                        elements[u] = elements[v];
                        elements[v] = p;
                    }
                }
            }

            return new Matrix(elements, m.rows_, m.columns_);
        }


        /*! \relates Matrix */
        public static double determinant(Matrix m)
        {
            int numColumns = m.columns_;
            double[] elements = m.GetData();

            int i, j, k, nis = 0, js = 0, l, u, v;
            double f, det, q, d;

            f = 1.0;
            det = 1.0;

            for (k = 0; k <= numColumns - 2; k++)
            {
                q = 0.0;
                for (i = k; i <= numColumns - 1; i++)
                {
                    for (j = k; j <= numColumns - 1; j++)
                    {
                        l = i * numColumns + j;
                        d = Math.Abs(elements[l]);
                        if (d > q)
                        {
                            q = d;
                            nis = i;
                            js = j;
                        }
                    }
                }

                if (q == 0.0)
                {
                    det = 0.0;
                    return (det);
                }

                if (nis != k)
                {
                    f = -f;
                    for (j = k; j <= numColumns - 1; j++)
                    {
                        u = k * numColumns + j;
                        v = nis * numColumns + j;
                        d = elements[u];
                        elements[u] = elements[v];
                        elements[v] = d;
                    }
                }

                if (js != k)
                {
                    f = -f;
                    for (i = k; i <= numColumns - 1; i++)
                    {
                        u = i * numColumns + js;
                        v = i * numColumns + k;
                        d = elements[u];
                        elements[u] = elements[v];
                        elements[v] = d;
                    }
                }

                l = k * numColumns + k;
                det = det * elements[l];
                for (i = k + 1; i <= numColumns - 1; i++)
                {
                    d = elements[i * numColumns + k] / elements[l];
                    for (j = k + 1; j <= numColumns - 1; j++)
                    {
                        u = i * numColumns + j;
                        elements[u] = elements[u] - d * elements[k * numColumns + j];
                    }
                }
            }

            det = f * det * elements[numColumns * numColumns - 1];

            return (det);
        }

        public void fill(double value) 
        {
            for (int i = 0; i < data_.Length; i++)
                data_[i] = value;
        }

        public void swap(int i1, int j1, int i2, int j2) 
        {
            double t = this[i2, j2];
            this[i2, j2] = this[i1, j1];
            this[i1, j1] = t;
        }
        
        public override string ToString()
        {
            string str = "";
            str = str + "Matrix size: " + this.rows().ToString() + " * " 
                + this.columns().ToString() + ", Value: \n";
            for (int i = 0; i < this.rows(); i++)
            {
                for(int j = 0; j< this.columns(); j++)
                {
                    str = str + this[i * columns_ + j].ToString("F4") + " ";
                }
                str = str + "\n";
            }
            return str;
        }

        /* add by K.L. Huang 20210505 */
        static public Matrix corrToCov(Matrix corrMtrx, double[] variance)
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

        public Matrix choleskyDecomp()
        {
            if (this.rows() != this.columns())
            {
                Console.WriteLine("not square");
                int theShort = Math.Min(this.rows_, this.columns_);
                return new Matrix(theShort, theShort);
            }
            Matrix theAns = new Matrix(this.rows_, this.columns_);

            //the [1, 1] data
            theAns.data_[0] = Math.Sqrt(this.data_[0]);
            //the first row
            for (int colIdx = 1; colIdx < this.columns_; colIdx++)
            {
                theAns[0, colIdx] = this.data_[colIdx] / theAns[0, 0];
            }

            for (int rowIdx = 1; rowIdx < this.rows_; rowIdx++)
            {
                for(int colIdx = rowIdx; colIdx < this.columns_; colIdx++)
                {
                    Console.WriteLine($"{rowIdx}, {colIdx}");
                    theAns[rowIdx, colIdx] = this[rowIdx, colIdx];
                    if (rowIdx == colIdx)
                    {
                        
                        for(int k = 0; k < colIdx; k++)
                        {
                            //Console.WriteLine(theAns[k, colIdx]);
                            theAns[rowIdx, colIdx] -= theAns[k, colIdx] * theAns[k, colIdx];
                        }
                        theAns[rowIdx, colIdx] = Math.Sqrt(theAns[rowIdx, colIdx]);
                    }else
                    {
                        for(int k =0; k < colIdx; k++)
                        {
                            theAns[rowIdx, colIdx] -= theAns[k, rowIdx] * theAns[k, colIdx];
                        }
                        theAns[rowIdx, colIdx] /=theAns[rowIdx, colIdx];
                    }
                }
            }

            return theAns;
        }
    }
}
