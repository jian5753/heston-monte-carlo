using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace DFinNR
{
    public static class FiniteDifference
    {        
        //////////////////// finite differences //////////////////
        public static double option_price_call_european_finite_diff_explicit(
            double S, double K, double r, double sigma, double time,
            int no_S_steps, int no_t_steps)
        {
            double sigma_sqr = Math.Pow(sigma, 2);

            int M = no_S_steps;
            if ((no_S_steps % 2) == 1)
            {
                ++M;
            } // need no_S_steps to be even:

            double delta_S = 2.0 * S / M;
            Vector S_values = new Vector(M + 1);
            for (int m = 0; m <= M; m++)
            {
                S_values[m] = m * delta_S;
            };

            int N = no_t_steps;
            double delta_t = time / N;

            Vector a = new Vector(M);
            Vector b = new Vector(M);
            Vector c = new Vector(M);

            double r1 = 1.0 / (1.0 + r * delta_t);
            double r2 = delta_t / (1.0 + r * delta_t);

            for ( int j = 1; j < M; j++)
            {
                a[j] = r2 * 0.5 * j * (-r + sigma_sqr * j);
                b[j] = r1 * (1.0 - sigma_sqr * j * j * delta_t);
                c[j] = r2 * 0.5 * j * (r + sigma_sqr * j);
            };

            Vector f_next = new Vector(M + 1);
            for (int m = 0; m <= M; ++m)
            {
                f_next[m] = Math.Max(0.0, S_values[m] - K);
            };

            Vector f = new Vector(M + 1);
            for (int t = N - 1; t >= 0; --t)
            {
                f[0] = 0;
                for (int m = 1; m < M; ++m)
                {
                    f[m] = a[m] * f_next[m - 1] + b[m] * f_next[m] + c[m] * f_next[m + 1];
                };

                f[M] = 0;
                for (int m = 0; m <= M; ++m)
                {
                    f_next[m] = f[m];
                };
            };

            return f[M / 2];
        }

        public static double option_price_put_european_finite_diff_explicit(
            double S, double X, double r, double sigma, double time, 
            int no_S_steps, int no_t_steps)
        {
            double sigma_sqr = Math.Pow(sigma, 2);

            int M = no_S_steps;
            if ((no_S_steps % 2) == 1)
            {
                ++M;
            }
            // need no_S_steps to be even:

            double delta_S = 2.0 * S / M;

            Vector S_values = new Vector(M + 1);
            for (int m = 0; m <= M; m++)
            {
                S_values[m] = m * delta_S;
            };

            int N = no_t_steps;
            double delta_t = time / N;

            Vector a = new Vector(M);
            Vector b = new Vector(M);
            Vector c = new Vector(M);

            double r1 = 1.0 / (1.0 + r * delta_t);
            double r2 = delta_t / (1.0 + r * delta_t);
            for (int j = 1; j < M; j++)
            {
                a[j] = r2 * 0.5 * j * (-r + sigma_sqr * j);
                b[j] = r1 * (1.0 - sigma_sqr * j * j * delta_t);
                c[j] = r2 * 0.5 * j * (r + sigma_sqr * j);
            };

            Vector f_next = new Vector(M + 1);
            for (int m = 0; m <= M; ++m)
            {
                f_next[m] = Math.Max(0.0, X - S_values[m]);
            };

            Vector f = new Vector(M + 1);
            for (int t = N - 1; t >= 0; --t)
            {
                f[0] = X;
                for (int m = 1; m < M; ++m)
                {
                    f[m] = a[m] * f_next[m - 1] + b[m] * f_next[m] + c[m] * f_next[m + 1];
                };

                f[M] = 0;
                for (int m = 0; m <= M; ++m)
                {
                    f_next[m] = f[m];
                };
            };

            return f[M / 2];
        }

        public static double option_price_call_american_finite_diff_explicit(
            double S, double K, double r, double sigma, double time,
            int no_S_steps, int no_t_steps)
        {
            double sigma_sqr = Math.Pow(sigma, 2);

            int M = no_S_steps + (no_S_steps % 2);
            // need no_S_steps to be even:

            double delta_S = 2.0 * S / M;
            Vector S_values = new Vector(M + 1, 0.0);
            for (int m = 1; m <= M; m++)
            {
                S_values[m] = m * delta_S;
            };

            int N = no_t_steps;
            double delta_t = time / N;

            Vector a = new Vector(M,0.0);
            Vector b = new Vector(M,0.0);
            Vector c = new Vector(M,0.0);

            double r1 = 1.0 / (1.0 + r * delta_t);
            double r2 = delta_t / (1.0 + r * delta_t);

            for (int j = 1; j < M; j++)
            {
                a[j] = r2 * 0.5 * j * (-r + sigma_sqr * j);
                b[j] = r1 * (1.0 - sigma_sqr * j * j * delta_t);
                c[j] = r2 * 0.5 * j * (r + sigma_sqr * j);
            };

            Vector f_next = new Vector(M + 1, 0.0);
            for (int m = 0; m <= M; ++m)
            {
                f_next[m] = Math.Max(0.0, S_values[m] - K);
            };

            Vector f = new Vector(M + 1, 0.0);
            for (int t = N - 1; t >= 0; --t)
            {
                f[0] = 0;
                for (int m = 1; m < M; ++m)
                {
                    f[m] = a[m] * f_next[m - 1] + b[m] * f_next[m] + c[m] * f_next[m + 1];
                    f[m] = Math.Max(f[m], S_values[m] - K);  // check for exercise
                };

                f[M] = S_values[M] - K;
                for (int m = 0; m <= M; ++m)
                {
                    f_next[m] = f[m];
                };
            };
            double C = f[M / 2];
            return C;

        }

        public static double option_price_put_american_finite_diff_explicit(
            double S, double K, double r, double sigma, double time,
            int no_S_steps, int no_t_steps)
        {
            double sigma_sqr = sigma * sigma;

            int M = no_S_steps + (no_S_steps % 2);
            // need no_S_steps to be even:

            double delta_S = 2.0 * S / M;

            Vector S_values = new Vector(M + 1);
            for (int m = 0; m <= M; m++)
            {
                S_values[m] = m * delta_S;
            };

            int N = no_t_steps;
            double delta_t = time / N;

            Vector a = new Vector(M);
            Vector b = new Vector(M);
            Vector c = new Vector(M);

            double r1 = 1.0 / (1.0 + r * delta_t);
            double r2 = delta_t / (1.0 + r * delta_t);
            for (int j = 1; j < M; j++)
            {
                a[j] = r2 * 0.5 * j * (-r + sigma_sqr * j);
                b[j] = r1 * (1.0 - sigma_sqr * j * j * delta_t);
                c[j] = r2 * 0.5 * j * (r + sigma_sqr * j);
            };

            Vector f_next = new Vector(M + 1);
            for (int m = 0; m <= M; ++m)
            {
                f_next[m] = Math.Max(0.0, K - S_values[m]);
            };

            Vector f = new Vector(M + 1);
            for (int t = N - 1; t >= 0; --t)
            {
                f[0] = K;
                for (int m = 1; m < M; ++m)
                {
                    f[m] = a[m] * f_next[m - 1] + b[m] * f_next[m] + c[m] * f_next[m + 1];
                    f[m] = Math.Max(f[m], K - S_values[m]);  // check for exercise
                };

                f[M] = 0;
                for (int m = 0; m <= M; ++m)
                {
                    f_next[m] = f[m];
                };
            };

            return f[M / 2];
        }

        public static double option_price_call_european_finite_diff_implicit(
            double S, double K, double r, double sigma, double time,
            int no_S_steps, int no_t_steps)
        {
            double sigma_sqr = sigma * sigma;

            int M = no_S_steps + (no_S_steps % 2);
            // need no_S_steps to be even:

            double delta_S = 2.0 * S / M;
            Vector S_values = new Vector(M + 1, 0.0);
            for (int m = 0; m <= M; m++)
            {
                S_values[m] = m * delta_S;
            };
            int N = no_t_steps;
            double delta_t = time / N;

            Matrix A = new Matrix(M + 1, M + 1, 0.0);
            //A = 0.0;

            A[0, 0] = 1.0;
            for (int j = 1; j < M; ++j)
            {
                A[j, j - 1] = 0.5 * j * delta_t * (r - sigma_sqr * j);    // a[j]
                A[j, j] = 1.0 + delta_t * (r + sigma_sqr * j * j);        // b[j];
                A[j, j + 1] = 0.5 * j * delta_t * (-r - sigma_sqr * j);   // c[j];
            };
            A[M, M] = 1.0;

            Vector B = new Vector(M + 1);
            for (int m = 0; m <= M; ++m)
            {
                B[m] = Math.Max(0.0, S_values[m] - K);
            };

            Vector F = new Vector(Matrix.inverse(A) * B);
            for (int t = N - 1; t > 0; --t)
            {
                B = F;
                F = Matrix.inverse(A) * B;
            };

            return F[M / 2];
        }

        public static double option_price_put_european_finite_diff_implicit(
            double S, double K, double r, double sigma, double time,
            int no_S_steps, int no_t_steps)
        {
            double sigma_sqr = sigma * sigma;
            int M = no_S_steps + (no_S_steps % 2); 
            // need no_S_steps to be even:

            double delta_S = 2.0 * S / M;
            Vector S_values = new Vector(M + 1);
            for (int m = 0; m <= M; m++)
            {
                S_values[m] = m * delta_S;
            };

            int N = no_t_steps;
            double delta_t = time / N;

            Matrix A = new Matrix(M + 1, M + 1, 0.0); //A = 0.0;
            A[0, 0] = 1.0;
            for (int j = 1; j < M; ++j)
            {
                A[j, j - 1] = 0.5 * j * delta_t * (r - sigma_sqr * j);    // a[j]
                A[j, j] = 1.0 + delta_t * (r + sigma_sqr * j * j);        // b[j];
                A[j, j + 1] = 0.5 * j * delta_t * (-r - sigma_sqr * j);   // c[j];
            };
            A[M, M] = 1.0;

            Vector B = new Vector(M + 1);
            for (int m = 0; m <= M; ++m)
            {
                B[m] = Math.Max(0.0, K - S_values[m]);
            };

            Vector F = new Vector(Matrix.inverse(A) * B);
            for (int t = N - 1; t > 0; --t)
            {
                B = F;
                F = Matrix.inverse(A) * B;
            };

            return F[M / 2];
        }

        public static double option_price_call_american_finite_diff_implicit(
            double S, double K, double r, double sigma, double time,
            int no_S_steps, int no_t_steps)
        {
            double sigma_sqr = Math.Pow(sigma, 2);
            int M = no_S_steps + (no_S_steps % 2);
            // need no_S_steps to be even:

            double delta_S = 2.0 * S / M;

            Vector S_values = new Vector(M + 1);
            for (int m = 0; m <= M; m++)
            {
                S_values[m] = m * delta_S;
            };

            int N = no_t_steps;
            double delta_t = time / N;

            Matrix A = new Matrix(M + 1, M + 1, 0.0); //A = 0.0;
            A[0, 0] = 1.0;
            for (int j = 1; j < M; ++j)
            {
                A[j, j - 1] = 0.5 * j * delta_t * (r - sigma_sqr * j);    // a[j]
                A[j, j] = 1.0 + delta_t * (r + sigma_sqr * j * j);  // b[j];
                A[j, j + 1] = 0.5 * j * delta_t * (-r - sigma_sqr * j);   // c[j];
            };
            A[M, M] = 1.0;

            Vector B = new Vector(M + 1);
            for (int m = 0; m <= M; ++m)
            {
                B[m] = Math.Max(0.0, S_values[m] - K);
            };

            Vector F = new Vector(Matrix.inverse(A) * B);
            for (int t = N - 1; t > 0; --t)
            {
                B = F;
                F = Matrix.inverse(A) * B;
                for (int m = 1; m < M; ++m)
                {   // now check for exercise
                    F[m] = Math.Max(F[m], S_values[m] - K);
                };
            };

            return F[M / 2];
        }

        public static double option_price_put_american_finite_diff_implicit(
            double S, double K, double r, double sigma, double time,
            int no_S_steps, int no_t_steps)
        {
            double sigma_sqr = sigma * sigma;
            int M = no_S_steps + (no_S_steps % 2);
            // need no_S_steps to be even:

            double delta_S = 2.0 * S / M;
            Vector S_values = new Vector(M + 1);
            for (int m = 0; m <= M; m++) { S_values[m] = m * delta_S; };
            int N = no_t_steps;
            double delta_t = time / N;

            Matrix A = new Matrix(M + 1, M+1, 0.0); //A = 0.0;
            A[0, 0] = 1.0;
            for (int j = 1; j < M; ++j)
            {
                A[j, j - 1] = 0.5 * j * delta_t * (r - sigma_sqr * j);    // a[j]
                A[j, j] = 1.0 + delta_t * (r + sigma_sqr * j * j);        // b[j];
                A[j, j + 1] = 0.5 * j * delta_t * (-r - sigma_sqr * j);   // c[j];
            };
            A[M, M] = 1.0;

            Vector B = new Vector(M + 1);
            for (int m = 0; m <= M; ++m)
            {
                B[m] = Math.Max(0.0, K - S_values[m]);
            };

            Vector F = new Vector(Matrix.inverse(A) * B);
            for (int t = N - 1; t > 0; --t)
            {
                B = F;
                F = Matrix.inverse(A) * B;
                for (int m = 1; m < M; ++m)
                {   // now check for exercise
                    F[m] = Math.Max(F[m], K - S_values[m]);
                };
            };

            return F[M / 2];
        }
    }

}
