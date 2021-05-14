using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace DFinNR
{
    public static class GeneralBSMOption
    {
        // Extensions of the Black Scholes model //////////////
        public static double option_price_european_call_payout(double S,
            double K, double r, double q, double sigma, double time)
        {
            double sigma_sqr = Math.Pow(sigma, 2);
            double time_sqrt = Math.Sqrt(time);
            double d1 = (Math.Log(S / K) + (r - q + 0.5 * sigma_sqr) * time)
                / (sigma * time_sqrt);
            double d2 = d1 - (sigma * time_sqrt);
            double call_price = S * Math.Exp(-q * time) * DStat.NormDist(d1)
                - K * Math.Exp(-r * time) * DStat.NormDist(d2);

            return call_price;
        }

        public static double option_price_european_put_payout(double S, 
            double K, double r, double q, double sigma, double time)
        {
            double sigma_sqr = Math.Pow(sigma, 2);
            double time_sqrt = Math.Sqrt(time);
            double d1 = (Math.Log(S / K) + (r - q + 0.5 * sigma_sqr) * time) 
                / (sigma * time_sqrt);
            double d2 = d1 - (sigma * time_sqrt);
            double put_price = K * Math.Exp(-r * time) * DStat.NormDist(-d2) 
                - S * Math.Exp(-q * time) * DStat.NormDist(-d1);

            return put_price;
        }

        public static double option_price_european_call_dividends(double S, 
            double K, double r, double sigma, double time_to_maturity,
            Vector dividend_times, Vector dividend_amounts)
        {
            double adjusted_S = S;

            for (int i = 0; i < dividend_times.size(); i++)
            {
                if (dividend_times[i] <= time_to_maturity)
                {
                    adjusted_S -= dividend_amounts[i] 
                        * Math.Exp(-r * dividend_times[i]);
                };
            };

            return BSMOption.option_price_call_black_scholes(adjusted_S,
                K, r, sigma, time_to_maturity);
        }

        public static double option_price_european_put_dividends(double S,
            double K, double r, double sigma, double time_to_maturity,
            Vector dividend_times, Vector dividend_amounts)
        {
            // reduce the current stock price by the amount of dividends. 
            double adjusted_S = S;
            for (int i = 0; i < dividend_times.size(); i++)
            {
                if (dividend_times[i] <= time_to_maturity)
                {
                    adjusted_S -= dividend_amounts[i] 
                        * Math.Exp(-r * dividend_times[i]);
                };
            };

            return BSMOption.option_price_put_black_scholes(
                adjusted_S, K, r, sigma, time_to_maturity);
        }

        public static double option_price_american_call_one_dividend(double S, 
            double K, double r, double sigma, double tau, double D1, double tau1)
        {
            if (D1 <= K * (1.0 - Math.Exp(-r * (tau - tau1)))) // check for no exercise
                return BSMOption.option_price_call_black_scholes(
                    S - Math.Exp(-r * tau1) * D1, K, r, sigma, tau);

             double ACCURACY = 1e-6;  // decrease this for more accuracy
            double sigma_sqr = sigma * sigma;
            double tau_sqrt = Math.Sqrt(tau);
            double tau1_sqrt = Math.Sqrt(tau1);
            double rho = -Math.Sqrt(tau1 / tau);

            double S_bar = 0;  // first find the S_bar that solves c=S_bar+D1-K 
            double S_low = 0;    // the simplest: binomial search
            double S_high = S;  // start by finding a very high S above S_bar

            double c = BSMOption.option_price_call_black_scholes(
                S_high, K, r, sigma, tau - tau1);

            double test = c - S_high - D1 + K;
            while ((test > 0.0) && (S_high <= 1e10))
            {
                S_high *= 2.0;
                c = BSMOption.option_price_call_black_scholes(
                    S_high, K, r, sigma, tau - tau1);
                test = c - S_high - D1 + K;
            };

            if (S_high > 1e10)
            { 
                // early exercise never optimal, find BS value
                return BSMOption.option_price_call_black_scholes(
                    S - D1 * Math.Exp(-r * tau1), K, r, sigma, tau);
            };

            S_bar = 0.5 * S_high;  // now find S_bar that solves c=S_bar-D+K
            c = BSMOption.option_price_call_black_scholes(
                S_bar, K, r, sigma, tau - tau1);

            test = c - S_bar - D1 + K;
            while ((Math.Abs(test) > ACCURACY) && ((S_high - S_low) > ACCURACY))
            {
                if (test < 0.0) { S_high = S_bar; }
                else { S_low = S_bar; };
                S_bar = 0.5 * (S_high + S_low);
                c = BSMOption.option_price_call_black_scholes(S_bar, K, r, sigma, tau - tau1);
                test = c - S_bar - D1 + K;
            };

            double a1 = (Math.Log((S - D1 * Math.Exp(-r * tau1)) / K) 
                + (r + 0.5 * sigma_sqr) * tau) / (sigma * tau_sqrt);
            double a2 = a1 - sigma * tau_sqrt;
            double b1 = (Math.Log((S - D1 * Math.Exp(-r * tau1)) / S_bar)
                + (r + 0.5 * sigma_sqr) * tau1) / (sigma * tau1_sqrt);
            double b2 = b1 - sigma * tau1_sqrt;
            double C = (S - D1 * Math.Exp(-r * tau1)) * DStat.NormDist(b1) 
                + (S - D1 * Math.Exp(-r * tau1)) * DStat.NormDist(a1, -b1, rho)
               - (K * Math.Exp(-r * tau)) * DStat.NormDist(a2, -b2, rho) 
               - (K - D1) * Math.Exp(-r * tau1) * DStat.NormDist(b2);

            return C;
        }

        public static double futures_option_price_call_european_black(double F, 
            double K, double r, double sigma, double time)
        {
            double sigma_sqr = sigma * sigma;
            double time_sqrt = Math.Sqrt(time);
            double d1 = (Math.Log(F / K) + 0.5 * sigma_sqr * time) / (sigma * time_sqrt);
            double d2 = d1 - sigma * time_sqrt;

            return Math.Exp(-r * time) * (F * DStat.NormDist(d1) 
                - K * DStat.NormDist(d2));
        }

        public static double futures_option_price_put_european_black(double F, 
            double K, double r, double sigma, double time)
        {
            double sigma_sqr = sigma * sigma;
            double time_sqrt = Math.Sqrt(time);
            double d1 = (Math.Log(F / K) + 0.5 * sigma_sqr * time) / (sigma * time_sqrt);
            double d2 = d1 - sigma * time_sqrt;

            return Math.Exp(-r * time) * (K * DStat.NormDist(-d2) 
                - F * DStat.NormDist(-d1));
        }

        public static double currency_option_price_call_european(double S, 
            double K, double r, double r_f, double sigma, double time)
        {
            double sigma_sqr = sigma * sigma;
            double time_sqrt = Math.Sqrt(time);
            double d1 = (Math.Log(S / K) + (r - r_f + (0.5 * sigma_sqr)) * time)
                / (sigma * time_sqrt);
            double d2 = d1 - sigma * time_sqrt;

            return S * Math.Exp(-r_f * time) * DStat.NormDist(d1) 
                - K * Math.Exp(-r * time) * DStat.NormDist(d2);
        }

        public static double currency_option_price_put_european(double S, 
            double K, double r, double r_f, double sigma, double time)
        {
            double sigma_sqr = sigma * sigma;
            double time_sqrt = Math.Sqrt(time);
            double d1 = (Math.Log(S / K) + (r - r_f + (0.5 * sigma_sqr)) * time)
                / (sigma * time_sqrt);
            double d2 = d1 - sigma * time_sqrt;

            return K * Math.Exp(-r * time) * DStat.NormDist(-d2)
                - S * Math.Exp(-r_f * time) * DStat.NormDist(-d1);
        }

        public static double option_price_american_perpetual_call(double S,
            double K, double r, double q, double sigma)
        {
            double sigma_sqr = Math.Pow(sigma, 2);
            double h1 = 0.5 - ((r - q) / sigma_sqr);
            h1 += Math.Sqrt(Math.Pow(((r - q) / sigma_sqr - 0.5), 2) 
                + 2.0 * r / sigma_sqr);
            double pric = (K / (h1 - 1.0)) * Math.Pow(((h1 - 1.0) / h1) 
                * (S / K), h1);

            return pric;
        }

        public static double option_price_american_perpetual_put(double S, 
            double K, double r, double q, double sigma)
        {
            double sigma_sqr = Math.Pow(sigma, 2);
            double h2 = 0.5 - ((r - q) / sigma_sqr) 
                - Math.Sqrt(Math.Pow(((r - q) / sigma_sqr - 0.5), 2) 
                + 2.0 * r / sigma_sqr);

            double pric = (K / (1.0 - h2)) * Math.Pow(((h2 - 1.0) / h2) 
                * (S / K), h2);

            return pric;
        }
    }

        /*    Heston OPtion Price
        //////////////////// Heston Function ////////////////////////////////
         
        public struct heston_parms
        {
            public double K, x, r, v, tau;
            public double kappa, theta, rho, sigma, lambda;
            public int j;
        };

        public static double heston_integrand_j(double phi, heston_parms p)
        {            
            heston_parms parms = p;
            double K = (parms.K);
            double x = (parms.x);
            double v = (parms.v);
            double r = (parms.r);

            double kappa = (parms.kappa);
            double theta = (parms.theta);
            double rho = (parms.rho);
            double sigma = (parms.sigma);
            double lambda = (parms.lambda);
            double tau = (parms.tau);

            double j = (parms.j);

            double sigma_sqr = Math.Pow(sigma, 2);
            double uj;
            double bj;

            if (j == 1)
            {
                uj = 0.5;
                bj = kappa + lambda - rho * sigma;
            }
            else
            {
                uj = -0.5;
                bj = kappa + lambda;
            };

            System.Numerics.Complex i = new System.Numerics.Complex(0, 1);
            double a = kappa * theta;

            System.Numerics.Complex d, g, C, D, f1, F = new System.Numerics.Complex();
            d = System.Numerics.Complex.Sqrt(System.Numerics.Complex.Pow(
                rho * sigma * phi * i - bj, 2) 
                - sigma_sqr * (2 * uj * phi * i - Math.Pow(phi, 2)));

            g = (bj - rho * sigma * phi * i + d) / (bj - rho * sigma * phi * i - d);
            C = r * phi * i * tau + (a / sigma_sqr) * ((bj - rho * sigma * phi * i + d) * tau
                - 2.0 * System.Numerics.Complex.Log((1.0 - g * System.Numerics.Complex.Exp(d * tau)) 
                / (1.0 - g)));

            D = (bj - rho * sigma * phi * i + d) / sigma_sqr 
                * ((1.0 - System.Numerics.Complex.Exp(d * tau)) 
                / (1.0 - g * System.Numerics.Complex.Exp(d * tau)));

            f1 = System.Numerics.Complex.Exp(C + D * v + i * phi * x);
            F = System.Numerics.Complex.Exp(-phi * i * Math.Log(K)) * f1 / (i * phi);

            return F.Real;
        }


        public static double heston_call_option_price(double S, double K, double r, double v, double tau,
                 double rho, double kappa, double lambda, double theta, double sigma)
        {
            return 0;
        }
    
        */

}
