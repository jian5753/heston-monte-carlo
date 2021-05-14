using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace DFinNR
{   
    public static class BSMOption
    {
        /// Black Scholes formula //////////////////////////////////////////
        public static double option_price_call_black_scholes(double S, double K,
            double r, double sigma, double time)
        {
            double time_sqrt = Math.Sqrt(time);
            double d1 = (Math.Log(S / K) + r * time) / (sigma * time_sqrt) + 0.5 * sigma * time_sqrt;
            double d2 = d1 - (sigma * time_sqrt);
            return S * DStat.NormDist(d1) - K * Math.Exp(-r * time) * DStat.NormDist(d2);
        }

        public static double option_price_put_black_scholes(double S, double K,
            double r, double sigma, double time)
        {
            double time_sqrt = Math.Sqrt(time);
            double d1 = (Math.Log(S / K) + r * time) / (sigma * time_sqrt) + 0.5 * sigma * time_sqrt;
            double d2 = d1 - (sigma * time_sqrt);
            return K * Math.Exp(-r * time) * DStat.NormDist(-d2) - S * DStat.NormDist(-d1);
        }

        public static double option_price_delta_call_black_scholes(double S, 
            double K, double r, double sigma, double time)
        {
            double time_sqrt = Math.Sqrt(time);
            double d1 = (Math.Log(S / K) + r * time) / (sigma * time_sqrt) 
                + 0.5 * sigma * time_sqrt;
            double delta = DStat.NormDist(d1);
            return delta;
        }

        public static double option_price_delta_put_black_scholes(double S, 
            double K, double r, double sigma, double time)
        {
            double time_sqrt = Math.Sqrt(time);
            double d1 = (Math.Log(S / K) + r * time) / (sigma * time_sqrt) 
                + 0.5 * sigma * time_sqrt;
            double delta = -DStat.NormDist(-d1);
            return delta;
        }

        public static void option_price_partials_call_black_scholes(double S, 
            double K, double r, double sigma, double time, out double Delta, 
            out double Gamma, out double Theta, out double Vega, out double Rho)
        {
            double time_sqrt = Math.Sqrt(time);
            double d1 = (Math.Log(S / K) + r * time) / (sigma * time_sqrt) 
                + 0.5 * sigma * time_sqrt;
            double d2 = d1 - (sigma * time_sqrt);
            Delta = DStat.NormDist(d1);
            Gamma = DStat.NormDist(d1) / (S * sigma * time_sqrt);
            Theta = -(S * sigma * DStat.NormDist(d1)) / (2 * time_sqrt) 
                - r * K * Math.Exp(-r * time) * DStat.NormDist(d2);
            Vega = S * time_sqrt * DStat.NormDist(d1);
            Rho = K * time * Math.Exp(-r * time) * DStat.NormDist(d2);
        }

        public static void option_price_partials_put_black_scholes(double S, 
            double K, double r, double sigma, double time, out double Delta,
            out double Gamma, out double Theta, out double Vega, out double Rho)
        {
            double time_sqrt = Math.Sqrt(time);
            double d1 = (Math.Log(S / K) + r * time) / (sigma * time_sqrt) 
                + 0.5 * sigma * time_sqrt;
            double d2 = d1 - (sigma * time_sqrt);
            Delta = -DStat.NormDist(-d1);
            Gamma = DStat.NormDist(d1) / (S * sigma * time_sqrt);
            Theta = -(S * sigma * DStat.NormDist(d1)) / (2 * time_sqrt) 
                + r * K * Math.Exp(-r * time) * DStat.NormDist(-d2);
            Vega = S * time_sqrt * DStat.NormDist(d1);
            Rho = -K * time * Math.Exp(-r * time) * DStat.NormDist(-d2);
        }

        public static double option_price_implied_volatility_call_black_scholes_bisections(
            double S, double K, double r, double time, double option_price)
        {
            if (option_price < 0.99 * (S - K * Math.Exp(-time * r)))
            {  // check for arbitrage violations. Option price is too low if this happens
                return 0.0; 
            };

            // simple binomial search for the implied volatility.
            // relies on the value of the option increasing in volatility
            const double ACCURACY = 1.0e-5; // make this smaller for higher accuracy
            const int MAX_ITERATIONS = 100;
            const double HIGH_VALUE = 1e10;
            const double ERROR = -1e40;

            // want to bracket sigma. first find a maximum sigma by finding a sigma
            // with a estimated price higher than the actual price.
            double sigma_low = 1e-5;
            double sigma_high = 0.3;
            double price = option_price_call_black_scholes(S, K, r, sigma_high, time);
            while (price < option_price)
            {
                sigma_high = 2.0 * sigma_high; // keep doubling.
                price = option_price_call_black_scholes(S, K, r, sigma_high, time);
                if (sigma_high > HIGH_VALUE) return ERROR; // panic, something wrong.
            };
            for (int i = 0; i < MAX_ITERATIONS; i++)
            {
                double sigma = (sigma_low + sigma_high) * 0.5;
                price = option_price_call_black_scholes(S, K, r, sigma, time);
                double test = (price - option_price);
                if (Math.Abs(test) < ACCURACY) { return sigma; };
                if (test < 0.0) { sigma_low = sigma; }
                else { sigma_high = sigma; }
            };
            return ERROR;
        }
        
        public static double option_price_implied_volatility_put_black_scholes_bisections(
            double S, double K, double r, double time, double option_price)
        {
            if (option_price < 0.99 * (S - K * Math.Exp(-time * r)))
            {  // check for arbitrage violations. Option price is too low if this happens
                return 0.0;
            };

            // simple binomial search for the implied volatility.
            // relies on the value of the option increasing in volatility
            const double ACCURACY = 1.0e-5; // make this smaller for higher accuracy
            const int MAX_ITERATIONS = 100;
            const double HIGH_VALUE = 1e10;
            const double ERROR = -1e40;

            // want to bracket sigma. first find a maximum sigma by finding a sigma
            // with a estimated price higher than the actual price.
            double sigma_low = 1e-5;
            double sigma_high = 0.3;
            double price = option_price_put_black_scholes(S, K, r, sigma_high, time);
            while (price < option_price)
            {
                sigma_high = 2.0 * sigma_high; // keep doubling.
                price = option_price_put_black_scholes(S, K, r, sigma_high, time);
                if (sigma_high > HIGH_VALUE) return ERROR; // panic, something wrong.
            };
            for (int i = 0; i < MAX_ITERATIONS; i++)
            {
                double sigma = (sigma_low + sigma_high) * 0.5;
                price = option_price_call_black_scholes(S, K, r, sigma, time);
                double test = (price - option_price);
                if (Math.Abs(test) < ACCURACY) { return sigma; };
                if (test < 0.0) { sigma_low = sigma; }
                else { sigma_high = sigma; }
            };
            return ERROR;
        }

        public static double option_price_implied_volatility_call_black_scholes_newton(
            double S, double K, double r, double time, double option_price)
        {
            if (option_price < 0.99 * (S - K * Math.Exp(-time * r)))
            {  // check for arbitrage violations. Option price is too low if this happens
                return 0.0;
            };

            const int MAX_ITERATIONS = 100;
            const double ACCURACY = 1.0e-5;
            double t_sqrt = Math.Sqrt(time);

            double sigma = (option_price / S) / (0.398 * t_sqrt);    // find initial value
            for (int i = 0; i < MAX_ITERATIONS; i++)
            {
                double price = option_price_call_black_scholes(S, K, r, sigma, time);
                double diff = option_price - price;
                if (Math.Abs(diff) < ACCURACY) return sigma;
                double d1 = (Math.Log(S / K) + r * time) / (sigma * t_sqrt) + 0.5 * sigma * t_sqrt;
                double vega = S * t_sqrt * DStat.NormDist(d1);
                sigma = sigma + diff / vega;
            };
            return -99e10;  // something screwy happened, should throw exception
        }
        
        public static double option_price_implied_volatility_put_black_scholes_newton(
            double S, double K, double r, double time, double option_price)
        {
            if (option_price < 0.99 * (S - K * Math.Exp(-time * r)))
            {  // check for arbitrage violations. Option price is too low if this happens
                return 0.0;
            };

            const int MAX_ITERATIONS = 100;
            const double ACCURACY = 1.0e-5;
            double t_sqrt = Math.Sqrt(time);

            double sigma = (option_price / S) / (0.398 * t_sqrt);    // find initial value
            for (int i = 0; i < MAX_ITERATIONS; i++)
            {
                double price = option_price_put_black_scholes(S, K, r, sigma, time);
                double diff = option_price - price;
                if (Math.Abs(diff) < ACCURACY) return sigma;
                double d1 = (Math.Log(S / K) + r * time) / (sigma * t_sqrt) + 0.5 * sigma * t_sqrt;
                double vega = S * t_sqrt * DStat.NormDist(d1);
                sigma = sigma + diff / vega;
            };
            return -99e10;  // something screwy happened, should throw exception
        }
    }
}
