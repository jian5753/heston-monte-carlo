using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace DFinNR
{
    public static class PayOff
    {
        /////////////////////////////
        // payoffs of various options, to be used as function arguments in above simulations
        public static double payoff_call(double S, double K)
        {
            return Math.Max(0.0, S - K);
        }

        public static double payoff_put(double S, double K)
        {
            return Math.Max(0.0, K - S);
        }

        public static double payoff_cash_or_nothing_call(double S, double K)
        {
            if (S >= K) return 1;
            return 0;
        }

        public static double payoff_asset_or_nothing_call(double S, double K)
        {
            if (S >= K) return S;
            return 0;
        }

        public static double payoff_binary_call(double S, double K)
        {
            if (S >= K) return 1;
            return 0;
        }

        public static double payoff_binary_put(double S, double K)
        {
            if (S <= K) return 1;
            return 0;
        }

        /////////////////////////////
        // payoffs of various options, to be used as function arguments in above simulations
        public static double payoff_arithmetric_average_call(Vector prices, double K)
        {
            double sum = prices.Sum(); //prices.begin(), prices.end(), 0.0);
            double avg = sum / prices.size();

            return Math.Max(0.0, avg - K);
        }

        public static double payoff_geometric_average_call(Vector prices, double K)
        {
            double logsum = Math.Log(prices[0]);
            for (int i = 1; i < prices.size(); ++i)
            {
                logsum += Math.Log(prices[i]);
            };

            double avg = Math.Exp(logsum / prices.size());

            return Math.Max(0.0, avg - K);
        }

        public static double payoff_lookback_call(Vector prices, double unused_variable)
        {
            double m = prices.Min(); //*min_element(prices.begin(), prices.end());
            return prices.Last() - m; // always positive or zero
        }

        public static double payoff_lookback_put(Vector prices, double unused_variable)
        {
            double m = prices.Max(); //*max_element(prices.begin(), prices.end());
            return m - prices.Last(); // max is always larger or equal.
        }
    }

    public static class Simulation
    {
        public static double simulate_lognormal_random_variable(
            double S, double r, double sigma, double time)
        {
            double R = (r - 0.5 * Math.Pow(sigma, 2)) * time;
            double SD = sigma * Math.Sqrt(time);

            Random RV = new Random(1234);
            return S * Math.Exp(R + SD * DStat.N_Inv(RV.NextDouble()));// random_normal()
        }

        public static Vector simulate_lognormally_distributed_sequence(double S,
            double r, double sigma, double time, int no_steps)
        {
            Vector prices = new Vector(no_steps);
            double delta_t = time / no_steps;
            double R = (r - 0.5 * Math.Pow(sigma, 2)) * delta_t;
            double SD = sigma * Math.Sqrt(delta_t);
            double S_t = S;                       // initialize at current price
            Random RV = new Random(1234);
            for (int i = 0; i < no_steps; ++i)
            {
                S_t = S_t * Math.Exp(R + SD * DStat.N_Inv(RV.NextDouble())); //random_normal()
                prices[i] = S_t;
            };
            return prices;

        }

        ///////////////////////// simulated option prices //////////////////////////////////////
        // Payoff only function of terminal price
        public static double option_price_call_european_simulated(double S,
            double K, double r, double sigma, double time, int no_sims)
        {
            double R = (r - 0.5 * Math.Pow(sigma, 2)) * time;
            double SD = sigma * Math.Sqrt(time);
            double sum_payoffs = 0.0;
            Random RV = new Random(1234);
            for (int n = 1; n <= no_sims; n++)
            {
                double S_T = S * Math.Exp(R + SD * DStat.N_Inv(RV.NextDouble()));
                sum_payoffs += Math.Max(0.0, S_T - K);
            };
            return Math.Exp(-r * time) * (sum_payoffs / no_sims);
        }

        public static double option_price_put_european_simulated(double S,
            double K, double r, double sigma, double time, int no_sims)
        {
            double sigma_sqr = sigma * sigma;
            double R = (r - 0.5 * sigma_sqr) * time;
            double SD = sigma * Math.Sqrt(time);
            double sum_payoffs = 0.0;
            Random RV = new Random(1234);
            for (int n = 1; n <= no_sims; n++)
            {
                double S_T = S * Math.Exp(R + SD * DStat.N_Inv(RV.NextDouble()));
                sum_payoffs += Math.Max(0.0, K - S_T);
            };
            return Math.Exp(-r * time) * (sum_payoffs / no_sims);
        }

        public static double option_price_delta_call_european_simulated(double S,
            double K, double r, double sigma, double time, int no_sims)
        {
            double R = (r - 0.5 * Math.Pow(sigma, 2)) * time;
            double SD = sigma * Math.Sqrt(time);
            double sum_payoffs = 0.0;
            double sum_payoffs_q = 0.0;
            double q = S * 0.01;
            Random RV = new Random(1234);
            for (int n = 1; n <= no_sims; n++)
            {
                double Z = DStat.N_Inv(RV.NextDouble());
                double S_T = S * Math.Exp(R + SD * Z);
                sum_payoffs += Math.Max(0.0, S_T - K);
                double S_T_q = (S + q) * Math.Exp(R + SD * Z);
                sum_payoffs_q += Math.Max(0.0, S_T_q - K);
            };
            double c = Math.Exp(-r * time) * (sum_payoffs / no_sims);
            double c_q = Math.Exp(-r * time) * (sum_payoffs_q / no_sims);
            return (c_q - c) / q;
        }

        public static double option_price_delta_put_european_simulated(double S,
            double K, double r, double sigma, double time, int no_sims)
        {
            double sigma_sqr = sigma * sigma;
            double R = (r - 0.5 * sigma_sqr) * time;
            double SD = sigma * Math.Sqrt(time);
            double sum_payoffs = 0.0;
            double sum_payoffs_2 = 0.0;
            double q = S * 0.01;
            Random RV = new Random(1234);
            for (int n = 1; n <= no_sims; n++)
            {
                double Z = DStat.N_Inv(RV.NextDouble());
                double S_T = S * Math.Exp(R + SD * Z);
                double S_T_2 = (S + q) * Math.Exp(R + SD * Z);
                sum_payoffs += Math.Max(0.0, K - S_T);
                sum_payoffs_2 += Math.Max(0.0, K - S_T_2);
            };
            double p = Math.Exp(-r * time) * (sum_payoffs / no_sims);
            double p2 = Math.Exp(-r * time) * (sum_payoffs_2 / no_sims);
            return (p2 - p) / q;
        }

        public static double derivative_price_simulate_european_option_generic(
            double S, double K, double r, double sigma, double time,
            Func<Vector, double, double> PayOffFunc, int no_sims)
        {
            double sum_payoffs = 0;

            for (int n = 0; n < no_sims; n++)
            {
                Vector prices = simulate_lognormally_distributed_sequence(S, r, sigma, time, no_sims);
                sum_payoffs += PayOffFunc(prices, K);
            };

            return Math.Exp(-r * time) * (sum_payoffs / no_sims);
        }

        public static double derivative_price_simulate_european_option_generic_with_control_variate(
            double S, double K,	double r, double sigma, double time, 
            Func<double, double, double> PayOffFunc, int no_sims)
        {
            double c_bs = BSMOption.option_price_call_black_scholes(S, S, r, sigma, time);
            // price an at the money Black Scholes call

            double sum_payoffs = 0;
            double sum_payoffs_bs = 0;
            for (int n = 0; n<no_sims; n++)
            {
	            double S_T = simulate_lognormal_random_variable(S, r, sigma, time);
                sum_payoffs += PayOffFunc(S_T, K);
                sum_payoffs_bs += PayOff.payoff_call(S_T, S); 
                // simulate at the money Black Scholes price
            };

            double c_sim = Math.Exp(-r * time) * (sum_payoffs / no_sims);
            double c_bs_sim = Math.Exp(-r * time) * (sum_payoffs_bs / no_sims);
            c_sim += (c_bs-c_bs_sim);

            return c_sim;
        }

        public static double derivative_price_simulate_european_option_generic_with_control_variate(
            double S, double K, double r, double sigma, double time,
            Func<Vector, double, double> PayOffFunc, int no_sims)
        {
            double c_bs = BSMOption.option_price_call_black_scholes(S, S, r, sigma, time);
            // price an at the money Black Scholes call

            double sum_payoffs = 0;
            double sum_payoffs_bs = 0;
            for (int n = 0; n < no_sims; n++)
            {
                Vector prices = simulate_lognormally_distributed_sequence(S, r, sigma, time, no_sims);
                double S1 = prices.Last();
                sum_payoffs += PayOffFunc(prices, K);
                sum_payoffs_bs += PayOff.payoff_call(S1, S); 
                // simulate at the money Black Scholes price
            };

            double c_sim = Math.Exp(-r * time) * (sum_payoffs / no_sims);
            double c_bs_sim = Math.Exp(-r * time) * (sum_payoffs_bs / no_sims);
            c_sim += (c_bs - c_bs_sim);

            return c_sim;
        }

        public static double derivative_price_simulate_european_option_generic_with_antithetic_variate(
            double S, double K, double r, double sigma, double time,
            Func<double, double, double> PayOffFunc, int no_sims)
        {
            double R = (r - 0.5 * Math.Pow(sigma, 2)) * time;
            double SD = sigma * Math.Sqrt(time);
            double sum_payoffs = 0;
            Random RV = new Random(1234);

            for (int n = 0; n < no_sims; n++)
            {
                double err = DStat.N_Inv(RV.NextDouble());

                double S1 = S * Math.Exp(R + SD * err);
                sum_payoffs += PayOffFunc(S1, K);

                double S2 = S * Math.Exp(R + SD * (-err));
                sum_payoffs += PayOffFunc(S2, K);
            };

            return Math.Exp(-r * time) * (sum_payoffs / (2 * no_sims));
        }
    }

}
