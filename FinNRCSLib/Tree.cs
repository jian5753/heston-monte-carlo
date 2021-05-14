using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace DFinNR
{
    public static class Binomial
    {
        // multiple periode binomial 
        public static List<List<double>> binomial_tree(double S0, double u,
            double d, int no_steps)
        {
            List<List<double>> tree = new List<List<double>>();
            for (int i = 0; i < no_steps; ++i)
            {
                List<double> S = new List<double>();
                for (int j = 0; j < i; ++j)
                {
                    S.Add(S0 * Math.Pow(u, j) * Math.Pow(d, i - j - 1));
                };
                tree.Add(S);
            };

            return tree;
        }

        // Binomial option pricing 
        // one periode binomial 
        public static double option_price_call_european_binomial_single_period(
            double S, double K, double r, double u, double d)
        {
            double p_up = (Math.Exp(r) - d) / (u - d);
            double p_down = 1.0 - p_up;
            double c_u = Math.Max(0.0, (u * S - K));
            double c_d = Math.Max(0.0, (d * S - K));
            double call_price = Math.Exp(-r) * (p_up * c_u + p_down * c_d);

            return call_price;
        }

        // multiple periode binomial 
        public static double option_price_call_european_binomial_multi_period_given_ud(
            double S, double K, double r, double u, double d, int no_periods)
        {
            double Rinv = Math.Exp(-r);   // inverse of interest rate
            double uu = u * u;
            double p_up = (Math.Exp(r) - d) / (u - d);
            double p_down = 1.0 - p_up;
            Vector prices = new Vector(no_periods + 1);     // price of underlying

            prices[0] = S * Math.Pow(d, no_periods);        // fill in the endnodes.

            for (int i = 1; i <= no_periods; ++i)
                prices[i] = uu * prices[i - 1];

            Vector call_values = new Vector(no_periods + 1);    // value of corresponding call 
            for (int i = 0; i <= no_periods; ++i)
            {
                // call payoffs at maturity
                call_values[i] = Math.Max(0.0, (prices[i] - K));
            };

            for (int step = no_periods - 1; step >= 0; --step)
            {
                for (int i = 0; i <= step; ++i)
                {
                    call_values[i] = (p_up * call_values[i + 1]
                        + p_down * call_values[i]) * Rinv;
                };
            };

            return call_values[0];
        }

        /////////////////////////////////////
        // generic binomial trees
        public static double option_price_generic_binomial(double S, double K,
            Func<double, double, double> PayOffFunc, double r,
            double sigma, double t, int steps)
        {
            double R = Math.Exp(r * (t / steps));            // interest rate for each step
            double Rinv = 1.0 / R;                    // inverse of interest rate
            double u = Math.Exp(sigma * Math.Sqrt(t / steps));    // up movement
            double d = 1.0 / u;
            double p_up = (R - d) / (u - d);
            double p_down = 1.0 - p_up;

            Vector prices = new Vector(steps + 1);       // price of underlying
            prices[0] = S * Math.Pow(d, steps);  // fill in the endnodes.
            double uu = u * u;
            for (int i = 1; i <= steps; ++i) prices[i] = uu * prices[i - 1];

            Vector values = new Vector(steps + 1);       // value of corresponding call 
            for (int i = 0; i <= steps; ++i)
                values[i] = PayOffFunc(prices[i], K); // payoffs at maturity

            for (int step = steps - 1; step >= 0; --step)
            {
                for (int i = 0; i <= step; ++i)
                {
                    values[i] = (p_up * values[i + 1] + p_down * values[i]) * Rinv; // value by not exercising
                    prices[i] = d * prices[i + 1];
                    values[i] = Math.Max(values[i], PayOffFunc(prices[i], K));       // check for exercise
                };
            };
            return values[0];
        }

        public static double option_price_delta_generic_binomial(double S,
            double K, Func<double, double, double> PayOffFunc,
            double r, double sigma, double t, int no_steps)
        {
            double R = Math.Exp(r * (t / no_steps));
            double Rinv = 1.0 / R;
            double u = Math.Exp(sigma * Math.Sqrt(t / no_steps));
            double d = 1.0 / u;
            double uu = u * u;
            double pUp = (R - d) / (u - d);
            double pDown = 1.0 - pUp;

            Vector prices = new Vector(no_steps + 1);
            prices[0] = S * Math.Pow(d, no_steps);
            for (int i = 1; i <= no_steps; ++i) prices[i] = uu * prices[i - 1];

            Vector values = new Vector(no_steps + 1);
            for (int i = 0; i <= no_steps; ++i) values[i] = PayOffFunc(prices[i], K);

            for (int CurrStep = no_steps - 1; CurrStep >= 1; --CurrStep)
            {
                for (int i = 0; i <= CurrStep; ++i)
                {
                    prices[i] = d * prices[i + 1];
                    values[i] = (pDown * values[i] + pUp * values[i + 1]) * Rinv;
                    values[i] = Math.Max(values[i], PayOffFunc(prices[i], K));
                };
            };
            double delta = (values[1] - values[0]) / (S * u - S * d);

            return delta;
        }

        // binomial option approximation ////////////////
        public static double option_price_call_european_binomial(double S,
            double K, double r, double sigma, double t, int steps)
        {
            double R = Math.Exp(r * (t / steps));            // interest rate for each step
            double Rinv = 1.0 / R;                    // inverse of interest rate
            double u = Math.Exp(sigma * Math.Sqrt(t / steps));    // up movement
            double uu = u * u;
            double d = 1.0 / u;
            double p_up = (R - d) / (u - d);
            double p_down = 1.0 - p_up;
            Vector prices = new Vector(steps + 1);       // price of underlying
            prices[0] = S * Math.Pow(d, steps);  // fill in the endnodes.
            for (int i = 1; i <= steps; ++i) prices[i] = uu * prices[i - 1];
            Vector call_values = new Vector(steps + 1);       // value of corresponding call 
            for (int i = 0; i <= steps; ++i) call_values[i] = Math.Max(0.0, (prices[i] - K)); // call payoffs at maturity
            for (int step = steps - 1; step >= 0; --step)
            {
                for (int i = 0; i <= step; ++i)
                {
                    call_values[i] = (p_up * call_values[i + 1] + p_down * call_values[i]) * Rinv;
                };
            };
            return call_values[0];
        }

        public static double option_price_put_european_binomial(double S,
            double K, double r, double sigma, double t, int steps)
        {
            double R = Math.Exp(r * (t / steps));            // interest rate for each step
            double Rinv = 1.0 / R;                    // inverse of interest rate
            double u = Math.Exp(sigma * Math.Sqrt(t / steps));    // up movement
            double uu = u * u;
            double d = 1.0 / u;
            double p_up = (R - d) / (u - d);
            double p_down = 1.0 - p_up;
            Vector prices = new Vector(steps + 1);       // price of underlying
            prices[0] = S * Math.Pow(d, steps);  // fill in the endnodes.
            for (int i = 1; i <= steps; ++i) prices[i] = uu * prices[i - 1];
            Vector put_values = new Vector(steps + 1);       // value of corresponding put 
            for (int i = 0; i <= steps; ++i) put_values[i] = Math.Max(0.0, (K - prices[i])); // put payoffs at maturity
            for (int step = steps - 1; step >= 0; --step)
            {
                for (int i = 0; i <= step; ++i)
                {
                    put_values[i] = (p_up * put_values[i + 1] + p_down * put_values[i]) * Rinv;
                };
            };
            return put_values[0];
        }

        public static double option_price_call_american_binomial(double S,
            double K, double r, double sigma, double t, int steps)
        {
            double R = Math.Exp(r * (t / steps));
            double Rinv = 1.0 / R;
            double u = Math.Exp(sigma * Math.Sqrt(t / steps));
            double d = 1.0 / u;
            double p_up = (R - d) / (u - d);
            double p_down = 1.0 - p_up;

            Vector prices = new Vector(steps + 1);       // price of underlying
            prices[0] = S * Math.Pow(d, steps);  // fill in the endnodes.
            double uu = u * u;
            for (int i = 1; i <= steps; ++i) prices[i] = uu * prices[i - 1];

            Vector call_values = new Vector(steps + 1);       // value of corresponding call 
            for (int i = 0; i <= steps; ++i) call_values[i] = Math.Max(0.0, (prices[i] - K)); // call payoffs at maturity

            for (int step = steps - 1; step >= 0; --step)
            {
                for (int i = 0; i <= step; ++i)
                {
                    call_values[i] = (p_up * call_values[i + 1] + p_down * call_values[i]) * Rinv;
                    prices[i] = d * prices[i + 1];
                    call_values[i] = Math.Max(call_values[i], prices[i] - K);       // check for exercise
                };
            };
            return call_values[0];
        }

        public static double option_price_put_american_binomial(double S,
            double K, double r, double sigma, double t, int steps)
        {
            double R = Math.Exp(r * (t / steps));            // interest rate for each step
            double Rinv = 1.0 / R;                    // inverse of interest rate
            double u = Math.Exp(sigma * Math.Sqrt(t / steps));    // up movement
            double uu = u * u;
            double d = 1.0 / u;
            double p_up = (R - d) / (u - d);
            double p_down = 1.0 - p_up;
            Vector prices = new Vector(steps + 1);       // price of underlying
            prices[0] = S * Math.Pow(d, steps);
            for (int i = 1; i <= steps; ++i) prices[i] = uu * prices[i - 1];

            Vector put_values = new Vector(steps + 1);       // value of corresponding put 
            for (int i = 0; i <= steps; ++i) put_values[i] = Math.Max(0.0, (K - prices[i])); // put payoffs at maturity

            for (int step = steps - 1; step >= 0; --step)
            {
                for (int i = 0; i <= step; ++i)
                {
                    put_values[i] = (p_up * put_values[i + 1] + p_down * put_values[i]) * Rinv;
                    prices[i] = d * prices[i + 1];
                    put_values[i] = Math.Max(put_values[i], (K - prices[i]));    // check for exercise
                };
            };
            return put_values[0];
        }

        public static double option_price_call_american_binomial(double S,
            double K, double r, double y, double sigma, double t, int steps)
        {
            double R = Math.Exp(r * (t / steps));            // interest rate for each step
            double Rinv = 1.0 / R;                    // inverse of interest rate
            double u = Math.Exp(sigma * Math.Sqrt(t / steps));    // up movement
            double uu = u * u;
            double d = 1.0 / u;
            double p_up = (Math.Exp((r - y) * (t / steps)) - d) / (u - d);
            double p_down = 1.0 - p_up;
            Vector prices = new Vector(steps + 1);       // price of underlying
            prices[0] = S * Math.Pow(d, steps);
            for (int i = 1; i <= steps; ++i) prices[i] = uu * prices[i - 1]; // fill in the endnodes.

            Vector call_values = new Vector(steps + 1);       // value of corresponding call 
            for (int i = 0; i <= steps; ++i) call_values[i] = Math.Max(0.0, (prices[i] - K)); // call payoffs at maturity

            for (int step = steps - 1; step >= 0; --step)
            {
                for (int i = 0; i <= step; ++i)
                {
                    call_values[i] = (p_up * call_values[i + 1] + p_down * call_values[i]) * Rinv;
                    prices[i] = d * prices[i + 1];
                    call_values[i] = Math.Max(call_values[i], prices[i] - K);       // check for exercise
                };
            };
            return call_values[0];
        }

        public static double option_price_put_american_binomial(double S,
            double K, double r, double y, double sigma, double t, int steps)
        {
            double R = Math.Exp(r * (t / steps));            // interest rate for each step
            double Rinv = 1.0 / R;                    // inverse of interest rate
            double u = Math.Exp(sigma * Math.Sqrt(t / steps));    // up movement
            double uu = u * u;
            double d = 1.0 / u;
            double p_up = (Math.Exp((r - y) * (t / steps)) - d) / (u - d);
            double p_down = 1.0 - p_up;
            Vector prices = new Vector(steps + 1);       // price of underlying
            Vector put_values = new Vector(steps + 1);       // value of corresponding put 

            prices[0] = S * Math.Pow(d, steps);  // fill in the endnodes.
            for (int i = 1; i <= steps; ++i) prices[i] = uu * prices[i - 1];
            for (int i = 0; i <= steps; ++i) put_values[i] = Math.Max(0.0, (K - prices[i])); // put payoffs at maturity
            for (int step = steps - 1; step >= 0; --step)
            {
                for (int i = 0; i <= step; ++i)
                {
                    put_values[i] = (p_up * put_values[i + 1] + p_down * put_values[i]) * Rinv;
                    prices[i] = d * prices[i + 1];
                    put_values[i] = Math.Max(put_values[i], (K - prices[i]));       // check for exercise
                };
            };
            return put_values[0];
        }

        public static double option_price_call_american_discrete_dividends_binomial(
            double S, double K, double r, double sigma, double t,
            int steps, Vector dividend_times, Vector dividend_amounts)
        {
            int no_dividends = dividend_times.size();
            if (no_dividends == 0) return option_price_call_american_binomial(S, K, r, sigma, t, steps);// just do regular
            int steps_before_dividend = (int)(dividend_times[0] / t * steps);
            double R = Math.Exp(r * (t / steps));
            double Rinv = 1.0 / R;
            double u = Math.Exp(sigma * Math.Sqrt(t / steps));
            double d = 1.0 / u;
            double pUp = (R - d) / (u - d);
            double pDown = 1.0 - pUp;
            double dividend_amount = dividend_amounts[0];
            Vector tmp_dividend_times = new Vector(no_dividends - 1);  // temporaries with 
            Vector tmp_dividend_amounts = new Vector(no_dividends - 1);  // one less dividend
            for (int i = 0; i < (no_dividends - 1); ++i)
            {
                tmp_dividend_amounts[i] = dividend_amounts[i + 1];
                tmp_dividend_times[i] = dividend_times[i + 1] - dividend_times[0];
            };
            Vector prices = new Vector(steps_before_dividend + 1);
            Vector call_values = new Vector(steps_before_dividend + 1);
            prices[0] = S * Math.Pow(d, steps_before_dividend);
            for (int i = 1; i <= steps_before_dividend; ++i) prices[i] = u * u * prices[i - 1];
            for (int i = 0; i <= steps_before_dividend; ++i)
            {
                double value_alive
                    = option_price_call_american_discrete_dividends_binomial(prices[i] - dividend_amount, K, r, sigma,
                                                 t - dividend_times[0],// time after first dividend
                                                 steps - steps_before_dividend,
                                                 tmp_dividend_times,
                                                 tmp_dividend_amounts);
                call_values[i] = Math.Max(value_alive, (prices[i] - K));  // compare to exercising now
            };
            for (int step = steps_before_dividend - 1; step >= 0; --step)
            {
                for (int i = 0; i <= step; ++i)
                {
                    prices[i] = d * prices[i + 1];
                    call_values[i] = (pDown * call_values[i] + pUp * call_values[i + 1]) * Rinv;
                    call_values[i] = Math.Max(call_values[i], prices[i] - K);
                };
            };
            return call_values[0];
        }

        public static double option_price_put_american_discrete_dividends_binomial(
            double S, double K, double r, double sigma, double t,
            int steps, Vector dividend_times, Vector dividend_amounts)
        {
            int no_dividends = dividend_times.size();
            if (no_dividends == 0)               // just take the regular binomial 
                return option_price_put_american_binomial(S, K, r, sigma, t, steps);
            int steps_before_dividend = (int)(dividend_times[0] / t * steps);

            double R = Math.Exp(r * (t / steps));
            double Rinv = 1.0 / R;
            double u = Math.Exp(sigma * Math.Sqrt(t / steps));
            double uu = u * u;
            double d = 1.0 / u;

            double pUp = (R - d) / (u - d);
            double pDown = 1.0 - pUp;
            double dividend_amount = dividend_amounts[0];
            Vector tmp_dividend_times = new Vector(no_dividends - 1);  // temporaries with 
            Vector tmp_dividend_amounts = new Vector(no_dividends - 1);  // one less dividend
            for (int i = 0; i < no_dividends - 1; ++i)
            {
                tmp_dividend_amounts[i] = dividend_amounts[i + 1];
                tmp_dividend_times[i] = dividend_times[i + 1] - dividend_times[0];
            };
            Vector prices = new Vector(steps_before_dividend + 1);
            Vector put_values = new Vector(steps_before_dividend + 1);

            prices[0] = S * Math.Pow(d, steps_before_dividend);
            for (int i = 1; i <= steps_before_dividend; ++i) prices[i] = uu * prices[i - 1];
            for (int i = 0; i <= steps_before_dividend; ++i)
            {
                double value_alive
                    = option_price_put_american_discrete_dividends_binomial(
                                                prices[i] - dividend_amount, K, r, sigma,
                                                t - dividend_times[0],               // time after first dividend
                                                steps - steps_before_dividend,
                                                tmp_dividend_times, tmp_dividend_amounts);
                // what is the value of keeping the option alive?  Found recursively, 
                // with one less dividend, the stock price is current value 
                // less the dividend.
                put_values[i] = Math.Max(value_alive, (K - prices[i]));  // compare to exercising now
            };
            for (int step = steps_before_dividend - 1; step >= 0; --step)
            {
                for (int i = 0; i <= step; ++i)
                {
                    prices[i] = d * prices[i + 1];
                    put_values[i] = (pDown * put_values[i] + pUp * put_values[i + 1]) * Rinv;
                    put_values[i] = Math.Max(put_values[i], K - prices[i]);         // check for exercise
                };
            };
            return put_values[0];
        }

        public static double option_price_call_american_proportional_dividends_binomial(
            double S, double K, double r, double sigma, double time, int no_steps,
            Vector dividend_times, Vector dividend_yields)
        {
            // note that the last dividend date should be before the expiry date, problems if dividend at terminal node
            int no_dividends = dividend_times.size();
            if (no_dividends == 0)
            {
                return option_price_call_american_binomial(S, K, r, sigma, time, no_steps); // price w/o dividends
            };
            double delta_t = time / no_steps;
            double R = Math.Exp(r * delta_t);
            double Rinv = 1.0 / R;
            double u = Math.Exp(sigma * Math.Sqrt(delta_t));
            double uu = u * u;
            double d = 1.0 / u;
            double pUp = (R - d) / (u - d);
            double pDown = 1.0 - pUp;
            List<int> dividend_steps = new List<int>(no_dividends); // when dividends are paid
            for (int i = 0; i < no_dividends; ++i)
            {
                dividend_steps[i] = (int)(dividend_times[i] / time * no_steps);
            };
            Vector prices = new Vector(no_steps + 1);
            Vector call_prices = new Vector(no_steps + 1);
            prices[0] = S * Math.Pow(d, no_steps); // adjust downward terminal prices by dividends
            for (int i = 0; i < no_dividends; ++i) { prices[0] *= (1.0 - dividend_yields[i]); };
            for (int i = 1; i <= no_steps; ++i) { prices[i] = uu * prices[i - 1]; };
            for (int i = 0; i <= no_steps; ++i) call_prices[i] = Math.Max(0.0, (prices[i] - K));

            for (int step = no_steps - 1; step >= 0; --step)
            {
                for (int i = 0; i < no_dividends; ++i)
                {   // check whether dividend paid
                    if (step == dividend_steps[i])
                    {
                        for (int j = 0; j <= (step + 1); ++j)
                        {
                            prices[j] *= (1.0 / (1.0 - dividend_yields[i]));
                        };
                    };
                };
                for (int i = 0; i <= step; ++i)
                {
                    call_prices[i] = (pDown * call_prices[i] + pUp * call_prices[i + 1]) * Rinv;
                    prices[i] = d * prices[i + 1];
                    call_prices[i] = Math.Max(call_prices[i], prices[i] - K);         // check for exercise
                };
            };
            return call_prices[0];
        }

        public static double option_price_put_american_proportional_dividends_binomial(
            double S, double K, double r, double sigma, double time, int no_steps,
            Vector dividend_times, Vector dividend_yields)
        {
            // when one assume a dividend yield, the binomial tree recombines 
            // note that the last dividend date should be before the expiry date
            int no_dividends = dividend_times.size();
            if (no_dividends == 0)               // just take the regular binomial 
                return option_price_put_american_binomial(S, K, r, sigma, time, no_steps);
            double R = Math.Exp(r * (time / no_steps));
            double Rinv = 1.0 / R;
            double u = Math.Exp(sigma * Math.Sqrt(time / no_steps));
            double uu = u * u;
            double d = 1.0 / u;
            double pUp = (R - d) / (u - d);
            double pDown = 1.0 - pUp;

            List<int> dividend_steps = new List<int>(no_dividends); // when dividends are paid
            for (int i = 0; i < no_dividends; ++i)
            {
                dividend_steps[i] = (int)(dividend_times[i] / time * no_steps);
            };

            Vector prices = new Vector(no_steps + 1);
            Vector put_prices = new Vector(no_steps + 1);
            prices[0] = S * Math.Pow(d, no_steps);
            for (int i = 0; i < no_dividends; ++i) { prices[0] *= (1.0 - dividend_yields[i]); };
            for (int i = 1; i <= no_steps; ++i) prices[i] = uu * prices[i - 1]; // terminal tree nodes
            for (int i = 0; i <= no_steps; ++i) put_prices[i] = Math.Max(0.0, (K - prices[i]));
            for (int step = no_steps - 1; step >= 0; --step)
            {
                for (int i = 0; i < no_dividends; ++i)
                {   // check whether dividend paid
                    if (step == dividend_steps[i])
                    {
                        for (int j = 0; j <= (step + 1); ++j)
                        {
                            prices[j] *= (1.0 / (1.0 - dividend_yields[i]));
                        };
                    };
                };
                for (int i = 0; i <= step; ++i)
                {
                    prices[i] = d * prices[i + 1];
                    put_prices[i] = (pDown * put_prices[i] + pUp * put_prices[i + 1]) * Rinv;
                    put_prices[i] = Math.Max(put_prices[i], K - prices[i]);         // check for exercise
                };
            };
            return put_prices[0];
        }

        public static double option_price_delta_american_call_binomial(double S,
            double K, double r, double sigma, double t, int no_steps)
        {
            double R = Math.Exp(r * (t / no_steps));
            double Rinv = 1.0 / R;
            double u = Math.Exp(sigma * Math.Sqrt(t / no_steps));
            double d = 1.0 / u;
            double uu = u * u;
            double pUp = (R - d) / (u - d);
            double pDown = 1.0 - pUp;

            Vector prices = new Vector(no_steps + 1);
            prices[0] = S * Math.Pow(d, no_steps);
            for (int i = 1; i <= no_steps; ++i) prices[i] = uu * prices[i - 1];

            Vector call_values = new Vector(no_steps + 1);
            for (int i = 0; i <= no_steps; ++i) call_values[i] = Math.Max(0.0, (prices[i] - K));

            for (int CurrStep = no_steps - 1; CurrStep >= 1; --CurrStep)
            {
                for (int i = 0; i <= CurrStep; ++i)
                {
                    prices[i] = d * prices[i + 1];
                    call_values[i] = (pDown * call_values[i] + pUp * call_values[i + 1]) * Rinv;
                    call_values[i] = Math.Max(call_values[i], prices[i] - K);        // check for exercise
                };
            };
            double delta = (call_values[1] - call_values[0]) / (S * u - S * d);
            return delta;
        }

        public static double option_price_delta_american_put_binomial(double S,
            double K, double r, double sigma, double t, int no_steps)
        {
            Vector prices = new Vector(no_steps + 1);
            Vector put_values = new Vector(no_steps + 1);
            double R = Math.Exp(r * (t / no_steps));
            double Rinv = 1.0 / R;
            double u = Math.Exp(sigma * Math.Sqrt(t / no_steps));
            double d = 1.0 / u;
            double uu = u * u;
            double pUp = (R - d) / (u - d);
            double pDown = 1.0 - pUp;
            prices[0] = S * Math.Pow(d, no_steps);
            int i;
            for (i = 1; i <= no_steps; ++i) prices[i] = uu * prices[i - 1];
            for (i = 0; i <= no_steps; ++i) put_values[i] = Math.Max(0.0, (K - prices[i]));
            for (int CurrStep = no_steps - 1; CurrStep >= 1; --CurrStep)
            {
                for (i = 0; i <= CurrStep; ++i)
                {
                    prices[i] = d * prices[i + 1];
                    put_values[i] = (pDown * put_values[i] + pUp * put_values[i + 1]) * Rinv;
                    put_values[i] = Math.Max(put_values[i], K - prices[i]);        // check for exercise
                };
            };
            double delta = (put_values[1] - put_values[0]) / (S * u - S * d);
            return delta;
        }

        public static void option_price_partials_american_call_binomial(double S,
            double K, double r, double sigma, double time, int no_steps,
            double delta, double gamma, double theta, double vega, double rho)
        {
            Vector prices = new Vector(no_steps + 1);
            Vector call_values = new Vector(no_steps + 1);
            double delta_t = (time / no_steps);
            double R = Math.Exp(r * delta_t);
            double Rinv = 1.0 / R;
            double u = Math.Exp(sigma * Math.Sqrt(delta_t));
            double d = 1.0 / u;
            double uu = u * u;
            double pUp = (R - d) / (u - d);
            double pDown = 1.0 - pUp;
            prices[0] = S * Math.Pow(d, no_steps);
            for (int i = 1; i <= no_steps; ++i) prices[i] = uu * prices[i - 1];
            for (int i = 0; i <= no_steps; ++i) call_values[i] = Math.Max(0.0, (prices[i] - K));
            for (int CurrStep = no_steps - 1; CurrStep >= 2; --CurrStep)
            {
                for (int i = 0; i <= CurrStep; ++i)
                {
                    prices[i] = d * prices[i + 1];
                    call_values[i] = (pDown * call_values[i] + pUp * call_values[i + 1]) * Rinv;
                    call_values[i] = Math.Max(call_values[i], prices[i] - K);        // check for exercise
                };
            };
            double f22 = call_values[2];
            double f21 = call_values[1];
            double f20 = call_values[0];
            for (int i = 0; i <= 1; i++)
            {
                prices[i] = d * prices[i + 1];
                call_values[i] = (pDown * call_values[i] + pUp * call_values[i + 1]) * Rinv;
                call_values[i] = Math.Max(call_values[i], prices[i] - K);        // check for exercise 
            };
            double f11 = call_values[1];
            double f10 = call_values[0];
            prices[0] = d * prices[1];
            call_values[0] = (pDown * call_values[0] + pUp * call_values[1]) * Rinv;
            call_values[0] = Math.Max(call_values[0], S - K);        // check for exercise on first date
            double f00 = call_values[0];
            delta = (f11 - f10) / (S * u - S * d);
            double h = 0.5 * S * (uu - d * d);
            gamma = ((f22 - f21) / (S * (uu - 1)) - (f21 - f20) / (S * (1 - d * d))) / h;
            theta = (f21 - f00) / (2 * delta_t);
            double diff = 0.02;
            double tmp_sigma = sigma + diff;
            double tmp_prices = option_price_call_american_binomial(S, K, r, tmp_sigma, time, no_steps);
            vega = (tmp_prices - f00) / diff;
            diff = 0.05;
            double tmp_r = r + diff;
            tmp_prices = option_price_call_american_binomial(S, K, tmp_r, sigma, time, no_steps);
            rho = (tmp_prices - f00) / diff;
        }

        public static void option_price_partials_american_put_binomial(double S,
            double K, double r, double sigma, double time, int no_steps,
            double delta, double gamma, double theta, double vega, double rho)
        {
            Vector prices = new Vector(no_steps + 1);
            Vector put_values = new Vector(no_steps + 1);
            double delta_t = (time / no_steps);
            double R = Math.Exp(r * delta_t);
            double Rinv = 1.0 / R;
            double u = Math.Exp(sigma * Math.Sqrt(delta_t));
            double d = 1.0 / u;
            double uu = u * u;
            double pUp = (R - d) / (u - d);
            double pDown = 1.0 - pUp;
            prices[0] = S * Math.Pow(d, no_steps);
            int i;
            for (i = 1; i <= no_steps; ++i) prices[i] = uu * prices[i - 1];
            for (i = 0; i <= no_steps; ++i) put_values[i] = Math.Max(0.0, (K - prices[i]));
            for (int CurrStep = no_steps - 1; CurrStep >= 2; --CurrStep)
            {
                for (i = 0; i <= CurrStep; ++i)
                {
                    prices[i] = d * prices[i + 1];
                    put_values[i] = (pDown * put_values[i] + pUp * put_values[i + 1]) * Rinv;
                    put_values[i] = Math.Max(put_values[i], K - prices[i]); // check for exercise
                };
            };
            double f22 = put_values[2];
            double f21 = put_values[1];
            double f20 = put_values[0];
            for (i = 0; i <= 1; i++)
            {
                prices[i] = d * prices[i + 1];
                put_values[i] = (pDown * put_values[i] + pUp * put_values[i + 1]) * Rinv;
                put_values[i] = Math.Max(put_values[i], K - prices[i]); // check for exercise
            };
            double f11 = put_values[1];
            double f10 = put_values[0];
            prices[0] = d * prices[1];
            put_values[0] = (pDown * put_values[0] + pUp * put_values[1]) * Rinv;
            put_values[0] = Math.Max(put_values[0], K - prices[i]); // check for exercise
            double f00 = put_values[0];
            delta = (f11 - f10) / (S * (u - d));
            double h = 0.5 * S * (uu - d * d);
            gamma = ((f22 - f21) / (S * (uu - 1.0)) - (f21 - f20) / (S * (1.0 - d * d))) / h;
            theta = (f21 - f00) / (2 * delta_t);
            double diff = 0.02;
            double tmp_sigma = sigma + diff;
            double tmp_prices = option_price_put_american_binomial(S, K, r, tmp_sigma, time, no_steps);
            vega = (tmp_prices - f00) / diff;
            diff = 0.05;
            double tmp_r = r + diff;
            tmp_prices = option_price_put_american_binomial(S, K, tmp_r, sigma, time, no_steps);
            rho = (tmp_prices - f00) / diff;
        }

        public static double futures_option_price_call_american_binomial(double F,
            double K, double r, double sigma, double time, int no_steps)
        {
            Vector futures_prices = new Vector(no_steps + 1);
            Vector call_values = new Vector(no_steps + 1);
            double t_delta = time / no_steps;
            double Rinv = Math.Exp(-r * (t_delta));
            double u = Math.Exp(sigma * Math.Sqrt(t_delta));
            double d = 1.0 / u;
            double uu = u * u;
            double pUp = (1 - d) / (u - d);   // note how probability is calculated
            double pDown = 1.0 - pUp;
            futures_prices[0] = F * Math.Pow(d, no_steps);
            int i;
            for (i = 1; i <= no_steps; ++i) futures_prices[i] = uu * futures_prices[i - 1]; // terminal tree nodes
            for (i = 0; i <= no_steps; ++i) call_values[i] = Math.Max(0.0, (futures_prices[i] - K));
            for (int step = no_steps - 1; step >= 0; --step)
            {
                for (i = 0; i <= step; ++i)
                {
                    futures_prices[i] = d * futures_prices[i + 1];
                    call_values[i] = (pDown * call_values[i] + pUp * call_values[i + 1]) * Rinv;
                    call_values[i] = Math.Max(call_values[i], futures_prices[i] - K); // check for exercise
                };
            };
            return call_values[0];
        }

        public static double futures_option_price_put_american_binomial(double F,
            double K, double r, double sigma, double time, int no_steps)
        {
            Vector futures_prices = new Vector(no_steps + 1);
            Vector put_values = new Vector(no_steps + 1);
            double t_delta = time / no_steps;
            double Rinv = Math.Exp(-r * (t_delta));
            double u = Math.Exp(sigma * Math.Sqrt(t_delta));
            double d = 1.0 / u;
            double uu = u * u;
            double uInv = 1.0 / u;
            double pUp = (1 - d) / (u - d);
            double pDown = 1.0 - pUp;
            futures_prices[0] = F * Math.Pow(d, no_steps);
            int i;
            for (i = 1; i <= no_steps; ++i)
                futures_prices[i] = uu * futures_prices[i - 1]; // terminal tree nodes
            for (i = 0; i <= no_steps; ++i) put_values[i] = Math.Max(0.0, (K - futures_prices[i]));
            for (int step = no_steps - 1; step >= 0; --step)
            {
                for (i = 0; i <= step; ++i)
                {
                    futures_prices[i] = uInv * futures_prices[i + 1];
                    put_values[i] = (pDown * put_values[i] + pUp * put_values[i + 1]) * Rinv;
                    put_values[i] = Math.Max(put_values[i], K - futures_prices[i]); // check for exercise
                };
            };
            return put_values[0];
        }

        public static double currency_option_price_call_american_binomial(double S,
            double K, double r, double r_f, double sigma, double time, int no_steps)
        {
            Vector exchange_rates = new Vector(no_steps + 1);
            Vector call_values = new Vector(no_steps + 1);
            double t_delta = time / no_steps;
            double Rinv = Math.Exp(-r * (t_delta));
            double u = Math.Exp(sigma * Math.Sqrt(t_delta));
            double d = 1.0 / u;
            double uu = u * u;
            double pUp = (Math.Exp((r - r_f) * t_delta) - d) / (u - d); // adjust for foreign int.rate 
            double pDown = 1.0 - pUp;
            exchange_rates[0] = S * Math.Pow(d, no_steps);
            int i;
            for (i = 1; i <= no_steps; ++i)
            {
                exchange_rates[i] = uu * exchange_rates[i - 1]; // terminal tree nodes
            }
            for (i = 0; i <= no_steps; ++i) call_values[i] = Math.Max(0.0, (exchange_rates[i] - K));
            for (int step = no_steps - 1; step >= 0; --step)
            {
                for (i = 0; i <= step; ++i)
                {
                    exchange_rates[i] = d * exchange_rates[i + 1];
                    call_values[i] = (pDown * call_values[i] + pUp * call_values[i + 1]) * Rinv;
                    call_values[i] = Math.Max(call_values[i], exchange_rates[i] - K); // check for exercise
                };
            };
            return call_values[0];
        }

        public static double currency_option_price_put_american_binomial(double S,
            double K, double r, double r_f, double sigma, double time, int no_steps)
        {
            Vector exchange_rates = new Vector(no_steps + 1);
            Vector put_values = new Vector(no_steps + 1);
            double t_delta = time / no_steps;
            double Rinv = Math.Exp(-r * (t_delta));
            double u = Math.Exp(sigma * Math.Sqrt(t_delta));
            double d = 1.0 / u;
            double uu = u * u;
            double pUp = (Math.Exp((r - r_f) * t_delta) - d) / (u - d); // adjust for foreign int.rate
            double pDown = 1.0 - pUp;
            exchange_rates[0] = S * Math.Pow(d, no_steps);
            int i;
            for (i = 1; i <= no_steps; ++i)
                exchange_rates[i] = uu * exchange_rates[i - 1]; // terminal tree nodes
            for (i = 0; i <= no_steps; ++i) put_values[i] = Math.Max(0.0, (K - exchange_rates[i]));
            for (int step = no_steps - 1; step >= 0; --step)
            {
                for (i = 0; i <= step; ++i)
                {
                    exchange_rates[i] = d * exchange_rates[i + 1];
                    put_values[i] = (pDown * put_values[i] + pUp * put_values[i + 1]) * Rinv;
                    put_values[i] = Math.Max(put_values[i], K - exchange_rates[i]); // check for exercise
                };
            };
            return put_values[0];
        }
    }

    public static class Trinomial
    {
        public static double option_price_call_american_trinomial(double S,
            double K, double r, double q, double sigma, double t, int steps)
        {
            double delta_t = t / steps;
            double Rinv = Math.Exp(-r * (delta_t));
            double sigma_sqr = Math.Pow(sigma, 2);

            double u = Math.Exp(sigma * Math.Sqrt(3.0 * delta_t));
            double d = 1.0 / u;
            double p_u = 1.0 / 6.0 + Math.Sqrt(delta_t / (12.0 * sigma_sqr)) * (r - q - 0.5 * sigma_sqr);
            double p_m = 2.0 / 3.0;
            double p_d = 1.0 / 6.0 - Math.Sqrt(delta_t / (12.0 * sigma_sqr)) * (r - q - 0.5 * sigma_sqr);

            List<Vector> Stree = new List<Vector>();       // price of underlying in a tree
            Vector Svec = new Vector();
            Svec.Add(S);
            for (int step = 1; step <= steps; ++step)
            {
                Stree.Add(Svec);
                Svec.Insert(0, Svec[0] * d); // use the fact that only the extreme values change. 
                Svec.Add(Svec[Svec.size() - 1] * u);
            };

            int m = Svec.size();
            Vector values_next = new Vector(m);       // value of option next step
            for (int i = 0; i < m; ++i) values_next[i] = Math.Max(0.0, Svec[i] - K); // call payoffs at maturity
            Vector values = new Vector();
            for (int step = steps - 1; step >= 0; --step)
            {
                m = Stree[step].size();
                values = new Vector(m);       // value of option
                for (int i = 0; i < m; ++i)
                {
                    values[i] = (p_u * values_next[i + 2] + p_m * values_next[i + 1] + p_d * values_next[i]) * Rinv;
                    values[i] = Math.Max(values[i], Stree[step][i] - K);       // check for exercise
                };
                values_next = values;
            };

            return values[0];
        }

        public static double option_price_put_american_trinomial(double S,
            double K, double r, double q, double sigma, double t, int steps)
        {
            double delta_t = t / steps;
            double Rinv = Math.Exp(-r * (delta_t));
            double sigma_sqr = Math.Pow(sigma, 2);

            double u = Math.Exp(sigma * Math.Sqrt(3.0 * delta_t));
            double d = 1.0 / u;
            double p_u = 1.0 / 6.0 + Math.Sqrt(delta_t / (12.0 * sigma_sqr)) * (r - q - 0.5 * sigma_sqr);
            double p_m = 2.0 / 3.0;
            double p_d = 1.0 / 6.0 - Math.Sqrt(delta_t / (12.0 * sigma_sqr)) * (r - q - 0.5 * sigma_sqr);

            List<Vector> Stree = new List<Vector>();       // price of underlying in a tree
            Vector Svec = new Vector();
            Svec.Add(S);
            for (int step = 1; step <= steps; ++step)
            {
                Stree.Add(Svec);
                Svec.Insert(0, Svec[0] * d); // use the fact that only the extreme values change. 
                Svec.Add(Svec[Svec.size() - 1] * u);
            };

            int m = Svec.size();
            Vector values_next = new Vector(m);       // value of option next step
            for (int i = 0; i < m; ++i) values_next[i] = Math.Max(0.0, K - Svec[i]); // call payoffs at maturity
            Vector values = new Vector();
            for (int step = steps - 1; step >= 0; --step)
            {
                m = Stree[step].size();
                values = new Vector(m);       // value of option
                for (int i = 0; i < m; ++i)
                {
                    values[i] = (p_u * values_next[i + 2] + p_m * values_next[i + 1] + p_d * values_next[i]) * Rinv;
                    values[i] = Math.Max(values[i], K - Stree[step][i]);       // check for exercise
                };

                values_next = values;
            };

            return values[0];
        }
    }

}
