using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace DFinNR
{
    public static class Futures
    {
        //// Futures pricing
        public static double futures_price(double S, double r, double time_to_maturity)
        {
            return Math.Exp(r * time_to_maturity) * S;
        }
    }

    public static class Warrant
    {
        /// warrant price
        public static double warrant_price_adjusted_black_scholes(double S,
            double K, double r, double sigma, double time, double m, double n)
        {
            double time_sqrt = Math.Sqrt(time);
            double w = (n / (n + m)) * BSMOption.option_price_call_black_scholes(S, K, r, sigma, time);
            double g = w - (n / (n + m)) * BSMOption.option_price_call_black_scholes(S + (m / n) * w, K, r, sigma, time);
            while (Math.Abs(g) > double.Epsilon)
            {
                double d1 = (Math.Log((S + (m / n)) / K) + r * time) / (sigma * time_sqrt) + 0.5 * sigma * time_sqrt;
                double gprime = 1 - (m / n) * DStat.NormDist(d1);
                w = w - g / gprime;
                g = w - (n / (n + m)) * BSMOption.option_price_call_black_scholes(S + (m / n) * w, K, r, sigma, time);
            };
            return w;


        }

        public static double warrant_price_adjusted_black_scholes(double S,
            double K, double r, double q, double sigma, double time, double m, double n)
        {
            double time_sqrt = Math.Sqrt(time);
            double w = (n / (n + m)) * GeneralBSMOption.option_price_european_call_payout(
                S, K, r, q, sigma, time);
            double g = w - (n / (n + m)) * GeneralBSMOption.
                option_price_european_call_payout(S + (m / n) * w, K, r, q, sigma, time);
            while (Math.Abs(g) > double.Epsilon)
            {
                double d1 = (Math.Log((S + (m / n)) / K) + (r - q) * time)
                    / (sigma * time_sqrt) + 0.5 * sigma * time_sqrt;
                double gprime = 1 - (m / n) * DStat.NormDist(d1);
                w = w - g / gprime;
                g = w - (n / (n + m)) * GeneralBSMOption.option_price_european_call_payout(
                    S + (m / n) * w, K, r, q, sigma, time);
            };
            return w;
        }
    }

    public static class Bond
    {
        // ////////////////////////////////////////////////////////////
        // fixed income derivatives, GBM assumption on bond price
        public static double bond_option_price_call_zero_black_scholes(double B,
            double K, double r, double sigma, double time)
        {
            double time_sqrt = Math.Sqrt(time);
            double d1 = (Math.Log(B / K) + r * time) / (sigma * time_sqrt)
                + 0.5 * sigma * time_sqrt;
            double d2 = d1 - (sigma * time_sqrt);
            double c = B * DStat.NormDist(d1) - K * Math.Exp(-r * time)
                * DStat.NormDist(d2);

            return c;
        }

        public static double bond_option_price_put_zero_black_scholes(double B,
            double K, double r, double sigma, double time)
        {
            double time_sqrt = Math.Sqrt(time);
            double d1 = (Math.Log(B / K) + r * time) / (sigma * time_sqrt)
                + 0.5 * sigma * time_sqrt;
            double d2 = d1 - (sigma * time_sqrt);
            double p = K * Math.Exp(-r * time) * DStat.NormDist(-d2)
                - B * DStat.NormDist(-d1);

            return p;
        }

        public static double bond_option_price_call_coupon_bond_black_scholes(double B,
            double K, double r, double sigma, double time,
            Vector coupon_times, Vector coupon_amounts)
        {
            double adjusted_B = B;

            for (int i = 0; i < coupon_times.size(); i++)
            { // subtract present value of coupons
                if (coupon_times[i] <= time)
                { // coupon paid befor option expiry
                    adjusted_B -= coupon_amounts[i] * Math.Exp(-r * coupon_times[i]);
                };
            };

            return bond_option_price_call_zero_black_scholes(adjusted_B, K, r, sigma, time);
        }

        public static double bond_option_price_put_coupon_bond_black_scholes(double B,
            double K, double r, double sigma, double time,
            Vector coupon_times, Vector coupon_amounts)
        {
            double adjusted_B = B;
            for (int i = 0; i < coupon_times.size(); i++)
            {
                if (coupon_times[i] <= time)
                {
                    adjusted_B -= coupon_amounts[i] * Math.Exp(-r * coupon_times[i]);
                };
            };

            return bond_option_price_put_zero_black_scholes(adjusted_B, K, r, sigma, time);
        }

        public static double bond_option_price_call_american_binomial(double B,
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
            Vector put_values = new Vector(steps + 1);       // value of corresponding put 

            prices[0] = B * Math.Pow(d, steps);  // fill in the endnodes.
            for (int i = 1; i <= steps; ++i)
                prices[i] = uu * prices[i - 1];
            for (int i = 0; i <= steps; ++i)
                put_values[i] = Math.Max(0.0, (prices[i] - K)); // put payoffs at maturity
            for (int step = steps - 1; step >= 0; --step)
            {
                for (int i = 0; i <= step; ++i)
                {
                    put_values[i] = (p_up * put_values[i + 1]
                        + p_down * put_values[i]) * Rinv;
                    prices[i] = d * prices[i + 1];
                    put_values[i] = Math.Max(put_values[i], (prices[i] - K));  // check for exercise
                };
            };

            return put_values[0];
        }

        public static double bond_option_price_put_american_binomial(double B,
            double K, double r, double sigma, double t, int steps)
        {
            double R = Math.Exp(r * (t / steps));     // interest rate for each step
            double Rinv = 1.0 / R;                    // inverse of interest rate
            double u = Math.Exp(sigma * Math.Sqrt(t / steps));    // up movement
            double uu = u * u;
            double d = 1.0 / u;
            double p_up = (R - d) / (u - d);
            double p_down = 1.0 - p_up;
            Vector prices = new Vector(steps + 1);       // price of underlying
            Vector put_values = new Vector(steps + 1);       // value of corresponding put 

            prices[0] = B * Math.Pow(d, steps);  // fill in the endnodes.
            for (int i = 1; i <= steps; ++i)
                prices[i] = uu * prices[i - 1];

            for (int i = 0; i <= steps; ++i)
                put_values[i] = Math.Max(0.0, (K - prices[i])); // put payoffs at maturity

            for (int step = steps - 1; step >= 0; --step)
            {
                for (int i = 0; i <= step; ++i)
                {
                    put_values[i] = (p_up * put_values[i + 1]
                        + p_down * put_values[i]) * Rinv;
                    prices[i] = d * prices[i + 1];
                    put_values[i] = Math.Max(put_values[i], (K - prices[i]));       // check for exercise
                };
            };
            return put_values[0];
        }

        // term structure derivatives, analytical solutions
        public static double bond_option_price_call_zero_vasicek(double X, double r,
            double option_time_to_maturity, double bond_time_to_maturity,
            double a, double b, double sigma)
        {
            double T_t = option_time_to_maturity;
            double s_t = bond_time_to_maturity;
            double T_s = s_t - T_t;
            double v_t_T;
            double sigma_P;
            if (a == 0.0)
            {
                v_t_T = sigma * Math.Sqrt(T_t);
                sigma_P = sigma * T_s * Math.Sqrt(T_t);
            }
            else
            {
                v_t_T = Math.Sqrt(sigma * sigma * (1 - Math.Exp(-2 * a * T_t))
                    / (2 * a));
                double B_T_s = (1 - Math.Exp(-a * T_s)) / a;
                sigma_P = v_t_T * B_T_s;
            };
            double h = (1.0 / sigma_P) * Math.Log(
                term_structure_utils.term_structure_discount_factor_vasicek(s_t, r, a, b, sigma)
                / (term_structure_utils.term_structure_discount_factor_vasicek(T_t, r, a, b, sigma) * X))
                + sigma_P / 2.0;
            double c = term_structure_utils.term_structure_discount_factor_vasicek(s_t, r, a, b, sigma)
                * DStat.NormDist(h)
                - X * term_structure_utils.term_structure_discount_factor_vasicek(T_t, r, a, b, sigma)
                * DStat.NormDist(h - sigma_P);

            return c;
        }

        public static double bond_option_price_put_zero_vasicek(double X, double r,
            double option_time_to_maturity, double bond_time_to_maturity,
            double a, double b, double sigma)
        {
            double s_t = bond_time_to_maturity;
            double T_t = option_time_to_maturity;
            double T_s = s_t - T_t;
            double v_t_T;
            double sigma_P;
            if (a == 0.0)
            {
                v_t_T = sigma * Math.Sqrt(T_t);
                sigma_P = sigma * T_s * Math.Sqrt(T_t);
            }
            else
            {
                v_t_T = Math.Sqrt(sigma * sigma * (1 - Math.Exp(-2 * a * T_t)) / (2 * a));
                double B_T_s = (1 - Math.Exp(-a * T_s)) / a;
                sigma_P = v_t_T * B_T_s;
            };
            double h = (1.0 / sigma_P) *
                Math.Log(term_structure_utils.term_structure_discount_factor_vasicek(s_t, r, a, b, sigma)
                / (term_structure_utils.term_structure_discount_factor_vasicek(T_t, r, a, b, sigma) * X))
                + sigma_P / 2.0;
            double p = term_structure_utils.term_structure_discount_factor_vasicek(T_t, r, a, b, sigma)
                * DStat.NormDist(-h + sigma_P)
                - term_structure_utils.term_structure_discount_factor_vasicek(s_t, r, a, b, sigma)
                * DStat.NormDist(-h);

            return p;
        }

        // //////////////////////////////////////////////////////////////////////////////
        // binomial term structure models
        // bond option, rendlemann bartter  (binomial)
        public static double bond_option_price_call_zero_american_rendleman_bartter(
            double K, double option_maturity, double S, double M, double interest,
            double bond_maturity, double maturity_payment, int no_steps)
        {
            double delta_t = bond_maturity / no_steps;

            double u = Math.Exp(S * Math.Sqrt(delta_t));
            double d = 1 / u;
            double p_up = (Math.Exp(M * delta_t) - d) / (u - d);
            double p_down = 1.0 - p_up;

            Vector r = new Vector(no_steps + 1);
            r[0] = interest * Math.Pow(d, no_steps);
            double uu = u * u;
            for (int i = 1; i <= no_steps; ++i)
                r[i] = r[i - 1] * uu;

            Vector P = new Vector(no_steps + 1);
            for (int i = 0; i <= no_steps; ++i)
                P[i] = maturity_payment;

            int no_call_steps = (int)(no_steps * option_maturity / bond_maturity);

            for (int curr_step = no_steps; curr_step > no_call_steps; --curr_step)
            {
                for (int i = 0; i < curr_step; i++)
                {
                    r[i] = r[i] * u;
                    P[i] = Math.Exp(-r[i] * delta_t) * (p_down * P[i] + p_up * P[i + 1]);
                };
            };
            Vector C = new Vector(no_call_steps + 1);
            for (int i = 0; i <= no_call_steps; ++i)
            {
                C[i] = Math.Max(0.0, P[i] - K);
            };

            for (int curr_step = no_call_steps; curr_step >= 0; --curr_step)
            {
                for (int i = 0; i < curr_step; i++)
                {
                    r[i] = r[i] * u;
                    P[i] = Math.Exp(-r[i] * delta_t) * (p_down * P[i] + p_up * P[i + 1]);
                    C[i] = Math.Max(P[i] - K, Math.Exp(-r[i] * delta_t) * (p_up * C[i + 1] + p_down * C[i]));
                };
            };
            return C[0];
        }

        public static List<Vector> interest_rate_trees_gbm_build(double r0, double u, double d, int n)
        {
            List<Vector> tree = new List<Vector>();
            Vector r = new Vector(1);
            r[0] = r0;
            tree.Add(r);

            for (int i = 1; i <= n; ++i)
            {
                double rtop = r[r.size() - 1] * u;
                for (int j = 0; j < i; ++j)
                {
                    r[j] = d * r[j];
                };
                r.Add(rtop);
                tree.Add(r);
            };

            return tree;
        }

        public static double interest_rate_trees_gbm_value_of_cashflows(Vector cflow, List<Vector> r_tree, double q)
        {
            int n = cflow.size();
            List<Vector> values = new List<Vector>(n);

            Vector value = new Vector(n);
            for (int i = 0; i < n; i++)
            {
                value[i] = cflow[n - 1];
            };

            values[n - 1] = value;

            for (int t = n - 1; t > 0; --t)
            {
                value = new Vector(t, 0.0);
                for (int i = 0; i < t; ++i)
                {
                    value[i] = cflow[t - 1] + Math.Exp(-r_tree[t - 1][i]) * (q * values[t][i] + (1 - q) * values[t][i + 1]);
                };

                values[t - 1] = value;
            };

            return values[0][0];
        }

        public static double interest_rate_trees_gbm_value_of_callable_bond(Vector cflows, List<Vector> r_tree,
            double q, int first_call_time, double call_price)
        {
            int n = cflows.size();
            List<Vector> values = new List<Vector>(n);

            Vector value = new Vector(n);
            for (int i = 0; i < n; i++)
            {
                value[i] = cflows[n - 1];
            };
            values[n - 1] = value;

            for (int t = n - 1; t > 0; --t)
            {
                value = new Vector(t, 0.0);
                for (int i = 0; i < t; ++i)
                {
                    value[i] = cflows[t - 1] + Math.Exp(-r_tree[t - 1][i]) * (q * values[t][i] + (1 - q) * values[t][i + 1]);
                    if (t >= first_call_time)
                    {
                        value[i] = Math.Min(value[i], call_price);
                    };
                };
                values[t - 1] = value;
            };

            return values[0][0];
        }

        public static List<time_contingent_cash_flows>
            build_time_series_of_bond_time_contingent_cash_flows(
            Vector initial_times, Vector initial_cflows)
        {
            List<time_contingent_cash_flows> vec_cf = new List<time_contingent_cash_flows>();
            Vector times = new Vector(initial_times);
            Vector cflows = new Vector(initial_cflows);
            while (times.size() > 0)
            {
                vec_cf.Add(new time_contingent_cash_flows(times, cflows));
                Vector tmp_times = new Vector();
                Vector tmp_cflows = new Vector();
                for (int i = 0; i < times.size(); ++i)
                {
                    if (times[i] - 1.0 >= 0.0)
                    {
                        tmp_times.Add(times[i] - 1);
                        tmp_cflows.Add(cflows[i]);
                    };
                };
                times = tmp_times; cflows = tmp_cflows;
            };

            return vec_cf;
        }

        public static double price_european_call_option_on_bond_using_ho_lee(
            term_structure_class initial, double delta, double pi,
            Vector underlying_bond_cflow_times, Vector underlying_bond_cflows,
            double K, double time_to_maturity)
        {
            int T = (int)(time_to_maturity + 0.0001);
            List<List<term_structure_class_ho_lee>> hl_tree
                = term_structure_utils.term_structure_ho_lee_build_term_structure_tree(initial, T + 1, delta, pi);
            List<time_contingent_cash_flows> vec_cf
                = build_time_series_of_bond_time_contingent_cash_flows(underlying_bond_cflow_times, underlying_bond_cflows);

            Vector values = new Vector(T + 1);
            for (int i = 0; i <= T; ++i)
            {
                values[i] = Math.Max(0.0, term_structure_utils.bonds_price(vec_cf[T + 1].times, vec_cf[T + 1].cash_flows, hl_tree[T + 1][i]) - K);
            };

            for (int t = T; t >= 0; --t)
            {
                Vector values_this = new Vector(t + 1);
                for (int i = 0; i <= t; ++i)
                {
                    values_this[i] = (pi * values[i + 1] + (1.0 - pi) * values[i]) * hl_tree[t][i].d(1);
                };

                values = values_this;
            };

            return values[0];
        }
    }

    public static class Option
    {
        /////////////////// alternative stochastic processes ////////////////
        public static double option_price_call_merton_jump_diffusion(double S,
            double K, double r, double sigma, double time_to_maturity,
            double lambda, double kappa, double delta)
        {
            int MAXN = 50;
            double tau = time_to_maturity;
            double sigma_sqr = sigma * sigma;
            double delta_sqr = delta * delta;
            double lambdaprime = lambda * (1 + kappa);
            double gamma = Math.Log(1 + kappa);
            double c = Math.Exp(-lambdaprime * tau) *
                BSMOption.option_price_call_black_scholes(S, K, r - lambda * kappa, sigma, tau);

            double log_n = 0;
            for (int n = 1; n <= MAXN; ++n)
            {
                log_n += Math.Log((double)n);
                double sigma_n = Math.Sqrt(sigma_sqr + n * delta_sqr / tau);
                double r_n = r - lambda * kappa + n * gamma / tau;
                c += Math.Exp(-lambdaprime * tau + n * Math.Log(lambdaprime * tau) - log_n) *
                    BSMOption.option_price_call_black_scholes(S, K, r_n, sigma_n, tau);
            };

            return c;
        }

        ////////////// path dependent and other exotic options ////////////////////////////////
        public static double option_price_call_bermudan_binomial(double S,
            double K, double r, double q, double sigma, double time,
            Vector potential_exercise_times, int steps)
        {
            double delta_t = time / steps;
            double R = Math.Exp(r * delta_t);
            double Rinv = 1.0 / R;
            double u = Math.Exp(sigma * Math.Sqrt(delta_t));
            double uu = u * u;
            double d = 1.0 / u;
            double p_up = (Math.Exp((r - q) * (delta_t)) - d) / (u - d);
            double p_down = 1.0 - p_up;
            Vector prices = new Vector(steps + 1);
            Vector call_values = new Vector(steps + 1);

            List<int> potential_exercise_steps = new List<int>();
            // create list of steps at which exercise may happen
            for (int i = 0; i < potential_exercise_times.size(); ++i)
            {
                double t = potential_exercise_times[i];
                if ((t > 0.0) && (t < time))
                {
                    potential_exercise_steps.Add((int)(t / delta_t));
                };
            };
            prices[0] = S * Math.Pow(d, steps);  // fill in the endnodes.
            for (int i = 1; i <= steps; ++i)
                prices[i] = uu * prices[i - 1];

            for (int i = 0; i <= steps; ++i)
                call_values[i] = Math.Max(0.0, (prices[i] - K));

            for (int step = steps - 1; step >= 0; --step)
            {
                bool check_exercise_this_step = false;

                for (int j = 0; j < potential_exercise_steps.Count; ++j)
                {
                    if (step == potential_exercise_steps[j]) { check_exercise_this_step = true; };
                };

                for (int i = 0; i <= step; ++i)
                {
                    call_values[i] = (p_up * call_values[i + 1]
                        + p_down * call_values[i]) * Rinv;
                    prices[i] = d * prices[i + 1];

                    if (check_exercise_this_step)
                        call_values[i] = Math.Max(call_values[i], prices[i] - K);
                };
            };

            return call_values[0];
        }

        public static double option_price_put_bermudan_binomial(double S,
            double K, double r, double q, double sigma, double time,
            Vector potential_exercise_times, int steps)
        {
            double delta_t = time / steps;
            double R = Math.Exp(r * delta_t);
            double Rinv = 1.0 / R;
            double u = Math.Exp(sigma * Math.Sqrt(delta_t));
            double uu = u * u;
            double d = 1.0 / u;
            double p_up = (Math.Exp((r - q) * delta_t) - d) / (u - d);
            double p_down = 1.0 - p_up;
            Vector prices = new Vector(steps + 1);
            Vector put_values = new Vector(steps + 1);

            List<int> potential_exercise_steps = new List<int>();
            // create list of steps at which exercise may happen
            for (int i = 0; i < potential_exercise_times.size(); ++i)
            {
                double t = potential_exercise_times[i];
                if ((t > 0.0) && (t < time))
                {
                    potential_exercise_steps.Add((int)(t / delta_t));
                };
            };

            prices[0] = S * Math.Pow(d, steps);  // fill in the endnodes.
            for (int i = 1; i <= steps; ++i)
                prices[i] = uu * prices[i - 1];

            for (int i = 0; i <= steps; ++i)
                put_values[i] = Math.Max(0.0, (K - prices[i])); // put payoffs at maturity

            for (int step = steps - 1; step >= 0; --step)
            {
                bool check_exercise_this_step = false;

                for (int j = 0; j < potential_exercise_steps.Count; ++j)
                {
                    if (step == potential_exercise_steps[j]) { check_exercise_this_step = true; };
                };

                for (int i = 0; i <= step; ++i)
                {
                    put_values[i] = (p_up * put_values[i + 1]
                        + p_down * put_values[i]) * Rinv;
                    prices[i] = d * prices[i + 1];

                    if (check_exercise_this_step)
                        put_values[i] = Math.Max(put_values[i], K - prices[i]);
                };
            };

            return put_values[0];
        }

        public static double option_price_european_lookback_call(double S,
            double Smin, double r, double q, double sigma, double time)
        {
            if (r == q)
                return 0;
            double sigma_sqr = sigma * sigma;
            double time_sqrt = Math.Sqrt(time);

            double a1 = (Math.Log(S / Smin) + (r - q + sigma_sqr / 2.0) * time) / (sigma * time_sqrt);
            double a2 = a1 - sigma * time_sqrt;
            double a3 = (Math.Log(S / Smin) + (-r + q + sigma_sqr / 2.0) * time) / (sigma * time_sqrt);
            double Y1 = 2.0 * (r - q - sigma_sqr / 2.0) * Math.Log(S / Smin) / sigma_sqr;

            return S * Math.Exp(-q * time) * DStat.NormDist(a1)
                - S * Math.Exp(-q * time) * (sigma_sqr / (2.0 * (r - q))) * DStat.NormDist(-a1)
                - Smin * Math.Exp(-r * time) * (DStat.NormDist(a2)
                - (sigma_sqr / (2 * (r - q))) * Math.Exp(Y1) * DStat.NormDist(-a3));
        }

        public static double option_price_european_lookback_put(double S,
            double Smax, double r, double q, double sigma, double time)
        {
            if (r == q)
                return 0;
            double sigma_sqr = sigma * sigma;
            double time_sqrt = Math.Sqrt(time);

            double b1 = (Math.Log(S / Smax) + (-r + q + sigma_sqr / 2.0) * time) / (sigma * time_sqrt);
            double b2 = b1 - sigma * time_sqrt;
            double b3 = (Math.Log(S / Smax) + (r - q - sigma_sqr / 2.0) * time) / (sigma * time_sqrt);

            double Y2 = (2.0 * (r - q - sigma_sqr / 2.0) * Math.Log(Smax / S)) / sigma_sqr;

            double p = Smax * Math.Exp(-r * time) * (DStat.NormDist(b1)
                - (sigma_sqr / (2 * (r - q))) * Math.Exp(Y2) * DStat.NormDist(-b3))
            + S * Math.Exp(-q * time) * (sigma_sqr / (2.0 * (r - q))) * DStat.NormDist(-b2)
            - S * Math.Exp(-q * time) * DStat.NormDist(b2);

            return p;
        }

        public static double option_price_asian_geometric_average_price_call(double S,
            double K, double r, double q, double sigma, double time)
        {
            double sigma_sqr = Math.Pow(sigma, 2);
            double adj_div_yield = 0.5 * (r + q + sigma_sqr / 6.0);
            double adj_sigma = sigma / Math.Sqrt(3.0);
            double adj_sigma_sqr = Math.Pow(adj_sigma, 2);
            double time_sqrt = Math.Sqrt(time);

            double d1 = (Math.Log(S / K) + (r - adj_div_yield + 0.5
                * adj_sigma_sqr) * time) / (adj_sigma * time_sqrt);
            double d2 = d1 - (adj_sigma * time_sqrt);

            double call_price = S * Math.Exp(-adj_div_yield * time)
                * DStat.NormDist(d1) - K * Math.Exp(-r * time) * DStat.NormDist(d2);

            return call_price;
        }

        /////////// approximated option prices ////////////////////////
        public static double option_price_american_put_approximated_johnson(double S,
            double X, double r, double sigma, double time)
        {
            double sigma_sqr = Math.Pow(sigma, 2);
            double a0 = 3.9649;
            double a1 = 0.032325;
            double b0 = 1.040803;
            double b1 = 0.00963;

            double gamma = 2 * r / sigma_sqr;
            double m = (sigma_sqr * time) / (b0 * sigma_sqr * time + b1);
            double Sc = X * Math.Pow(((gamma) / (1 + gamma)), m);
            double l = (Math.Log(S / Sc)) / (Math.Log(X / Sc));

            double alpha = Math.Pow(((r * time) / (a0 * r * time + a1)), l);

            double P = alpha * BSMOption.option_price_put_black_scholes(S, X * Math.Exp(r * time), r, sigma, time)
                + (1 - alpha) * BSMOption.option_price_put_black_scholes(S, X, r, sigma, time);

            double p = BSMOption.option_price_put_black_scholes(S, X, r, sigma, time);  // for safety use the Black Scholes as lower bound

            return Math.Max(p, P);
        }

        public static double option_price_american_call_approximated_baw(double S,
            double X, double r, double b, double sigma, double time)
        {
            double ACCURACY = 1.0e-6;
            double sigma_sqr = sigma * sigma;
            double time_sqrt = Math.Sqrt(time);
            double nn = 2.0 * b / sigma_sqr;
            double m = 2.0 * r / sigma_sqr;
            double K = 1.0 - Math.Exp(-r * time);
            double q2 = (-(nn - 1) + Math.Sqrt(Math.Pow((nn - 1), 2.0) + (4 * m / K))) * 0.5;

            double q2_inf = 0.5 * (-(nn - 1) + Math.Sqrt(Math.Pow((nn - 1), 2.0) + 4.0 * m));    // seed value from paper
            double S_star_inf = X / (1.0 - 1.0 / q2_inf);
            double h2 = -(b * time + 2.0 * sigma * time_sqrt) * (X / (S_star_inf - X));
            double S_seed = X + (S_star_inf - X) * (1.0 - Math.Exp(h2));

            int no_iterations = 0; // iterate on S to find S_star, using Newton steps
            double Si = S_seed;
            double g = 1;
            double gprime = 1.0;
            double c = 0;
            while ((Math.Abs(g) > ACCURACY)
               && (Math.Abs(gprime) > ACCURACY) // to avoid exploding Newton's  
               && (no_iterations++ < 500)
               && (Si > 0.0))
            {
                c = GeneralBSMOption.option_price_european_call_payout(Si, X, r, b, sigma, time);
                double d1 = (Math.Log(Si / X) + (b + 0.5 * sigma_sqr) * time) / (sigma * time_sqrt);
                g = (1.0 - 1.0 / q2) * Si - X - c + (1.0 / q2) * Si * Math.Exp((b - r) * time) * DStat.NormDist(d1);
                gprime = (1.0 - 1.0 / q2) * (1.0 - Math.Exp((b - r) * time) * DStat.NormDist(d1))
                    + (1.0 / q2) * Math.Exp((b - r) * time) * DStat.NormProb(d1) * (1.0 / (sigma * time_sqrt));
                Si = Si - (g / gprime);
            };
            double S_star = 0;
            if (Math.Abs(g) > ACCURACY) { S_star = S_seed; } // did not converge
            else { S_star = Si; };
            double C = 0;
            c = GeneralBSMOption.option_price_european_call_payout(S, X, r, b, sigma, time);
            if (S >= S_star)
            {
                C = S - X;
            }
            else
            {
                double d1 = (Math.Log(S_star / X) + (b + 0.5 * sigma_sqr) * time) / (sigma * time_sqrt);
                double A2 = (1.0 - Math.Exp((b - r) * time) * DStat.NormDist(d1)) * (S_star / q2);
                C = c + A2 * Math.Pow((S / S_star), q2);
            };
            return Math.Max(C, c); // know value will never be less than BS value
        }

        public static double option_price_american_put_approximated_baw(double S,
            double X, double r, double b, double sigma, double time)
        {
            double ACCURACY = 1.0e-6;
            double sigma_sqr = sigma * sigma;
            double time_sqrt = Math.Sqrt(time);
            double M = 2.0 * r / sigma_sqr;
            double NN = 2.0 * b / sigma_sqr;
            double K = 1.0 - Math.Exp(-r * time);
            double q1 = 0.5 * (-(NN - 1) - Math.Sqrt(Math.Pow((NN - 1), 2.0) + (4.0 * M / K)));

            // now find initial S value 
            double q1_inf = 0.5 * (-(NN - 1) - Math.Sqrt(Math.Pow((NN - 1), 2.0) + 4.0 * M));
            double S_star_star_inf = X / (1.0 - 1.0 / q1_inf);
            double h1 = (b * time - 2 * sigma * time_sqrt) * (X / (X - S_star_star_inf));
            double S_seed = S_star_star_inf + (X - S_star_star_inf) * Math.Exp(h1);

            double Si = S_seed;  // now do Newton iterations to solve for S**
            int no_iterations = 0;
            double g = 1;
            double p = 0;
            double gprime = 1;
            while ((Math.Abs(g) > ACCURACY)
               && (Math.Abs(gprime) > ACCURACY) // to avoid non-convergence
               && (no_iterations++ < 500)
               && Si > 0.0)
            {
                p = GeneralBSMOption.option_price_european_put_payout(Si, X, r, b, sigma, time);
                double d1 = (Math.Log(Si / X) + (b + 0.5 * sigma_sqr) * time) / (sigma * time_sqrt);
                g = X - Si - p + (1 - Math.Exp((b - r) * time) * DStat.NormDist(-d1)) * Si / q1;
                gprime = (1.0 / q1 - 1.0) * (1.0 - Math.Exp((b - r) * time) * DStat.NormDist(-d1))
                    + (1.0 / q1) * Math.Exp((b - r) * time) * (1.0 / (sigma * time_sqrt)) * DStat.NormProb(-d1);
                Si = Si - (g / gprime);
            };
            double S_star_star = Si;
            if (g > ACCURACY)
            {
                S_star_star = S_seed;
            };  // if not found g**
            double P = 0;
            p = GeneralBSMOption.option_price_european_put_payout(S, X, r, b, sigma, time);
            if (S > S_star_star)
            {
                double d1 = (Math.Log(S_star_star / X)
                         + (b + 0.5 * sigma_sqr) * time) / (sigma * time_sqrt);
                double A1 = -(S_star_star / q1) * (1 - Math.Exp((b - r) * time) * DStat.NormDist(-d1));
                P = p + A1 * Math.Pow((S / S_star_star), q1);
            }
            else
            {
                P = X - S;
            };
            return Math.Max(P, p);  // should not be lower than Black Scholes value
        }
    }

    /////////// approximated option prices ////////////////////////
    /* Geske Johnson Approximation
    public static double d1(double S, double X, double r, double sigma, 
        double tau)
    {
        return (Math.Log(S / X) + (r + 0.5 * Math.Pow(sigma, 2)) * tau) 
            / (sigma * Math.Sqrt(tau));
    }

    public static double d2(double S, double X, double r, double sigma, 
        double tau)
    {
        return d1(S, X, r, sigma, tau) - sigma * Math.Sqrt(tau);
    }

    public static double calcP2(double S, double X, double r, double sigma,
        double time, double t2, double S2_bar, double rho12)
    {
        double P2 = X * Math.Exp(-r * t2) 
            * DStat.NormDist(-d2(S, S2_bar, r, sigma, t2));

        P2 -= S * DStat.NormDist(-d1(S, S2_bar, r, sigma, t2));

        P2 += X * Math.Exp(-r * time) * DStat.NormDist(
            d2(S, S2_bar, r, sigma, t2), -d2(S, X, r, sigma, time), -rho12);

        P2 -= S * DStat.NormDist(d1(S, S2_bar, r, sigma, t2), 
            -d1(S, X, r, sigma, time), -rho12);

        return P2;
    }

    public static double option_price_american_put_approximated_geske_johnson(
        double S, double X, double r, double sigma, double time)
    {
        double ACCURACY = 1e-6;

        double P1 = BSMOption.option_price_put_black_scholes(S, X, r, sigma, time);
        double rho12 = 1.0 / Math.Sqrt(2.0);
        double rho13 = 1.0 / Math.Sqrt(3.0);
        double rho23 = Math.Sqrt(2.0 / 3.0);
        double t2 = time / 2.0;
        double t23 = time * 2.0 / 3.0;
        double t3 = time / 3.0;
        double Si = S;
        double S2_bar = S;
        double g = 1;
        double gprime = 1;
        while (Math.Abs(g) > ACCURACY)
        {
            g = Si - X + BSMOption.option_price_put_black_scholes(Si, X, r, sigma, t2);
            gprime = 1.0 + BSMOption.option_price_delta_put_black_scholes(Si, X, r, sigma, t2);
            S2_bar = Si;
            Si = Si - g / gprime;
        };

        double P2 = calcP2(S, X, r, sigma, time, t2, S2_bar, rho12);
        P2 = Math.Max(P1, P2); // for safety, use one less step as lower bound

        double S23_bar = S2_bar;
        g = 1;
        while (Math.Abs(g) > ACCURACY)
        {
            g = Si - X + BSMOption.option_price_put_black_scholes(Si, X, r, sigma, t23);
            gprime = 1.0 + BSMOption.option_price_delta_put_black_scholes(Si, X, r, sigma, t23);
            S23_bar = Si;
            Si = Si - g / gprime;
        };

        double S3_bar = S23_bar;
        g = 1;
        while (Math.Abs(g) > ACCURACY)
        {
            g = Si - X + BSMOption.option_price_put_black_scholes(Si, X, r, sigma, t3);
            gprime = 1.0 + BSMOption.option_price_delta_put_black_scholes(Si, X, r, sigma, t3);
            S3_bar = Si;
            Si = Si - g / gprime;
        };

        double P3 = X * Math.Exp(-r * t3) * DStat.NormDist(-d2(S, S3_bar, r, sigma, t3));
        P3 -= S * DStat.NormDist(-d1(S, S3_bar, r, sigma, t3));
        P3 += X * Math.Exp(-r * time) * DStat.NormDist(d2(S, S3_bar, r, sigma, t3),
            -d2(S, S23_bar, r, sigma, t23), -rho12);
        P3 -= S * DStat.NormDist(d1(S, S3_bar, r, sigma, t3),
            -d1(S, S23_bar, r, sigma, t23), -rho12);
        P3 += X * Math.Exp(-r * t23) * N3(d1(S, S3_bar, r, sigma, t3),
            d1(S, S23_bar, r, sigma, t23), -d1(S, X, r, sigma, time), 
            rho12, -rho13, -rho23);
        P3 -= S * N3(d2(S, S3_bar, r, sigma, t3), d2(S, S23_bar, r, sigma, t23), 
            -d2(S, X, r, sigma, time), rho12, -rho13, -rho23);
        P3 = Math.Max(P2, P3); // for safety, use one less step as lower bound

        return P3 + 3.5 * (P3 - P2) - 0.5 * (P2 - P1);
    }
    */

    //public static double option_price_american_call_approximated_bjerksund_stensland(
    //    double S, double X, double r, double q, double sigma, double time) { }

    //public static double option_price_american_put_approximated_bjerksund_stensland(
    //    double S, double X, double r, double q, double sigma, double time ) { }        

}
