using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace DFinNR
{
    public static class term_structure_utils
    {
        public static double term_structure_yield_from_discount_factor(double d_t, double t)
        {
            return (-Math.Log(d_t) / t);
        }

        public static double term_structure_discount_factor_from_yield(double r, double t)
        {
            return Math.Exp(-r * t);
        }

        public static double term_structure_forward_rate_from_discount_factors(double d_t1,
            double d_t2, double time)
        {
            return (Math.Log(d_t1 / d_t2)) / time;
        }

        public static double term_structure_forward_rate_from_yields(double r_t1, double r_t2,
            double t1, double t2)
        {
            return r_t2 * t2 / (t2 - t1) - r_t1 * t1 / (t2 - t1);
        }

        public static double term_structure_yield_linearly_interpolated(double time,
            Vector obs_times, Vector obs_yields)
        {
            // assume the yields are in increasing time to maturity order.
            int no_obs = obs_times.size();
            if (no_obs < 1)
                return 0;

            double t_min = obs_times[0];
            if (time <= t_min)
                return obs_yields[0];  // earlier than lowest obs.

            double t_max = obs_times[no_obs - 1];
            if (time >= t_max)
                return obs_yields[no_obs - 1]; // later than latest obs

            int t = 1;  // find which two observations we are between
            while ((t < no_obs) && (time > obs_times[t]))
            { ++t; };

            double lambda = (obs_times[t] - time) / (obs_times[t] - obs_times[t - 1]);
            // by ordering assumption, time is  between t-1,t
            double r = obs_yields[t - 1] * lambda + obs_yields[t] * (1.0 - lambda);

            return r;
        }

        public static double bonds_price(Vector cashflow_times, Vector cashflows,
            term_structure_class d)
        {
            double p = 0;
            for (int i = 0; i < cashflow_times.size(); i++)
            {
                p += d.d(cashflow_times[i]) * cashflows[i];
            };
            return p;
        }

        public static double bonds_duration(Vector cashflow_times, Vector cashflow_amounts,
            term_structure_class d)
        {
            double S = 0;
            double D1 = 0;
            for (int i = 0; i < cashflow_times.size(); i++)
            {
                S += cashflow_amounts[i] * d.d(cashflow_times[i]);
                D1 += cashflow_times[i] * cashflow_amounts[i] * d.d(cashflow_times[i]);
            };
            return D1 / S;
        }

        public static double bonds_convexity(Vector cashflow_times, Vector cashflow_amounts,
            term_structure_class d)
        {
            double B = 0;
            double Cx = 0;
            for (int i = 0; i < cashflow_times.size(); i++)
            {
                B += cashflow_amounts[i] * d.d(cashflow_times[i]);
                Cx += Math.Pow(cashflow_times[i], 2) * cashflow_amounts[i] * d.d(cashflow_times[i]);
            };
            return Cx / B;
        }


        ////////////////////////////////////////////////////////////////////////////////
        // term structure models formulas for calculation

        public static double term_structure_yield_nelson_siegel(double t,
            double beta0, double beta1, double beta2, double lambda)
        {
            if (t == 0.0) return beta0;
            double tl = t / lambda;
            double r = beta0 + (beta1 + beta2) * ((1 - Math.Exp(-tl)) / tl) + beta2 * Math.Exp(-tl);
            return r;
        }

        public static double term_structure_yield_svensson( double t,
            double beta0,  double beta1,  double beta2,  double beta3,
            double tau1,   double tau2 )
        {
            if (t == 0.0) return beta0;
            double r = beta0;
            r += beta1 * ((1 - Math.Exp(-t / tau1)) / (t / tau1));
            r += beta2 * (((1 - Math.Exp(-t / tau1)) / (t / tau1)) - Math.Exp(-t / tau1));
            r += beta3 * (((1 - Math.Exp(-t / tau2)) / (t / tau2)) - Math.Exp(-t / tau2));
            return r;
        }

        public static double term_structure_discount_factor_cubic_spline( double t,
            double b1, double c1, double d1, Vector f, Vector knots)
        {
            double d = 1.0 + b1 * t + c1 * (Math.Pow(t, 2)) + d1 * (Math.Pow(t, 3));
            for (int i = 0; i < knots.size(); i++)
            {
                if (t >= knots[i]) { d += f[i] * (Math.Pow((t - knots[i]), 3)); }
                else { break; };
            };
            return d;
        }

        public static double term_structure_discount_factor_cir( double t, 
            double r,  double kappa,   double lambda,  double theta,  double sigma)
        {
            double sigma_sqr = Math.Pow(sigma, 2);
            double gamma = Math.Sqrt(Math.Pow((kappa + lambda), 2) + 2.0 * sigma_sqr);
            double denum = (gamma + kappa + lambda) * (Math.Exp(gamma * t) - 1) + 2 * gamma;
            double p = 2 * kappa * theta / sigma_sqr;
            double enum1 = 2 * gamma * Math.Exp(0.5 * (kappa + lambda + gamma) * t);
            double A = Math.Pow((enum1 / denum), p);
            double B = (2 * (Math.Exp(gamma * t) - 1)) / denum;
            double dfact = A * Math.Exp(-B * r);
            return dfact;
        }

        public static double term_structure_discount_factor_vasicek( double time,
            double r,  double a, double b,  double sigma)
        {
            double A, B;
            double sigma_sqr = sigma * sigma;
            double aa = a * a;
            if (a == 0.0)
            {
                B = time;
                A = Math.Exp(sigma_sqr * Math.Pow(time, 3)) / 6.0;
            }
            else
            {
                B = (1.0 - Math.Exp(-a * time)) / a;
                A = Math.Exp(((B - time) * (aa * b - 0.5 * sigma_sqr)) / aa 
                    - ((sigma_sqr * B * B) / (4 * a)));
            };
            double dfact = A * Math.Exp(-B * r);
            return dfact;
        }

        public static List<List<term_structure_class_ho_lee>>
            term_structure_ho_lee_build_term_structure_tree(
            term_structure_class initial, int no_steps, double delta, double pi)
        {
            List<List<term_structure_class_ho_lee>> hl_tree = 
                new List<List<term_structure_class_ho_lee>>();
            for (int t = 0; t < 5; ++t)
            {
                hl_tree.Add(new List<term_structure_class_ho_lee>());
                for (int j = 0; j <= t; ++j)
                {
                    term_structure_class_ho_lee hl = 
                        new term_structure_class_ho_lee(initial, t, j, delta, pi);
                    hl_tree[t].Add(hl);
                };
            };
            return hl_tree;
        }
        
    }

    public class term_structure_class
    {
        public term_structure_class() { }

        public virtual double r(double t) // yield on zero coupon bond
        {
            return term_structure_utils.term_structure_yield_from_discount_factor(d(t), t);
        }

        public virtual double d(double t) // discount factor/price of zero coupon bond
        {
            return term_structure_utils.term_structure_discount_factor_from_yield(r(t), t);
        }

        public virtual double f(double t1, double t2) // forward rate
        {
            double d1 = d(t1);
            double d2 = d(t2);
            return term_structure_utils.term_structure_forward_rate_from_discount_factors(
                d1, d2, t2 - t1);
        }
    }

    public class term_structure_class_flat : term_structure_class
    {
        private double R_;   // interest rate

        public term_structure_class_flat(double r)
        {
            R_ = r;
        }

        public override double r(double T)
        {
            if (T >= 0) return R_;
            return 0;
        }

        public void set_int_rate(double r)
        {
            R_ = r;
        }
    }

    public class term_structure_class_interpolated : term_structure_class
    {
        private Vector times_;     // use to keep a list of yields
        private Vector yields_;
        private void clear()
        {
            times_.Clear();
            yields_.Clear();
        }

        public term_structure_class_interpolated() : base()
        {
            clear();
        }

        public term_structure_class_interpolated(Vector in_times, Vector in_yields)
        {
            clear();
            if (in_times.size() != in_yields.size()) return;
            times_ = new Vector(in_times.size());
            yields_ = new Vector(in_yields.size());
            for (int i = 0; i < in_times.size(); i++)
            {
                times_[i] = in_times[i];
                yields_[i] = in_yields[i];
            };
        }

        public term_structure_class_interpolated(term_structure_class_interpolated term)
        {
            times_ = new Vector(term.no_observations());
            yields_ = new Vector(term.no_observations());
            for (int i = 0; i < term.no_observations(); i++)
            {
                times_[i] = term.times_[i];
                yields_[i] = term.yields_[i];
            };
        }

        public int no_observations()
        {
            return times_.size();
        }

        public override double r(double T)
        {
            return term_structure_utils.term_structure_yield_linearly_interpolated(T, times_, yields_);
        }

        public void set_interpolated_observations(Vector in_times, Vector in_yields)
        {
            clear();
            if (in_times.size() != in_yields.size()) return;
            times_ = new Vector(in_times.size());
            yields_ = new Vector(in_yields.size());
            for (int i = 0; i < in_times.size(); i++)
            {
                times_[i] = in_times[i];
                yields_[i] = in_yields[i];
            };
        }
    }

    public class term_structure_class_nelson_siegel : term_structure_class
    {
        private double beta0_, beta1_, beta2_, lambda_;
        public term_structure_class_nelson_siegel(double beta0, double beta1,
            double beta2, double lambda)
        {
            beta0_ = beta0;
            beta1_ = beta1;
            beta2_ = beta2;
            lambda_ = lambda;
        }

        public override double r(double T)
        {
            if (T <= 0.0)
                return beta0_;

            return term_structure_utils.term_structure_yield_nelson_siegel(
                T, beta0_, beta1_, beta2_, lambda_);
        }
    }

    public class term_structure_class_svensson : term_structure_class
    {
        private double beta0_, beta1_, beta2_, beta3_;
        private double tau1_, tau2_;

        public term_structure_class_svensson(double beta0, double beta1, 
            double beta2, double beta3, double tau1, double tau2)
        {
            beta0_ = beta0;
            beta1_ = beta1;
            beta2_ = beta2;
            beta3_ = beta3;
            tau1_ = tau1;
            tau2_ = tau2;
        }

        public override double r(double T)
        {
            if (T <= 0.0)
                return beta0_;
            return term_structure_utils.term_structure_yield_svensson(
                T, beta0_, beta1_, beta2_, beta3_, tau1_, tau2_);
        }
    }

    public class term_structure_class_cubic_spline : term_structure_class
    {
        private double b_, c_, d_;
        private Vector f_ = new Vector();
        private Vector knots_ = new Vector();

        public term_structure_class_cubic_spline(double b, double c,
            double d, Vector f, Vector knots)
        {
            b_ = b;
            c_ = c;
            d_ = d;
            f_.Clear();
            knots_.Clear();
            if (f.size() != knots.size())
            {
                Utils.QL_Require(false, "Size not matched!");
            };

            for (int i = 0; i < f.size(); ++i)
            {
                f_.Add(f[i]);
                knots_.Add(knots[i]);
            };

        }
        public override double d(double T) // discount factor
        {
            return term_structure_utils.
                term_structure_discount_factor_cubic_spline(T, b_, c_, d_, f_, knots_);
        }
    }

    public class term_structure_class_cir : term_structure_class
    {
        private double r_, kappa_, lambda_, theta_, sigma_;

        public term_structure_class_cir(double r, double k, double l,
            double th, double sigma)
        {
            r_ = r;
            kappa_ = k;
            lambda_ = l;
            theta_ = th;
            sigma_ = sigma;
        }
        public override double d(double T)  // discount factor
        {
            return term_structure_utils.term_structure_discount_factor_cir(
                T, r_, kappa_, lambda_, theta_, sigma_);
        }
    }

    public class term_structure_class_vasicek : term_structure_class
    {
        private double r_, a_, b_, sigma_;

        public term_structure_class_vasicek(double r, double a, 
            double b, double sigma)
        {
            r_ = r;
            a_ = a;
            b_ = b;
            sigma_ = sigma;
        }
        public override double d(double T)
        {
            return term_structure_utils.term_structure_discount_factor_vasicek(
                T, r_, a_, b_, sigma_);
        }
    }

    public class term_structure_class_ho_lee : term_structure_class
    {
        private term_structure_class initial_term_;
        private int n_;
        private int i_;
        private double delta_;
        private double pi_;

        public term_structure_class_ho_lee(term_structure_class fitted_term,
            int n, int i, double delta, double pi)
        {
            initial_term_ = fitted_term;
            n_ = n;
            i_ = i;
            delta_ = delta;
            pi_ = pi;
        }

        private double hT(double T, double delta, double pi)
        {
            return (1.0 / (pi + (1 - pi) * Math.Pow(delta, T)));
        }

        public override double d(double T)
        {
            double d = initial_term_.d(T + n_) / initial_term_.d(n_);

            for (int j = 1; j < n_; ++j)
            {
                d *= hT(T + (n_ - j), delta_, pi_) / hT(n_ - j, delta_, pi_);
            };

            d *= hT(T, delta_, pi_) * Math.Pow(delta_, T * (n_ - i_));

            return d;
        }
    }

    public class time_contingent_cash_flows
    {
        public Vector times;
        public Vector cash_flows;
        public time_contingent_cash_flows(Vector in_times, Vector in_cflows)
        {
	        times=in_times;	cash_flows=in_cflows;
        }
        public int no_cflows() { return times.size();
    }
}

}
