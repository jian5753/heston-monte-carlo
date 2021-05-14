using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace DFinNR
{
    public static class PresentValue
    {
        ///////// present value ////////////////////////////////////
        // discrete coumpounding
        /////////////////////////////////
        // discrete, annual compounding

        public static double cash_flow_pv_discrete(Vector cflow_times, 
            Vector cflow_amounts, double r)
        {
            double PV = 0.0;
            for (int t = 0; t < cflow_times.size(); t++)
            {
                PV += cflow_amounts[t] / Math.Pow(1.0 + r, cflow_times[t]);
            };
            return PV;
        }

        public static double cash_flow_irr_discrete(Vector cflow_times, 
            Vector cflow_amounts)
        {
            // simple minded irr function.  Will find one root (if it exists.)
            // adapted from routine in Numerical Recipes in C.
            if (cflow_times.size() != cflow_amounts.size())
                Utils.QL_Require(false, "Times size != Cash Flow size");

            const double ACCURACY = 1.0e-5;
            const int MAX_ITERATIONS = 50;
            double x1 = 0.0;
            double x2 = 0.2;

            // create an initial bracket, with a root somewhere between bot,top
            double f1 = cash_flow_pv_discrete(cflow_times, cflow_amounts, x1);
            double f2 = cash_flow_pv_discrete(cflow_times, cflow_amounts, x2);
            int i;
            for (i = 0; i < MAX_ITERATIONS; i++)
            {
                if ((f1 * f2) < 0.0) { break; }; // 
                if (Math.Abs(f1) < Math.Abs(f2))
                {
                    f1 = cash_flow_pv_discrete(cflow_times, cflow_amounts, x1 += 1.6 * (x1 - x2));
                }
                else
                {
                    f2 = cash_flow_pv_discrete(cflow_times, cflow_amounts, x2 += 1.6 * (x2 - x1));
                };
            };
            if (f2 * f1 > 0.0)
            {
                Utils.QL_Require(false, "Error");
            };
            double f = cash_flow_pv_discrete(cflow_times, cflow_amounts, x1);
            double rtb;
            double dx = 0;
            if (f < 0.0)
            {
                rtb = x1;
                dx = x2 - x1;
            }
            else
            {
                rtb = x2;
                dx = x1 - x2;
            };
            for (i = 0; i < MAX_ITERATIONS; i++)
            {
                dx *= 0.5;
                double x_mid = rtb + dx;
                double f_mid = cash_flow_pv_discrete(cflow_times, cflow_amounts, x_mid);
                if (f_mid <= 0.0) { rtb = x_mid; }
                if ((Math.Abs(f_mid) < ACCURACY) || (Math.Abs(dx) < ACCURACY)) return x_mid;
            };

            Utils.QL_Require(false, "Error");
            return 0;
        }

        public static int sgn(double r)
        {
            if (r >= 0)
                return 1;
            else
                return -1;
        }

        public static bool cash_flow_unique_irr(Vector cflow_times, 
            Vector cflow_amounts)
        {
            int sign_changes = 0;     // first check Descartes rule
            for (int t = 1; t < cflow_times.size(); ++t)
            {
                if (sgn(cflow_amounts[t - 1]) != sgn(cflow_amounts[t])) sign_changes++;
            };
            if (sign_changes == 0) return false;  // can not find any irr
            if (sign_changes == 1) return true;

            double A = cflow_amounts[0]; // check the aggregate cash flows, due to Norstrom
            sign_changes = 0;
            for (int t = 1; t < cflow_times.size(); ++t)
            {
                if (sgn(A) != sgn(A += cflow_amounts[t])) sign_changes++;
            };
            if (sign_changes <= 1) return true;
            return false;
        }


        public static double bonds_price_discrete(Vector times, 
            Vector cashflows, double r)
        {
            double p = 0;
            for (int i = 0; i < times.size(); i++)
            {
                p += cashflows[i] / (Math.Pow((1 + r), times[i]));
            };
            return p;
        }

        public static double bonds_yield_to_maturity_discrete(Vector times, 
            Vector cashflows, double bondprice)
        {
            const double ACCURACY = 1e-5;
            const int MAX_ITERATIONS = 200;
            double bot = 0, top = 1.0;
            while (bonds_price_discrete(times, cashflows, top) > bondprice) { top = top * 2; };
            double r = 0.5 * (top + bot);
            for (int i = 0; i < MAX_ITERATIONS; i++)
            {
                double diff = bonds_price_discrete(times, cashflows, r) - bondprice;
                if (Math.Abs(diff) < ACCURACY) return r;
                if (diff > 0.0) { bot = r; }
                else { top = r; };
                r = 0.5 * (top + bot);
            };
            return r;
        }

        public static double bonds_duration_discrete(Vector times, 
            Vector cashflows, double r)
        {
            double B = 0;
            double D = 0;
            for (int i = 0; i < times.size(); i++)
            {
                D += times[i] * cashflows[i] / Math.Pow(1 + r, times[i]);
                B += cashflows[i] / Math.Pow(1 + r, times[i]);
            };
            return D / B;
        }

        public static double bonds_duration_macaulay_discrete(Vector times, 
            Vector cashflows, double bond_price)
        {
            double y = bonds_yield_to_maturity_discrete(times, cashflows, bond_price);
            return bonds_duration_discrete(times, cashflows, y); // use YTM in duration calculation
        }

        public static double bonds_duration_modified_discrete(Vector times, 
            Vector cashflows, double bond_price)
        {
            double y = bonds_yield_to_maturity_discrete(times, cashflows, bond_price);
            double D = bonds_duration_discrete(times, cashflows, y);
            return D / (1 + y);
        }

        public static double bonds_convexity_discrete(Vector times, 
            Vector cashflows, double r)
        {
            double Cx = 0;
            for (int i = 0; i < times.size(); i++)
            {
                Cx += cashflows[i] * times[i] * (times[i] + 1) / (Math.Pow((1 + r), times[i]));
            };
            double B = bonds_price_discrete(times, cashflows, r);
            return (Cx / (Math.Pow(1 + r, 2))) / B;
        }

        /////////////////////////////////
        // continous compounding. 
        public static double cash_flow_pv(Vector cflow_times, 
            Vector cflow_amounts, double r)
        {
            double PV = 0.0;
            for (int t = 0; t < cflow_times.size(); t++)
            {
                PV += cflow_amounts[t] * Math.Exp(-r * cflow_times[t]);
            };
            return PV;
        }

        public static double cash_flow_irr(Vector cflow_times, 
            Vector cflow_amounts)
        {
            // simple minded irr function.  Will find one root (if it exists.)
            // adapted from routine in Numerical Recipes in C.
            if (cflow_times.size() != cflow_amounts.size())            
                Utils.QL_Require(false, "Error");

            const double ACCURACY = 1.0e-5;
            const int MAX_ITERATIONS = 50;
            double x1 = 0.0;
            double x2 = 0.2;

            // create an initial bracket, with a root somewhere between bot,top
            double f1 = cash_flow_pv(cflow_times, cflow_amounts, x1);
            double f2 = cash_flow_pv(cflow_times, cflow_amounts, x2);
            int i;
            for (i = 0; i < MAX_ITERATIONS; i++)
            {
                if ((f1 * f2) < 0.0) { break; }; // 
                if (Math.Abs(f1) < Math.Abs(f2))
                {
                    f1 = cash_flow_pv(cflow_times, cflow_amounts, x1 += 1.6 * (x1 - x2));
                }
                else
                {
                    f2 = cash_flow_pv(cflow_times, cflow_amounts, x2 += 1.6 * (x2 - x1));
                };
            };
            if (f2 * f1 > 0.0)
                Utils.QL_Require(false, "Error");

            double f = cash_flow_pv(cflow_times, cflow_amounts, x1);
            double rtb;
            double dx = 0;
            if (f < 0.0)
            {
                rtb = x1;
                dx = x2 - x1;
            }
            else
            {
                rtb = x2;
                dx = x1 - x2;
            };
            for (i = 0; i < MAX_ITERATIONS; i++)
            {
                dx *= 0.5;
                double x_mid = rtb + dx;
                double f_mid = cash_flow_pv(cflow_times, cflow_amounts, x_mid);
                if (f_mid <= 0.0) { rtb = x_mid; }
                if ((Math.Abs(f_mid) < ACCURACY) || (Math.Abs(dx) < ACCURACY)) return x_mid;
            };

            Utils.QL_Require(false, "Error");
            return 1.0;   // error.
        }

        public static double bonds_price(Vector cashflow_times, 
            Vector cashflows, double r)
        {
            double p = 0;
            for (int i = 0; i < cashflow_times.size(); i++)
            {
                p += Math.Exp(-r * cashflow_times[i]) * cashflows[i];
            };
            return p;
        }

        public static double bonds_price(Vector coupon_times, Vector coupon_amounts,
               Vector principal_times, Vector principal_amounts, double r)
        {
            double p = 0;
            for (int i = 0; i < coupon_times.size(); i++)
            {
                p += Math.Exp(-r * coupon_times[i]) * coupon_amounts[i];
            };
            for (int i = 0; i < principal_times.size(); i++)
            {
                p += Math.Exp(-r * principal_times[i]) * principal_amounts[i];
            };
            return p;
        }

        public static double bonds_duration(Vector cashflow_times, 
            Vector cashflows, double r)
        {
            double S = 0;
            double D1 = 0;
            for (int i = 0; i < cashflow_times.size(); i++)
            {
                S += cashflows[i] * Math.Exp(-r * cashflow_times[i]);
                D1 += cashflow_times[i] * cashflows[i] * Math.Exp(-r * cashflow_times[i]);
            };
            return D1 / S;
        }

        public static double bonds_yield_to_maturity(Vector cashflow_times, 
            Vector cashflow_amounts, double bondprice)
        {
            const double ACCURACY = 1e-5;
            const int MAX_ITERATIONS = 200;
            double bot = 0, top = 1.0;
            while (bonds_price(cashflow_times, cashflow_amounts, top) > bondprice)
            {
                top = top * 2;
            };
            double r = 0.5 * (top + bot);
            for (int i = 0; i < MAX_ITERATIONS; i++)
            {
                double diff = bonds_price(cashflow_times, cashflow_amounts, r) - bondprice;
                if (Math.Abs(diff) < ACCURACY) return r;
                if (diff > 0.0) { bot = r; }
                else { top = r; };
                r = 0.5 * (top + bot);
            };
            return r;
        }

        public static double bonds_duration_macaulay(Vector cashflow_times, 
            Vector cashflows, double bond_price)
        {
            double y = bonds_yield_to_maturity(cashflow_times, cashflows, bond_price);
            return bonds_duration(cashflow_times, cashflows, y); // use YTM in duration 
        }

        public static double bonds_convexity(Vector times, Vector cashflows, 
            double r)
        {
            double C = 0;
            for (int i = 0; i < times.size(); i++)
            {
                C += cashflows[i] * Math.Pow(times[i], 2) * Math.Exp(-r * times[i]);
            };
            double B = bonds_price(times, cashflows, r);
            return C / B;
        }
    }

}
