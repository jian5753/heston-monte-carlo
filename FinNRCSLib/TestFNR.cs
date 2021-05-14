using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace QuantLibNet.FinNR
{
    public class TestFNR
    {
        public static Dictionary<String, Double> Output = new Dictionary<string, double>();

        public static void test_implicit_finite_differences_using_newmat()
        {
            double S = 50.0;
            double K = 50.0;
            double r = 0.1;
            double sigma = 0.4;
            double time = 0.5;
            int no_S_steps = 20;
            int no_t_steps = 20;

            double bs_price = BSMOption.option_price_put_black_scholes(S, K, r, sigma, time);
            Console.WriteLine(" black scholes put price = " + bs_price.ToString());

            double fd_exp_price = FiniteDifference.option_price_put_european_finite_diff_explicit(
                    S, K, r, sigma, time, no_S_steps, no_t_steps);
            Console.WriteLine(" explicit Euro put price = " + fd_exp_price.ToString());

            double fd_imp_price = FiniteDifference.option_price_put_european_finite_diff_implicit(
                    S, K, r, sigma, time, no_S_steps, no_t_steps);
            Console.WriteLine(" implicit Euro put price = " + fd_imp_price.ToString());

            Output.Add("BS_Price", bs_price);
            Output.Add("FD_Exp_Price", fd_exp_price);
            Output.Add("FD_Imp_Price", fd_imp_price);
        }

        public static void examples_finite_diffs_using_newmat()
        {
            Console.WriteLine("----------------------------");
            Console.WriteLine("Finite Differences examples using newmat ");
            Console.WriteLine("----------------------------");
            test_implicit_finite_differences_using_newmat();
        }
    }
}
