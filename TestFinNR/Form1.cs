using System;
using System.Collections.Generic;
using System.ComponentModel;
using System.Data;
using System.Drawing;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows.Forms;

using DFinNR;

namespace TestFinNR
{
    public partial class Form1 : Form
    {
        public Form1()
        {
            InitializeComponent();
        }

        private void button1_Click(object sender, EventArgs e)
        {
            Application.Exit();
        }

        private void button2_Click(object sender, EventArgs e)
        {
            double S = double.Parse(textBox1.Text);
            double K = double.Parse(textBox2.Text);
            double time = double.Parse(textBox3.Text);
            double sigma = double.Parse(textBox4.Text);
            double r = double.Parse(textBox5.Text);

            int no_S_steps = int.Parse(textBox6.Text);
            int no_t_steps = int.Parse(textBox7.Text);

            double bs_price = BSMOption.option_price_call_black_scholes(S, K, r, sigma, time);
            textBox8.Text = bs_price.ToString("F6");

            double fd_exp_price = FiniteDifference.option_price_call_european_finite_diff_explicit(
                    S, K, r, sigma, time, no_S_steps, no_t_steps);
            textBox9.Text = fd_exp_price.ToString("F6");

            double fd_imp_price = FiniteDifference.option_price_call_european_finite_diff_implicit(
                    S, K, r, sigma, time, no_S_steps, no_t_steps);
            textBox10.Text = fd_imp_price.ToString("F6");
        }
    }
}
