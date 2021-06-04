using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Threading;
using DFinNR;

namespace parrallelTest
{
    class NormRv
    {
        public double data;

        public NormRv()
        {
            Random rv = new Random();
            data = DStat.N_Inv(rv.NextDouble());
        }

        public NormRv(Random rv)
        {
            data = DStat.N_Inv(rv.NextDouble());
        }
    }
    class Program
    {
        
        static void Main(string[] args)
        {
            Random rv = new Random();
            double[,,] data = new double[5, 4, 3];
            Parallel.For(0, 5, i =>
            {
                Parallel.For(0, 4, j =>
                {
                    NormRv nr = new NormRv();
                    data[i, j, 0] = i;
                    data[i, j, 1] = j;
                    data[i, j, 2] = nr.data;
                    Console.WriteLine($"{i}, {j}, {nr.data}");
                    Thread.Sleep(500);
                });
            });
            data.ToString();
            Console.ReadKey();
        }
    }
}
