using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Xaml;
using System.Windows;
using OxyPlot;
using System.Windows.Controls;
using System.Windows.Forms;

namespace Aberdeen
{
    class Program
    {
        
        static Random rand = new Random(); // random number class

        enum DiscretisationScheme
        {
            Euler,
            Milstein
        }

        static double CaculateUnderlyingValueAtExpiry(double S0, double mu, double sigma, double timeToExpiry, int numberOfSteps, DiscretisationScheme method)
        {
            // Calculate value of underlying S at expiry, given S0
            // S0 - initial price
            // method - Euler or Milstein

            double dt = timeToExpiry / numberOfSteps;
            double S = S0;

            for (int i = 0; i < numberOfSteps; i++)
            {
                double r1 = rand.NextDouble(); // uniform random number (0, 1)
                double r2 = rand.NextDouble();
                double Z = Math.Sqrt(-2 * Math.Log(r1)) * Math.Sin(2 * Math.PI * r2); // standard normal distribution (0, 1)
                
                S += mu * S * dt + sigma * S * Z * Math.Sqrt(dt);

                if (method == DiscretisationScheme.Milstein)
                {
                    S += 0.5 * sigma * sigma * dt * (Z * Z - 1);
                }
            }
            
            return S;
        }

        static double CalculateCallOptionValue(double S0, double mu, double sigma,
            double timeToExpiry, int numberOfSteps, double strike, DiscretisationScheme method)
        {
            // Calculate European call option value

            double expectedUnderlyingValue = CaculateUnderlyingValueAtExpiry(S0, mu, sigma, timeToExpiry, numberOfSteps, method);
            if (expectedUnderlyingValue > strike)
            {
                return expectedUnderlyingValue - strike;
            }

            return 0;
        }

        static void CalculateCallOptionValue(double S0, double mu, double sigma,
            double timeToExpiry, int numberOfSteps, int numberOfSimulations, double strike, DiscretisationScheme method,
            out double optionValue, out double optionStandardDeviation)
        {
            // Simulation to calculate European call option value
            double[] optionValues = new double[numberOfSimulations];
            for (int i = 0; i < numberOfSimulations; i++)
            {
                optionValues[i] = CalculateCallOptionValue(S0, mu, sigma, timeToExpiry, numberOfSteps, strike, method);
            }

            double average = optionValues.Average();
            double sumOfSquaresOfDifferences = optionValues.Select(val => (val - average) * (val - average)).Sum();
            optionValue = average;
            optionStandardDeviation = Math.Sqrt(sumOfSquaresOfDifferences / numberOfSimulations);
        }

        static void Plot(double[] x, double[] y)
        {
            // Plot (x, y) on a chart

            OxyPlot.WindowsForms.Plot plot = new OxyPlot.WindowsForms.Plot();
            plot.Dock = DockStyle.Fill;
            PlotModel plotModel = new PlotModel{ Title = "Call option price" };
            plot.Model = plotModel;
            plotModel.PlotType = PlotType.XY;

            OxyPlot.Series.LineSeries line = new OxyPlot.Series.LineSeries();

            for (int i = 0; i < x.Length; i++)
            {
                line.Points.Add(new DataPoint(x[i], y[i]));
            }
            plotModel.Series.Add(line);

            Form form = new Form();
            form.FormBorderStyle = FormBorderStyle.None;
            PictureBox pb = new PictureBox();
            pb.ImageLocation = "price.png";
            pb.SizeMode = PictureBoxSizeMode.StretchImage;
            pb.SizeMode = PictureBoxSizeMode.AutoSize;
            form.Controls.Add(plot);
            form.ShowDialog();
        }
        

        [STAThread]
        static void Main(string[] args)
        {
            int numberOfPoints = 100; // number of points to plot

            double[] S = new double[numberOfPoints];
            double[] P = new double[numberOfPoints];
            double[] errors = new double[numberOfPoints];

            double riskFreeRate = 0.05; // Risk free rate
            double strike = 100; // Strike price
            double timeToExpiry = 1; // Time to expiry (years)

            // Underlying
            double mu = 0.2; // Expected growth
            double sigma = 0.1; // standard deviation
            
            // Simulation parameters
            int numberOfSteps = 100;
            int numberOfSimulations = 100;


            Console.WriteLine("Underlying price", "Option price", "Standard deviation (error)");
            for (int i = 0; i < numberOfPoints; i++)
            {
                S[i] = strike / numberOfPoints * 2 * i;
                CalculateCallOptionValue(S[i], mu, sigma, timeToExpiry, numberOfSteps, numberOfSimulations, strike, DiscretisationScheme.Euler, out P[i], out errors[i]);
                P[i] = Math.Exp(-riskFreeRate * timeToExpiry) * P[i]; // discount to today's value
                Console.WriteLine(String.Format("{0}, {1}, {2}", S[i], P[i], errors[i]));
            }
            
            Plot(S, P);
        }
    }
}
