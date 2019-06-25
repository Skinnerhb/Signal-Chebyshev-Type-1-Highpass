using System;
using System.Collections.Generic;
using System.Linq;
using System.Threading.Tasks;
using System.Windows.Forms;
using System.Data;
using System.ComponentModel;
using System.Drawing;
using System.Text;
using System.Numerics;
using MathNet.Numerics;
using System.Windows.Forms.DataVisualization.Charting;
using MathNet.Numerics.IntegralTransforms;
using System.Diagnostics;

namespace EERI_414_Practical_HB_Skinner
{
    public partial class Form1 : Form
    {
        
        //Variables
        static double Ap = 0.25;
        static double As = 40;
        static double Ft = 150000;
        static double Fs = 36000;
        static double Fp = 64000;
        static int SampleSize = 150000;
        static int SampleSizer = 150000;

        static double wp = (2 * Math.PI * Fp) / Ft;
        static double ws = (2 * Math.PI * Fs) / Ft;


        static double Wpp = Math.Tan(wp / 2);
        static double Wss = Math.Tan(ws / 2);

        static double Wp = 1;
        static double Ws = Wpp / Wss;

        static double e = Math.Sqrt(Math.Pow(10, Ap / 10) - 1);

        static double k = Wp / Ws;
        static double d = e / Math.Sqrt(Math.Pow(10, As / 10) - 1);
        static double N = Math.Ceiling((Math.Log((1 / d) + Math.Sqrt((1 / d) - 1))) / (Math.Log((1 / k) + Math.Sqrt((1 / k) - 1))));
        static int N1 = Convert.ToInt32(N);

        static double B = Math.Pow((Math.Sqrt(1 + Math.Pow(e, 2)) + 1) / (e), 1 / N);
        static double R = Wp * ((Math.Pow(B, 2) + 1) / (2 * B));
        static double r = Wp * ((Math.Pow(B, 2) - 1) / (2 * B));

        static void SineWave(double input, double[] output)
        {
            double[] time = new double[SampleSize];
            for (int i = 0; i < SampleSize; i++)
            {
                if (i == 0)
                {

                    time[i] = 0;
                    output[i] = Math.Sin((input) * time[i]);
                }
                else
                {
                    time[i] = time[i - 1] + 1;
                    output[i] = Math.Sin((input) * time[i]);
                }

            }
        }

        static void Th(Complex[] ss, int n, double q, double Q)
        {
            double[] A = new double[n];
            double[] x = new double[n];
            double[] y = new double[n];
            
            for (int i = 0; i < n; i++)
            {
                A[i] = (Math.PI / 2) + ((((2 * i) + 1) * Math.PI) / (2 * n));
                x[i] = q * Math.Cos(A[i]);
                y[i] = Q * Math.Sin(A[i]);
                ss[i] = new Complex(x[i], y[i]);

            }
            
        }

        static void LP(Complex[] ss, double[] input, Complex[] output, double[] mag, double[] phas)
        {
            Complex a = new Complex(0, 0);
            Complex b = new Complex(0, 0);
            Complex c = new Complex(0, 0);
            Complex d = new Complex(0, 0);

            Complex j = new Complex((double)0, (double)1);

            a = ss[0] + ss[1] + ss[2] + ss[3];
            b = ((ss[0] + ss[1]) * (ss[2] + ss[3])) + (ss[0] * ss[1]);
            c = ((ss[0] + ss[1]) * (ss[2] * ss[3])) + ((ss[2] + ss[3]) * (ss[0] * ss[1]));
            d = (ss[0] * ss[1] * ss[2] * ss[3]);

            double b0 = d.Real;
            double k1 = b0 / Math.Sqrt(1 + Math.Pow(e, 2));

            for (int i = 0; i < SampleSize; i++)
            {
                output[i] = k1 / (Complex.Pow((input[i] * j), 4) - (a.Real * (Complex.Pow((input[i] * j), 3))) + (b.Real * (Complex.Pow((input[i] * j), 2))) - (c.Real * (input[i] * j)) + d.Real);
                mag[i] = 20 * Math.Log10(output[i].Magnitude);
                phas[i] = output[i].Phase * (180 / Math.PI);
            }


        }

        public void HP(Complex[] ss, double[] input, Complex[] output, double[] mag, double[] phas)
        {
            Complex a = new Complex(0, 0);
            Complex b = new Complex(0, 0);
            Complex c = new Complex(0, 0);
            Complex d = new Complex(0, 0);

            Complex j = new Complex((double)0, (double)1);

            a = ss[0] + ss[1] + ss[2] + ss[3];
            b = ((ss[0] + ss[1]) * (ss[2] + ss[3])) + (ss[0] * ss[1]);
            c = ((ss[0] + ss[1]) * (ss[2] * ss[3])) + ((ss[2] + ss[3]) * (ss[0] * ss[1]));
            d = (ss[0] * ss[1] * ss[2] * ss[3]);

            double b0 = d.Real;
            double k1 = b0 / Math.Sqrt(1 + Math.Pow(e, 2));

            for (int i = 1; i < SampleSize; i++)
            {
                output[i] = k1 / (Complex.Pow((Wpp / (input[i] * j)), 4) - (a.Real * (Complex.Pow((Wpp / (input[i] * j)), 3))) + (b.Real * (Complex.Pow((Wpp / (input[i] * j)), 2))) - (c.Real * (Wpp / (input[i] * j))) + d.Real);
                //k1 / (Complex.Pow((Wpp / (input[i] * j)), 4) - (a.Real * (Complex.Pow((Wpp / (input[i] * j)), 3))) + (b.Real * (Complex.Pow((Wpp / (input[i] * j)), 2))) - (c.Real * (Wpp/(input[i] * j))) + d.Real);
                //(0.9717 * Complex.Pow(input[i] * j, 4)) / ((Complex.Pow(input[i] * j, 4)) + (11.18 * Complex.Pow(input[i] * j, 3)) + (70.61 * Complex.Pow(input[i] * j, 2)) + (212.8 * (input[i] * j)) + 625.3);
                mag[i] = 20 * Math.Log10(output[i].Magnitude);
                phas[i] = output[i].Phase * (180 / Math.PI);
            }

        }

        static void BI(Complex[] ss, double[] input, Complex[] output, double[] mag, double[] phas)
        {
            Complex a = new Complex(0, 0);
            Complex b = new Complex(0, 0);
            Complex c = new Complex(0, 0);
            Complex d = new Complex(0, 0);

            Complex j = new Complex((double)0, (double)1);

            a = ss[0] + ss[1] + ss[2] + ss[3];
            b = ((ss[0] + ss[1]) * (ss[2] + ss[3])) + (ss[0] * ss[1]);
            c = ((ss[0] + ss[1]) * (ss[2] * ss[3])) + ((ss[2] + ss[3]) * (ss[0] * ss[1]));
            d = (ss[0] * ss[1] * ss[2] * ss[3]);
            

            double b0 = d.Real;
            double k1 = b0 / Math.Sqrt(1 + Math.Pow(e, 2));

            for (int i = 0; i < SampleSize; i++)
            {
                output[i] = k1 / (Complex.Pow(Wpp / ((1 - (Complex.Exp(input[i] * j))) / (1 + (Complex.Exp(input[i] * j)))), 4) - (a.Real * (Complex.Pow(Wpp / ((1 - (Complex.Exp(input[i] * j))) / (1 + (Complex.Exp(input[i] * j)))), 3))) + (b.Real * (Complex.Pow(Wpp / ((1 - (Complex.Exp(input[i] * j))) / (1 + (Complex.Exp(input[i] * j)))), 2))) - (c.Real * (Wpp / ((1 - (Complex.Exp(input[i] * j))) / (1 + (Complex.Exp(input[i] * j)))))) + d.Real);
                mag[i] = 20 * Math.Log10(output[i].Magnitude);
                phas[i] = output[i].Phase * (180 / Math.PI);
            }
        }
        

        static void Grey(Complex[] inputs, double[] z,double[] output)
        {
            double[] p = { 0.001055, -0.004221, 0.006331, -0.004221, 0.001055 };
            double[] d = { 1, 3.15, 3.927, 2.274, 0.5135 };

            double d1 = (d[0] - (d[4]*d[3] * d[2] * d[1])) / (1 - (Math.Pow(d[4], 2)));
            double d2 = (d[1] - (d[4]*d[3] * d[2] * d[0])) / (1 - (Math.Pow(d[4], 2)));
            double d3 = (d[2] - (d[4]*d[3] * d[0] * d[1])) / (1 - (Math.Pow(d[4], 2)));

            double d11 = (d1 - (d3 * d2)) / (1 - (Math.Pow(d3, 2)));
            double d22 = (d2 - (d3 * d1)) / (1 - (Math.Pow(d3, 2)));

            double d111 = (d11) / (1 - (Math.Pow(d22, 2)));

            double a1 = p[4];
            double a2 = p[3] - (a1 * d[0]);
            double a3 = p[2] - (a1 * d[1]) - (a2 * d1);
            double a4 = p[1] - (a1 * d[2]) - (a2 * d2) - (a3 * d11);
            double a5 = p[0] - (a1 * d[3]) - (a2 * d3) - (a3 * d22) - (a4 * d111);

            double[] X1 = new double[SampleSize];
            double[] X2 = new double[SampleSize];
            double[] X3 = new double[SampleSize];
            double[] X4 = new double[SampleSize];
            double[] X5 = new double[SampleSize];
            double[] X6 = new double[SampleSize];
            double[] X7 = new double[SampleSize];
            double[] X8 = new double[SampleSize];
            double[] X9 = new double[SampleSize];
            double[] X10 = new double[SampleSize];
            double[] X11 = new double[SampleSize];
            double[] X12 = new double[SampleSize];


            for (int n = 1; n < SampleSize; n++)
            {
                X2[n] = (-d[4]*z[n]*X8[n]) + inputs[n].Real;
                X3[n] = (-d3 * z[n]*X7[n]) + X2[n];
                X4[n] = (-d22 * z[n]*X6[n]) + X3[n];
                X5[n] = (-d111 * z[n]*X5[n]) + X4[n];
                X6[n] = (d111 * X5[n]) + (X5[n]*z[n]);
                X7[n] = (d22 * X4[n]) + (X6[n] * z[n]);
                X8[n] = (d3 * X3[n]) + (X7[n] * z[n]);
                X9[n] = (d[4] * X2[n]) + (X8[n] * z[n]);
                X10[n] = (a1 * X9[n]) + (a2 * X8[n]);
                X11[n] = X10[n] + (a3 * X7[n]);
                X12[n] = X11[n] + (a4 * X6[n]);

                output[n] = X12[n];
            }
            
        }

        public Form1()
        {
            InitializeComponent();

            textBox1.Text = 0.ToString();
            textBox2.Text = 150000.ToString();
        }

        private void chart1_Click(object sender, EventArgs e)
        {

        }

        private void textBox1_TextChanged(object sender, EventArgs e)
        {

        }

        private void textBox2_TextChanged(object sender, EventArgs e)
        {

        }

        private void button1_Click(object sender, EventArgs e)
        {
            chart1.Titles.Clear();
            chart2.Titles.Clear();
            chart3.Titles.Clear();
            var input1 = textBox1.Text;
            var input2 = textBox2.Text;
            double inputmin = 0;
            double min = 0;
            double inputmax = 0; ;
            double max = 0;
            double dinput1 = 0;
            double dinput2 = 0;
            Complex[] output = new Complex[SampleSize];
            double[] phas = new double[SampleSize];
            double[] mag = new double[SampleSize];
            double[] input = new double[SampleSize];
            double[] inputt = new double[SampleSize];
            double[] inputs = new double[SampleSize];
            chart1.Series[0].Points.Clear();
            chart1.Titles.Add("Magnitude");
            chart1.ChartAreas[0].AxisX.Title = "Frequency (Hz)";
            chart1.ChartAreas[0].AxisY.Title = "Magnitude (dB)";
            chart2.Series[0].Points.Clear();
            chart2.Titles.Add("Phase");
            chart2.ChartAreas[0].AxisX.Title = "Frequency (Hz)";
            chart2.ChartAreas[0].AxisY.Title = "Phase";
            chart3.Series[0].Points.Clear();
            chart3.Titles.Add("Transfer Function");
            chart3.ChartAreas[0].AxisX.Title = "Frequency (Hz)";
            chart3.ChartAreas[0].AxisY.Title = "Amplitude";

            if (double.TryParse(input1, out inputmin))
            {
                min = 2 * Math.PI * inputmin / Ft;
                dinput1 = inputmin;
            }
            else
            {
                MessageBox.Show("Invalid value");
            }

            if (double.TryParse(input2, out inputmax))
            {
                max = 2 * Math.PI * inputmax / Ft;
                dinput2 = inputmax;
            }
            else
            {
                MessageBox.Show("Invalid value");
            }


            for (int i = 0; i < SampleSize; i++)
            {
                if (i == 0)
                {
                    input[i] = min;
                    inputs[i] = dinput1;
                }
                else
                {
                    input[i] = input[i - 1] + (max / (SampleSize - 1));
                    inputs[i] = inputs[i - 1] + (dinput2 / (SampleSize - 1));
                }

            }

            Complex[] ss = new Complex[SampleSize];
            Th(ss, N1, r, R);
            LP(ss, input, output, mag, phas);

            for (int i = 0; i < SampleSize; i++)
            {

                chart1.Series["Magnitude"].Points.AddXY(inputs[i], mag[i]);
                chart2.Series["Phase"].Points.AddXY(inputs[i], phas[i]);
                chart3.Series["Transfer"].Points.AddXY(inputs[i], output[i].Real);

            }
            chart1.ChartAreas[0].AxisX.Maximum = 50000;
            chart1.ChartAreas[0].AxisX.Minimum = 0;
            chart1.ChartAreas[0].AxisX.Interval = 10000;

            chart2.ChartAreas[0].AxisX.Maximum = 50000;
            chart2.ChartAreas[0].AxisX.Minimum = 0;
            chart2.ChartAreas[0].AxisX.Interval = 10000;

            chart3.ChartAreas[0].AxisX.Maximum = 50000;
            chart3.ChartAreas[0].AxisX.Minimum = 0;
            chart3.ChartAreas[0].AxisX.Interval = 10000;
        }

        

        private void chart2_Click(object sender, EventArgs e)
        {

        }

        private void button3_Click(object sender, EventArgs e)
        {
            chart1.Titles.Clear();
            chart2.Titles.Clear();
            chart3.Titles.Clear();
            var input1 = textBox1.Text;
            var input2 = textBox2.Text;
            double inputmin = 0;
            double min = 0;
            double inputmax = 0; ;
            double max = 0;
            double dinput1 = 0;
            double dinput2 = 0;
            Complex[] output = new Complex[SampleSize];
            double[] phas = new double[SampleSize];
            double[] mag = new double[SampleSize];
            double[] input = new double[SampleSize];
            double[] inputt = new double[SampleSize];
            double[] inputs = new double[SampleSize];
            chart1.Series[0].Points.Clear();
            chart1.Titles.Add("Magnitude");
            chart1.ChartAreas[0].AxisX.Title = "Frequency (Hz)";
            chart1.ChartAreas[0].AxisY.Title = "Magnitude (dB)";
            chart2.Series[0].Points.Clear();
            chart2.Titles.Add("Phase");
            chart2.ChartAreas[0].AxisX.Title = "Frequency (Hz)";
            chart2.ChartAreas[0].AxisY.Title = "Phase";
            chart3.Series[0].Points.Clear();
            chart3.Titles.Add("Transfer Function");
            chart3.ChartAreas[0].AxisX.Title = "Frequency (Hz)";
            chart3.ChartAreas[0].AxisY.Title = "Amplitude";


            if (double.TryParse(input1, out inputmin))
            {
                min = 2 * Math.PI * inputmin / Ft;
                dinput1 = inputmin;
            }
            else
            {
                MessageBox.Show("Invalid value");
            }

            if (double.TryParse(input2, out inputmax))
            {
                max = 2 * Math.PI * inputmax / Ft;
                dinput2 = inputmax;
            }
            else
            {
                MessageBox.Show("Invalid value");
            }


            for (int i = 0; i < SampleSize; i++)
            {
                if (i == 0)
                {
                    input[i] = min;
                    inputs[i] = dinput1;
                }
                else
                {
                    input[i] = input[i - 1] + (max / (SampleSize - 1));
                    inputs[i] = inputs[i - 1] + (dinput2 / (SampleSize - 1));
                }

            }

            Complex[] ss = new Complex[SampleSize];

            Th(ss, N1, r, R);

            BI(ss, input, output, mag, phas);

            //int minn = Convert.ToInt32(dinput1);
            //int maxx = Convert.ToInt32(dinput2);

            for (int i = 0; i < SampleSize; i++)
            {
                

                chart1.Series["Magnitude"].Points.AddXY(inputs[i], mag[i]);
                chart2.Series["Phase"].Points.AddXY(inputs[i], phas[i]);
                chart3.Series["Transfer"].Points.AddXY(inputs[i], output[i].Real);

            }

            chart1.ChartAreas[0].AxisX.Maximum = 75000;
            chart1.ChartAreas[0].AxisX.Minimum = 0;
            chart1.ChartAreas[0].AxisX.Interval = 10000;

            chart1.ChartAreas[0].AxisY.Minimum = -200;

            chart2.ChartAreas[0].AxisX.Maximum = 75000;
            chart2.ChartAreas[0].AxisX.Minimum = 0;
            chart2.ChartAreas[0].AxisX.Interval = 10000;

            chart3.ChartAreas[0].AxisX.Maximum = 75000;
            chart3.ChartAreas[0].AxisX.Minimum = 0;
            chart3.ChartAreas[0].AxisX.Interval = 10000;
        }

        private void button4_Click(object sender, EventArgs e)
        {
            chart1.Titles.Clear();
            chart2.Titles.Clear();
            chart3.Titles.Clear();
            double[] inputt1 = new double[SampleSizer];
            double[] inputt2 = new double[SampleSizer];
            double[] inputt3 = new double[SampleSizer];
            double[] inputt4 = new double[SampleSizer];
            double[] output = new double[SampleSizer];
            double[] inputs = new double[SampleSizer];
            chart1.Series[0].Points.Clear();
            chart1.Titles.Add("Direct Form Realization");
            chart1.ChartAreas[0].AxisX.Title = "Time (s)";
            chart1.ChartAreas[0].AxisY.Title = "Amplitude";
            chart2.Series[0].Points.Clear();
            chart2.Titles.Add("Generated Signal");
            chart2.ChartAreas[0].AxisX.Title = "Time (s)";
            chart2.ChartAreas[0].AxisY.Title = "Amplitude";
            chart3.Series[0].Points.Clear();
            chart3.Titles.Add("Frequency Spectrum");
            chart3.ChartAreas[0].AxisX.Title = "Time (s)";
            chart3.ChartAreas[0].AxisY.Title = "Amplitude";

            
            double[] GSignal = new double[SampleSizer];

            SineWave(2 * Math.PI * 20000 / Ft, inputt1);
            SineWave(2 * Math.PI * 40000 / Ft, inputt2);
            SineWave(2 * Math.PI * 55000 / Ft, inputt3);
            SineWave(2 * Math.PI * 70000 / Ft, inputt4);

            Complex[] scignal = new Complex[SampleSizer];

            for (int m = 0; m < SampleSizer; m++)
            {
                GSignal[m] = inputt1[m] + inputt2[m] + inputt3[m] + inputt4[m];
                scignal[m] = new Complex(GSignal[m], (double)0);
            }

            double[] z = { 0.001055, -0.004221, 0.006331, -0.004221, 0.001055 };
            double[] p = { -1, -3.15, -3.927, -2.274, -0.5135 };

            double[] w = new double[SampleSizer];

            for (int n = 1; n < SampleSizer-5; n++)
            {
                w[n + 4] = -1*scignal[n + 4].Real + (p[1] * w[n + 3]) + (p[2] * w[n + 2]) + (p[3] * w[n + 1]) + (p[4] * w[n]);
                output[n + 4] = (z[0] * w[n + 4]) +(z[1] * w[n + 3]) + (z[2] * w[n + 2]) + (z[3] * w[n + 1]) + (z[4] * w[n]);

                chart1.Series["Magnitude"].Points.AddXY(n, output[n+4]);
                chart2.Series["Phase"].Points.AddXY(n, scignal[n+4].Real);
            }

            Fourier.Forward(scignal);

            for (int q = 1; q < SampleSizer/2; q++)
            {
                chart3.Series["Transfer"].Points.AddXY(inputs[q]*2, scignal[q].Magnitude);
            }

            chart1.ChartAreas[0].AxisX.Maximum = 1000;
            chart1.ChartAreas[0].AxisX.Minimum = 0;
            chart1.ChartAreas[0].AxisX.Interval = 1000;

            chart1.ChartAreas[0].AxisY.Maximum = 2;
            chart1.ChartAreas[0].AxisY.Minimum = -2;
            chart1.ChartAreas[0].AxisY.Interval = 1;

            chart2.ChartAreas[0].AxisX.Maximum = 1000;
            chart2.ChartAreas[0].AxisX.Minimum = 0;
            chart2.ChartAreas[0].AxisX.Interval = 1000;

            chart2.ChartAreas[0].AxisY.Maximum = 5;
            chart2.ChartAreas[0].AxisY.Minimum = -5;
            chart2.ChartAreas[0].AxisY.Interval = 1;

            chart3.ChartAreas[0].AxisX.Maximum = 75000;
            chart3.ChartAreas[0].AxisX.Minimum = 0;
            chart3.ChartAreas[0].AxisX.Interval = 5000;
            
            chart3.ChartAreas[0].AxisY.Minimum = 0;
            chart3.ChartAreas[0].AxisY.Interval = 10;
        }
        

        private void button6_Click(object sender, EventArgs e)
        {
            chart1.Titles.Clear();
            chart2.Titles.Clear();
            chart3.Titles.Clear();
            var input1 = textBox1.Text;
            var input2 = textBox2.Text;
            double inputmin = 0;
            double min = 0;
            double inputmax = 0; ;
            double max = 0;
            double dinput1 = 0;
            double dinput2 = 0;
            Complex[] output = new Complex[SampleSize];
            double[] phas = new double[SampleSize];
            double[] mag = new double[SampleSize];
            double[] input = new double[SampleSize];
            double[] inputt = new double[SampleSize];
            double[] inputs = new double[SampleSize];
            chart1.Series[0].Points.Clear();
            chart1.Titles.Add("Magnitude");
            chart1.ChartAreas[0].AxisX.Title = "Frequency (Hz)";
            chart1.ChartAreas[0].AxisY.Title = "Magnitude (dB)";
            chart2.Series[0].Points.Clear();
            chart2.Titles.Add("Phase");
            chart2.ChartAreas[0].AxisX.Title = "Frequency (Hz)";
            chart2.ChartAreas[0].AxisY.Title = "Phase";
            chart3.Series[0].Points.Clear();
            chart3.Titles.Add("Transfer Function");
            chart3.ChartAreas[0].AxisX.Title = "Frequency (Hz)";
            chart3.ChartAreas[0].AxisY.Title = "Amplitude";

            if (double.TryParse(input1, out inputmin))
            {
                min = 2 * Math.PI * inputmin / Ft;
                dinput1 = inputmin;
            }
            else
            {
                MessageBox.Show("Invalid value");
            }

            if (double.TryParse(input2, out inputmax))
            {
                max = 2 * Math.PI * inputmax / Ft;
                dinput2 = inputmax;
            }
            else
            {
                MessageBox.Show("Invalid value");
            }


            for (int i = 0; i < SampleSize; i++)
            {
                if (i == 0)
                {
                    input[i] = min;
                    inputs[i] = dinput1;
                }
                else
                {
                    input[i] = input[i - 1] + (max / (SampleSize - 1));
                    inputs[i] = inputs[i - 1] + (dinput2 / (SampleSize - 1));
                }

            }
            Complex[] ss = new Complex[SampleSize];
            Th(ss, N1, r, R);
            HP(ss, input, output, mag, phas);
            int minn = Convert.ToInt32(dinput1);
            int maxx = Convert.ToInt32(dinput2);
            //Complex[] tf = new Complex[max+1];
            //double[] mag = new double[max+1];
            //double[] phase = new double[max+1];
            //double[] inputs = new double[max+1];

           for (int i = 0; i < SampleSize; i++)
            {
                //if (i == min)
                //{
                    //inputs[i] = dinput1;
                //}
                //else
                //{
                    //inputs[i] = inputs[i - 1] + i;
                //}

                //double w = 2 * inputs[i] * Math.PI / Ft;
                //Complex s = new Complex(0, Math.Tan(w / 2));

                //tf[i] = (0.9717 * Complex.Pow(s, 4)) / ((Complex.Pow(s, 4)) + (11.18 * Complex.Pow(s, 3)) + (70.61 * Complex.Pow(s, 2)) + (212.8 * (s)) + 625.3);
                //mag[i] = 20 * Math.Log10(tf[i].Magnitude);
                //phase[i] = tf[i].Phase * (180/Math.PI);
                
                chart1.Series["Magnitude"].Points.AddXY(inputs[i], mag[i]);
                chart2.Series["Phase"].Points.AddXY(inputs[i], phas[i]);
                chart3.Series["Transfer"].Points.AddXY(inputs[i], output[i].Real);

            }
            chart1.ChartAreas[0].AxisX.Maximum = 110000;
            chart1.ChartAreas[0].AxisX.Minimum = 0;
            chart1.ChartAreas[0].AxisX.Interval = 10000;

            chart2.ChartAreas[0].AxisX.Maximum = 110000;
            chart2.ChartAreas[0].AxisX.Minimum = 0;
            chart2.ChartAreas[0].AxisX.Interval = 10000;

            chart3.ChartAreas[0].AxisX.Maximum = 110000;
            chart3.ChartAreas[0].AxisX.Minimum = 0;
            chart3.ChartAreas[0].AxisX.Interval = 10000;
        }

        private void button2_Click(object sender, EventArgs e)
        {
            chart1.Titles.Clear();
            chart2.Titles.Clear();
            chart3.Titles.Clear();
            double[] inputt1 = new double[SampleSizer];
            double[] inputt2 = new double[SampleSizer];
            double[] inputt3 = new double[SampleSizer];
            double[] inputt4 = new double[SampleSizer];
            double[] output = new double[SampleSizer];
            double[] inputs = new double[SampleSizer];
            chart1.Series[0].Points.Clear();
            chart1.Titles.Add("Direct Form Realization");
            chart1.ChartAreas[0].AxisX.Title = "Frequency (Hz)";
            chart1.ChartAreas[0].AxisY.Title = "Amplitude";
            chart2.Series[0].Points.Clear();
            chart2.Titles.Add("Generated Signal");
            chart2.ChartAreas[0].AxisX.Title = "Frequency (Hz)";
            chart2.ChartAreas[0].AxisY.Title = "Amplitude";
            chart3.Series[0].Points.Clear();
            chart3.Titles.Add("Time Domain");
            chart3.ChartAreas[0].AxisX.Title = "Time (s)";
            chart3.ChartAreas[0].AxisY.Title = "Amplitude";


            double[] GSignal = new double[SampleSizer];
            int freq = 0;
            double[] Toets = new double[SampleSizer];

            for (int baie = 0; baie < 500; baie++)
            {
                SineWave(2 * Math.PI * freq / Ft, inputt1);
                freq += 150;
                for (int nm = 0; nm < SampleSizer; nm++)
                {
                    Toets[nm] = Toets[nm] +inputt1[nm];
                }
            }
           

            Complex[] scignal = new Complex[SampleSizer];

            for (int m = 0; m < SampleSizer; m++)
            {
                //GSignal[m] = inputt1[m] + inputt2[m] + inputt3[m] + inputt4[m];
                scignal[m] = new Complex(Toets[m], (double)0);
            }

            

            double[] z = { 0.001055, -0.004221, 0.006331, -0.004221, 0.001055 };
            double[] p = { -1, -3.15, -3.927, -2.274, -0.5135 };

            double[] w = new double[SampleSizer];

            for (int n = 1; n < SampleSizer-5; n++)
            {
                w[n+4] = p[0]*scignal[n+4].Real + (p[1] * w[n +3]) + (p[2] * w[n + 2]) + (p[3] * w[n +1]) + (p[4] * w[n]);
                output[n+4] = (z[0] * w[n+4]) + (z[1] * w[n +3]) + (z[2] * w[n + 2]) + (z[3] * w[n +1]) + (z[4] * w[n]);
            }

            Complex[] soignal = new Complex[SampleSizer];

            for (int v = 0; v < output.Length; v++)
            {
                soignal[v] = new Complex(output[v], (double)0);
            }

            Fourier.Forward(scignal);
            Fourier.Forward(soignal);

            for (int b = 1; b < SampleSizer; b++)
            {
                chart3.Series["Transfer"].Points.AddXY(b, Toets[b]);
            }

            for (int q = 1; q < SampleSizer; q++)
            {
                chart1.Series["Magnitude"].Points.AddXY(inputs[q] , soignal[q].Magnitude);
            }

            for (int c = 1; c < SampleSizer; c++)
            {
                chart2.Series["Phase"].Points.AddXY(inputs[c], scignal[c].Magnitude);
            }

            chart1.ChartAreas[0].AxisX.Maximum = 75000;
            chart1.ChartAreas[0].AxisX.Minimum = 0;
            chart1.ChartAreas[0].AxisX.Interval = 5000;
            chart1.ChartAreas[0].AxisY.Maximum = 200;
            chart1.ChartAreas[0].AxisY.Minimum = 0;
            chart1.ChartAreas[0].AxisY.Interval = 10;

            chart2.ChartAreas[0].AxisX.Maximum = 75000;
            chart2.ChartAreas[0].AxisX.Minimum = 0;
            chart2.ChartAreas[0].AxisX.Interval = 5000;
            chart2.ChartAreas[0].AxisY.Maximum = 200;
            chart2.ChartAreas[0].AxisY.Minimum = 0;
            chart2.ChartAreas[0].AxisY.Interval = 10;
        }
    }   
}
