using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows;
using System.Windows.Controls;
using System.Windows.Data;
using System.Windows.Documents;
using System.Windows.Input;
using System.Windows.Media;
using System.Windows.Media.Imaging;
using System.Windows.Navigation;
using System.Windows.Shapes;

namespace Find_min_genaration_algorithm
{
    /// <summary>
    /// Interaction logic for MainWindow.xaml
    /// </summary>
    public partial class MainWindow : Window
    {
        public MainWindow()
        {
            InitializeComponent();
        }
        Random rand = new Random();

        indv SelfAdaptationMethodES()
        {
            var pop = new List<indv>();
            int f = f1.IsChecked.HasValue ? 1 : (f2.IsChecked.HasValue ? 2 : 3);
            for (int i = 0; i < int.Parse(numPop.Text); i++)
                pop.Add(new indv(int.Parse(nd.Text), f, rand, double.Parse(a.Text), double.Parse(b.Text), double.Parse(c.Text), double.Parse(initSig.Text)));
            indv minF = pop[1];

            var par = new List<indv>();
            for (int i = 0; i < int.Parse(numGen.Text); i++)
            {
                par.Clear();
                for (int j = 0; j < int.Parse(numPop.Text); j++)
                    par.Add(pop[rand.Next(0, pop.Count)]);
                var childs = new List<indv>();
                for (int k = 0; k < int.Parse(numChild.Text) / 2; k++)
                {
                    var p1 = par[rand.Next(0, par.Count)];
                    var p2 = par[rand.Next(0, par.Count)];
                    if (rand.NextDouble() < 0.1)
                    {
                        var ch = simpleCrossover(p1, p2);
                        p1 = ch.Item1;
                        p2 = ch.Item2;
                    }
                    if (rand.NextDouble() <= 0.91)
                    {
                        var mch = mutation(p1);
                        if (mch.fitness > p1.fitness) p1 = mch;
                    }
                    if (rand.NextDouble() <= 0.91)
                    {
                        var mch = mutation(p2);
                        if (mch.fitness > p2.fitness) p2 = mch;
                    }
                    if (p1.fitness > minF.fitness) minF = p1;
                    if (p2.fitness > minF.fitness) minF = p2;
                    childs.Add(p1);
                    childs.Add(p2);
                }
                if (rouletW.IsChecked.Value)
                {
                    pop = rouletweelSelection(childs, int.Parse(numPop.Text));
                }
                else if (sus.IsChecked.Value)
                {
                    pop = SusSelection(childs, int.Parse(numPop.Text));
                }
                else pop = TournamentSelection(childs, int.Parse(numPop.Text));
            }

            return minF;
        }

        indv oneFifthMethodES()
        {
            var pop = new List<indv>();
            int f = f1.IsChecked.HasValue ? 1 : (f2.IsChecked.HasValue ? 2 : 3);
            for (int i = 0; i < int.Parse(numPop.Text); i++)
                pop.Add(new indv(int.Parse(nd.Text), f, rand, double.Parse(a.Text), double.Parse(b.Text), double.Parse(c.Text)));
            indv minF = pop[1];
            double sigma = rand.NextDouble() * double.Parse(initSig.Text);
            var par = new List<indv>();
            for (int i = 0; i < int.Parse(numGen.Text); i++)
            {
                par.Clear();
                for (int j = 0; j < int.Parse(numPop.Text); j++)
                    par.Add(pop[rand.Next(0, pop.Count)]);
                var childs = new List<indv>();
                int success = 0, im = 0;
                for (int k = 0; k < int.Parse(numChild.Text) / 2; k++)
                {
                    var p1 = par[rand.Next(0, par.Count)];
                    var p2 = par[rand.Next(0, par.Count)];
                    if (rand.NextDouble() < 0.1)
                    {
                        var ch = simpleCrossover(p1, p2);
                        p1 = ch.Item1;
                        p2 = ch.Item2;
                    }
                    if (rand.NextDouble() <= 0.91)
                    {
                        im++;
                        var mch = mutation(p1, sigma);
                        if (mch.fitness > p1.fitness)
                        {
                            p1 = mch;
                            success++;
                        }
                    }
                    if (rand.NextDouble() <= 0.91)
                    {
                        im++;
                        var mch = mutation(p2, sigma);
                        if (mch.fitness > p2.fitness)
                        {
                            p2 = mch;
                            success++;
                        }
                    }
                    if (p1.fitness > minF.fitness) minF = p1;
                    if (p2.fitness > minF.fitness) minF = p2;
                    childs.Add(p1);
                    childs.Add(p2);
                }
                if (rouletW.IsChecked.Value)
                {
                    pop = rouletweelSelection(childs, int.Parse(numPop.Text));
                }
                else if (sus.IsChecked.Value)
                {
                    pop = SusSelection(childs, int.Parse(numPop.Text));
                }
                else pop = TournamentSelection(childs, int.Parse(numPop.Text));

                if (success / im > 0.2) sigma /= .8;
                else if (success / im < 0.2) sigma *= .8;
            }
            return minF;
        }

        (indv, indv) simpleCrossover(indv p1, indv p2)
        {
            int sp = rand.Next(1, p1.jen.Length - 1);
            double[] j1 = new double[p1.jen.Length], j2 = new double[p1.jen.Length], s1 = null, s2 = null;
            for (int i = 0; i < j1.Length; i++)
            {
                if (i < sp)
                {
                    j1[i] = p1.jen[i];
                    j2[i] = p2.jen[i];
                }
                else j1[i] = j2[i] = (p1.jen[i] + p2.jen[i]) / 2;
            }
            if (p1.sig != null)
            {
                s1 = new double[p1.sig.Length];
                s2 = new double[p2.sig.Length];
                for (int i = 0; i < s1.Length; i++)
                {
                    if (i < sp)
                    {
                        s1[i] = p1.sig[i];
                        s2[i] = p2.sig[i];
                    }
                    else s1[i] = s2[i] = (p1.sig[i] + p2.sig[i]) / 2;
                }
            }
            int f = f1.IsChecked.HasValue ? 1 : (f2.IsChecked.HasValue ? 2 : 3);
            indv ch1 = new indv(j1, s1, f, double.Parse(a.Text), double.Parse(b.Text), double.Parse(c.Text));
            indv ch2 = new indv(j2, s2, f, double.Parse(a.Text), double.Parse(b.Text), double.Parse(c.Text));
            return (ch1, ch2);
        }

        indv mutation(indv ch, double sig)
        {
            double[] j = new double[ch.jen.Length];
            for (int i = 0; i < j.Length; i++)
                j[i] = ch.jen[i] + (rand.NextDouble() * 2 * sig - sig);

            indv mch = new indv(j, null, f1.IsChecked.HasValue ? 1 : (f2.IsChecked.HasValue ? 2 : 3)
                , double.Parse(a.Text), double.Parse(b.Text), double.Parse(c.Text));
            return mch;
        }

        indv mutation(indv ch)
        {
            double t = 1.0 / Math.Sqrt(2 * Math.Sqrt(ch.jen.Length)),
                tp = 1.0 / Math.Sqrt(2 * ch.jen.Length);

            double[] j = new double[ch.jen.Length], s = new double[ch.sig.Length];
            for (int i = 0; i < s.Length; i++)
            {
                s[i] = ch.sig[i] * Math.Exp(2 * t * rand.NextDouble() - t + 2 * tp * rand.NextDouble() - tp);
                s[i] = s[i] < .1 ? .1 : s[i];
            }

            for (int i = 0; i < j.Length; i++)
                j[i] = ch.jen[i] + (rand.NextDouble() * 2 * s[i] - s[i]);

            indv mch = new indv(j, s, f1.IsChecked.HasValue ? 1 : (f2.IsChecked.HasValue ? 2 : 3)
                , double.Parse(a.Text), double.Parse(b.Text), double.Parse(c.Text));
            return mch;
        }

        class indv
        {
            public double[] jen;
            public double[] sig;
            public double fitness;
            public indv(double[] j, double[] s, int f, double a = 0, double b = 0, double c = 0)
            {
                jen = j;
                sig = s;
                switch (f)
                {
                    case 1:
                        fitness = 1.0 / (f1(jen) + 1);
                        break;
                    case 2:
                        fitness = 1.0 / (f2(jen, a, b, c) + 1);
                        break;
                    case 3:
                        fitness = 1.0 / (f3(jen) + 1);
                        break;
                    default: break;
                }
            }
            public indv(int d, int f, Random rand, double a = 0, double b = 0, double c = 0, double s = 0)
            {
                jen = new double[d];
                for (int i = 0; i < jen.Length; i++) jen[i] = rand.NextDouble() * 20 - 10;
                if (s != 0)
                {
                    sig = new double[d];
                    for (int i = 0; i < jen.Length; i++) sig[i] = rand.NextDouble() * s;
                }
                else sig = null;
                switch (f)
                {
                    case 1:
                        fitness = 1.0 / (f1(jen) + 1);
                        break;
                    case 2:
                        fitness = 1.0 / (f2(jen, a, b, c) + 1);
                        break;
                    case 3:
                        fitness = 1.0 / (f3(jen) + 1);
                        break;
                    default: break;
                }
            }

            double f1(double[] x)
            {
                double res = 10 * x.Length;
                foreach (var xi in x)
                    res += xi * xi - 10 * Math.Cos(2 * Math.PI * xi);
                return res;
            }
            double f2(double[] x, double a, double b, double c)
            {
                double res = 0, temp = 0;
                foreach (var xi in x)
                    temp += xi * xi;
                res = -a * Math.Exp(-b * Math.Sqrt(temp / x.Length));
                temp = 0;
                foreach (var xi in x)
                    temp += Math.Cos(c * xi);
                res -= Math.Exp(temp / x.Length);
                res += a + Math.Exp(1);
                return res;
            }
            double f3(double[] x)
            {
                double res = 10 * x.Length;
                foreach (var xi in x)
                    if (xi > -5.12 && xi < 5.12) res += 10 * xi * xi;
                    else res += xi * xi - 10 * Math.Cos(2 * Math.PI * xi);
                return res;
            }
        }

        List<indv> rouletweelSelection(List<indv> pop, int num)
        {
            var sel = new List<indv>();
            double lengthRoulet = 0;
            foreach (var item in pop) lengthRoulet += item.fitness;
            while (num-- > 0)
            {
                double pointer = rand.NextDouble() * lengthRoulet;
                foreach (var item in pop)
                {
                    if (pointer <= item.fitness)
                    {
                        sel.Add(item);
                        break;
                    }
                    pointer -= item.fitness;
                }
            }
            return sel;
        }
        List<indv> SusSelection(List<indv> pop, int num)
        {
            var sel = new List<indv>();
            List<indv> popSorted = (from indv in pop orderby indv.fitness descending select indv).ToList();
            double lengthRoulet = 0;
            foreach (var item in popSorted) lengthRoulet += item.fitness;
            double pointer = rand.NextDouble() * lengthRoulet / num;
            double iter = 0;
            foreach (var item in popSorted)
            {
                iter += item.fitness;
                if (pointer < iter)
                {
                    sel.Add(item);
                    if (sel.Count == num) break;
                    pointer += lengthRoulet / num;
                }
            }
            return sel;
        }
        List<indv> TournamentSelection(List<indv> pop, int num)
        {
            var sel = new List<indv>();
            indv maxOftour;
            while (num-- > 0)
            {
                maxOftour = null;
                for (int i = 0; i < int.Parse(torK.Text); i++)
                {
                    var t = pop[rand.Next(0, pop.Count)];
                    if ((maxOftour?.fitness ?? 0) < t.fitness)
                        maxOftour = t;
                }
                sel.Add(maxOftour);
            }
            return sel;
        }

        private void Button_Click(object sender, RoutedEventArgs e)
        {
            var watch = new System.Diagnostics.Stopwatch();
            watch.Start();

            indv minf;
            if (SAdaptation.IsChecked.Value) minf = SelfAdaptationMethodES();
            else minf = oneFifthMethodES();
            double f = 1.0 / minf.fitness - 1;
            string xi = "(";
            foreach (var item in minf.jen) xi += ((float)item).ToString() + ',';
            xi = xi.Substring(0, xi.Length - 1) + ')';

            watch.Stop();
            MessageBox.Show($"min(f)={f:F8}" + '\n' + $"(x0,...,xi)={xi}" + '\n' + $"Execution Time: {watch.ElapsedMilliseconds} ms");
        }
    }
}
