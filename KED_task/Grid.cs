using System.Collections;
using System.Collections.Specialized;
using System.Drawing;

namespace KED_task;

public class FiniteElement
{
    public int[] Nodes { get; set; } = new int[4];
    public double Sigma { get; set; }
    public int this[int i] { get => Nodes[i]; set => Nodes[i] = value; }
}

public class Grid
{
    public List<PointRZ> Nodes { get; init; }
    public List<FirstCondition> Boundary { get; init; }
    public FiniteElement[] Elements { get; init; }
    public double VELRadius { get => _rPoints is null ? 0.001 : _rPoints[1]; }
    public double r0 { get => _rPoints is null ? 0.001 : _rPoints[0]; }
    public double[] Time { get; init; }

    private double[] _rPoints;
    private int[] _rSteps;
    private double[] _rCoefs;

    private double[] _zPoints;
    private int[] _zSteps;
    private double[] _zCoefs;
    private double[] _sigmas;

    public int NumR { get; init; }
    public int NumZ { get; init; }

    private double current;
    public double Current { get => current; }

    public Grid(string spaceGridPath, string timeGridPath)
    {
        string[] data;
        int rSplits, zSplits;

        // чтение данных из файла сетки
        using (StreamReader sr = new(spaceGridPath))
        {
            data = sr.ReadLine()!.Split(" ").Where(str => str != "").ToArray();
            current = Convert.ToDouble(data[0]);


            data = sr.ReadLine()!.Split(" ").Where(str => str != "").ToArray();
            rSplits = Convert.ToInt32(data[1]);

            _rPoints = new double[rSplits + 1];
            _rPoints[0] = Convert.ToDouble(data[0]);
            _rSteps = new int[rSplits];
            _rCoefs = new double[rSplits];

            data = sr.ReadLine()!.Split(" ").Where(str => str != "").ToArray();
            for (int i = 0; i < rSplits; i++)
            {
                _rPoints[i + 1] = Convert.ToDouble(data[i]);
            }

            data = sr.ReadLine()!.Split(" ").Where(str => str != "").ToArray();
            for (int i = 0; i < rSplits; i++)
            {
                _rSteps[i] = Convert.ToInt32(data[i]);
            }

            data = sr.ReadLine()!.Split(" ").Where(str => str != "").ToArray();
            for (int i = 0; i < rSplits; i++)
            {
                _rCoefs[i] = Convert.ToDouble(data[i]);
            }

            data = sr.ReadLine()!.Split(" ").Where(str => str != "").ToArray();
            zSplits = Convert.ToInt32(data[1]);

            _zPoints = new double[zSplits + 1];
            _zPoints[0] = Convert.ToDouble(data[0]);
            _zSteps = new int[zSplits];
            _zCoefs = new double[zSplits];
            _sigmas = new double[zSplits];

            data = sr.ReadLine()!.Split(" ").Where(str => str != "").ToArray();
            for (int i = 0; i < zSplits; i++)
            {
                _zPoints[i + 1] = Convert.ToDouble(data[i]);
            }

            data = sr.ReadLine()!.Split(" ").Where(str => str != "").ToArray();
            for (int i = 0; i < zSplits; i++)
            {
                _zSteps[i] = Convert.ToInt32(data[i]);
            }

            data = sr.ReadLine()!.Split(" ").Where(str => str != "").ToArray();
            for (int i = 0; i < zSplits; i++)
            {
                _zCoefs[i] = Convert.ToDouble(data[i]);
            }

            data = sr.ReadLine()!.Split(" ").Where(str => str != "").ToArray();
            for (int i = 0; i < zSplits; i++)
            {
                _sigmas[i] = Convert.ToDouble(data[i]);
            }
        }

        using (StreamReader sr = new(timeGridPath))
        {
            data = sr.ReadLine().Split(" ").ToArray();
            var steps = Convert.ToInt32(data[1]);
            Time = new double[steps + 1];
            Time[0] = Convert.ToDouble(data[0]);
            Time[^1] = Convert.ToDouble(data[2]);
            var step = (Time[^1] - Time[0]) / steps;

            for (int i = 1; i < steps; i++)
            {
                Time[i] = Time[0] + step * i;
            }
        }

        // Генерация сетки и глобальная нумерация узлов
        var rUniq = new double[_rSteps.Sum() + 1];
        var zUniq = new double[_zSteps.Sum() + 1];
        NumR = rUniq.Length;
        NumZ = zUniq.Length;

        rUniq[0] = _rPoints[0];
        int ridx = 1;
        for (int i = 0; i < rSplits; i++)
        {
            double sumCoef = 0;
            for (int j = 0; j < _rSteps[i]; j++)
            {
                sumCoef += Math.Pow(_rCoefs[i], j);
            }

            double r_step = (_rPoints[i + 1] - _rPoints[i]) / sumCoef;

            for (int j = 0; j < _rSteps[i] - 1; j++)
            {
                rUniq[ridx] = rUniq[ridx - 1] + r_step;
                r_step *= _rCoefs[i];
                ridx++;
            }

            rUniq[ridx++] = _rPoints[i + 1];
        }

        zUniq[0] = _zPoints[0];
        int zidx = 1;
        for (int i = 0; i < zSplits; i++)
        {
            double sumCoef = 0;
            for (int j = 0; j < _zSteps[i]; j++)
            {
                sumCoef += Math.Pow(_zCoefs[i], j);
            }

            double z_step = (_zPoints[i + 1] - _zPoints[i]) / sumCoef;

            for (int j = 0; j < _zSteps[i] - 1; j++)
            {
                zUniq[zidx] = zUniq[zidx - 1] + z_step;
                z_step *= _zCoefs[i];
                zidx++;
            }

            zUniq[zidx++] = _zPoints[i + 1];
        }

        Nodes = new List<PointRZ>();

        for (int i = 0; i < zUniq.Length; i++)
            for (int j = 0; j < rUniq.Length; j++)
                Nodes.Add(new(rUniq[j], zUniq[i]));


        // генерация элементов
        Elements = new FiniteElement[(rUniq.Length - 1) * (zUniq.Length - 1)];

        int elidx = 0;
        zidx = 0;
        for (int i = 0; i < _zSteps.Length; i++)
        {
            for (int j = 0; j < _zSteps[i]; j++)
            {
                for (int rIdx = 0; rIdx < rUniq.Length - 1; rIdx++)
                {
                    Elements[elidx] = new();
                    Elements[elidx].Nodes[0] = (zidx + 1) * rUniq.Length + rIdx;
                    Elements[elidx].Nodes[1] = (zidx + 1) * rUniq.Length + rIdx + 1;
                    Elements[elidx].Nodes[2] = zidx * rUniq.Length + rIdx;
                    Elements[elidx].Nodes[3] = zidx * rUniq.Length + rIdx + 1;
                    Elements[elidx].Sigma = _sigmas[i];
                    elidx++;
                }
                zidx++;
            }
        }

        // указываем узлы с краевые узлы
        Boundary = new List<FirstCondition>();

        // нижняя граница
        for (int i = 1; i <= rUniq.Length; i++)
            Boundary.Add(new FirstCondition(Nodes[Nodes.Count - i], Nodes.Count - i, 0));

        // боковые границы
        for (int i = 1; i < zUniq.Length - 1; i++)
        {
            Boundary.Add(new FirstCondition(Nodes[rUniq.Length * i], rUniq.Length * i, 0));
            Boundary.Add(new FirstCondition(Nodes[rUniq.Length * (i + 1) - 1], rUniq.Length * (i + 1) - 1, 0));
        }

        // верхняя граница
        for (int i = _rSteps[0] + 1; i < rUniq.Length; i++)
            Boundary.Add(new FirstCondition(Nodes[i], i, 0));
        for (int i = 0; i <= _rSteps[0]; i++)
        {
            Boundary.Add(new FirstCondition(Nodes[i], i, 1));
        }
    }
}
