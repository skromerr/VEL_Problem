using System.Drawing;
using static System.Net.Mime.MediaTypeNames;

namespace KED_task;

public class FEM
{
    private const double mu0 = 4.0 * Math.PI * 1E-07;
    private Grid grid;
    private SparseMatrix globalMatrix;
    private Vector globalVector;
    private Solver slae;
    private Vector localVector;
    private Matrix stiffnessMatrix;
    private Matrix massMatrix;
    private Matrix massApproxMatrix;
    private Vector[] layers;
    private IBasis2D Basis;
    (int maxIters, double eps) slaeParametres = (10000, 1e-14);
    private List<PointRZ> receivers;

    public FEM(Grid grid)
    {
        this.grid = grid;

        Basis = new BilinearBasis();

        stiffnessMatrix = new(Basis.Size);
        massMatrix = new(Basis.Size);
        massApproxMatrix = new(Basis.Size);
        localVector = new(Basis.Size);
        slae = new Solver(slaeParametres.maxIters, slaeParametres.eps);

        layers = new Vector[2].Select(_ => new Vector(grid.Nodes.Count)).ToArray();

        receivers = new List<PointRZ>();
    }

    public void SetSlaeParametres(int maxiters, double eps)
    {
        slaeParametres = (maxiters, eps);
        slae.SetParametres(maxiters, eps);
    }

    public void SetReceivers(string path)
    {
        using StreamReader sr = new(path);

        string[] data;

        data = sr.ReadLine().Split(" ").ToArray();
        int numRec = Convert.ToInt32(data[0]);

        for (int i = 0; i < numRec; i++)
        {
            data = sr.ReadLine().Split(" ").Where(str => str != "").ToArray();
            receivers.Add(new(Convert.ToDouble(data[0]), Convert.ToDouble(data[1])));
        }
    }

    public void Compute()
    {
        BuildPortrait();

        PrepareStartConditions();

        for (int itime = 1; itime < grid.Time.Length; itime++)
        {
            AssemblySLAE(itime);
            AccountFirstConditions();

            slae.SetSLAE(globalVector, globalMatrix);
            slae.CGM();

            Vector.Copy(layers[1], layers[0]);
            Vector.Copy(slae.solution, layers[1]);

            PrintValueAtReceivers(itime);
            CheckResult(itime);
        }
    }

    public void PrepareStartConditions()
    {
        AssemblySLAE(0);
        AccountFirstConditions(grid.Current);

        slae.SetSLAE(globalVector, globalMatrix);
        slae.CGM();

        Vector.Copy(layers[1], layers[0]);
        Vector.Copy(slae.solution, layers[1]);

        PrintValueAtReceivers(0, false);
    }

    public void AccountFirstConditions(double cur = 0)
    {
        foreach (var fc in grid.Boundary)
        {
            globalMatrix.Di[fc.NodeNumber] = 1;
            var value = fc.Value(fc.Point.R, cur);
            globalVector[fc.NodeNumber] = value;
            for (int i = globalMatrix.Ig[fc.NodeNumber]; i < globalMatrix.Ig[fc.NodeNumber + 1]; i++)
            {
                globalVector[globalMatrix.Jg[i]] -= value * globalMatrix.Gg[i];
                globalMatrix.Gg[i] = 0;
            }
            for (int i = fc.NodeNumber + 1; i < globalMatrix.Size; i++)
            {
                for (int j = globalMatrix.Ig[i]; j < globalMatrix.Ig[i + 1]; j++)
                {
                    if (globalMatrix.Jg[j] == fc.NodeNumber)
                    {
                        globalVector[i] -= value * globalMatrix.Gg[j];
                        globalMatrix.Gg[j] = 0;
                    }
                }
            }
        }
    }
        public void AccountFirstConditionsWithBigNumber()
    {
        foreach (var fc in grid.Boundary)
        {
            globalMatrix.Di[fc.NodeNumber] = 1e30;
            globalVector[fc.NodeNumber] = 1e30 * fc.Value(grid.r0);
        }
    }

    public void BuildPortrait()
    {
        var list = new HashSet<int>[grid.Nodes.Count].Select(_ => new HashSet<int>()).ToList();
        foreach (var element in grid.Elements)
        {
            foreach (var pos in element.Nodes)
            {
                foreach (var node in element.Nodes)
                {
                    if (pos > node)
                    {
                        list[pos].Add(node);
                    }
                }
            }
        }

        int count = list.Sum(childList => childList.Count);

        globalMatrix = new(grid.Nodes.Count, count);
        globalVector = new(grid.Nodes.Count);

        globalMatrix.Ig[0] = 0;

        for (int i = 0; i < list.Count; i++)
            globalMatrix.Ig[i + 1] = globalMatrix.Ig[i] + list[i].Count;

        int k = 0;

        foreach (var childList in list)
        {
            foreach (var value in childList.Order())
            {
                globalMatrix.Jg[k++] = value;
            }
        }
    }
        private void AddElement(int i, int j, double value)
    {
        if (i == j)
        {
            globalMatrix.Di[i] += value;
            return;
        }

        for (int icol = globalMatrix.Ig[i]; icol < globalMatrix.Ig[i + 1]; icol++)
        {
            if (globalMatrix.Jg[icol] == j)
            {
                globalMatrix.Gg[icol] += value;
                return;
            }
        }
    }

    private void AssemblySLAE(int itime)
    {
        double t01 = itime > 0 ? grid.Time[itime] - grid.Time[itime - 1] : 1;
        double t02 = itime > 1 ? grid.Time[itime] - grid.Time[itime - 2] : 1;
        double t12 = itime > 1 ? grid.Time[itime - 1] - grid.Time[itime - 2] : 1;

        globalMatrix.Clear();
        globalVector.Fill(0);

        for (int ielem = 0; ielem < grid.Elements.Length; ielem++)
        {
            AssemblyLocalMatrixes(ielem);

            stiffnessMatrix = 1 / (grid.Elements[ielem].Sigma) * stiffnessMatrix;
            massMatrix = 1 / (grid.Elements[ielem].Sigma) * massMatrix;
            massApproxMatrix = mu0 * (itime > 0 ? 1.0 / t01 + (itime > 1 ? 1.0 / t02 : 0) : 0) * massApproxMatrix;

            stiffnessMatrix += massMatrix + massApproxMatrix;

            for (int i = 0; i < Basis.Size; i++)
                for (int j = 0; j < Basis.Size; j++)
                    AddElement(grid.Elements[ielem].Nodes[i], grid.Elements[ielem].Nodes[j], stiffnessMatrix[i, j]);

            AssemblyLocallVector(ielem, itime, t01, t02, t12);

            AddElementToVector(ielem);

            stiffnessMatrix.Clear();
            localVector.Fill(0);
        }
    }

    void AssemblyLocallVector(int ielem, int itime, double t01, double t02, double t12)
    {
        for (int i = 0; i < Basis.Size; i++)
        {
            localVector[i] = 0;
        }

        Vector qj1 = new(Basis.Size);
        Vector qj2 = new(Basis.Size);
        if (itime == 1)
        {
            for (int i = 0; i < qj1.Length; i++)
            {
                for (int j = 0; j < qj1.Length; j++)
                    qj1[i] += massApproxMatrix[i, j] * layers[1][grid.Elements[ielem].Nodes[j]];
            }

            localVector += 1.0 / t01 * mu0 * qj1;
        }
        else if (itime > 1)
        {
            for (int i = 0; i < qj1.Length; i++)
            {
                for (int j = 0; j < qj1.Length; j++)
                {
                    qj2[i] += massApproxMatrix[i, j] * layers[0][grid.Elements[ielem].Nodes[j]];

                    qj1[i] += massApproxMatrix[i, j] * layers[1][grid.Elements[ielem].Nodes[j]];
                }
            }

            localVector += t02 / (t12 * t01) * mu0 * qj1 - t01 / (t02 * t12) * mu0 * qj2;
        }
    }

    private void AddElementToVector(int ielem)
    {
        for (int i = 0; i < Basis.Size; i++)
        {
            globalVector[grid.Elements[ielem].Nodes[i]] += localVector[i];
        }
    }

    private void AssemblyLocalMatrixes(int ielem)
    {
        var elem = new Rectangle(grid.Nodes[grid.Elements[ielem].Nodes[0]], grid.Nodes[grid.Elements[ielem].Nodes[3]]);
        Basis.SetElem(elem);

        for (int i = 0; i < Basis.Size; i++)
            for (int j = 0; j <= i; j++)
            {
                double stifFunc(PointRZ point)
                    => (Basis.GetDPsi(i, VarType.R, point) * Basis.GetDPsi(j, VarType.R, point) +
                    Basis.GetDPsi(i, VarType.Z, point) * Basis.GetDPsi(j, VarType.Z, point) * point.R);

                double massFunc(PointRZ point)
                    => Basis.GetPsi(i, point) * Basis.GetPsi(j, point) / point.R;

                double massApproxFunc(PointRZ point)
                    => Basis.GetPsi(i, point) * Basis.GetPsi(j, point) * point.R;

                stiffnessMatrix[i, j] = stiffnessMatrix[j, i] = Integration.Gauss2D(stifFunc, elem);

                massMatrix[i, j] = massMatrix[j, i] = Integration.Gauss2D(massFunc, elem);

                massApproxMatrix[i, j] = massApproxMatrix[j, i] = Integration.Gauss2D(massApproxFunc, elem);
            }
    }
    
    public int FindElement(PointRZ point)
    {
        int numR, numZ;
        if (point.R < grid.Nodes[0].R || point.Z > grid.Nodes[0].Z) return -1;

        for (numR = 1; numR < grid.NumR; numR++)
        {
            if (point.R <= grid.Nodes[numR].R) break;
        }

        if (numR == grid.NumR) return -1;

        for (numZ = 1; numZ < grid.NumZ; numZ++)
        {
            if (point.Z >= grid.Nodes[numZ * grid.NumR].Z) 
                return (numZ - 1) * (grid.NumR - 1) + numR - 1;
        }
        return -1;
    }

    public double ValueAtPoint(PointRZ point)
    {
        double res = 0;
        int ielem = FindElement(point);
        if (ielem == -1) return 0;
        var elem = new Rectangle(grid.Nodes[grid.Elements[ielem].Nodes[0]], grid.Nodes[grid.Elements[ielem].Nodes[3]]);
        Basis.SetElem(elem);

        for (int i = 0; i < Basis.Size; i++)
            res += layers[1][grid.Elements[ielem][i]] * Basis.GetPsi(i, point);
        return res;
    }

    public double ValueAtPointAtPreviousLayer(PointRZ point)
    {
        {
            double res = 0;
            int ielem = FindElement(point);
            if (ielem == -1) return 0;
            var elem = new Rectangle(grid.Nodes[grid.Elements[ielem].Nodes[0]], grid.Nodes[grid.Elements[ielem].Nodes[3]]);
            Basis.SetElem(elem);

            for (int i = 0; i < Basis.Size; i++)
                res += layers[0][grid.Elements[ielem][i]] * Basis.GetPsi(i, point);
            return res;
        }
    }

    private void PrintValueAtReceivers(int itime, bool append = true)
    {
        using StreamWriter sw = new("..\\..\\..\\output\\ReceiversResults.txt", append);
        sw.Write($"{grid.Time[itime]:E4}");
        for (int i = 0; i < receivers.Count; i++)
            sw.Write($" {ValueAtPoint(receivers[i]):E7}");
        sw.WriteLine();

        sw.Close();
    }

    public void PrintSolution()
    {
        for (int i = 0; i < layers[1].Length; i++)
        {
            Console.WriteLine($"x = {grid.Nodes[i].R:e3}, y = {grid.Nodes[i].Z:e3}, Hphi = {layers[1][i]}");
        }
    }

    private double rotE(PointRZ point)
    {
        int ielem = FindElement(point);
        double sigma = grid.Elements[ielem].Sigma;

        double res;
        //double dPhiDz(PointRZ point) => Derivative(ValueAtPoint, point, 1);
        //double dPhiDr(PointRZ point) => Derivative(ValueAtPoint, point, 0);

        res = SecondDerivative(ValueAtPoint, point, 1) + SecondDerivative(ValueAtPoint, point, 0) + Derivative(ValueAtPoint, point, 0) / point.R - ValueAtPoint(point) / (point.R * point.R);
        return - res / sigma;
    }

    private double minusdBdT(PointRZ point, int itime)
    {
        double res;

        res = TimeDerivative(point, itime);

        return -mu0 * res;
    }

    public static double Derivative(Func<PointRZ, double> func, PointRZ point, int varType)
    {
        double res;
        double h = 1e-7;
        switch (varType)
        {
            case 0:
                res = (func(point + (Math.Abs(point.R) * h, 0)) - func(point)) / (Math.Abs(point.R) * h);
                break;
            default:
                res = (func(point + (0, Math.Abs(point.Z) * h)) - func(point)) / (Math.Abs(point.Z) * h);
                break;
        }

        return res;
    }

    public static double SecondDerivative(Func<PointRZ, double> func, PointRZ point, int varType)
    {
        double res;
        double h = 1e-7;
        switch (varType)
        {
            case 0:
                res = (func(point + (Math.Abs(point.R) * h, 0)) - 2 * func(point) + func(point - (Math.Abs(point.R) * h, 0))) / (Math.Abs(point.R) * h * Math.Abs(point.R) * h);
                break;
            default:
                res = (func(point + (0, Math.Abs(point.Z) * h)) - 2 * func(point) + func(point - (0, Math.Abs(point.Z) * h))) / (Math.Abs(point.Z) * h * Math.Abs(point.Z) * h);
                break;
        }

        return res;
    }

    private double TimeDerivative(PointRZ point, int itime)
    {
        double res;

        res = ValueAtPoint(point) - ValueAtPointAtPreviousLayer(point);
        res /= (grid.Time[itime] - grid.Time[itime - 1]);

        return res;
    }

    private void CheckResult(int itime)
    {
        for (int i = 0; i < receivers.Count; i++)
            Console.WriteLine($"t = {grid.Time[itime]:E3}:   " +
                $"rotE_phi = {rotE(receivers[i]):E7};      " +
                $"-dB/dt = {minusdBdT(receivers[i], itime):E7}");
    }
}

