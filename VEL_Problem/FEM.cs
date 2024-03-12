using System.Drawing;
using static System.Net.Mime.MediaTypeNames;

namespace VEL_Problem;

public class FEM
{
    private Grid grid;
    private SparseMatrix globalMatrix;
    private Vector globalVector;
    private Solver slae;
    private Vector localVector;
    private Matrix stiffnessMatrix;
    private Matrix massMatrix;
    private Vector[] layers;
    private FirstCondition[] firstConditions;
    private IBasis2D Basis;
    (int maxIters, double eps) slaeParametres = (10000, 1e-14);

    public FEM(Grid grid)
    {
        grid = grid;

        Basis = new BilinearBasis();


        stiffnessMatrix = new(Basis.Size);
        massMatrix = new(Basis.Size);
        localVector = new(Basis.Size);
        slae = new Solver(slaeParametres.maxIters, slaeParametres.eps);


        layers = new Vector[2].Select(_ => new Vector(grid.Nodes.Count)).ToArray();
    }

    public void SetSlaeParametres(int maxiters, double eps)
    {
        slaeParametres = (maxiters, eps);
        slae.SetParametres(maxiters, eps);
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
        }

    }

    public void PrepareStartConditions()
    {
        // начальное условие
        double U(PointRZ point, double time)
            => point.R * time;

        // первый слой
        for (int i = 0; i < layers[0].Length; i++)
        {
            layers[1][i] = U(grid.Nodes[i], grid.Time[0]);
        }
    }


    public void AccountFirstConditions()
    {
        foreach (var fc in firstConditions)
        {
            globalMatrix.Di[fc.NodeNumber] = 1;
            globalVector[fc.NodeNumber] = 0;
            for (int i = globalMatrix.Ig[fc.NodeNumber]; i < globalMatrix.Ig[fc.NodeNumber + 1]; i++)
            {
                globalMatrix.Gg[i] = 0;
            }
            for (int i = fc.NodeNumber + 1; i < globalMatrix.Size; i++)
            {
                for (int j = globalMatrix.Ig[i]; j < globalMatrix.Ig[i + 1]; j++)
                {
                    if (globalMatrix.Jg[j] == fc.NodeNumber)
                    {
                        globalMatrix.Gg[j] = 0;
                    }
                }
            }
        }
    }

    public void AccountFirstConditionsWithBigNumber()
    {
        foreach (var fc in firstConditions)
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
        double t01 = grid.Time[itime] - grid.Time[itime - 1];
        double t02 = itime > 1 ? grid.Time[itime] - grid.Time[itime - 2] : 1;
        double t12 = itime > 1 ? grid.Time[itime - 1] - grid.Time[itime - 2] : 1;

        globalMatrix.Clear();
        globalVector.Fill(0);

        for (int ielem = 0; ielem < grid.Elements.Length; ielem++)
        {
            AssemblyLocalMatrixes(ielem);

            stiffnessMatrix = 1 / (grid.Elements[ielem].Sigma) * stiffnessMatrix;

            for (int i = 0; i < Basis.Size; i++)
                for (int j = 0; j < Basis.Size; j++)
                    AddElement(grid.Elements[ielem].Nodes[i], grid.Elements[ielem].Nodes[j], stiffnessMatrix[i, j]);

            AssemblyLocallVector(ielem, itime, t01, t02, t12);

            AddElementToVector(ielem);

            stiffnessMatrix.Clear();
            localVector.Fill(0);
        }
        //globalVector = mu0 * globalVector;
    }

    void AssemblyLocallVector(int ielem, int itime, double t01, double t02, double t12)
    {

        for (int i = 0; i < Basis.Size; i++)
        {
            localVector[i] = 0;
        }

        Vector qj1 = new(6);
        Vector qj2 = new(6);
        if (itime == 1)
        {
            for (int i = 0; i < qj1.Length; i++)
            {
                for (int j = 0; j < qj1.Length; j++)
                    qj1[i] += massMatrix[i, j] * layers[1][grid.Elements[ielem].Nodes[j]];
            }

            localVector += 1.0 / t01 * qj1;
        }
        else
        {
            for (int i = 0; i < qj1.Length; i++)
            {
                for (int j = 0; j < qj1.Length; j++)
                {
                    qj2[i] += massMatrix[i, j] * layers[0][grid.Elements[ielem].Nodes[j]];

                    qj1[i] += massMatrix[i, j] * layers[1][grid.Elements[ielem].Nodes[j]];
                }
            }

            localVector += t02 / (t12 * t01) * qj1 - t01 / (t02 * t12) * qj2;
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
                    => Basis.GetPsi(i, point) * Basis.GetPsi(j, point) * point.R;

                stiffnessMatrix[i, j] = stiffnessMatrix[j, i] = Integration.Gauss2D(stifFunc, elem);

            }
    }


   

    private int FindElement(PointRZ point)
    {
        for (int ielem = 0; ielem < grid.Elements.Length; ielem++)
        {
            
        }
        return -1;
    }



    public void PrintSolution()
    {
        for (int i = 0; i < layers[1].Length; i++)
        {
            Console.WriteLine($"x = {grid.Nodes[i].R:e3}, y = {grid.Nodes[i].Z:e3}, Hphi = {layers[1][i]}");
        }
    }

    public static Func<PointRZ, double> Mult(Func<PointRZ, double> fst, Func<PointRZ, double> scnd)
        => (point) => fst(point) * scnd(point);
    public static Func<PointRZ, double> Sum(Func<PointRZ, double> fst, Func<PointRZ, double> scnd)
        => (point) => fst(point) + scnd(point);
}

