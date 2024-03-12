using System;
using System.Collections.Generic;
using System.Drawing;
using System.Linq;
using System.Reflection;
using System.Text;
using System.Threading.Tasks;

namespace VEL_Problem;

public enum VarType { R = 0, Z = 1 }

public interface IBasis2D
{
    int Size { get; }

    double GetPsi(int number, PointRZ point);

    double GetDPsi(int number, VarType varType, PointRZ point);
}

public class BilinearBasis : IBasis2D
{
    public int Size => 4;

    public double GetPsi(int number, PointRZ point)
        => number switch
        {
            0 => (1.0 - point.R) * (1.0 - point.Z),
            1 => point.R * (1.0 - point.Z),
            2 => (1.0 - point.R) * point.Z,
            3 => point.R * point.Z,
            _ => throw new ArgumentOutOfRangeException(nameof(number), number, "Not expected function number!")
        };

    public double GetDPsi(int number, VarType varType, PointRZ point)
        => varType switch
        {
            VarType.R => number switch
            {
                0 => point.Z - 1.0,
                1 => 1.0 - point.Z,
                2 => -point.Z,
                3 => point.Z,
                _ => throw new ArgumentOutOfRangeException(nameof(number), number, "Not expected function number!")
            },

            VarType.Z => number switch
            {
                0 => point.R - 1.0,
                1 => -point.R,
                2 => 1.0 - point.R,
                3 => point.R,
                _ => throw new ArgumentOutOfRangeException(nameof(number), number, "Not expected function number!")
            },
            _ => throw new ArgumentOutOfRangeException(nameof(varType), varType, "Not expected var type!")
        };
}