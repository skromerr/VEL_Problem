using System;
using System.Collections.Generic;
using System.Drawing;
using System.Linq;
using System.Reflection;
using System.Text;
using System.Threading.Tasks;
using System.Xml;

namespace VEL_Problem;

public enum VarType { R = 0, Z = 1 }

public interface IBasis2D
{
    int Size { get; }

    void SetElem(Rectangle elem);
    double GetPsi(int number, PointRZ point);

    double GetDPsi(int number, VarType varType, PointRZ point);
}

public class BilinearBasis : IBasis2D
{
    public int Size => 4;

    private Rectangle cur_elem;
    private double hx => cur_elem.TopRight.R - cur_elem.TopLeft.R;
    private double hy => cur_elem.TopLeft.Z - cur_elem.BottomLeft.Z;

    public void SetElem(Rectangle elem)
    {
        cur_elem = elem;
    }

    private double X(int i, double x)
        => i switch
        {
            0 => (cur_elem.TopRight.R - x) / hx,
            _ => (x - cur_elem.TopLeft.R) / hx,
        };

    private double Y(int i, double y)
    => i switch
    {
        0 => (cur_elem.TopRight.Z - y) / hy,
        _ => (y - cur_elem.BottomLeft.Z) / hy,
    };

    public double GetPsi(int number, PointRZ point)
        => number switch
        {
            0 => X(0, point.R) * Y(0, point.Z),
            1 => X(1, point.R) * Y(0, point.Z),
            2 => X(0, point.R) * Y(1, point.Z),
            3 => X(1, point.R) * Y(1, point.Z),
            _ => throw new ArgumentOutOfRangeException(nameof(number), number, "Not expected function number!")
        };

    public double GetDPsi(int number, VarType varType, PointRZ point)
        => varType switch
        {
            VarType.R => number switch
            {
                0 => -Y(0, point.Z) / hx,
                1 => Y(0, point.Z) / hx,
                2 => -Y(1, point.Z) / hx,
                3 => Y(1, point.Z) / hx,
                _ => throw new ArgumentOutOfRangeException(nameof(number), number, "Not expected function number!")
            },

            VarType.Z => number switch
            {
                0 => -X(0, point.R) / hy,
                1 => -X(1, point.R) / hy,
                2 => X(0, point.R) / hy,
                3 => X(1, point.R) / hy,
                _ => throw new ArgumentOutOfRangeException(nameof(number), number, "Not expected function number!")
            },
            _ => throw new ArgumentOutOfRangeException(nameof(varType), varType, "Not expected var type!")
        };
}