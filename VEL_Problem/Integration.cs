using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;

namespace CED_Problem;

public record struct Rectangle(PointRZ BottomLeft, PointRZ TopRight)
{
    public PointRZ BottomRight => new(TopRight.R, BottomLeft.Z);
    public PointRZ TopLeft => new(BottomLeft.R, TopRight.Z);

    public bool IsIn(PointRZ point)
    {
        if (point.R < TopLeft.R) return false;
        if (point.R > TopRight.R) return false;
        if (point.Z > TopLeft.Z) return false;
        if (point.Z < BottomLeft.Z ) return false;

        return true;
    }
}

public static class Integration
{

    public static double Gauss1D(Func<double, double> func, double start, double end)
    {
        double result = 0;
        double h = end - start;

        (double point, double weight)[] q = GaussOrder5Quadrature().ToArray();

        for (int i = 0; i < q.Length; i++)
        {
            double newpoint = (h * q[i].point + h) / 2.0 + start;
            result += q[i].weight * func(newpoint);
        }

        return result * h / 2.0;
    }

    public static double Gauss2D(Func<PointRZ, double> func, Rectangle area)
    {
        double result = 0;
        double hx = area.BottomRight.R - area.BottomLeft.R;
        double hy = area.TopLeft.Z - area.BottomLeft.Z;

        (double point, double weight)[] q = GaussOrder5Quadrature().ToArray();

        for (int i = 0; i < q.Length; i++)
            for (int j = 0; j < q.Length; j++)
            {
                PointRZ point = new((hx * q[i].point + hx) / 2.0 + area.TopLeft.R, (hy * q[j].point + hy) / 2.0 + area.BottomLeft.Z);
                result += q[i].weight * q[j].weight * func(point);
            }

        return result * hx * hy / 4.0;
    }

    public static IEnumerable<(double, double)> GaussOrder5Quadrature()
    {
        const int n = 3;
        double[] points =
        {
            0,
            -Math.Sqrt(3.0 / 5.0),
            Math.Sqrt(3.0 / 5.0)
        };
        double[] weights =
        {
            8.0 / 9.0,
            5.0 / 9.0,
            5.0 / 9.0
        };

        for (int i = 0; i < n; i++)
        {
            yield return new(points[i], weights[i]);
        }
    }
}
