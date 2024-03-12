namespace VEL_Problem;

public record struct PointRZ(double R = 0, double Z = 0)
{
    public static PointRZ operator +(PointRZ a, PointRZ b) => new(a.R + b.R, a.Z + b.Z);

    public static PointRZ operator -(PointRZ a, PointRZ b) => new(a.R - b.R, a.Z - b.Z);

    public static PointRZ operator *(double coef, PointRZ a) => new(coef * a.R, coef * a.Z);

    public static PointRZ operator *(PointRZ a, double coef) => coef * a;

    public static PointRZ operator /(PointRZ a, double coef) => new(a.R / coef, a.R / coef);

    public override string ToString() => R.ToString() + ' ' + Z.ToString();
}
