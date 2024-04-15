namespace KED_task;

public class FirstCondition
{
    public PointRZ Point { get; }
    public int NodeNumber { get; }
    public int Type { get; }

    public FirstCondition(PointRZ node, int nodeNumber, int type)
    {
        Point = node;
        NodeNumber = nodeNumber;
        Type = type;
    }

    public double Value(double r, double cur = 1.0)
        => Type switch
        {
            0 => 0.0,
            _ => cur / (2.0 * Math.PI * r)
        };
}

public class SecondCondition
{
   public int[] Edge { get; }
   public int ElemNumber { get; }
   public int EdgeType { get; }   // 0 - bottom, 1 - right
                                       // 2 - top, 3 - left

   public SecondCondition(int elemNumber, int edgeType, int[] edge)
   {
      ElemNumber = elemNumber;
      EdgeType = edgeType;
      Edge = edge;
   }
}
