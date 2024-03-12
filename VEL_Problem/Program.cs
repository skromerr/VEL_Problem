using VEL_Problem;

Thread.CurrentThread.CurrentCulture = new System.Globalization.CultureInfo("en-US");

string gridPath = "..\\..\\..\\grid.txt";
string timePath = "..\\..\\..\\time.txt";

Grid grid = new(gridPath, timePath);