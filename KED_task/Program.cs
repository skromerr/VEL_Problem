using KED_task;

Thread.CurrentThread.CurrentCulture = new System.Globalization.CultureInfo("en-US");

string gridPath = "..\\..\\..\\input\\grid.txt";
string timePath = "..\\..\\..\\input\\time.txt";
string receiversPath = "..\\..\\..\\input\\receivers.txt";

Grid grid = new(gridPath, timePath);
FEM fem = new(grid);
fem.SetReceivers(receiversPath);
fem.Compute();