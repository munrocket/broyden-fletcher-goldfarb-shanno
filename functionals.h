typedef double(*double_in_double)(double);

double differential(double_in_double f, double x);
double simpson(double_in_double f, double left, double right, int k);
double golden_section(double_in_double q, double A, double B, double epsilon);
double max(double, double);
double min(double, double);
