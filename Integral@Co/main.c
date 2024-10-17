#include "math.h"
#include <stdio.h>

#define GNUPLOT_NAME "gnuplot -persist"

typedef double (*dfun)(double);

double fun1(double x) {
	return 3 / ((x - 1) * (x - 1) + 1);
}

double dfun1(double x) {
	return -(6 * x - 6) / (x * x * x * x - 4 * x * x * x + 8 * x * x - 8 * x + 4);
}

double fun2(double x) {
	return sqrt(x + 0.5);
}

double dfun2(double x) {
	return 1 / (sqrt(2) * sqrt(2 * x + 1));
}

double fun3(double x) {
	return exp(-x);
}

double dfun3(double x) {
	return -exp(-x);
}


dfun derivative(double (*fun)(double)) {
	if (fun == fun1) return dfun1;
	else if (fun == fun2) return dfun2;
	else if (fun == fun3) return dfun3;
}


//via Newton method
double root(double (*f)(double), double (*g)(double), double a, double b, double eps) {
	double x = (a + b) / 2;
	double xl = b;
	double p = 1;
	while (fabs(xl - x) >= eps) {
		p = f(x) - g(x);
		double df = (derivative(f))(x)-(derivative(g))(x);
		xl = x;
		x = x - p / df;
	}
	return x;
}

//via rectangle method
double integral(double (*f)(double), double a, double b, double n, double eps) {

	double h, S, s1, x;
	h = (b - a) / n;
	S = 0;
	for (int i = 0; i < n - 1; i++)
	{
		x = a + i * h;
		S += f(x);
	}
	S *= h;
	
	do {
		s1 = S;     //второе приближение
		n = 2 * n;  //увеличение числа шагов в два раза, 
		//т.е. уменьшение значения шага в два раза
		h = (b - a) / n;
		S = 0;
		for (int i = 0; i < n - 1; i++)
		{
			x = a + i * h;
			S += f(x);
		}
		S *= h;
	} while (fabs(s1 - S) > eps);
	
	return S;
}

//variant 9
void main() {
	double eps1 = 1e-5; double eps2 = eps1*10;
	double x1 = root(fun1, fun3, -0.25, -0.2, eps1);
	double x2 = root(fun3, fun2, 0.1, 0.2, eps1);
	double x3 = root(fun1, fun2, 1.8, 2, eps1);
	printf("root of fun1 and fun3 is %f\n", x1);
	printf("root of fun3 and fun2 is %f\n", x2);
	printf("root of fun1 and fun2 is %f\n", x3);
	double y1 = integral(fun1, x1, x3, 10000, eps2);
	double y2 = integral(fun2, x2, x3, 10000, eps2);
	double y3 = integral(fun3, x1, x2, 10000, eps2);
	printf("integral of fun1 in area [%f;%f] is %f\n", x1, x3, y1);
	printf("integral of fun2 in area [%f;%f] is %f\n", x2, x3, y2);
	printf("integral of fun3 in area [%f;%f] is %f\n", x1, x2, y3);
	double res = y1 - y2 - y3;
	printf("Result is %f", res);

	FILE* pipe = _popen(GNUPLOT_NAME, "w");
	if (pipe != NULL)
	{
		fprintf(pipe, "set xlabel \"X\"\n");
		fprintf(pipe, "set ylabel \"Y\"\n");
		fprintf(pipe, "set yrange[0:4]\n");
		fprintf(pipe, "set xrange[-3:5]\n");
		fprintf(pipe, "set grid\n");
		fprintf(pipe, "set title \"Variant 9\"\n");
		fprintf(pipe, "set label 1 \"\" at %f,%f point pointtype 1\n", x1, fun1(x1));
		fprintf(pipe, "set label 2 \"\" at %f,%f point pointtype 1\n", x2, fun2(x2));
		fprintf(pipe, "set label 3 \"\" at %f,%f point pointtype 1\n", x3, fun1(x3));
		fprintf(pipe, "f(x)=3 / ((x - 1) * (x - 1) + 1); g(x)=sqrt(x + 0.5); h(x)=exp(-x);\n");
		fprintf(pipe, "max(a,b)=(a>b)?a:b\n");
		fprintf(pipe, "lowerb(x)=max(g(x),h(x))\n");
		fprintf(pipe, "upperb(x)=f(x)\n");
		fprintf(pipe, "plot f(x) lc rgb \"red\", g(x) lc rgb \"green\", h(x) lc rgb \"blue\" \
				, '+' using 1:(lowerb($1)):(upperb($1)) ti \"internal area\" with filledcurves below \
				fc \"red\" fs pattern 1\n");
		fprintf(pipe, "set term pngcairo\n");
		fprintf(pipe, "set output \"myFile.png\"\n");
		fprintf(pipe, "replot\n");
		fflush(pipe);
		_pclose(pipe);
	}
	else
		printf("\nCould not open pipe");
}
