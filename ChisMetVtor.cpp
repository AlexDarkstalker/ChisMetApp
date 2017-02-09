//3.17 Решение системы уравнений методом минимальных ошибок
#include "stdafx.h"
#include "MinErr.h"


using namespace std;
#define N 100; //размер

const double eps = 1.e-10; //DBL_EPSILON
const double pi = acos(-1);

double MatrixNorm(double **f, int n)
{
	int i, j;
	double s, norm = 0.;

	for (i = 0; i < n; i++)
	{
		for (s = 0., j = 0; j < n; j++) s += fabs(f[i][j]);
		if (s > norm) norm = s;
	}

	return norm;
}

double calculNormNevSPEC(double **ar, int n, double **f)
{
	double h = 1. / n;
	double hhInverse = (double)(n*n);
	double norm = 0.;
	double r;

	for (int i = 1; i < n - 1; i++)
	for (int j = 1; j < n - 1; j++)
	{
		r = (4.*ar[i][j] - ar[i][j - 1] - ar[i][j + 1] - ar[i - 1][j] - ar[i + 1][j])*hhInverse - f[i][j];
		if (fabs(r) > norm)
			norm = fabs(r);
	}

	return norm;
}

//Правая часть уравнения
double fun(double x, double y)
{
	//return -x*x + x - y*y + y
	return 5.*pi*pi*sin(pi*x)*sin(2.*pi*y);
}

//точное решение
double u(double x, double y)
{
	//return 0.5*(x*x*y*y - x*y*y - x*x*y + x*y);
	return sin(pi*x)*sin(2.*pi*y);
}

int _tmain(int argc, _TCHAR* argv[])
{
	setlocale(LC_ALL, "Russian");
	int n = N;

	double h = 1. / (double)n;
	double xi, yj;


	double **phi = new double*[n + 1];
	for (int i = 0; i < n + 1; i++)
		phi[i] = new double[n + 1];


	for (int i = 0; i < n + 1; i++)
	for (int j = 0; j < n + 1; j++)
		phi[i][j] = 0.;


	double **phiCopy = new double*[n + 1];
	for (int i = 0; i < n + 1; i++)
		phiCopy[i] = new double[n + 1];


	for (int i = 0; i < n + 1; i++)
	for (int j = 0; j < n + 1; j++)
		phiCopy[i][j] = phi[i][j];


	double **f = new double*[n + 1];
	for (int i = 0; i < n + 1; i++)
		f[i] = new double[n + 1];


	h = 1. / n;


	for (int i = 0; i < n + 1; i++)
	{
		xi = i * h;
		for (int j = 0; j < n + 1; j++)
		{
			yj = j * h;
			f[i][j] = fun(xi, yj);
		}
	}

	//Вызов функции 


	double **fCopy = new double*[n + 1];
	for (int i = 0; i < n + 1; i++)
		fCopy[i] = new double[n + 1];

	for (int i = 0; i < n + 1; i++)
	for (int j = 0; j < n + 1; j++)
		fCopy[i][j] = f[i][j];

	int code = MethodMinErr(phi, n, eps, f);

	double normaphi = calculNormNevSPEC(phi, n, fCopy);

	double normaf = MatrixNorm(fCopy, n);

	double z_c;
	double norm = 0.;
	for (int i = 1; i < n; i++)
	{
		xi = (double)i*h;
		for (int j = 1; j < n; j++)
		{
			yj = (double)j*h;
			z_c = fabs(phi[i][j] - u(xi, yj));
			if (z_c > norm)
			{
				norm = z_c;
			}
		}
	}

	cout << "Размерность = : " << n << endl;
	cout << "eps  " << eps << endl;
	cout << "Количество итераций: " << code << endl;
	cout << "Норма ошибки в решении : " << norm << endl;
	cout << "Отн невязка: " << normaphi / normaf << endl;
	//system("pause");

	return 0;
}

