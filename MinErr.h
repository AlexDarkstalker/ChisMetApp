#include "stdafx.h"

//умножение матрицы на вектор в схеме крест для уравнения Лапласса
void MultVektSpec(double **phi, int n, double **R)
//n - размерность, , R - результат. На краях квадрата должны стоять нули
{
	for (int i = 1; i < n; i++)
	{
		for (int j = 1; j < n; j++)
			R[i][j] = 4.*phi[i][j] - phi[i - 1][j] - phi[i + 1][j] - phi[i][j - 1] - phi[i][j + 1];
	}
}

double ScalMultVect(double** a, double**b, int n)
{
	double sum = 0.;
	for (int i = 1; i < n; i++)
	{
		for (int j = 1; j < n; j++)
			sum += a[i][j] * b[i][j];
	}
	return sum;
}

void ZeroContur(double** a, int n)
{
	for (int i = 0; i < n; i++)
		a[i][0] = 0.;

	for (int i = 0; i < n; i++)
		a[n][i] = 0.;

	for (int i = 0; i < n; i++)
		a[n - i][n] = 0.;

	for (int i = 0; i < n; i++)
		a[0][n - i] = 0.;
}


//Метод минимальных ошибок
int MethodMinErr(double **phi, int n, double eps, double **f)
//  phi - матрица, n - размерность, eps - точность, f - вектор свободных членов
{



	/*double ** phi_b = new double*[n + 1]; //вектор предыдущих значений
	for (int i = 0; i < n + 1; i++)
	phi_b[i] = new double[n + 1];


	for (int i = 0; i <= n; i++)
	for (int j = 0; j <= n; j++)
	{
	phi_b[i][j] = phi[i][j];
	}*/

	double h = 1. / (double)n;
	double hh = h*h;
	double Tau = 0.;
	double feps;
	double norm_f;
	double norm_r;
	int j = 0;
	double s;
	int count = 0;//  количество итераций

	double **r = new double*[n + 1];
	for (int i = 0; i < n + 1; i++)
		r[i] = new double[n + 1];

	ZeroContur(r, n);

	double **Ar = new double*[n + 1];
	for (int i = 0; i < n + 1; i++)
		Ar[i] = new double[n + 1];

	ZeroContur(Ar, n);
	norm_f = 0.;
	for (int i = 1; i < n; i++)
	{
		for (s = 0., j = 1; j < n; j++)
			s += fabs(f[i][j]);
		if (s > norm_f)
			norm_f = s;
	}

	feps = eps * norm_f;

	for (int i = 1; i < n; i++)
	for (int j = 1; j < n; j++)
		f[i][j] *= hh;

	do
	{
		MultVektSpec(phi, n, r);
		for (int i = 1; i < n; i++)
		{
			for (int j = 1; j < n; j++)
				r[i][j] -= f[i][j];//Вектор невязок
		}

		MultVektSpec(r, n, Ar);

		Tau = ScalMultVect(r, r, n) / ScalMultVect(Ar, r, n);

		for (int i = 1; i < n; i++)
		for (int j = 1; j < n; j++)
			phi[i][j] -= Tau*r[i][j];

		norm_r = 0.;
		for (int i = 1; i < n; i++)
		{
			for (s = 0., j = 1; j < n; j++)
				s += fabs(r[i][j]);
			if (s > norm_r)
				norm_r = s;
		}

		count++;
		if (count > 1.e+5)
			return -1;
	} while (norm_r > feps);
	return count;
}