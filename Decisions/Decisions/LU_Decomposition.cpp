#include <iostream>

using namespace std;


void Print(double * A, int n)
{
	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
			cout << A[i*n + j] << "\t";
		cout << endl;
	}

	cout << endl;
}

void multiplying(double * A, double * B, double * C, int n)
{
	for (int i = 0; i < n; i++)
		for (int j = 0; j < n; j++)
			for (int k = 0; k < n; k++)
				C[i*n + j] += A[i*n + k] * B[k*n + j];
}


void LU_Decomposition1(double * A, double * L, double * U, int n)
{
	L[0] = 1;
	U[0] = A[0];

	for (int i = 0; i < n; i++)
		for(int j = 0; j < n; j++)
		{
			U[i*n + j] = A[i*n + j];
			for (int k = 0; k < i; k++)
				U[i*n + j] -= L[i*n + k] * U[k*n + j];
			if (i > j) U[i*n + j] = 0.0;
			
			L[i*n + j] = A[i*n + j];
			for (int k = 0; k < j; k++)
				L[i*n + j] -= L[i*n + k] * U[k*n + j];
			if (i >= j) L[i*n + j] /= U[j*n + j];
		}
}

void LU_Decomposition2(double * A, double * L, double * U, int n)
{

	for (int i = 0; i < n; i++)
		for (int j = 0; j < n; j++)
		{
			U[i*n + j] = 0.0;
			L[i*n + j] = 0.0;
		}

	for (int i = 0; i < n - 1; i++)
		for (int j = i + 1; j < n; j++)
		{
			A[j*n + i] = A[j*n + i] / A[i*n + i];
			for (int k = i + 1; k < n; k++)
				A[j*n + k] = A[j*n + k] - A[j*n + i] * A[i*n + k];
		}

	for (int i = 0; i < n; i++)
		for (int j = 0; j < n; j++)
		{
			if (i <= j) U[i*n + j] = A[i*n + j];
			if (i > j) L[i*n + j] = A[i*n + j];
			if (i == j) L[i*n + j] = 1.0;
		}
}

void LU_Decomposition3(double * A, double * L, double * U, int n)
{
	int size = 10;

	for (int blok = 0; blok < n / size; blok++)
	{
		L[blok*size*(n + 1)] = 1;
		U[blok*size*(n + 1)] = A[blok*size*(n + 1)];
		for (int i = blok * size; i < (blok + 1)*size; i++)
			for (int j = blok * size; j < (blok + 1)*size; j++)
			{
				U[i*n + j] = A[i*n + j];
				for (int k = blok * size; k < i; k++)
					U[i*n + j] -= L[i*n + k] * U[k*n + j];
				if (i > j) U[i*n + j] = 0.0;

				L[i*n + j] = A[i*n + j];
				for (int k = blok * size; k < j; k++)
					L[i*n + j] -= L[i*n + k] * U[k*n + j];
				if (i >= j) L[i*n + j] /= U[j*n + j];
			}

		for (int i = (blok + 1) * size; i < n; i++)
			for (int j = blok * size; j < (blok + 1)*size; j++)
			{
				L[i*n + j] = A[i*n + j];
				for (int k = blok * size; k < j; k++)
					L[i*n + j] -= L[i*n + k] * U[k*n + j];
				if (i >= j) L[i*n + j] /= U[j*n + j];

				U[j*n + i] = A[j*n + i];
				for (int k = blok * size; k < j; k++)
					U[j*n + i] -= L[j*n + k] * U[k*n + i];
				if (i < j) U[j*n + i] /= L[i*n + i];
			}

		double SUM = 0.0;
		for (int i = (blok + 1) * size; i < n; i++)
			for (int j = (blok + 1) * size; j < n; j++)
			{
				SUM = 0.0;
				for (int k = blok * size; k < (blok + 1)*size; k++)
					SUM += L[i*n + k] * U[k*n + j];
				A[i*n + j] -= SUM;
			}
	}

	int blok = n / size;

	for (int i = blok * size; i < n; i++)
		for (int j = blok * size; j < n; j++)
		{
			U[i*n + j] = A[i*n + j];
			for (int k = blok * size; k < i; k++)
				U[i*n + j] -= L[i*n + k] * U[k*n + j];
			if (i > j) U[i*n + j] = 0.0;

			L[i*n + j] = A[i*n + j];
			for (int k = blok * size; k < j; k++)
				L[i*n + j] -= L[i*n + k] * U[k*n + j];
			if (i >= j) L[i*n + j] /= U[j*n + j];
		}

	for (int i = 0; i < n; i++)
		for (int j = 0; j < i; j++)
		{
			U[i*n + j] = 0.0;
			L[j*n + i] = 0.0;
		}
}

void LU_Decomposition(double * A, double * L, double * U, int n)
{
	int size =10;

	for (int blok = 0; blok < n / size; blok++)
	{
		L[blok*size*(n + 1)] = 1;
		U[blok*size*(n + 1)] = A[blok*size*(n + 1)];
		for (int i = blok * size; i < (blok + 1)*size; i++)
			for (int j = blok * size; j < (blok + 1)*size; j++)
			{
				U[i*n + j] = A[i*n + j];
				for (int k = blok * size; k < i; k++)
					U[i*n + j] -= L[i*n + k] * U[k*n + j];
				if (i > j) U[i*n + j] = 0.0;

				L[i*n + j] = A[i*n + j];
				for (int k = blok * size; k < j; k++)
					L[i*n + j] -= L[i*n + k] * U[k*n + j];
				if (i >= j) L[i*n + j] /= U[j*n + j];
			}

#pragma omp parallel for
		for (int i = (blok + 1) * size; i < n; i++)
			for (int j = blok * size; j < (blok + 1)*size; j++)
			{
				L[i*n + j] = A[i*n + j];
				for (int k = blok * size; k < j; k++)
					L[i*n + j] -= L[i*n + k] * U[k*n + j];
				if (i >= j) L[i*n + j] /= U[j*n + j];

				U[j*n + i] = A[j*n + i];
				for (int k = blok * size; k < j; k++)
					U[j*n + i] -= L[j*n + k] * U[k*n + i];
				if (i < j) U[j*n + i] /= L[i*n + i];
			}

		double SUM = 0.0;
		int start = (blok + 1) * size;
#pragma omp parallel for private(SUM)
		for(int ii = 0; ii < (n - start) / size; ii++)
			for (int jj = 0; jj < (n - start) / size; jj++)
				for (int i = ii * size + start; i < (ii + 1) * size + start; i++)
					for (int j = jj * size + start; j < (jj + 1) * size + start; j++)
					{
						SUM = 0.0;
						for (int k = blok * size; k < (blok + 1)*size; k++)
							SUM += L[i*n + k] * U[k*n + j];
						A[i*n + j] -= SUM;
					}

		for (int i = start; i < n; i++)
			for (int j = ((n - start) / size) * size + start; j < n; j++)
			{
				SUM = 0.0;
				for (int k = blok * size; k < (blok + 1)*size; k++)
					SUM += L[i*n + k] * U[k*n + j];
				A[i*n + j] -= SUM;
			}

		for (int i = (n - start) / size * size + start; i < n; i++)
			for (int j = start; j < ((n - start) / size) * size + start; j++)
			{
				SUM = 0.0;
				for (int k = blok * size; k < (blok + 1)*size; k++)
					SUM += L[i*n + k] * U[k*n + j];
				A[i*n + j] -= SUM;
			}
	}

	int blok = n / size;

	for (int i = blok * size; i < n; i++)
		for (int j = blok * size; j < n; j++)
		{
			U[i*n + j] = A[i*n + j];
			for (int k = blok * size; k < i; k++)
				U[i*n + j] -= L[i*n + k] * U[k*n + j];
			if (i > j) U[i*n + j] = 0.0;

			L[i*n + j] = A[i*n + j];
			for (int k = blok * size; k < j; k++)
				L[i*n + j] -= L[i*n + k] * U[k*n + j];
			if (i >= j) L[i*n + j] /= U[j*n + j];
		}

	for (int i = 0; i < n; i++)
		for (int j = 0; j < i; j++)
		{
			U[i*n + j] = 0.0;
			L[j*n + i] = 0.0;
		}
}


void main()
{
	int n = 105;

	double *A1 = new double[n*n];
	double *L1 = new double[n*n];
	double *U1 = new double[n*n];

	double *A2 = new double[n*n];
	double *L2 = new double[n*n];
	double *U2 = new double[n*n];

	for (int i = 0; i < n; i++)
		for (int j = 0; j < n; j++)
		{
			A1[i*n + j] = rand() % 100 + 1;
			A2[i*n + j] = A1[i*n + j];
		}

	LU_Decomposition(A1, L1, U1, n);
	LU_Decomposition2(A2, L2, U2, n);

	double Uerr = 0.0;
	double Lerr = 0.0;
	for (int i = 0; i < n; i++)
		for (int j = 0; j < n; j++)
		{
			Uerr += abs(U1[i*n + j] - U2[i*n + j]);
			Lerr += abs(L1[i*n + j] - L2[i*n + j]);
		}

	cout << "n = " << n << "\tUerr = " << Uerr << "\tLerr = " << Lerr << endl;

	cin.get();
}