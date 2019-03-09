void LU_Decomposition(double * A, double * L, double * U, int n)
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
		for (int ii = 0; ii < (n - start) / size; ii++)
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