/*

Copyright (c) 2015, Douglas Santry
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, is permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
   list of conditions and the following disclaimer.
2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

*/

#include <unistd.h>
#include <stdlib.h>
#include <time.h>
#include <assert.h>
#include <float.h>

#include <GMRES.h>
#include <francis.h>

int __DIM = 1000;

void run ();
void unitary (Krylov_t &);

int __m;

int main (int argc, char *argv[])
{
	long seed = time (0);

	if (argc > 1)
		seed = atol (argv[1]);

	__m = atoi (argv[2]);

//	printf ("Using seed %ld\n", seed);

	srand (seed);

	run ();

	return 0;
}

void run ()
{
	Ms_t A (__DIM, __DIM);
	Md_t b (__DIM, 1);

	/*
	 * Generate a diagonaly dominant matrix that is fairly sparse
	  and asymetric.
	 *
	 */
	for (int i = 0; i < __DIM; ++i)
	{
		int band = 0;
		int start;
		int end;

		while (band == 0)
			band = rand () % 3;

		(i - band < 0 ? start = 0 : start = i - band);
		(i + band >= __DIM ? end = __DIM - 1 : end = i + band);

		for (int j = start; j <= end; ++j)
		{
			double sample = 0.0;

			while (sample == 0.0)
				sample = (double) (rand () % (100));

			if (i == j)
				A[i][j] = 4 * sample;
			else
				A[i][j] = sample;
		}
	}

	b.randomly_fill (5);

	Md_t K = A.Copy ();
	Md_t Spectrum = K;
	Md_t __b = b;
	Md_t b_save = b;
	Md_t x = K.solveQR (__b);

#if 0
	EigenFrancis_t Aleph;
	int N = Aleph.CalcEigenValues (Spectrum);
	Aleph.SortEigenValues ();
	printf ("%d eigen values ∆ = %f\n", 
		N, 
		Aleph.ef_EigenValues[N - 1].real - 
		Aleph.ef_EigenValues[0].real);

	for (int i = 0; i < N; ++i)
		printf ("%f\t%fi\n", 
			Aleph.ef_EigenValues[i].real, 
			Aleph.ef_EigenValues[i].imag);

	printf ("%f\t%fi\n", 
		Aleph.ef_EigenValues[0].real, 
		Aleph.ef_EigenValues[0].imag);

	printf ("%f\t%fi\n", 
		Aleph.ef_EigenValues[N - 1].real, 
		Aleph.ef_EigenValues[N - 1].imag);
#endif

	__b = A * x;
	assert (__b.equal_eps (b, 1e-10));

	/*
	 * Ok, we have confirmed that there is a solution.  We can test
	 * GMRES now.
	 *
	 */

//	printf ("STARTING RUN\n");

	double residual;
	bool solved;
	Md_t _x;

	/*
	 * Increase projected matrix size until we find a solution.
	 *
	 */
	for (int m = 100; m < 500; m += 50)
	{
		GMRES_t Z (__m, A, b);
		Z.SetTolerance (0.5); // Euclidean magnitude of residual vector
		solved = Z.Solve (_x, residual);

		printf ("%f\t%f\t%d\n", 
			residual, 
			(x - _x).vec_magnitude (),
			__m);

		if (solved)
			break;
	}

//	printf ("Error = %f\n", (x - _x).vec_magnitude ());
}

