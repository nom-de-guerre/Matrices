/*

Copyright (c) 2020, Douglas Santry
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

#include <sys/types.h>
#include <sys/uio.h>
#include <unistd.h>

#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <float.h>

#include <utility>

#include <francis.h>

#define __DIM 1000

void run (void);

int main (int argc, char *argv[])
{
	unsigned seed = time (0);

	if (argc == 2)
		seed = atoi (argv[1]);

	printf ("Using seed %d\n", seed);

	srand (seed);

	run ();
}

void run (void)
{
	Md_t A (__DIM, __DIM);

	for (int i = 0; i < __DIM; ++i)
		for (int j = 0; j < __DIM; ++j )
			A(i, j) = rand () % __DIM;

	Md_t _A = A;
	_A.copy ();

	EigenFrancis_t FR;

	int N = FR.CalcEigenValuesGeneral (A);

	printf ("Finished Processing: %d iterations\n", FR.N_Iterations ());

	assert (N == FR.ef_N);

	int real = 0;
	for (int i = 0; i < FR.ef_N; ++i) 
	{
		double residual;

		if (FR.ef_EigenValues[i].imag)
			continue;

		++real;
		A = _A;

		Md_t u (A.rows (), 1);
		for (int j = A.rows () - 1; j >= 0; --j)
			u(j, 0) = 1;
		u.vec_norm ();
	
		printf ("Looking for eigenvector associated with %.4f (%d)\t",
					FR.ef_EigenValues[i].real,
					i);

		if (!FR.FindEigenVectorReal (FR.ef_EigenValues[i].real, A, u)) {

			printf (" - did not converge\n");
			continue;
		}

		Md_t r = (A * u - FR.ef_EigenValues[i].real * u);
		residual = r.vec_magnitude ();
		printf (" - residual %e\n", residual);
		assert (fabs (residual) < 1e-7);
	}

	printf ("Finished processing %d eigenvalues, %d are real\n", N, real);
}

