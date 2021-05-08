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

#include <GMRES.h> // defines typedef Matrix_t<double> Md_t

int __DIM = 1000;

void run ();
void unitary (Krylov_t &);

int KrylovDim = 200;

int main (int argc, char *argv[])
{
	long seed = time (0);
	char opt;

	while (true)
	{
		opt = getopt (argc, argv, "s:m:");
		if (opt == -1)
			break;

		switch (opt)
		{
		case 's':

			seed = atol (optarg);
			break;

		case 'm':

			KrylovDim = atoi (optarg);
			break;

		default:

			printf ("usage: %s [-s seed] [-m size of Km]\n", argv[0]);
			exit (-1);
		}
	}

	printf ("Using seed %ld\n", seed);

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
	 * and asymetric.
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

	Md_t QR = A.Copy ();
	Md_t __b = b;
	Md_t b_save = b;

	printf ("Solving with QR\n");
	Md_t x = QR.solveQR (__b);

	__b = A * x;
	assert (__b.equal_eps (b, 1e-10));

	/*
	 * Ok, we have confirmed that there is a solution.  We can test
	 * GMRES now.
	 *
	 */

	printf ("Solving with GMRES\n");

	double residual;
	bool solved;
	Md_t _x;

	GMRES_t K (KrylovDim, A, b, 1000);
	K.SetTolerance (0.5); // Euclidean magnitude of residual vector
	solved = K.Solve (_x, residual);

	printf ("Restarts = %d\tResidual = %f\tError = %f\t%d\n", 
		K.GetIterations (),
		residual, 
		(x - _x).vec_magnitude (),
		KrylovDim);
}

