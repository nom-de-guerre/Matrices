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

#include <ConjugateGradient.h> // defines typedef Matrix_t<double> Md_t

int __DIM = 1000;

void run ();
Md_t BuildPDM (void);

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
	Md_t A = BuildPDM ();

	Md_t b (__DIM, 1);
	b.randomly_fill (__DIM);

	Md_t QR = A;
	QR.copy ();
	Md_t __b = b;
	Md_t b_save = b;

	printf ("Solving with QR\n");
	Md_t x = QR.solveQR (__b);

	__b = A * x;
	assert (__b.equal_eps (b, 1e-8));

	ConjugateGrad_t CG (A, b);
	CG.Compute ();

	Md_t _x = CG.Answer ();

	printf ("METRICS\t%e\t%e\t%e\n",
		(A * _x - b).vec_magnitude (),
		(A * _x - b).vec_magnitude () / b.vec_magnitude (),
		(_x - x).vec_magnitude () / x.vec_magnitude ());
}

/*
 * Construct a positive definite matrix.
 *
 */
Md_t BuildPDM (void)
{
	Md_t A (__DIM, __DIM);

	for (int i = 0; i < __DIM; ++i)
	{
		for (int j = 0; j < __DIM; ++j)
			A (i, j) = (double) rand () / RAND_MAX;
	}

	A = 0.5 * (A * A.transpose ());
	for (int i = 0; i < __DIM; ++i)
		A (i, i) *= __DIM;

	return A;
}

