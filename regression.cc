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

#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <assert.h>

#include <matrix.h>

typedef Matrix_t<double> Md_t;

#define N_POINTS 10
#define UNIVERSE (10 * N_POINTS)

double M55 [] = 
{
13 , 
5 , 
54 , 
98 , 
67 , 
66 , 
25 , 
70 , 
77 , 
38 , 
9 , 
71 , 
80 , 
97 , 
45 , 
72 , 
84 , 
61 , 
77 , 
60 , 
55 , 
29 , 
38 , 
87 , 
63
};

double M53 [] =
{
42 , 
85 , 
94 , 
56 , 
85 , 
26 , 
68 , 
2 , 
15 , 
72 , 
12 , 
76 , 
45 , 
5 , 
73
};

double AnswerProd [] =
{
14569 , 
3149 , 
14501 , 
16186 , 
8989 , 
16530 , 
18803 , 
8349 , 
14549 , 
20120 , 
14606 , 
20099 , 
15617 , 
8575 , 
17705
};

double scaler [] =
{
6.5 , 
2.5 , 
27 , 
49 , 
33.5 , 
33 , 
12.5 , 
35 , 
38.5 , 
19 , 
4.5 , 
35.5 , 
40 , 
48.5 , 
22.5 , 
36 , 
42 , 
30.5 , 
38.5 , 
30 , 
27.5 , 
14.5 , 
19 , 
43.5 , 
31.5
};


void VerifyMultiplication (void);
void VerifyQR (void);
void VerifySymetricSolution (void);
void VerifyQRSolution (void);

int main (void)
{
	srand (time (0));

	VerifyMultiplication ();
	VerifyQR ();
	VerifySymetricSolution ();
	VerifyQRSolution ();

	printf ("\nSUCCESS: All tests passed.\n");

	return 0;
}

void VerifyMultiplication (void)
{
	Md_t A (5, 5, M55);
	Md_t B (5, 3, M53);
	Md_t Answer (5, 3, AnswerProd);
	Md_t scalerMult (5, 5, scaler);

	Md_t G = A * B;
	assert (Answer == G);

	G = 0.5 * A;
	assert (scalerMult == G);

	G = G + G;
	assert (G == A);

	assert (A.rcount () == 1);
	assert (B.rcount () == 1);
	assert (Answer.rcount () == 1);
	assert (scalerMult.rcount () == 1);
	assert (G.rcount () == 1);

	printf ("Basic Operations:\t\t\tPassed.\n");
}

void VerifyQR (void)
{
	Md_t A (N_POINTS, N_POINTS);

	for (int i = 0; i < N_POINTS; ++i)
		for (int j = 0; j < N_POINTS; ++j)
			A (i, j) = rand () % UNIVERSE;

	Md_t A_saved = A; // A will contain R after the call
	Md_t Q (N_POINTS, N_POINTS);

	A.QR (Q);

	Md_t I (N_POINTS, N_POINTS, 1);

	Md_t Qt = Q.transpose ();
	Md_t G = Q * Qt;
	assert (G.equal_eps (I, 0.1));

	Md_t R = A;
	A = Q * R;
	assert (A.equal_eps (A_saved, 0.1));

	assert (A.rcount () == 1);
	assert (Q.rcount () == 1);
	assert (Qt.rcount () == 1);
	assert (I.rcount () == 1);
	assert (G.rcount () == 1);

	printf ("QR:\t\t\t\t\tPassed.\n");
}

void VerifySymetricSolution (void)
{
	Md_t A (N_POINTS, N_POINTS);
	Md_t b (N_POINTS, 1);
	Md_t x (N_POINTS, 1);

	/*
	 * Construct a positive definite matrix.
	 *
	 */
	for (int i = 0; i < N_POINTS; ++i)
	{
		for (int j = 0; j < N_POINTS; ++j)
			A (i, j) = rand () % UNIVERSE;

		b (i, 0) = rand () % UNIVERSE;
	}

	A = 0.5 * (A * A.transpose ());
	for (int i = 0; i < N_POINTS; ++i)
		A (i, i) *= N_POINTS;

	Md_t A_saved = A;
	bool rc = A.SolveSymmetric (b, x);
	assert (rc);

	Md_t G = A * x;

	assert (G.equal_eps (b, 0.1));

	printf ("Cholesky solve for x:\t\t\tPassed.\n");
}

void VerifyQRSolution (void)
{
	Md_t A (N_POINTS, N_POINTS);
	Md_t b (N_POINTS, 1);

	for (int i = 0; i < N_POINTS; ++i)
	{
		for (int j = 0; j < N_POINTS; ++j)
			A (i, j) = rand () % UNIVERSE;

		b (i, 0) = rand () % UNIVERSE;
	}

	Md_t A_save = A;
	Md_t b_save = b;
	Md_t x = A.solveQR (b);
	Md_t verify = A_save * x;

	assert (verify.equal_eps (b_save, 0.1));

	printf ("QR Solve for x:\t\t\t\tPassed.\n");
}

