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

#ifndef __DJS_ARNOLDI__H__
#define __DJS_ARNOLDI__H__

#include <matrix.h>
#include <SparseMatrix.h>

typedef Matrix_t<double> Md_t;
typedef SparseMatrix::SparseMatrix_t Ms_t;

/*
 * Computes a Krylov subspace, Kn = { b, An, ..., A^(n-1)b }, with
 * Arnoldi iteration:
 *
 * AVm = Vm+1Hm (Saad, Iterative Methods for Sparse Linear Systems)
 *
 * It is a struct as to be useful it needs to be inherited by a class
 * that uses it (e.g. GMRES or IRAM).
 *
 */

struct Krylov_t 
{
	Ms_t		&k_A;
	Md_t		k_b;
	Md_t		k_x0;
	Md_t		k_e1;

	Md_t		k_H;			// Hessenberg matrix of projection
	Md_t		k_Q;			// Orthonormal basis of Kn

	int			k_n;			// maximum number of iterations
	int			k_i;			// iterations so far

	Krylov_t (Ms_t &A, Md_t &b, int n) :
		k_A (A),
		k_b (b)
	{
		Restart (b, n);
	}

	~Krylov_t (void)
	{
	}

	void Restart (Md_t &x0, int n)
	{
		k_n = n;
		k_i = 0;

		Md_t H (k_n + 1, k_n, 0.0, true);
		Md_t Q (k_A.rows (), k_n + 1);
		Md_t e1 (H.rows (), 1, 0.0);

		k_e1 = e1;

		k_H = H;
		k_Q = Q;

		Md_t v = k_Q.vec_view (0);
		v.pipe (x0);
		double bnorm = x0.vec_magnitude ();
		v /= bnorm;

		k_e1 (0, 0) = bnorm;
	}

	/*
	 * Verifies the basis of Kn is orthogonal.  If it breaks down
	 * then either restart, re-orthogonal or use Walker's Householder
	 * version to compute Kn.
	 *
	 */
	bool Orthogonal (void)
	{
		Md_t I (k_i + 2, k_i + 2, 1.0);
		Md_t Vm = k_Q.view (0, 0, k_A.rows (), k_i + 2);
		Md_t _I = transpose (Vm) * Vm;

		return I.equal_eps (_I, 1e-10);
	}

	int RunArnoldi (int runs);
};

int Krylov_t::RunArnoldi (int runs)
{
	Md_t v;
	Md_t qi;
	Md_t hi;
	Md_t __x;
	Md_t r;
	Md_t r0;
	int rows = k_Q.rows ();

	if (k_i + runs > k_n)
		runs = k_n - k_i;

	for (int i = 0; i < runs; ++i, ++k_i)
	{
		v = k_Q.vec_view (k_i + 1);
		qi = k_Q.vec_view (k_i);
		/*
		 * This matrix vector product is critical to performance.  It is
		 * the most expensive operation in the procedure.
		 *
		 */
		v.pipe (k_A * qi);

		hi = k_H.vec_view (k_i);

		/*
		 * Compute the projections for the Hessenberg in column i
		 *
		 */

		hi.pipe (transpose (k_Q) * v);

		/*
		 * vi -= ∑ Hj * qj
		 *
		 * Modified Gram-Schmidt process
		 *
		 */

		double *Qptr = k_Q.raw ();
		double *vptr;
		double alpha;

		for (int j = 0; j <= k_i; ++j)
		{
			vptr = v.raw ();
			alpha = k_H (j, k_i);

			for (int k = 0; k < rows; ++k)
				*vptr++ -= alpha * *Qptr++;
		}

		k_H (k_i + 1, k_i) = v.vec_magnitude ();
		if (k_H (k_i + 1, k_i) == 0)
			return k_i;

		v /= k_H (k_i + 1, k_i);
	}

	return k_i;
}

#endif // header inclusion

