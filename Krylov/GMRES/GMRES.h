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

#ifndef __DJS_GMRES_GIVENS__H__
#define __DJS_GMRES_GIVENS__H__

#include <Krylov.h>

/*
 * An example use of Krylov_t by implementing a naive Generalized Minimum 
 * Residual (GMRES).
 *
 * It is a simple and naive implementation intended to demonstrate how GMRES
 * works in code.
 *
 */

class GMRES_t : private Krylov_t
{
	int				gm_restarts;
	double			gm_residue;

	void init ()
	{
	}

	Md_t step (double &);
	double rotate (void);

public:

	GMRES_t (int n, Ms_t &A, Md_t &b) :
		Krylov_t (A, b, n),
		gm_restarts (10),
		gm_residue (0.5)
	{
		init ();
	}

	GMRES_t (int n, Ms_t &A, Md_t &b, int restarts) :
		Krylov_t (A, b, n),
		gm_restarts (restarts),
		gm_residue (0.5)
	{
		init ();
	}

	~GMRES_t ()
	{
	}

	void SetTolerance (double residue)
	{
		gm_residue = residue;
	}

	bool Solve (Md_t &, double &);
};

bool
GMRES_t::Solve (Md_t &xm, double &residue)
{
	xm = Md_t (k_b.rows (), 1, 0.0);
	bool rc = false;
	double last = DBL_MAX;
	int restarts = 0;
	Md_t _x;
	Md_t r;

//	for (int i = 0; i < gm_restarts; ++i)
	while (rc == false)
	{
		_x = step (residue);

		xm += _x;

		r = k_b - k_A * xm;

		if (residue <= gm_residue)
		{
			rc = true;
			break;
		}

		if (residue == last)
		{
			// GMRES has stalled.  See Saad for further suggestions.
			k_n += 50;
			printf ("Increasing m = %d\n", k_n);
		}

		last = residue;

		Restart (r, k_n);

		++restarts;

		if ((restarts % 10) == 0)
			printf ("|r| = %f\n", residue);
	}

	return rc;
}

/*
 * GMRES (Saad, chapter 6, Iterative Methods for Sparse Linear Systems)
 *
 * Find an aproximation, x*, to Ax = b based on AV = VH, where V spans
 * Kn (Krylov subspace generated from A).
 *
 *    i) Build Kn
 *   ii) find x* that is a member of Kn, so x* = Vy
 *  iii) Ax* = AVy = b, AV = VH, so VHy = b
 *   iv) Use least squares to minimize |VHy - b|, so Hy = Vtb
 *    v) b is orthogonal to V by construction (Arnoldi) so we have Hy = |b|e1
 *   vi) x* = Vy
 *
 * H is a Hessenberg matrix and much smaller than A, so this is very fast.
 *
 */

Md_t
GMRES_t::step (double &residue)
{
	int runs;

	runs = RunArnoldi (k_n);
	if (runs < k_n) // broke down
	{
		residue = nan (NULL);
		return Md_t ();
	}

	residue = rotate (); // H -> upper triangular (R)
	Md_t H = k_H.view (0, 0, k_i, k_i, false);
	Md_t Q = k_Q.view (0, 0, k_A.rows (), k_i);
	Md_t y = k_e1.view (0, 0, k_i, 1);
	y = H.find_x (y); // Hx = b', where H is upper triangular

	Md_t _x = Q * y;

	return _x;
}

/*
 * Use Givens rotations to transform Hessenberg matrix to upper triangular.
 * GVL (4th edition) Section 5.1.8
 *
 * |b|e1 is taken care of while computing the upper triangular.
 *
 */
double 
GMRES_t::rotate (void)
{
	int m = k_i;
	double tmp[2];
	double denom;
	double ci;
	double si;

	for (int i = 0; i < m; ++i)
	{
		denom = hypot (k_H(i, i), k_H(i + 1, i));
		ci = k_H(i, i) / denom;
		si = k_H(i + 1, i) / denom;

		tmp[0] = ci * k_e1(i, 0) + si * k_e1(i + 1, 0);
		tmp[1] = -si * k_e1(i, 0) + ci * k_e1(i + 1, 0);

		k_e1 (i, 0) = tmp[0];
		k_e1 (i + 1, 0) = tmp[1];

		k_H (i, i) = ci * k_H(i, i) + si * k_H(i + 1, i);
		k_H (i + 1, i) = 0.0;

		for (int j = i + 1; j < m; ++j) 
		{
			tmp[0] = ci * k_H(i, j) + si * k_H(i + 1, j);
			tmp[1] = -si * k_H(i, j) + ci * k_H(i + 1, j);
			k_H (i, j) = tmp[0];
			k_H (i + 1, j) = tmp[1];
		}
	}

	return fabs (k_e1(m, 0));
}

#endif // header inclusion

