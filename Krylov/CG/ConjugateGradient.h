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

#ifndef __DJS_CONJUGATE_GRAD__H__
#define __DJS_CONJUGATE_GRAD__H__

#define __DEBUG
#include <matrix.h>

typedef Matrix_t<double> Md_t;

struct ConjugateGrad_t
{
	Md_t			cg_A;
	Md_t			cg_b;
	Md_t			cg_x;
	Md_t			cg_r;
	Md_t			cg_p;
	double			cg_rho;
	double			cg_rhoMinus;

	double			cg_halt;
	int				cg_step;

public:

	ConjugateGrad_t (Md_t &A, Md_t &b) :
		cg_A (A),
		cg_b (b),
		cg_halt (MACH_EPS * cg_b.vec_magnitude ()),
		cg_step (0)
	{
		Reset ();
	}

	ConjugateGrad_t (Md_t &A, Md_t &b, double halt) :
		cg_A (A),
		cg_b (b),
		cg_halt (halt),
		cg_step (0)
	{
		Reset ();
	}

	~ConjugateGrad_t (void)
	{
	}

	void Reset (void)
	{
		cg_x = cg_b;
		cg_r = cg_b - cg_A * cg_x;
		cg_rho = cg_r.vec_dotT ();
	}

	void Compute (void)
	{
		int N = cg_A.rows () - 1;

		for (int i = 0; i < N && cg_halt < sqrt (cg_rho); ++i)
			Step ();
	}

	void Step (void)
	{
		++cg_step; 

		if (cg_step == 1)

			cg_p = cg_r;

		else {

			double tau = cg_rho / cg_rhoMinus;
			cg_p = cg_r + tau * cg_p;
		}

		Md_t w = cg_A * cg_p;
		double mu = cg_rho / cg_p.vec_dot (w);
		cg_x = cg_x + mu * cg_p;
		cg_r = cg_r - mu * w;
		cg_rhoMinus = cg_rho;
		cg_rho = cg_r.vec_dotT ();
	}

	Md_t Answer (void)
	{
		return cg_x;
	}
};

#endif // header inclusion

