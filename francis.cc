/*

Copyright (c) 2015, Douglas Santry
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

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

#include <math.h>
#include <float.h>

#include <utility>

#include <francis.h>

/*
 * Accepts an arbitrary matrix, it will be put into Hessenberg form
 *
 * This will destroy the argument, A.
 *
 */

int EigenFrancis_t::CalcEigenValues (Md_t &A)
{
	A.HessenbergSimilarity ();

	return CalcEigenValuesHessenberg (A);
}

/*
 * Calculate the Eigen values of A.  The eigen values are stored in the class.
 *
 * This will destroy the argument, A.
 *
 */

int EigenFrancis_t::CalcEigenValuesHessenberg (Md_t &A)
{
	Md_t Francis;
	int shift;
	int rows = A.rows ();

	ef_totalIterations = 0;
	if (ef_EigenValues)
		delete [] ef_EigenValues;

	ef_EigenValues = new conj_t [A.rows ()];
	ef_N = 0;

	A.set_WiP ();

	Francis = A;

	for (int i = rows; i > 0; i -= shift)
	{
		Md_t Ai = Francis;
		Ai.set_WiP ();

		// shift and iterate (QR) for Eigenvalue
		shift = IterateAndShift (Ai);
		if (shift < 0)
			break;	// didn't converge - we're done

		Francis = Ai;
	}

	return ef_N;
}

void ComplexEigen (Md_t &A, int row, conj_t &conj)
{
	/*
	 * use quadratic equation to solve sub-block 
	 * (roots of 2x2 characteristic polynomial)
	 *
	 */

	double b = -(A[row + 1][row + 1] + A[row][row]);
	double c = A[row + 1][row + 1] * A[row][row] -
					A[row][row + 1] * A[row + 1][row];

	conj.real = -b / 2; // real
	conj.imag = sqrt (fabs (b * b - 4 * c)) / 2; // imag (+/-)
}

int EigenFrancis_t::DetectConvergence (Md_t &H)
{
	int rows = H.rows ();
	int index = -1;

	for (int i = rows - 1; i > 0; --i) {

		// equation 7.5.4 in GVL4
		double diag = MACH_EPS * (fabs (H[i][i]) + fabs (H[i - 1][i - 1]));

		if (H[i][i - 1] == 0 || fabs (H[i][i - 1]) <= diag) {

			H[i][i - 1] = 0;
			return i;
		}
	}

	return index;
}

void EigenFrancis_t::FrancisStep (Md_t &A, double shift)
{
	Md_t e1 (A.rows (), 1, 0.0);
	int last = A.rows () - 1;
	double s = A[last - 1][last - 1] + A[last][last];
	double t = A[last - 1][last - 1] * A[last][last] -
				A[last - 1][last] * A[last][last - 1];

	// From GVL4
	e1[0][0] = A[0][0] * A[0][0] + A[0][1] * A[1][0] - s * A[0][0] + t;
	e1[1][0] = A[1][0] * (A[0][0] + A[1][1] - s);
	e1[2][0] = A[2][1] * A[1][0];

	ApplyBulge (A, e1);
	ChaseBulge (A);

	++ef_totalIterations;
}

int EigenFrancis_t::SchurSubMatrix (Md_t &A, int index, conj_t EigenValues[])
{
	// Matrix_t<>::CoW is expensive => cache values
	double a = A[index][index];
	double b = A[index][index + 1];
	double c = A[index + 1][index];
	double d = A[index + 1][index + 1];
	double tmp = a - d;
	double p = 0.5 * tmp;
	double bcmax = fmax (fabs(b), fabs(c));
	double bcmis = fmin (fabs(b), fabs(c)) * SIGN(b) * SIGN(c);
	double scale = fmax (fabs(p), bcmax);
	double z = (p / scale) * p + (bcmax / scale) * bcmis;
	double tau;
	double cs, sn;

	if (z >= 4.0 * MACH_EPS) {

		/* real eigenvalues, compute a and d */
		z = p + SIGN(p) * fabs(sqrt(scale) * sqrt(z));
		a = d + z;
		d -= (bcmax / z) * bcmis;

		/* compute b and the rotation matrix */
		tau = sqrt (c * c + z * z);
		cs = z / tau;
		sn = c / tau;
		b -= c;
		c = 0.0;

		EigenValues[0].real = d;
		EigenValues[1].real = a;
		EigenValues[0].imag = EigenValues[1].imag = 0;

		return 2;

	} else {

		ComplexEigen (A, index, EigenValues[0]);

		return 1;
	}
}

int EigenFrancis_t::IterateAndShift (Md_t &A)
{
	int rows = A.rows ();
	int last = rows - 1;
	bool working = true;
	int iterations = 0;
	int deflate = -10000;
	int pivot;

	if (rows == 2) {

		ef_N += SchurSubMatrix (A, 0, ef_EigenValues + ef_N);

		deflate = 2;
		working = false;

	} else if (rows == 1) {

		ef_EigenValues[ef_N].real = A[0][0];
		ef_EigenValues[ef_N].imag = 0;
		++ef_N;

		deflate = 1;
		working = false;
	}

	while (working)
	{
		++iterations;

		/*
		 * This choice of when to give up is entirely arbitrary. More
		 * work is required to make this better.
		 *
		 * Need a means of introducing an extraordinary shift when we
		 * get stuck.
		 *
		 */
		if (iterations == last + 30)
			return -1;

		pivot = DetectConvergence (A);

		if (pivot == -1) {

			double shift;

			if ((iterations % 10) == 0)
				shift = A[last][last]; // Rayleigh Quotient
			else {

				SchurSubMatrix (A, last - 1, ef_EigenValues + ef_N);
				shift = ef_EigenValues[ef_N].real;
			}

			FrancisStep (A, shift);

			continue;
		}

		if (pivot == last) {

			ef_EigenValues[ef_N].real = A[last][last];
			ef_EigenValues[ef_N].imag = 0;
			++ef_N;

			A = A.view (0, 0, last, last);

			deflate = 1;
			working = false;

			continue;

		} else if (pivot == last - 1) {

			ef_N += SchurSubMatrix (A, pivot, ef_EigenValues + ef_N);

			A = A.view (0, 0, last - 1, last - 1);
			deflate = 2;
			working = false;

			continue;

		} else if (pivot == 1) {

			ef_EigenValues[ef_N].real = A[0][0];
			ef_EigenValues[ef_N].imag = 0;
			++ef_N;

			A = A.view (1, 1, last, last);

			deflate = 1;
			working = false;

			continue;

		} else if (pivot == 2) {

			ef_N += SchurSubMatrix (A, 0, ef_EigenValues + ef_N);

			A = A.view (2, 2, rows - 2, rows - 2);

			deflate = 2;
			working = false;

			continue;

		} else {

			// de-couple

			Md_t UL = A.view (0, 0, pivot, pivot);

			Md_t LR = A.view (pivot, pivot, rows - pivot, rows - pivot);

			CalcEigenValues (LR);
			CalcEigenValues (UL);

			deflate = last + 1;
			working = false;

			continue;
		}
	}

	return deflate;
}

void EigenFrancis_t::makeHeap (int place, int var_N)
{
    int left = 2 * place + 1;
    int right = left + 1;
    int swap = place;

    if (left < var_N && ef_EigenValues[left] > ef_EigenValues[swap])
        swap = left;

    if (right < var_N && ef_EigenValues[right] > ef_EigenValues[swap])
        swap = right;

    if (swap != place) {

        conj_t tmp = ef_EigenValues[place];
        ef_EigenValues[place] = ef_EigenValues[swap];
        ef_EigenValues[swap] = tmp;

        makeHeap (swap, var_N);
    }
}

void EigenFrancis_t::SortEigenValues (void)
{
    int var_N = ef_N;

	if (ef_N <= 0)
		return;

    for (int i = ef_N - 1; i >= 0; --i)
        makeHeap (i, var_N);

    for (int i = ef_N - 1; i; --i) {

        conj_t tmp = ef_EigenValues[0];
        ef_EigenValues[0] = ef_EigenValues[i];
        ef_EigenValues[i] = tmp;
        --var_N;

        makeHeap (0, var_N);
    }
}


void EigenFrancis_t::ChaseBulge (Md_t &A)
{
	int runs = A.rows () - 2;

	for (int i = 0; i < runs; ++i)
		RawStep (A, i);
}

bool EigenFrancis_t::RawStep (Md_t &A, int step)
{
	int rows = A.rows ();
	int start = step + 1;
	double alpha;
	double beta;
	int halt = (start + 3 < rows ? start + 3 : rows);
	double * __restrict _A = A.raw ();
	double * __restrict w = new double [rows];
	double * __restrict Aprodw = new double [rows];
	int prows = A.prows ();

	beta = alpha = A[start][step] * A[start][step];
	w[start] = A[start][step];

	for (int i = start + 1, k = step * prows + i; i < halt; ++i, ++k) {

		w[i] = _A[k];
		alpha += w[i] * w[i];
	}

	beta = alpha - beta;
	alpha = sqrt (alpha);

	if (w[start] > 0.0)
		w[start] += alpha;
	else
		w[start] -= alpha;

	beta = 2 / (beta + w[start] * w[start]);

	// A' = (I - vvT)A
	for (int c = step; c < rows; ++c)
	{
		Aprodw[c] = 0;
		for (int i = start, k = c * prows + i; i < halt; ++i, ++k)
			Aprodw[c] += w[i] * _A[k];

		Aprodw[c] *= beta;
	}

	for (int r = start; r < halt; ++r)
		for (int c = step, k = c * prows + r; c < rows; ++c, k += prows)
			_A[k] -= w[r] * Aprodw[c];

	// A'' = A'(I - vvT)
	for (int r = 0; r <= halt && r < rows; ++r) {

		Aprodw[r] = 0.0;

		for (int c = start, k = c * prows + r; c < halt; ++c, k += prows)
			Aprodw[r] += _A[k] * w[c];

		Aprodw[r] *= beta;
	}

	for (int r = 0; r <= halt && r < rows; ++r)
		for (int c = start, k = c * prows + r; c < halt; ++c, k += prows)
			_A[r + c * prows] -= Aprodw[r] * w[c];

	delete [] Aprodw;
	delete [] w;

	return true;
}

bool EigenFrancis_t::ApplyBulge (Md_t &A, Md_t &x)
{
	int rows = A.rows ();
	double alpha;
	double beta;
	int halt = (3 < rows ? 3 : rows);
	double * __restrict _A = A.raw ();
	double * __restrict w = new double [rows];
	double * __restrict Aprodw = new double [rows];
	double * __restrict p = x.raw ();

	int prows = A.prows ();

	beta = alpha = p[0] * p[0];
	w[0] = p[0];

	for (int i = 1; i < halt; ++i) {

		w[i] = p[i];
		alpha += w[i] * w[i];
	}

	beta = alpha - beta;
	alpha = sqrt (alpha);

	if (w[0] > 0.0)
		w[0] += alpha;
	else
		w[0] -= alpha;

	beta = 2 / (beta + w[0] * w[0]);

	// A' = (I - vvT)A
	// Aprodw = vTA
	for (int c = 0; c < rows; ++c)
	{
		Aprodw[c] = 0;
		for (int i = 0, k = c * prows + i; i < halt; ++i, ++k)
			Aprodw[c] += w[i] * _A[k];

		Aprodw[c] *= beta;
	}

	// A -= v*Aprodw
	for (int r = 0; r < halt; ++r)
		for (int c = 0, k = c * prows + r; c < rows; ++c, k += prows)
			_A[k] -= w[r] * Aprodw[c];

	// A'' = A'(I - vvT)
	for (int r = 0; r <= halt && r < rows; ++r) {

		Aprodw[r] = 0.0;

		for (int c = 0, k = c * prows + r; c < halt; ++c, k += prows)
			Aprodw[r] += _A[k] * w[c];

		Aprodw[r] *= beta;
	}

	for (int r = 0; r <= halt && r < rows; ++r)
		for (int c = 0, k = c * prows + r; c < halt; ++c, k += prows)
			_A[r + c * prows] -= Aprodw[r] * w[c];

	delete [] Aprodw;
	delete [] w;

	return true;
}

/*
 * Implements inverse iteration.  If the inverse of A is _A, then
 * _Ax_i = Ax_(i + 1), --> x_i = Ax_(i + 1), factor with QR to find x_(i+1)
 *
 */
bool EigenFrancis_t::FindEigenVectorReal (double lambda, Md_t &A, Md_t &u)
{
    double halt = 10 * MACH_EPS * A.norm_inf ();
    if (isnan (halt))
        return false;

	int rows = A.rows ();
	int iterations = rows;
    Md_t _A = A;
    Md_t scratch_A;
	Md_t d;
	double r_inf;
    Md_t x;

	for (int i = 0; i < rows; ++i)
		_A[i][i] -= lambda;

    u.set_WiP ();

    // _Ax_n = x_n+1 -> Ax_n+1=x_n, factor with QR
    while (true) {

        scratch_A = _A;
        scratch_A.copy ();
        x = scratch_A.solve_b (u);
        u = x.vec_norm ();

		d = _A * u;

		r_inf = DBL_MIN;
		for (int i = 0; i < rows; ++i)
			if (r_inf < fabs (d[i][0]))
				r_inf = fabs (d[i][0]);

#if 0
       	Md_t Rayleigh = A * u;
       	Rayleigh = u.transpose () * Rayleigh;
       	double quotient = Rayleigh[0][0];
        double r = ((A * u - lambda * u).vec_magnitude ());
		double R = lambda - quotient;

		printf ("Halt conditions %e\t%e\t%e\n", r, halt, A.norm_inf ());
		printf ("Rayleigh %e %f\n", quotient / lambda, quotient);
#endif

        if (r_inf <= halt)
            break;

		--iterations;
		if (iterations < 0)
			return false;
    }

    return true;
}


