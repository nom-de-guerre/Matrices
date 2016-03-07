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

#ifndef __DJS_FRANCIS__H__
#define __DJS_FRANCIS__H__

#include <matrix.h>
typedef Matrix_t<double> Md_t;

/*
 * Type used to represent Eigen values
 *
 */

struct conj_t {

    double real;
    double imag;

	/*
	 * Used for sorting.
	 *
	 */
	double modulus (void) {

		return sqrt (real * real + imag * imag);
//		return real;
	}

	bool operator== (conj_t &X) {

		if (real == X.real && imag == X.imag)
			return true;

		return false;
	}

	bool operator< (conj_t &X) {

		if (modulus () < X.modulus ())
			return true;

		return false;
	}

	bool operator> (conj_t &X) {

		if (modulus () > X.modulus ())
			return true;

		return false;
	}
};

/*
 * Encapsulates the Francis algoritm.  Assumes the matrices are of type
 * double.
 *
 */

class EigenFrancis_t {

	int				ef_totalIterations;
	Md_t			ef_A;

	int IterateAndShift (Md_t &);
	int DetectConvergence (Md_t &);
	void FrancisStep (Md_t &, double);
	int SchurSubMatrix (Md_t &, int, conj_t []);
	void ChaseBulge (Md_t &);
	bool RawStep (Md_t &, int);
	bool ApplyBulge (Md_t &, Md_t &);

	void makeHeap (int, int);

public:

	conj_t			*ef_EigenValues;
	int				ef_N;

	EigenFrancis_t (void) :
		ef_EigenValues (0)
	{
	}

	~EigenFrancis_t (void) {

		if (ef_EigenValues)
			delete [] ef_EigenValues;
	}

	int N_Iterations (void) {

		return ef_totalIterations;
	}

	void display (const char *name = 0)
	{
		if (name)
			printf ("%s:\n", name);

		for (int i = 0; i < ef_N; ++i)
			printf ("<%f, %f>\n", ef_EigenValues[i].real,
								ef_EigenValues[i].imag);
	}

	int CalcEigenValuesHessenberg (Md_t &); // general square matrix
	int CalcEigenValues (Md_t &); // already in Hessenberg form

	/*
	 * This will find the eigen vector associated with the eigen value
	 *
	 */
	static bool FindEigenVectorReal (double, Md_t &, Md_t &);

	void SortEigenValues (void);
};


#endif // header inclusion

