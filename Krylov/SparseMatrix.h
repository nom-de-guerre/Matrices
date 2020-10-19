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

#ifndef __DJS_SPARSE_MATRIX__H__
#define __DJS_SPARSE_MATRIX__H__

#include <assert.h>

#include <matrix.h>
typedef Matrix_t<double> Md_t;

namespace SparseMatrix
{

#define __MAX_SMALL 25

struct row_t
{
	struct element_t
	{
		int column;
		double datum;
	};

	int					index;
	element_t 			columns[__MAX_SMALL];

	row_t (void) :
		index (0)
	{
	}

	double &operator[] (int column)
	{
		assert (index < __MAX_SMALL);

		columns[index].column = column;
		++index;
		return columns[index - 1].datum;
	}
};

class SparseMatrix_t
{
	int			mm_rows;
	int			mm_columns;
	row_t		*mm_matrix;

public:

	SparseMatrix_t (int rows, int columns) :
		mm_rows (rows),
		mm_columns (columns),
		mm_matrix (new row_t [rows])
	{
	}

	~SparseMatrix_t (void)
	{
		delete [] mm_matrix;
	}

	int rows (void)
	{
		return mm_rows;
	}

	int columns (void)
	{
		return mm_columns;
	}

	row_t &operator[] (int row)
	{
		return mm_matrix[row];
	}

	void display (const char *name = "")
	{
		double element;

		for (int i = 0; i < mm_rows; ++i)
		{
			row_t *p = mm_matrix + i;

			for (int j = 0, index = 0; j < mm_columns; ++j)
			{
				if (index < p->index && p->columns[index].column == j)
				{
					element = p->columns[index].datum;
					++index;

				} else
					element = 0.0;

				printf ("%.1f\t", element);
			}

			printf ("\n");
		}
	}

	Md_t Copy (void)
	{
		Md_t A (mm_rows, mm_columns);

		double element;

		for (int i = 0; i < mm_rows; ++i)
		{
			row_t *p = mm_matrix + i;

			for (int j = 0, index = 0; j < mm_columns; ++j)
			{
				if (index < p->index && p->columns[index].column == j)
				{
					element = p->columns[index].datum;
					++index;

				} else
					element = 0.0;

				A (i, j) = element;
			}
		}

		return A;
	}

	// optimized for matrix vector multiplication - incorrect for all else
	friend Md_t operator* (SparseMatrix_t &A, Md_t &v)
	{
		assert (v.columns () == 1);

		Md_t u (A.rows (), v.columns ());

		MatrixVectorProduct (A, v, u);

		return u;
	}

	friend void MatrixVectorProduct (SparseMatrix_t &A, Md_t &v, Md_t &u);

};

/*
 * Computes u = Av
 *
 */
void MatrixVectorProduct (SparseMatrix_t &A, Md_t &v, Md_t &u)
{
	double * __restrict x = v.raw ();
	double * __restrict y = u.raw ();
	int rows = A.rows ();

	for (int i = 0; i < rows; ++i)
	{
		*y = 0;
		for (int k = 0; k < A[i].index; ++k)
			*y += A[i].columns[k].datum * x[A[i].columns[k].column];

		++y;
	}
}

};

#endif // header inclusion

