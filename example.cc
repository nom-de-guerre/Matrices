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

#include <matrix.h>

typedef Matrix_t<double> Md_t;

#define __DIM 10

int main ()
{
	Md_t A (__DIM, __DIM);		// 10x10 matrix - undefined contents...

	// ...so fill it
	for (int i = 0, count = 0; i < __DIM; ++i)
		for (int j = 0; j < __DIM; ++j, ++count)
			A(i, j) = count;

	A.display ("A");

	Md_t B = A;					// B and A share memory - copy-on-write (CoW)
	B.display ("B - shares view and data");
	B(__DIM >> 1, __DIM >> 1) = 3.1415926; // This will trigger a CoW
	// second argument is number of decimal places.  default is 2
	B.display ("B now has own copy of memory (5, 5) = π", "3");
	A.display ("A");

	/*
	 * Create a view of a column in A.  The default for vec_view is
	 * write-in-place, so writing the vector will be reflected in A.
	 * It can be turned off by calling vec_view (column, false).  It
	 * can also be turned off/on, like any matrix, with set_CoW/set_WiP.
	 *
	 */
	Md_t x = A.vec_view (__DIM >> 1); // x and A share memory
	Md_t z = x;				// z, x and A share memory
	x.display ("x = A(, 5)");
	z.display ("z = x");
	z.set_CoW (); 			// set to copy-on-write (sharing until update)
	z(1, 0)	= 2.71828182;	// not reflected in A (z now has private memory)
	x(0, 0)	= 2.71828182;	// default is write in place, reflected in A
	z.display ("z CoW");

	A.display ("A: x modified me");
	z.display ("z: not me though!");

	x.viewOriginal (); // explode x to full underlying size
	x.display ("independent view of A");

	// create a matrix with an initial value of 0.5 - the last argument
	// means zero the matrix and initialize the diagonal
	Md_t I (__DIM, __DIM, 0.5, false);
	I *= 2;
	I.display ("Identity matrix");
	Md_t C = A * I;

	(C - A).display ("Should be all naughts");

	Md_t D = A;
	D.copy (); // force a CoW - so changes will not be reflected in A
	/*
	 * Create a window into E, from (2, 2) to (4, 4)
	 * default is write in place (changes relfected in E) - last arg
	 * should be false for a CoW window
	 *
	 */
	Md_t E = D.view (2, 2, 2, 2, true);
	E.display ("E");
	E(0, 0) = 0.99;
	D(3, 6) = 7;
	D.display ("(2, 2) is 0.99 now");

	E.set_CoW ();
	E(1, 1) = 0.99;
	E.display ("E has its own memory now");
	D.display ("(3, 3) untouched");
}


