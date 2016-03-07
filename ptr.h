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

#ifndef __MATRIX_REFERENCES__PTR_H__
#define __MATRIX_REFERENCES__PTR_H__

/*
 * This is a basic smart pointer implementation.  It is used by the matrix
 * code to implement CoW and views into common memory.
 *
 */

#define VASSERT(X) if (!(X)) {						\
	printf ("ASSERT %s:%d\n", __FILE__, __LINE__);	\
	*((long *) 0) = 0xdead;							\
}

#ifdef __DEBUG_PARANOID
extern int __total_objs;
#define VERIFY_CLEAN VASSERT (__total_objs == 0)
#else
#define VERIFY_CLEAN
#endif

template<class T> struct __ptr_t
{
	T		*_r_data;
	int		_r_ref;

	__ptr_t()
	{
		_r_data = 0;
		_r_ref = 0;
	}

	__ptr_t(T *ptr)
	{
		_r_data = ptr;
		_r_ref = 1;
#ifdef __DEBUG_PARANOID
	++__total_objs;
#endif
	}

	~__ptr_t() 
	{ 
		if (_r_data) 
			delete _r_data; 
#ifdef __DEBUG_PARANOID
	--__total_objs;
#endif
	}

	void 		pget()
	{
		_r_ref++;
	}

	bool		pput()
	{
		_r_ref--;

		if (_r_ref == 0)
			return true;

		return false;
	}	
};

template<typename T> class ptr_t
{
protected:

	__ptr_t<T>	*r_data;

	void pput ()
	{
		if (r_data == 0)
			return;

		if (r_data->pput ()) 
			delete r_data;
	}

public:

	ptr_t()
	{
		r_data = 0;
	}

	ptr_t(ptr_t<T> &X)
	{
		r_data = X.r_data;
		r_data->pget();
	}

	ptr_t(const ptr_t<T> &X)
	{
		r_data = X.r_data;
		r_data->pget();
	}

	ptr_t(T *X)
	{
		r_data = new __ptr_t<T> (X);
	}

	~ptr_t()
	{
		pput ();
	}

	T		*get()
	{
		return r_data->_r_data;
	}

	ptr_t<T>	&operator=(T *X)
	{
		pput ();

		r_data = new __ptr_t<T> (X);

		return *this;
	}

	ptr_t<T>	&operator=(ptr_t<T> &X)
	{
		pput ();

		r_data = X.r_data;
		if (r_data) r_data->pget ();

		return *this;
	}

	T		*operator->()
	{
		return (r_data ? r_data->_r_data : 0);
	}

	bool valid ()
	{
		return r_data != 0;
	}

	bool exclusive ()
	{
		return valid () && (r_data->_r_ref == 1);
	}
};

#endif // header inclusion

