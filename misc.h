#ifndef	__MATRIX_MISC__H__
#define	__MATRIX_MISC__H__

#include <algorithm>

template<typename T> void makeHeap (T A[], const int N, const int place)
{
	int left = 2 * place + 1;
	int right = left + 1;
	int swap = place;

	if (left < N && A[left] > A[place])
		swap = left;

	if (right < N && A[right] > A[swap])
		swap = right;

	if (swap == place)
		return;

	std::swap (A[place], A[swap]);

	makeHeap (A, N, swap);
}

template<typename T> void heapSort (T A[], const int N)
{
	for (int i = N - 1; i >= 0; --i)
		makeHeap (A, N, i);

	for (int i = N - 1; i >= 0; --i)
	{
		std::swap (A[0], A[i]);
		makeHeap (A, i, 0);
	}
}

template<typename T> bool search (T A[], const int N, const T key)
{
	int left = 0;
	int right = N - 1;
	int pivot;
	bool rc = false;

	while (left <= right)
	{
		pivot = left + (right - left) / 2;

		if (A[pivot] == key)
		{
			rc = true;
			break;
		}

		if (A[pivot] > key)
			right = pivot - 1;
		else
			left = pivot + 1;
	}

	return rc;
}

#endif // header

