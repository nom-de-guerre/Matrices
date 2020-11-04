#ifndef __DataSetN_H__

#include <float.h>
#include <math.h>
#include <assert.h>

typedef double * TrainingRow_t;

struct DataSet_t
{
	int						t_N;
	int						t_Nin;
	int						t_Nout;
	TrainingRow_t			t_data;

	DataSet_t (int N, int Nin, int Nout) :
		t_N (N),
		t_Nin (Nin),
		t_Nout (Nout),
		t_data (new double [N * Nin + N * Nout])
	{
		assert (t_Nout == 1);
	}

	DataSet_t (int N, int Nin, int Nout, double *datap) :
		t_N (N),
		t_Nin (Nin),
		t_Nout (Nout),
		t_data (datap)
	{
	}

	~DataSet_t (void)
	{
		delete [] t_data;
	}

	TrainingRow_t entry (const int index) const
	{
		return t_data + index * (t_Nout + t_Nin);
	}

	TrainingRow_t operator[] (int index) const
	{
		return entry (index);
	}

	double Result (const int index) const
	{
		return *(t_data + index * (t_Nout + t_Nin) + t_Nin);
	}

	int N (void)
	{
		return t_N;
	}

	double Max (const int feature)
	{
		double best = -DBL_MAX;
		int stride = t_Nin + t_Nout;
		double *column = t_data + feature;

		for (int i = 0; i < t_N; ++i, column += stride)
			if (best < *column)
				best = *column;

		return best;
	}

	double Min (const int feature)
	{
		double best = DBL_MAX;
		int stride = t_Nin + t_Nout;
		double *column = t_data + feature;

		for (int i = 0; i < t_N; ++i, column += stride)
			if (best > *column)
				best = *column;

		return best;
	}

	double Mean (const int feature)
	{
		double sum = 0;
		int stride = t_Nin + t_Nout;
		double *column = t_data + feature;

		for (int i = 0; i < t_N; ++i, column += stride)
			sum += *column;

		return sum / t_N;
	}

	double Variance (const int feature)
	{
		double sum = 0;
		double sumsq = 0;
		double var;
		int stride = t_Nin + t_Nout;
		double *column = t_data + feature;

		for (int i = 0; i < t_N; ++i, column += stride)
		{
			sum += *column;
			sumsq += *column * *column;
		}

		sum /= t_N;
		sum *= sum;

		var = sumsq / (t_N - 1) - sum;

		return var;
	}

	void Center (const int feature)
	{
		double centre = Mean (feature);
		int stride = t_Nin + t_Nout;
		double *column = t_data + feature;

		for (int i = 0; i < t_N; ++i, column += stride)
			*column -= centre;
	}

	void Zscore (const int feature)
	{
	}
};

#endif // header inclusion

