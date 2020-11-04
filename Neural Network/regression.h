#ifndef __NN_REGRESSION__H__
#define __NN_REGRESSION__H__

#include <NNLM.h>

class Regression_t : public NNet_t<Regression_t>
{

public:

	Regression_t (int *width, int levels) :
		NNet_t (width, levels)
	{
	}

	double bprop (const TrainingRow_t &, const int, int &);
	double f (double *);
	double error (DataSet_t const *);
	void Cycle (void) {}
	bool Test (DataSet_t const * const);
};

double Regression_t::f (double *x)
{
#if 0
	for (int i = 1; i <= n_Nout; ++i)
		result[i - 1] = output[i];
#endif

	/*
	 * x[0] is the bias - 1.0
	 *
	 */
	return x[1];
}

double Regression_t::bprop (const TrainingRow_t &x, const int row, int &dEj)
{
	double Result;
	double error = 0;

	Result = Compute (x);

	double y; 			// y = ak below
	double delta_k;
	double aj;
	double dAct;
	perceptron_t *outlayer = n_nn[n_levels - 1];

	dEj = 0;

	for (int output_i = 0; output_i < n_Nout; ++output_i)
	{
		y = Result; // s[output_i];

		// + 1 because output_i is the input to Compute (0)
		error += y - x[output_i + 1];

#ifdef __RELU
		if (outlayer->RELU ()
			dAct = SIGMOID_FN (outlayer->p_iprod);
		else
			dAct = DERIVATIVE_FN (y);
#else
		dAct = DERIVATIVE_FN (y);
#endif
 
		if (n_RPROP_pending)
			delta_k = (y - x[output_i + 1]) * dAct;
		else
			delta_k = dAct;

		aj = 1; // the bias

		outlayer[output_i].p_delta = delta_k;

		/*
		 * initiate the recurrence
		 *
		 */
		for (int i = 0; i < outlayer[output_i].p_N; ++i, ++dEj)
		{

		// output layer, first input is always the bias, and virtual
		// ∆wj = (ak - t)ak(1 - ak)aj = delta_k * aj

			if (i) // not the (virtual) bias
				aj = n_nn[0][i - 1].signal ();

			n_J(row, dEj) = delta_k * aj;

			outlayer[output_i].p_error[i] += delta_k * aj;
		}
	}

	return Result;
}

double Regression_t::error (DataSet_t const * tp)
{
	for (int i = 0; i < tp->t_N; ++i)
		n_Results (i, 0) -= (*tp)[i][1];

	return n_Results.vec_dotT ();
}

bool Regression_t::Test (DataSet_t const * const tp)
{
	if (n_error <= n_halt)
		return true;

	return false;
}

#endif // header inclusion

