#ifndef __NN_REGRESSION__H__
#define __NN_REGRESSION__H__

#include <math.h>

#include <NNLM.h>

#define NATURAL_NUMBER		2.718281828459045

class Softmax_t : public NNet_t<Softmax_t>
{
	double			*c_P;

public:

	Softmax_t (int *width, int levels) :
		NNet_t (width, levels)
	{
		c_P = new double [n_Nout];
		setRPROP (INT_MAX);
	}

	~Softmax_t (void)
	{
		delete [] c_P;
	}

	double bprop (const TrainingRow_t &, const int, int &);
	double f (double *);
	double error (DataSet_t const *);
	void Cycle (void);
	bool Test (DataSet_t const * const);

	int ComputeSoftmax (void);
};

double Softmax_t::f (double *x)
{
	assert (x[0] == 1); // the bias

	for (int i = 1; i <= n_Nout; ++i)
		c_P[i - 1] = x[i];

	return ComputeSoftmax ();
}

/*
 * Convert network outputs, c_P[], to softmax "probabilities"
 *
 */
int Softmax_t::ComputeSoftmax ()
{
	double denom = 0;
	double best = -100000000;
	double max = -10000000;
	int factor = -1;

	for (int i = 0; i < n_Nout; ++i)
	{
		if (c_P[i] > max)
			max = c_P[i];
	}

	for (int i = 0; i < n_Nout; ++i)
	{
		c_P[i] = exp (c_P[i] - max);
		denom += c_P[i];
	}

	for (int i = 0; i < n_Nout; ++i)
	{
		c_P[i] /= denom;

		if (c_P[i] > best)
		{
			best = c_P[i];
			factor = i;
		}
	}

	assert (factor > -1 && factor < n_Nout);

	return factor;
}

double Softmax_t::bprop (const TrainingRow_t &x, const int row, int &dEj)
{
	double loss;
	int answer = static_cast<int> (x[n_Nin]);

	Compute (x); // forces computation of Softmax Pi

	loss = -log (c_P[answer]);

#if 0
printf ("bprop\t%d\t%d\t%e\t%e\n", answer, response, c_P[answer], loss);
printf ("SM\t");
if (row == 4)
{
for (int i = 0; i < n_Nout; ++i)
	printf ("%f/%f\t", c_P[i], n_nn[n_levels - 1][i].signal ());
printf ("\n");
}
#endif
	double y; 			// y = ak below
	double delta_k;
	double aj;
	double dAct;
	double dL;

	perceptron_t *outlayer = n_nn[n_levels - 1];

	dEj = 0;

	for (int output_i = 0; output_i < n_Nout; ++output_i)
	{
		y = outlayer[output_i].signal ();

		dL = c_P[output_i];
		if (output_i == answer)
			dL -= 1;

		// + 1 because output_i is the input to Compute (0)

#ifdef __RELU
		if (outlayer->RELU ()
			dAct = SIGMOID_FN (outlayer->p_iprod);
		else
			dAct = DERIVATIVE_FN (y);
#else
		dAct = DERIVATIVE_FN (y);
#endif
 
		if (n_RPROP_pending)
			delta_k = dL * dAct;
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
		// ∂L/∂wj = dL · ak(1 - ak) · aj = delta_k * aj

			if (i) // not the (virtual) bias
				aj = n_nn[0][i - 1].signal ();

			n_J(row, dEj) = delta_k * aj;

			outlayer[output_i].p_error[i] += delta_k * aj;
		}
	}

	return loss;
}

double Softmax_t::error (DataSet_t const * tp)
{
	double loss = 0;
	for (int i = 0; i < n_Results.rows (); ++i)
		loss += n_Results (i, 0);

	loss /= tp->t_N;

	return loss;
}

void Softmax_t::Cycle (void)
{
}

bool Softmax_t::Test (DataSet_t const * const tp)
{
	int correct = 0;

	for (int i = 0; i < tp->t_N; ++i)
	{
		const TrainingRow_t p = (*tp)[i];

		if (Compute (p) == p[2])
			++correct;
	}

	if (correct == tp->t_N)
		return true;

	return false;
}

#endif // header inclusion

