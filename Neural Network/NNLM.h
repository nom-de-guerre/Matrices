#ifndef __LEVENBERG_MARQUARDT__NN_H__
#define __LEVENBERG_MARQUARDT__NN_H__

#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>

#include <data.h>

#include <matrix.h>

typedef Matrix_t<double> Md_t;

#define RECTIFIER(X) log (1 + exp (X))
#define SIGMOID_FN(X) (1 / (1 + exp (-X))) // derivative of rectifier

#ifdef __TANH_ACT_FN
#define ACTIVATION_FN(X) tanh(X)
#define DERIVATIVE_FN(Y) (1 - Y*Y)
#else
#define ACTIVATION_FN(X) SIGMOID_FN(X)
#define DERIVATIVE_FN(Y) (Y * (1 - Y))
#endif

// RPROP+ update parameters
#define DELTA0				1e-2
#define DELTA_MIN			1e-8
#define DELTA_MAX			50

#define ETA_PLUS			1.2
#define ETA_MINUS			0.5

// Levenberg-Marquardt parameters for dilating the spectrum
#define MU_INIT		0.01
#define MU_THETA    10

#ifndef SIGN
#define SIGN(X) (signbit (X) != 0 ? -1.0 : 1.0)
#endif

typedef double weight_t;
struct ModelSLP_t;

struct perceptron_t
{
	int					p_N;
	double				*p_weights;
	double				*p_Ei;
	double				*p_error;
	double				*p_deltaW;
	double				p_delta;
	double				p_iprod;
	double				p_output;
#ifdef __RELU
	bool				p_RELU;
#endif

	perceptron_t (void) :
		p_weights(NULL),
		p_Ei (NULL),
		p_deltaW (NULL)
#ifdef __RELU
	, p_RELU (false)
#endif
	{
	}

	perceptron_t (int N)
	{
		init (N);
	}

	~perceptron_t (void)
	{
		if (p_weights == NULL)
			return;

		delete [] p_weights;
		delete [] p_Ei;
		delete [] p_error;
		delete [] p_deltaW;
	}

#ifdef __RELU
	void setRELU (void) 
	{ 
		p_RELU = true; 
	}

	bool RELU (void) 
	{ 
		return p_RELU; 
	}
#endif

	void init (int N)
	{
		// N accounts for the bias, so + 1
		if (p_weights == NULL)
		{
			p_N = N;
			p_weights = new double [N];
			p_Ei = new double [N];
			p_error = new double [N];
			p_deltaW = new double [N];

		} else
			assert (p_N == N);

		// Glorot for sigmoid: W ~ U[-r, r], r = (6/(f_in + f_out))^0.5
		double r = sqrt (6 / ((double)  2 * N + 1));

		for (int i = 0; i < p_N; ++i)
		{
			double w = (double) rand () / RAND_MAX;

			p_weights[i] = 2 * r * w - r;
			p_Ei[i] = 0;
			p_error[i] = 0;
			p_deltaW[i] = DELTA0;
		}
	}

	void rprop (void)
	{
		for (int i = 0; i < p_N; ++i)
			rprop (i);
	}

	void rprop (int index);
	void bprop (perceptron_t *, 
			perceptron_t *, 
			const int, 
			const int, 
			const int, 
			const int, 
			Md_t &);

	double signal (double *xi)
	{
		p_iprod = 0;
		for (int i = 0; i < p_N; ++i)
			p_iprod += xi[i] * p_weights[i];

#ifdef __RELU
		if (p_RELU)
			p_output = RECTIFIER (p_iprod);
		else
			p_output = ACTIVATION_FN (p_iprod);
#else
		p_output = ACTIVATION_FN (p_iprod);
#endif

		return (p_output);
	}

	double signal (void)
	{
		return p_output;
	}

	void reset_weight (void)
	{
		double w = (double) rand () / RAND_MAX;
		double r = sqrt (6 / (p_N + 1));

		double error = 0;
		int weight = -1;

		for (int i = 0; i < p_N; ++i)
			if (fabs (p_error[i]) > error)
			{
				error = p_error[i];
				weight = i;
			}

if (weight < 0)
	assert (error == 0.0);
assert (weight < p_N);

//		int weight = rand () % p_N;

		p_weights[weight] = 2 * r * w - r;
		p_Ei[weight] = 0.0;
		p_error[weight] = 0.0;
		p_deltaW[weight] = DELTA0;
	}
};

template<typename T> class NNet_t
{
protected:

	int					n_steps;

	// morphology of the net
	int					n_Nin;
	int					n_Nout;
	int					n_levels;
	int					*n_width;		// array of lengths of n_nn
	perceptron_t		**n_nn;			// array of perceptron_t vectors

	double				*n_buffers[2];	// scratch space for compute ()

	int					n_RPROP_pending;
	int					n_modes;
	double				n_switch;		// RPROP+ threshold, switch to LM

	int					n_RPROP_pending_steps;
	int					n_LM_steps;

	double				n_mu;
	int					n_Nweights;
	Md_t				n_J;			// Jacobian for LM
	Md_t				n_Results;		// results

	double				n_halt;			// solution accuracy (sumsq, not derivs)
	double				n_error;

	int					n_me;

	enum e_tcodes { WORKING, FINISHED, STALLED };

	void PresentExamples (const DataSet_t * const, bool);
	double PresentExamplesLoss (const DataSet_t * const, bool);
	void UpdateWeights (Md_t &);
	bool Step (const DataSet_t * const training, double &);
	void Reset (const int);
	void RPROPStep (const DataSet_t * const training, double &);
	void LMStep (const DataSet_t * const training, double &);
	void InflectWeightSpace (void);

	double ComputeDerivative (const TrainingRow_t, const int);
	void Start (void);
	bool Halt (DataSet_t const * const);
	
public:

	/*
	 * levels is the number of layers, including the output.  width is an
	 * arrays specifying the width of each level.  e.g., { 1, 4, 1 }, is an SLP
	 * with a single input, 4 hidden and 1 output perceptron.
	 *
	 */
	NNet_t (int *width, int levels) :
		n_steps (0),
		n_Nin (width[0]),
		n_Nout (width[levels - 1]),
		n_levels (levels - 1), // no state for input
		n_RPROP_pending (0),
		n_modes (0),
		n_switch (0.8),
		n_RPROP_pending_steps (0),
		n_LM_steps (0),
		n_mu (MU_INIT),
		n_halt (9e-20),
		n_error (nan (NULL)),
		n_me (-1)
	{
		int max = 0;
		n_Nweights = 0;

		// width = # inputs, width 1, ..., # outputs

		for (int i = levels - 1; i > 0; --i)
			n_Nweights += width[i] * (width[i - 1] + 1); // + 1 for bias

		n_nn = new perceptron_t * [n_levels];
		n_width = new int [n_levels];

printf ("CONFIG:\t%d\t%d\t%d\n", n_Nin, n_Nout, n_Nweights);

		for (int i = 0; i <= n_levels; ++i)
		{
			if (width[i] > max)
				max = width[i];

			if (i == 0)
				continue; // ignore the inputs

			n_width[i - 1] = width[i];
			n_nn[i - 1] = new perceptron_t [width[i]];
		}

		n_buffers[0] = new double [max + 1]; // signal propigation
		n_buffers[0][0] = 1.0; // Bias
		n_buffers[1] = new double [max + 1]; // signal propigation
		n_buffers[1][0] = 1.0; // Bias
	}

	~NNet_t (void)
	{
		delete [] n_buffers[0];
		delete [] n_buffers[1];

		for (int i = 0; i < n_levels; ++i)
			delete [] n_nn[i];

		delete n_nn;
		delete [] n_width;
	}

	void setMSE (double mse)
	{
		n_halt = mse;
	}

	double whatIsHalt (void)
	{
		return n_halt;
	}

	void setRPROP (int RPROP)
	{
		n_RPROP_pending = RPROP;
	}

	bool TrainLifeSign (const DataSet_t * const, int, int);
	bool Train (const DataSet_t * const, int);
	bool UpdateTrain (const DataSet_t * const, int);
	void DisplayCurve (const DataSet_t * const);
	void DisplayWeights (void);
	double error (void)
	{
		return n_error;
	}

	double Compute (double *);
	inline double Loss (DataSet_t const *);
};

#include <NNLM.tcc>

#endif // header inclusion

