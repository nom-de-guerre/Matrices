#ifndef __LEVENBERG_MARQUARDT__NN_H__
#define __LEVENBERG_MARQUARDT__NN_H__

#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>

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
typedef double * TrainingRow_t;
struct ModelSLP_t;

struct training_t
{
	int						t_N;
	int						t_Nout;
	TrainingRow_t			t_data;

	training_t (int N, int Nout) :
		t_N (N),
		t_Nout (Nout),
		t_data (new double [N + N * Nout])
	{
	}

	~training_t (void)
	{
		delete [] t_data;
	}

	TrainingRow_t operator[] (int index) const
	{
		return t_data + index * (t_Nout + 1);
	}
};

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
		double r = sqrt (6 / ((double) N + 1));

		for (int i = 0; i < p_N; ++i)
		{
			double w = 10 * (double) rand () / RAND_MAX - 5;

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

class NNet_t
{
	int					n_steps;

	// morphology of the net
	int					n_levels;
	int					*n_width;		// array of lengths of n_nn
	perceptron_t		**n_nn;			// array of perceptron_t vectors

	double				*n_buffers[2];	// scratch space for compute ()

	int					n_RPROP;
	int					n_modes;
	double				n_switch;		// RPROP+ threshold, switch to LM

	int					n_RPROP_steps;
	int					n_LM_steps;

	double				n_mu;
	int					n_Nweights;
	Md_t				n_J;

	double				n_halt;			// solution accuracy (sumsq, not derivs)
	double				n_error;

	int					n_me;

	enum e_tcodes { WORKING, FINISHED, STALLED };

	double ComputeDerivative (const TrainingRow_t, const int);
	Md_t PresentExamples (const training_t * const, bool);
	void UpdateWeights (Md_t &);
	bool Step (const training_t * const training, double &);
	void Reset (const int);
	e_tcodes RPROPStep (const training_t * const training, double &);
	e_tcodes LMStep (const training_t * const training, double &);
	void InflectWeightSpace (void);

public:

	/*
	 * levels is the number of layers, including the output.  width is an
	 * arrays specifying the width of each level.  e.g., { 4, 1}, is an SLP
	 * with 4 hidden and 1 output perceptron.
	 *
	 */
	NNet_t (int *width, int levels) :
		n_steps (0),
		n_levels (levels),
		n_RPROP (5000),
		n_modes (0),
		n_switch (0.8),
		n_RPROP_steps (0),
		n_LM_steps (0),
		n_mu (MU_INIT),
		n_halt (9e-20),
		n_error (nan (NULL)),
		n_me (-1)
	{
		int max = 0;
		n_Nweights = 0;

		n_Nweights += width[0] * 2; // from input

		for (int i = n_levels - 1; i > 0; --i)
			n_Nweights += width[i] * (width[i - 1] + 1); // + 1 for bias

		n_nn = new perceptron_t * [n_levels];
		n_width = new int [n_levels];
		for (int i = 0; i < n_levels; ++i)
		{
			if (width[i] > max)
				max = width[i];

			n_width[i] = width[i];
			n_nn[i] = new perceptron_t [width[i]];
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
		n_RPROP = RPROP;
	}

	void Compute (double x, double *);
	bool TrainLifeSign (const training_t * const, int, int);
	bool Train (const training_t * const, int);
	bool UpdateTrain (const training_t * const, int);
	void DisplayCurve (const training_t * const);
	void DisplayWeights (void);
	double error (void)
	{
		return n_error;
	}

	ModelSLP_t *buildModel (void);
	void buildModel (ModelSLP_t *);
};

#endif // header inclusion

