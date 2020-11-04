#include <NNLM.h>

void perceptron_t::bprop (
	perceptron_t *input,
	perceptron_t *output,
	const int Noutputs,
	const int me,
	const int row,
	const int start,
	Md_t &M)
{
	double total_dE = 0;
	double d_act = DERIVATIVE_FN (p_output);

	/*
	 * 1. Compute the total derivative for ∆j
	 *
	 * Note: me + 1, account for bias
	 *
	 */
	for (int i = 0; i < Noutputs; ++i)
		total_dE += output[i].p_delta * output[i].p_weights[me + 1];

	p_delta = total_dE * d_act;

	p_error[0] += p_delta; // bias is virtual
	M (row, start) = p_delta;

	// 2. Update weight errors
	for (int i = 1; i < p_N; ++i)
	{
		M (row, start + i) = p_delta * input[i - 1].signal ();
		p_error[i] += p_delta * input[i - 1].signal ();
	}
}

void perceptron_t::rprop (int index)
{
	double delta;
	double backtrack;

	if (p_Ei[index] == 0.0 || p_error[index] == 0.0)
	{
		delta = -SIGN (p_error[index]) * p_deltaW[index];
		p_weights[index] += delta;

		p_Ei[index] = p_error[index];

	} else if (signbit (p_error[index]) == signbit (p_Ei[index])) {

		// (1)
		delta = p_deltaW[index] * ETA_PLUS;
		if (delta > DELTA_MAX)
			delta = DELTA_MAX;

		p_deltaW[index] = delta;

		// (2)
		delta *= -(SIGN (p_error[index]));
		// (3)
		p_weights[index] += delta;

		p_Ei[index] = p_error[index];

	} else {

		backtrack = p_deltaW[index] * SIGN (p_Ei[index]);

		// (1)
		delta = p_deltaW[index] * ETA_MINUS;
		if (delta < DELTA_MIN)
			delta = DELTA_MIN;

		p_deltaW[index] = delta;

		// (2)
		p_weights[index] += backtrack;

		// (3)
		p_Ei[index] = 0.0;
	}

	p_error[index] = 0.0;

	assert (p_deltaW[index] > 0);
}

template<typename T> void 
NNet_t<T>::UpdateWeights (Md_t &dw)
{
	int index = 0;

	perceptron_t *subject;

	for (int level = n_levels - 1; level >= 0; --level)
	{
		subject = n_nn[level];

		for (int i = 0; i < n_width[level]; ++i)
			for (int wi = 0; wi < subject[i].p_N; ++wi, ++index)
			{
				subject[i].p_error[wi] = 0.0;

				if (isnan (dw(index, 0)))
					throw ("Update degradation");
				subject[i].p_weights[wi] -= dw(index, 0);
			}
	}
}

template<typename T> void
NNet_t<T>::PresentExamples (
	const DataSet_t * const training, 
	bool ComputeDerivatives)
{
	for (int i = 0; i < training->t_N; ++i)
	{
		const TrainingRow_t p = (*training)[i];

		if (ComputeDerivatives)
			n_Results(i, 0) = ComputeDerivative (p, i);
		else
			n_Results(i, 0) = Compute (p);
	}
}

template<typename T> double
NNet_t<T>::PresentExamplesLoss (
	const DataSet_t * const training, 
	bool ComputeDerivatives)
{
	PresentExamples (training, ComputeDerivatives);
	double error;

	error = Loss (training);

	return error;
}

int cycles = 0;

template<typename T> bool 
NNet_t<T>::Step (const DataSet_t * const training, double &progress)
{
	double error;

	Start ();

	if (n_RPROP_pending > 0)
	{
		RPROPStep (training, error);
		/*
		 * We begin by letting RPROP+ shape the initial curve
		 * to the correct one.  Once the error is down to the 
		 * best RPROP+ can do we switch to using Levenberg-Marquardt 
		 * optimization.  After that n_RPROP_pending is used to push
		 * the solution if LM gets stuck. (LM degrades to gradient
		 * descent outside of the trust area).
		 *
		 */

		--n_RPROP_pending;

	} else
		LMStep (training, error);

	progress = n_error - error;

	n_error = error;

	if (Halt (training))
		return false;

	if (n_steps && (n_steps % 10000) == 0)
	{
		printf ("%s\t%e\t%e\t%d\t%d\t%d\t%d\t%e\n",
			(n_RPROP_pending ? "RPROP+" : "LM"),
			error,
			progress,
			n_steps,
			n_RPROP_pending_steps,
			n_steps - n_RPROP_pending_steps,
			n_LM_steps,
			n_mu);
	}

	if (fabs (progress) < 1e-300)
		Reset (training->t_N);

	return true;
}

template<typename T> void 
NNet_t<T>::RPROPStep (const DataSet_t * const training, double &sumsq)
{
	++n_RPROP_pending_steps;

	sumsq = PresentExamplesLoss (training, true);

	for (int i = 0; i < n_levels; ++i)
		for (int j = 0; j < n_width[i]; ++j)
			n_nn [i][j].rprop ();
}

template<typename T> void 
NNet_t<T>::LMStep (const DataSet_t * const training, double &sumsq)
{
	bool progress = false;
	double attempt;
	int N_attempts = 0;
	Md_t z;

	PresentExamples (training, true);
	sumsq = Loss (training); // will convert n_Results to g

	Md_t grad = transpose (n_J) * n_Results;
	Md_t H = transpose (n_J) * n_J;
	Md_t dW (H.rows (), 1);
	Md_t dW_s (H.rows (), 1);
	bool solved;

	while (!progress)
	{
		Md_t grad_ = grad; // it will be destroyed in QR
		Md_t H_ = H;

		++n_LM_steps;

		for (int i = 0; i < H.rows (); ++i)
			H_(i, i) += n_mu;

		/*
		 * Attempt a Cholesky factorization first.  If that fails then
		 * we resort to full QR (rare).
		 *
		 * Should we just back off lambda instead?
		 *
		 */
		solved = H_.SolveSymmetric (grad, dW);
		if (!solved)
			dW = H_.solveQR (grad_);

		UpdateWeights (dW);

		attempt = PresentExamplesLoss (training, false);
		if (!isnan (attempt) && attempt < sumsq)
		{
			n_mu /= MU_THETA;
			progress = true;

			continue;
		}

		dW *= -1;
		UpdateWeights (dW);
		n_mu *= MU_THETA;

		++N_attempts;

		if (N_attempts > 10)
		{
			n_RPROP_pending = 200;
			n_mu = MU_INIT;
			progress = true;
		}
	}

	sumsq = attempt;
}

template<typename T> void 
NNet_t<T>::Reset (const int N)
{
	for (int i = 0; i < n_levels; ++i)
		for (int j = 0; j < n_width[i]; ++j)
			n_nn[i][j].init ((i == 0 ? n_Nin + 1 : n_width[i - 1] + 1));

#ifdef __RELU
	n_nn[n_width[n_levels - 1], 0].setRELU ();
#endif

}

template<typename T> void 
NNet_t<T>::Start (void)
{
	return static_cast<T *> (this)->Cycle ();
}

template<typename T> bool
NNet_t<T>::Halt (DataSet_t const * const tp)
{
	return static_cast<T *> (this)->Test (tp);
}

template<typename T> double
NNet_t<T>::Loss (DataSet_t const *tp)
{
	return static_cast<T *> (this)->error (tp);
}

/*
 * ∂E     ∂E  ∂ak  ∂net
 * --   = --  ---  ----
 * ∂wkj   ∂ak ∂net ∂wkj
 *
 * y = ak
 *
 * ∆k = (ak - t)ak(1 - ak)
 * wkj = ∆k(aj)
 *
 */

template<typename T> double 
NNet_t<T>::ComputeDerivative (const TrainingRow_t x, const int row)
{
	int dEj;
	/*
	 * Initiate the recurrence by triggering the network
	 * specialization.
	 *
	 */
	double error = static_cast<T *>(this)->bprop (x, row, dEj);

	int ncols = n_J.columns ();
	perceptron_t input;
	perceptron_t *subject, *outputp, *inputp;
	input.p_output = x[0];

	for (int level = n_levels - 2; level >= 0; --level)
	{
		subject = n_nn[level];
		outputp = n_nn[level + 1];
		inputp = (level > 0 ? n_nn[level - 1] : &input);
		for (int i = 0; i < n_width[level]; ++i)
		{
			subject[i].bprop (
				inputp,
				outputp, 
				n_width[level + 1],
				i,
				row,
				dEj,
				n_J);

			dEj += subject[i].p_N;
		}

		assert (dEj <= ncols);
	}

	return error;
}

template<typename T> double
NNet_t<T>::Compute (double *x)
{
	double *input;
	double *output;

	assert (n_buffers[0][0] == 1);
	assert (n_buffers[1][0] == 1);

	for (int i = 0; i < n_Nin; ++i)
		n_buffers[1][i + 1] = x[i]; // input to the first layer

	for (int layers = 0; layers < n_levels; ++layers)
	{
		if (layers % 2) {

			input = n_buffers[0];
			output = n_buffers[1];

		} else {

			input = n_buffers[1];
			output = n_buffers[0];
		}

		perceptron_t *p = n_nn[layers];
		
		for (int i = 0; i < n_width[layers]; ++i)
			output[i + 1] = p[i].signal (input);
	}

	return static_cast<T *>(this)->f (output);
}

template<typename T> bool 
NNet_t<T>::TrainLifeSign (
	const DataSet_t * const training, 
	int me, 
	int maxIterations)
{
	n_me = me;

	return Train (training, maxIterations);
}

template<typename T> bool 
NNet_t<T>::Train (const DataSet_t * const training, int maxIterations)
{
	Reset (training->t_N);
	return UpdateTrain (training, maxIterations);
}

template<typename T> bool 
NNet_t<T>::UpdateTrain (const DataSet_t * const training, int maxIterations)
{
	bool solving = true;
	double progress = 1;

	n_J = Md_t (training->t_N, n_Nweights, 0.0);
	n_J.set_WiP ();
	n_Results = Md_t (training->t_N, 1);

	n_error = 1.0;

	for (n_steps = 0; 
//		(n_steps < maxIterations) && solving; 
		(n_steps < maxIterations || progress > 0) && solving; 
		++n_steps)
	{
		try {

			solving = Step (training, progress);

		} catch (const char *error) {

			printf ("%d\t%s\tstill %d steps to try.\n",
				n_me,
				error,
				maxIterations - n_steps);

			Reset (training->t_N);
		}
	}

	if (n_steps >= maxIterations)
		throw ("Exceeded Iterations");

	printf ("FINISHED TRAINING (%d)\t%d\t%e\tRPROP+ %d\tLM\t%d\t%d\n", 
		n_me,
		n_steps, 
		n_error, 
		n_RPROP_pending_steps, 
		n_steps - n_RPROP_pending_steps,
		n_LM_steps);
	if (n_RPROP_pending)
		printf ("SOLVED BY RPROP\n");

	return true;
}

template<typename T> void 
NNet_t<T>::DisplayWeights (void)
{
	for (int i = 0; i < n_levels; ++i)
		for (int j = 0; j < n_width[i]; ++j)
		{
			for (int wi = 0; wi < n_nn[i][j].p_N; ++wi)
				printf ("%f\t", n_nn[i][j].p_weights[wi]);

			printf ("\n");
		}
}

