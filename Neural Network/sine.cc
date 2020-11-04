#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>

#include <regression.h>

#define N_POINTS	32

#define PI 			3.141592653589793
#define PI_2		1.570796326794897
#define PI_DELTA	(PI_2 / N_POINTS)

DataSet_t *BuildTrainingSet (int);
void Run (int *);

int main (int argc, char *argv[])
{
	if (argc < 2)
	{
		printf ("Usage: LoW hidden-layers output-layers\n");
		exit (-1);
	}

	long seed = time (0);
	printf ("Seed %ld\n", seed);

	srand (seed);

	int N_layers = argc - 1;
	int *layers = new int [N_layers + 2];	// widths plus length prefix, inputs
	layers[0] = N_layers + 1;
	layers[1] = 1;							// one input
	for (int i = 0; i < N_layers; ++i)
		layers[i + 2] = atoi (argv[i + 1]);

	Run (layers);

	delete [] layers;
}

void Run (int *layers)
{
	double soln_MSE = 5e-7;

	DataSet_t *O = BuildTrainingSet (N_POINTS);
	Regression_t *Np = NULL;
	double guess;

	Np = new Regression_t (layers + 1, layers[0]);
	Np->setMSE (soln_MSE);
	Np->setRPROP (200);

	try {

		Np->TrainLifeSign (O, 1, 500000);

	} catch (const char *excep) {

		printf ("ERROR: %s\n", excep);
		exit (-1);
	}

	bool accept_soln = true;
	double error;

	for (int i = 0; i < N_POINTS; ++i)
	{
		guess = Np->Compute ((*O)[i]);

		error = (*O)[i][1] - guess;
		error *= error;

		if (error > soln_MSE)
			accept_soln = false;

		printf ("DJS_RESULT\t%1.8f\t%1.8f\t%1.8f\n",
			(*O)[i][0],
			(*O)[i][1],
			guess);
	}

#if 1
	O = BuildTrainingSet (64);
	for (int i = 0; i < O->t_N; ++i)
	{
		guess = Np->Compute ((*O)[i]);
		printf ("DJS_INFER\t%1.8f\t%1.8f\t%1.8f\n",
			(*O)[i][0],
			(*O)[i][1],
			guess);
	}
#endif

	if (accept_soln)
		printf (" *** Solution ACCEPTED.\n");
	else
		printf (" *** Solution REJECTED.\n");

}

DataSet_t *BuildTrainingSet (int N)
{
	DataSet_t *O = new DataSet_t (N, 1, 1);

	for (int i = 0; i < N; ++i)
	{
		double sample = (double) rand () / RAND_MAX;

		(*O)[i][0] = sample * PI_2;
		(*O)[i][1] = sin ((*O)[i][0]);
	}

	return O;
}

