#include <threads.h>
#include <NNLM.h>
#include <NN_misc.h>

static lock_t count_guard;
static volatile int N_working = 0;

static void StartWorking (void)
{
	count_guard.lock ();

	++N_working;

	count_guard.unlock ();
}

static void StopWorking (sleep_t *cond)
{
	count_guard.lock ();

	--N_working;

	if (N_working == 0)
		cond->wakeup ();

	count_guard.unlock ();
}

struct NN_args_t
{
	training_t		*data;
	sleep_t			*cond;
	NNet_t			**home;
	int				*layers;
	int				me;
	double			mse;

	void init (training_t *_data, 
				sleep_t *_cond, 
				NNet_t **_home, 
				int *_layers, 
				int _me,
				double _mse)
	{
		data = _data;
		cond = _cond;
		home = _home;
		layers = _layers;
		me = _me;
		mse = _mse;
	}
};

struct Trainer_t : public thread_t
{
	Trainer_t (NN_args_t *p) : thread_t ()
	{
		init (p);
	}

#if 0
	static void cleanup (void *netp)
	{
		delete (NN_t *) netp
	}
#endif

	void run (void *p)
	{
		NN_args_t *args = (NN_args_t *) p;
		NNet_t *Np;
		bool success = true;

		Np = new NNet_t (args->layers + 1, args->layers[0]);
		Np->setMSE (args->mse);

		*args->home = Np;

		StartWorking ();

#if 0
MacOS pthread_cancel does not work :-(

		pthread_setcanceltype (PTHREAD_CANCEL_ASYNCHRONOUS, NULL);
		pthread_cleanup_push (cleanup, Np);
#endif

		try {

			Np->TrainLifeSign (args->data, args->me, 50000);

		} catch (char const *error) {

			success = false;
			printf ("(%d) FAILED: %s\n", args->me, error);

		}

		StopWorking (args->cond);

		if (success)
			args->cond->wakeup ();

#if 0
		pthread_cleanup_pop (0);
#endif

		pthread_exit (NULL);
	}
};

#define N_THREADS	1
bool t_keep_working;

NNet_t *TrainNN (training_t *tp, int *layers, bool &found)
{
	NN_args_t Comms[N_THREADS];
	NNet_t *Results[N_THREADS];
	Trainer_t *threads[N_THREADS];
	sleep_t X;
	NNet_t *soln = NULL;;

	found = false;

	t_keep_working = true;

	for (int i = 0; i < N_THREADS; ++i)
	{
		Results[i] = NULL;
		Comms[i].init (tp, &X, Results + i, layers, i, 9e-19);
		threads[i] = new Trainer_t (Comms + i);
	}

	X.sleep ();
	t_keep_working = false; // MacOS pthread_cancel does not work

	soln = Results[0];

	for (int i = 0; i < N_THREADS; ++i)
	{
		threads[i]->join (NULL);
		delete threads[i];

		if (i == 0)
			continue;

		if (Results[i]->error () < soln->error ()) {

			delete soln;
			soln = Results[i];

		} else
			delete Results[i];
	}

	t_keep_working = true;

	printf ("Accepting %e\n", soln->error ());

	if (soln->error () <= soln->whatIsHalt ())
		found = true;

	return soln;
}

