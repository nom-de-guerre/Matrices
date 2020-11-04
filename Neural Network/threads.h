#ifndef __THREAD_MANAGER__H__
#define __THREAD_MANAGER__H__

#include <pthread.h>

class lock_t
{
	pthread_mutex_t	vl_mx;

public:

	lock_t ()
	{ 
		pthread_mutex_init (&vl_mx, 0); 
	}

	void lock()     
	{ 
		pthread_mutex_lock (&vl_mx); 
	}

	void unlock()   
	{ 
		pthread_mutex_unlock (&vl_mx); 
	}

	friend class sleep_t;
};

class sleep_t : public lock_t
{
	pthread_cond_t st_cond;

public:

	sleep_t () : lock_t ()
	{
		pthread_cond_init (&st_cond, 0);
	}

	void sleep()
	{
		lock();
		pthread_cond_wait (&st_cond, &vl_mx);
		unlock();
	}

	void wakeup()
	{
		pthread_cond_signal (&st_cond);
	}

	void wakeup_all()
	{
		pthread_cond_broadcast (&st_cond);
	}
};


class thread_t
{
	pthread_t		t_state;
	void			*t_arg;

	static void *run_thread (void *arg)
	{
		thread_t *pThread = (thread_t *) arg;
		pThread->run (pThread->t_arg);

		return 0;
	}

public:

	thread_t (void)
	{
	}

	virtual ~thread_t (void)
	{
		// please by quiet, Mr Clang
	}

	void init (void *arg)
	{
		t_arg = arg;
		pthread_create (&t_state, 0, run_thread, (void *) this);
	}

	int join (void **pResult)
	{
		return (pthread_join (t_state, pResult));
	}

	/*
	 * THIS DOES NOT WORK ON MACOS
	 *
	 */
	int kill (void)
	{
		return pthread_cancel (t_state);
	}

	virtual void run ( void * ) = 0;
};

#endif // header inclusion

