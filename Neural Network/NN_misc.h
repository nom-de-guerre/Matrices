#ifndef __NEURAL_SUPPORT__H__
#define __NEURAL_SUPPORT__H__

#ifndef __NN_MSE // mean SQUARED error
#define __NN_MSE 9e-12
#endif

NNet_t *TrainNN (training_t *, int *, bool &);

#endif // header inclusion

