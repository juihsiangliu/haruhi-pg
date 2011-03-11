#ifndef PARALLEL_H
#define PARALLEL_H

#include <pthread.h>

inline int block_low(int pid, int numprocs, int numelement);
inline int block_high(int pid, int numprocs, int numelement);
inline int block_size(int pid, int numprocs, int numelement);
inline int block_owner(int index, int numprocs,int numelement);

#endif
