#include "parallel.h"

inline int block_low(int pid, int numprocs, int numelement)
{
	return pid*numelement / numprocs;
}

inline int block_high(int pid, int numprocs, int numelement)
{
	return block_low(pid+1,numprocs,numelement) -1;
}

inline int block_size(int pid, int numprocs, int numelement)
{
	return block_low(pid+1,numprocs,numelement) - block_low(pid,numprocs,numelement);
}

inline int block_owner(int index, int numprocs,int numelement)
{
	return (numprocs * (index+1) -1) / numelement;
}
