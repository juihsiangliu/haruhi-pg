#ifndef HARUHI_PG_H
#define HARUHI_PG_H


#include <string.h>
#include <stdio.h>

#include "partition_double.h"
#include "mempool.h"
#include "sparsedoublematrix.h"
#include "solvesparsedoublematrix.h"
#include "parallel_lu_double.h"
#include "parallel_pcg_double.h"
#include "parse_spice.h"



struct InputPar
{
	enum OOCFlag oocFlag; // -ooc_flag [ooc|ic]
	enum SolveMethod solveMethod; // -solve_method [direct|iter]
	enum OrderMethod orderMethod; // -order_method [amd|metis]
	int threadNum; // -thead [1|2|4|8]
};


typedef struct InputPar InputPar;



#endif
