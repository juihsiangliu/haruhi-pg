/*
 * Copyright 2003, Regents of the University of Minnesota
 *
 * cmetis.c
 *
 * This file contains the driving routine for partitioning for 
 * sub-domain direct factorization.
 *
 * Started 3/25/03
 * George
 *
 * $Id: cmetis.c,v 1.5 2003/04/08 12:52:59 karypis Exp $
 *
 */

#include <metisbin.h>



/*************************************************************************
* Let the game begin
**************************************************************************/
int main(int argc, char *argv[])
{
  idxtype i, options[10], nnodes, nclean, naclean, ndirty, maxdepth;
  idxtype *part, *sflag;
  float lbvec[MAXNCON];
  GraphType graph;
  idxtype numflag = 0, wgtflag = 0, edgecut;
  ParamType params;
  void *cinfo;
  double TOTALTmr, METISTmr, IOTmr;

  parse_cmdline(&params, argc, argv);

  if (params.nparts < 2) {
    mprintf("The number of partitions should be greater than 1!\n");
    exit(0);
  }

  gk_clearcputimer(TOTALTmr);
  gk_clearcputimer(METISTmr);
  gk_clearcputimer(IOTmr);

  gk_startcputimer(TOTALTmr);
  gk_startcputimer(IOTmr);
  ReadGraph(&graph, params.filename, &wgtflag);
  if (graph.nvtxs <= 0) {
    mprintf("Empty graph. Nothing to do.\n");
    exit(0);
  }
  ReadCoordinates(&graph, params.xyzfilename);
  gk_stopcputimer(IOTmr);

  mprintf("**********************************************************************\n");
  mprintf("%s", METISTITLE);
  mprintf("Graph Information ---------------------------------------------------\n");
  mprintf("  Name: %s, #Vertices: %D, #Edges: %D, #Parts: %D\n", params.filename, graph.nvtxs, graph.nedges/2, params.nparts);
  if (graph.ncon > 1)
    mprintf("  Balancing Constraints: %D\n", graph.ncon);
  mprintf("\nRecursive Partitioning... -------------------------------------------\n");

  part  = idxmalloc(graph.nvtxs, "main: part");
  sflag = idxsmalloc(graph.nvtxs, 1, "main: sflag");
  options[0] = 0;

  gk_startcputimer(METISTmr);
  cinfo = METIS_PartGraphForContact(&graph.nvtxs, graph.xadj, graph.adjncy, graph.coords, sflag, 
                    &numflag, &(params.nparts), options, &edgecut, part);

  METIS_UpdateContactInfo(cinfo, &graph.nvtxs, graph.coords, sflag);
  gk_stopcputimer(METISTmr);

  ComputePartitionBalance(&graph, params.nparts, part, lbvec);

  graph.vwgt = idxsmalloc(graph.nvtxs, 1, "main: graph->vwgt");

  mprintf("  %D-way Edge-Cut: %7D, Volume: %7D, Balance: ", params.nparts, edgecut, 
                                                           ComputeVolume(&graph, part));
  for (i=0; i<graph.ncon; i++)
    mprintf("%5.2f ", lbvec[i]);
  mprintf("\n");

  gk_startcputimer(IOTmr);
  WritePartition(params.filename, part, graph.nvtxs, params.nparts); 
  gk_stopcputimer(IOTmr);
  gk_stopcputimer(TOTALTmr);

  mprintf("\nTiming Information --------------------------------------------------\n");
  mprintf("  I/O:          \t\t %7.3f\n", gk_getcputimer(IOTmr));
  mprintf("  Partitioning: \t\t %7.3f   (CMETIS time)\n", gk_getcputimer(METISTmr));
  mprintf("  Total:        \t\t %7.3f\n", gk_getcputimer(TOTALTmr));
  mprintf("**********************************************************************\n");

  gk_free((void **)&graph.xadj, &graph.adjncy, &graph.vwgt, &graph.adjwgt, &graph.coords, &part, &sflag, LTERM);
}  

