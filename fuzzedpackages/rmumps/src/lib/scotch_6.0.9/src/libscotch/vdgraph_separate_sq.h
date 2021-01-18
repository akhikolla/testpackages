/* Copyright 2007,2018 IPB, Universite de Bordeaux, INRIA & CNRS
**
** This file is part of the Scotch software package for static mapping,
** graph partitioning and sparse matrix ordering.
**
** This software is governed by the CeCILL-C license under French law
** and abiding by the rules of distribution of free software. You can
** use, modify and/or redistribute the software under the terms of the
** CeCILL-C license as circulated by CEA, CNRS and INRIA at the following
** URL: "http://www.cecill.info".
** 
** As a counterpart to the access to the source code and rights to copy,
** modify and redistribute granted by the license, users are provided
** only with a limited warranty and the software's author, the holder of
** the economic rights, and the successive licensors have only limited
** liability.
** 
** In this respect, the user's attention is drawn to the risks associated
** with loading, using, modifying and/or developing or reproducing the
** software by the user in light of its specific status of free software,
** that may mean that it is complicated to manipulate, and that also
** therefore means that it is reserved for developers and experienced
** professionals having in-depth computer knowledge. Users are therefore
** encouraged to load and test the software's suitability as regards
** their requirements in conditions enabling the security of their
** systems and/or data to be ensured and, more generally, to use and
** operate it in the same conditions as regards security.
** 
** The fact that you are presently reading this means that you have had
** knowledge of the CeCILL-C license and that you accept its terms.
*/
/************************************************************/
/**                                                        **/
/**   NAME       : vdgraph_separate_sq.h                   **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : These lines are the data declaration    **/
/**                for the sequential vertex separation    **/
/**                routine for distributed graphs.         **/
/**                                                        **/
/**   DATES      : # Version 5.0  : from : 15 feb 2006     **/
/**                                 to   : 01 mar 2006     **/
/**                # Version 6.0  : from : 06 jun 2018     **/
/**                                 to     06 jun 2018     **/
/**                                                        **/
/************************************************************/

/*
**  The type and structure definitions.
*/

/*+ This structure holds the method parameters. +*/

typedef struct VdgraphSeparateSqParam_ {
  Strat *                   strat;                /*+ Sequential separation strategy used +*/
} VdgraphSeparateSqParam;

/*
**  The function prototypes.
*/

#ifdef VDGRAPH_SEPARATE_SQ
static void                 vdgraphSeparateSqOpBest (const Gnum * const, Gnum * const, const int * const, const MPI_Datatype * const);
#endif /* VDGRAPH_SEPARATE_SQ */

int                         vdgraphSeparateSq   (Vdgraph * const, const VdgraphSeparateSqParam * const);
