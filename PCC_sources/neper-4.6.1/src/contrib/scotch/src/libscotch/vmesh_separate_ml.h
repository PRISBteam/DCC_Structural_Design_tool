/* Copyright 2004,2007,2018 IPB, Universite de Bordeaux, INRIA & CNRS
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
/**   NAME       : vmesh_separate_ml.h                     **/
/**                                                        **/
/**   AUTHOR     : Francois PELLEGRINI                     **/
/**                                                        **/
/**   FUNCTION   : These lines are the data declarations   **/
/**                for the multi-level node separation     **/
/**                routines.                               **/
/**                                                        **/
/**   DATES      : # Version 4.0  : from : 20 feb 2003     **/
/**                                 to   : 31 aug 2005     **/
/**                # Version 6.0  : from : 31 may 2018     **/
/**                                 to   : 31 may 2018     **/
/**                                                        **/
/************************************************************/

/*
**  The type and structure definitions.
*/

/*+ This structure holds the method parameters. +*/

typedef struct VmeshSeparateMlParam_ {
  INT                       vnodnbr;             /*+ Minimum number of node vertices +*/
  double                    coarrat;             /*+ Coarsening ratio                +*/
  MeshCoarsenType           coartype;            /*+ Element matching function type  +*/
  Strat *                   stratlow;            /*+ Strategy at lowest level        +*/
  Strat *                   stratasc;            /*+ Strategy at ascending levels    +*/
} VmeshSeparateMlParam;

/*
**  The function prototypes.
*/

#ifdef VMESH_SEPARATE_ML
static int                  vmeshSeparateMlCoarsen (const Vmesh * restrict const, Vmesh * restrict const, Gnum * restrict * const, const VmeshSeparateMlParam * const);
static int                  vmeshSeparateMlUncoarsen (Vmesh * const, const Vmesh * const, const Gnum * restrict const);
static int                  vmeshSeparateMl2    (Vmesh * const, const VmeshSeparateMlParam * const);
#endif /* VMESH_SEPARATE_ML */

int                         vmeshSeparateMl     (Vmesh * const, const VmeshSeparateMlParam * const);
