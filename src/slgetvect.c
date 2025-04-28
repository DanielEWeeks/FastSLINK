/* This file is part of a C version of SLINK a simulation program based
on the LINKAGE package. SLINK was originally described in: D. E.
Weeks, J. Ott, G. M. Lathrop (1990), SLINK: a general simulation
program for linkage analysis, American Journal of Human Genetics
47(1990), A204 (abstract). This version of SLINK is adapted to use the
faster pedigree traversal algorithms described in: R. W. Cottingham
Jr., R. M. Idury, A. A. Schaffer (1993), Faster sequential genetic
linkage computations, American Journal of Human Genetics 53(1993), pp.
252-263. */
/* This file is a piece of SLINK. It contains some low level
   probability routines
   8 December 1993 Sliced out of SLINK by A. A. Schaffer
  14 March 1994 A.A. Schaffer fixed non-ANSI problem with prototypes*/

#include "commondefs.h"
#include "sldefs.h"

/* Local variables for getvect: */
struct LOC_getvect {
  struct LOC_likelihood *LINK;
  thisperson *p;
  hapvector hap1, hap2;
} ;

/*static void getgene PP((int syste, double val, struct LOC_getvect *LINK));
static void ugetgene PP((int syste, double val, struct LOC_getvect *LINK)); */

/* void getgene(int syste, double val, struct LOC_getvect *LINK);
void ugetgene(int syste, double val, struct LOC_getvect *LINK); */

void getgene();
void ugetgene();

static double quanfun(phen, thislocus, i, j, mean, LINK)
phenotype *phen;
locusvalues *thislocus;
int i, j;
double *mean;
struct LOC_getvect *LINK;
{
  double val;
  int it, jt;   /*quanfun*/
  int FORLIM, FORLIM1;

  val = 1.0;
  if (phen->missing)
    return val;
  val = 0.0;
  FORLIM = thislocus->UU.U1.ntrait;
  for (it = 0; it < FORLIM; it++) {
    FORLIM1 = thislocus->UU.U1.ntrait;
    for (jt = 0; jt < FORLIM1; jt++) {
      if (i == j)
	val += thislocus->UU.U1.vmat[it]
	       [jt] * (phen->x[jt] - mean[jt]) *
	       (phen->x[it] - mean[it]);
      else
	val += thislocus->UU.U1.conmat * thislocus->UU.U1.vmat[it]
	       [jt] * (phen->x[jt] - mean[jt]) *
	       (phen->x[it] - mean[it]);
    }
  }
  val = thislocus->UU.U1.det * exp(-val * 0.5);
  if (i != j)
    val *= thislocus->UU.U1.contrait;
  return val;
}  /* quanfun */

/*quanfun*/

static void getval(syste, i, j, val, LINK)
int syste, i, j;
double *val;
struct LOC_getvect *LINK;
{
  /*getval*/
  locusvalues *WITH;
  phenotype *WITH1;

  WITH = thislocus[syste - 1];
  switch (WITH->which) {

  case quantitative:
    *val *= quanfun(LINK->p->phen[syste - 1], thislocus[syste - 1], i, j,
		    WITH->UU.U1.pm[i][j - 1], LINK);
    break;

  case affection:
    WITH1 = LINK->p->phen[syste - 1];
    *val *= WITH->UU.U0.pen[i][j - 1][WITH1->aff]
      [WITH1->liability - 1];
    break;
  }
}  /* getval */

#include "commongetvect.c"

