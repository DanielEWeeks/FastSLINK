/* This file is part of a C version of SLINK a simulation program based
on the LINKAGE package. SLINK was originally described in: D. E.
Weeks, J. Ott, G. M. Lathrop (1990), SLINK: a general simulation
program for linkage analysis, American Journal of Human Genetics
47(1990), A204 (abstract). This version of SLINK is adapted to use the
faster pedigree traversal algorithms described in: R. W. Cottingham
Jr., R. M. Idury, A. A. Schaffer (1993), Faster sequential genetic
linkage computations, American Journal of Human Genetics 53(1993), pp.
252-263.*/
/*This file contains many of the I/O routines*/
/* December 1993 Sliced out of slink.c by A. A. Schaffer*/
/* May 2008 Edited to allow > 32 alleles at a numbered allele locus*/

#include <stdio.h>
#include "commondefs.h"
#include "sldefs.h"

extern int P_eof(FILE *);

/* Local variables for recombination: */
struct LOC_recombination {
  int here;
  double p1, p2, p3, p4;
} ;

/* Local variables for recombine: */
struct LOC_recombine {
  struct LOC_recombination *LINK;
  double *theta;
  double *segprob;
  int there;
  int nhap;
  boolean thishet[maxlocus];
  hapvector hap1, hap2;
  hapvector thishap1[maxseg], thishap2[maxseg];
} ;

static void scramble(LINK)
struct LOC_recombine *LINK;
{
  int whichhap, start, length, i, j, k, FORLIM2;
  double recval, val;

  start = 0;
  do {
    start++;
  } while (!LINK->thishet[start - 1]);
  length = LINK->there - LINK->LINK->here;
  for (i = 1; i < length; i++) {
    memcpy(LINK->hap1, LINK->thishap1[i], sizeof(hapvector));
    for (j = 1; j <= length; j++) {
      val = 0.5;
      whichhap = 1;
      recval = LINK->theta[start - 1];
      FORLIM2 = nlocus;
      for (k = start; k < FORLIM2; k++) {
	if (!LINK->thishet[k])
	  recval = recval * (1.0 - LINK->theta[k]) +
		   LINK->theta[k] * (1.0 - recval);
	else {
	  if (whichhap == 1) {
	    if (LINK->thishap1[j - 1][k] == LINK->hap1[k])
	      val *= 1 - recval;
	    else {
	      val *= recval;
	      whichhap = 2;
	    }
	  } else {
	    if (LINK->thishap2[j - 1][k] == LINK->hap1[k])
	      val *= 1 - recval;
	    else {
	      val *= recval;
	      whichhap = 1;
	    }
	  }
	  recval = LINK->theta[k];
	}
      }
      LINK->there++;
      LINK->segprob[LINK->there - 1] = val;
    }
  }
}  /* scramble */


static void setrec(val, LINK)
double val;
struct LOC_recombine *LINK;
{
  LINK->nhap++;
  memcpy(LINK->thishap1[LINK->nhap - 1], LINK->hap1, sizeof(hapvector));
  memcpy(LINK->thishap2[LINK->nhap - 1], LINK->hap2, sizeof(hapvector));
  LINK->there++;
  LINK->segprob[LINK->there - 1] = val;
}  /* setrec */


static void dointer(LINK)
struct LOC_recombine *LINK;
{
  int i;
  boolean temphet[3];
  int FORLIM;

  FORLIM = nlocus;
  for (i = 0; i < FORLIM; i++)
    temphet[i] = LINK->thishet[i];

  if (temphet[0] && temphet[1] && !temphet[2]) {
    setrec(0.5 - 0.5 * LINK->theta[0], LINK);
    setrec(0.5 * LINK->theta[0], LINK);
    setrec(0.5 * LINK->theta[0], LINK);
    setrec(0.5 - 0.5 * LINK->theta[0], LINK);
    return;
  }
  if (temphet[2] && temphet[1] && !temphet[0]) {
    setrec(0.5 - 0.5 * LINK->theta[nlocus - 2], LINK);
    setrec(0.5 * LINK->theta[nlocus - 2], LINK);
    setrec(0.5 * LINK->theta[nlocus - 2], LINK);
    setrec(0.5 - 0.5 * LINK->theta[nlocus - 2], LINK);
    return;
  }
  if (temphet[2] && temphet[0] && !temphet[1]) {
    setrec(0.5 - 0.5 * LINK->theta[nlocus - 1], LINK);
    setrec(0.5 * LINK->theta[nlocus - 1], LINK);
    setrec(0.5 * LINK->theta[nlocus - 1], LINK);
    setrec(0.5 - 0.5 * LINK->theta[nlocus - 1], LINK);
    return;
  }
  if (!(temphet[0] && temphet[1] && temphet[2])) {
    setrec(0.5, LINK);
    return;
  }
  setrec(0.5 * LINK->LINK->p4, LINK);
  setrec(0.5 * LINK->LINK->p2, LINK);
  setrec(0.5 * LINK->LINK->p3, LINK);
  setrec(0.5 * LINK->LINK->p1, LINK);
  setrec(0.5 * LINK->LINK->p2, LINK);
  setrec(0.5 * LINK->LINK->p4, LINK);
  setrec(0.5 * LINK->LINK->p1, LINK);
  setrec(0.5 * LINK->LINK->p3, LINK);
  setrec(0.5 * LINK->LINK->p3, LINK);
  setrec(0.5 * LINK->LINK->p1, LINK);
  setrec(0.5 * LINK->LINK->p4, LINK);
  setrec(0.5 * LINK->LINK->p2, LINK);
  setrec(0.5 * LINK->LINK->p1, LINK);
  setrec(0.5 * LINK->LINK->p3, LINK);
  setrec(0.5 * LINK->LINK->p2, LINK);
  setrec(0.5 * LINK->LINK->p4, LINK);

  /*not informative*/
}  /* dointer */


static void nexthet(i, val, inphase, LINK)
int i;
double val;
boolean inphase;
struct LOC_recombine *LINK;
{
  double newval, recval;

  recval = LINK->theta[i - 1];
  do {
    i++;
    LINK->hap1[i - 1] = 0;
    LINK->hap2[i - 1] = 0;
    if (!LINK->thishet[i - 1] && i != nlocus)
      recval = recval * (1 - LINK->theta[i - 1]) +
	       (1 - recval) * LINK->theta[i - 1];
  } while (!(i == nlocus || LINK->thishet[i - 1]));
  if (i != nlocus) {
    if (inphase)
      newval = val * (1 - recval);
    else
      newval = val * recval;
    LINK->hap1[i - 1] = 1;
    LINK->hap2[i - 1] = 2;
    nexthet(i, newval, true, LINK);
    LINK->hap2[i - 1] = 1;
    LINK->hap1[i - 1] = 2;
    if (!inphase)
      newval = val * (1 - recval);
    else
      newval = val * recval;
    nexthet(i, newval, false, LINK);
    return;
  }
  if (!LINK->thishet[i - 1]) {
    setrec(val, LINK);
    return;
  }
  if (inphase)
    newval = val * (1 - recval);
  else
    newval = val * recval;
  LINK->hap1[i - 1] = 1;
  LINK->hap2[i - 1] = 2;
  setrec(newval, LINK);
  if (!inphase)
    newval = val * (1 - recval);
  else
    newval = val * recval;
  LINK->hap2[i - 1] = 1;
  LINK->hap1[i - 1] = 2;
  setrec(newval, LINK);
}  /* nexthet */


static void getrecprob(LINK)
struct LOC_recombine *LINK;
{
  int i;

  LINK->nhap = 0;
  LINK->there = LINK->LINK->here;
  i = 0;
  do {
    i++;
    if (LINK->thishet[i - 1]) {
      LINK->hap1[i - 1] = 1;
      LINK->hap2[i - 1] = 2;
    } else {
      LINK->hap1[i - 1] = 0;
      LINK->hap2[i - 1] = 0;
    }
  } while (!(LINK->thishet[i - 1] || i == nlocus));
  if (i == nlocus)
    setrec(0.5, LINK);
  else if (interfer)
    dointer(LINK);
  else
    nexthet(i, 0.5, true, LINK);
  if (LINK->nhap > 1 && !interfer)
    scramble(LINK);
  LINK->LINK->here = LINK->there;
}  /* getrecprob */


static void gethet(system, LINK)
int *system;
struct LOC_recombine *LINK;
{
  int newsystem;

  newsystem = *system + 1;
  LINK->thishet[*system - 1] = false;
  if (*system != nlocus)
    gethet(&newsystem, LINK);
  else
    getrecprob(LINK);
  LINK->thishet[*system - 1] = true;
  if (*system != nlocus)
    gethet(&newsystem, LINK);
  else
    getrecprob(LINK);
}  /* gethet */

static void recombine(theta_, segprob_, LINK)
double *theta_;
double *segprob_;
struct LOC_recombination *LINK;
{
  struct LOC_recombine V;
  int system;


  V.LINK = LINK;
  V.theta = theta_;
  V.segprob = segprob_;
  LINK->here = 0;
  system = 1;
  gethet(&system, &V);
}  /* recombine */


static void getfemaletheta(LINK)
struct LOC_recombination *LINK;
{
  double dist;
  int ntheta, i;

  if (interfer)
    ntheta = nlocus;
  else
    ntheta = nlocus - 1;
  for (i = 0; i < ntheta; i++) {
    dist = getdist(&maletheta->theta[i]) * distratio;
    femaletheta->theta[i] = invdist(&dist);
  }
}  /* getfemaletheta */


void recombination()
{
  struct LOC_recombination V;
  int i;
  thetarray oldtheta;
  thetavalues *WITH;
  int FORLIM;


  if (interfer) {
    WITH = maletheta;
    if (mapping)
      WITH->theta[nlocus - 1] = mapfunction(WITH->theta[nlocus - 2],
					    WITH->theta[0]);
    memcpy(oldtheta, WITH->theta, sizeof(thetarray));
    if (!mapping && !dolod) {
      FORLIM = nlocus;
      for (i = 0; i < FORLIM; i++)
	oldtheta[i] = 1 / (1 + exp(oldtheta[i]));
      WITH->theta[0] = oldtheta[0] + oldtheta[nlocus - 1];
      WITH->theta[nlocus - 2] = oldtheta[nlocus - 2] + oldtheta[nlocus - 1];
      WITH->theta[nlocus - 1] = oldtheta[0] + oldtheta[nlocus - 2];
      V.p1 = oldtheta[0];
      V.p2 = oldtheta[nlocus - 2];
      V.p3 = oldtheta[nlocus - 1];
      V.p4 = 1.0 - V.p1 - V.p2 - V.p3;
    } else {
      V.p1 = (WITH->theta[0] +
	      WITH->theta[nlocus - 1] - WITH->theta[nlocus - 2]) / 2.0;
      V.p2 = (WITH->theta[nlocus - 2] +
	      WITH->theta[nlocus - 1] - WITH->theta[0]) / 2.0;
      V.p3 = (WITH->theta[nlocus - 2] +
	      WITH->theta[0] - WITH->theta[nlocus - 1]) / 2.0;
      V.p4 = 1.0 - V.p1 - V.p2 - V.p3;
    }
    recombine(WITH->theta, WITH->segprob, &V);
  } else
    recombine(maletheta->theta, maletheta->segprob, &V);
  if (sexdif) {
    if (!readfemale) {
      if (interfer && !dolod) {
	WITH = maletheta;
	WITH->theta[0] = oldtheta[0] + oldtheta[nlocus - 1];
	WITH->theta[nlocus - 2] = oldtheta[nlocus - 2] + oldtheta[nlocus - 1];
	WITH->theta[nlocus - 1] = oldtheta[0] + oldtheta[nlocus - 2];
      }
      getfemaletheta(&V);
    }
    if (interfer) {
      WITH = femaletheta;
      if (mapping)
	WITH->theta[nlocus - 1] = mapfunction(WITH->theta[nlocus - 2],
					      WITH->theta[0]);
      if (readfemale && !mapping && !dolod) {
	memcpy(oldtheta, WITH->theta, sizeof(thetarray));
	FORLIM = nlocus;
	for (i = 0; i < FORLIM; i++)
	  oldtheta[i] = 1 / (1 + exp(oldtheta[i]));
	WITH->theta[0] = oldtheta[0] + oldtheta[nlocus - 1];
	WITH->theta[nlocus - 2] = oldtheta[nlocus - 2] + oldtheta[nlocus - 1];
	WITH->theta[nlocus - 1] = oldtheta[0] + oldtheta[nlocus - 2];
	V.p1 = oldtheta[0];
	V.p2 = oldtheta[nlocus - 2];
	V.p3 = oldtheta[nlocus - 1];
	V.p4 = 1.0 - V.p1 - V.p2 - V.p3;
      } else {
	V.p1 = (WITH->theta[0] +
		WITH->theta[nlocus - 1] - WITH->theta[nlocus - 2]) / 2.0;
	V.p2 = (WITH->theta[nlocus - 2] +
		WITH->theta[nlocus - 1] - WITH->theta[0]) / 2.0;
	V.p3 = (WITH->theta[nlocus - 2] +
		WITH->theta[0] - WITH->theta[nlocus - 1]) / 2.0;
	V.p4 = 1.0 - V.p1 - V.p2 - V.p3;
      }
      recombine(WITH->theta, WITH->segprob, &V);
    } else
      recombine(femaletheta->theta, femaletheta->segprob, &V);
  }
  if (firsttime) {
    /*Introduced by Ramana*/
    nuneed=V.here;
    /*if (V.here < maxneed)
      printf(" Maxneed may be reduced to %12d\n", V.here); */
  }
}  /* recombination */


/* Local variables for getlocations: */
struct LOC_getlocations {
  int ngene, nseg, here, there, start, nhet, thisseg;
  boolean rarepresent, riskhom, riskhet;
  hapvector hap1, hap2;
  boolean thishet[maxlocus];
} ;

static boolean checkrare(LINK)
struct LOC_getlocations *LINK;
{
  int i;
  boolean check;
  int FORLIM;
  locusvalues *WITH;

  check = false;
  FORLIM = nlocus;
  for (i = 0; i < FORLIM; i++) {
    if (nohom[i]) {
      WITH = thislocus[i];
      if (WITH->freq[LINK->hap1[i] - 1] < minfreq ||
	  WITH->freq[LINK->hap2[i] - 1] < minfreq)
	check = true;
    }
  }
  return check;
}  /* checkrare */


static void checkrisk(riskhet, riskhom, LINK)
boolean *riskhet, *riskhom;
struct LOC_getlocations *LINK;
{
  *riskhet = false;
  *riskhom = false;
  if (LINK->hap1[risksys - 1] == riskall && LINK->hap2[risksys - 1] == riskall)
    *riskhom = true;
  else if ((LINK->hap1[risksys - 1] != riskall &&
	    LINK->hap2[risksys - 1] == riskall) ||
	   (LINK->hap2[risksys - 1] != riskall &&
	    LINK->hap1[risksys - 1] == riskall))
    *riskhet = true;
}  /* checkrisk */


static int gethapn(hap, LINK)
int *hap;
struct LOC_getlocations *LINK;
{
  int i, n, FORLIM;

  n = 1;
  FORLIM = nlocus;
  for (i = 0; i < FORLIM; i++)
    n += increment[i] * (hap[i] - 1);
  return n;
}  /* gethapn */

/* Local variables for domalerisk: */
struct LOC_domalerisk {
  struct LOC_getlocations *LINK;
} ;

static void setrisk(LINK)
struct LOC_domalerisk *LINK;
{
  int n;

  n = gethapn(LINK->LINK->hap1, LINK->LINK);
  if (LINK->LINK->hap1[risksys - 1] == riskall)
    riskmale[n - 1] = true;
  else
    riskmale[n - 1] = false;
}  /* setrisk */


static void getriskhap(system, LINK)
int system;
struct LOC_domalerisk *LINK;
{
  int i;
  locusvalues *WITH;
  int FORLIM;

  WITH = thislocus[system - 1];
  FORLIM = WITH->nallele;
  for (i = 1; i <= FORLIM; i++) {
    LINK->LINK->hap1[system - 1] = i;
    if (system != nlocus)
      getriskhap(system + 1, LINK);
    else
      setrisk(LINK);
  }
}  /* getriskhap */


static void domalerisk(LINK)
struct LOC_getlocations *LINK;
{
  struct LOC_domalerisk V;


  V.LINK = LINK;
  getriskhap(1, &V);
}  /* domalerisk */

/* Local variables for domutation: */
struct LOC_domutation {
  struct LOC_getlocations *LINK;
} ;

static void setmutation(LINK)
struct LOC_domutation *LINK;
{
  int i, n;

  n = gethapn(LINK->LINK->hap1, LINK->LINK);
  if (LINK->LINK->hap1[mutsys - 1] == thislocus[mutsys - 1]->nallele)
    muthap[n - 1] = n;
  else {
    i = LINK->LINK->hap1[mutsys - 1];
    LINK->LINK->hap1[mutsys - 1] = thislocus[mutsys - 1]->nallele;
    muthap[n - 1] = gethapn(LINK->LINK->hap1, LINK->LINK);
    LINK->LINK->hap1[mutsys - 1] = i;
  }
  hapstore[mutsys - 1][muthap[n - 1] - 1] = LINK->LINK->hap1[mutsys - 1];
}  /* setmutation */


static void getmuthap(system, LINK)
int system;
struct LOC_domutation *LINK;
{
  int i;
  locusvalues *WITH;
  int FORLIM;

  WITH = thislocus[system - 1];
  FORLIM = WITH->nallele;
  for (i = 1; i <= FORLIM; i++) {
    LINK->LINK->hap1[system - 1] = i;
    if (system != nlocus)
      getmuthap(system + 1, LINK);
    else
      setmutation(LINK);
  }
}  /* getmuthap */


static void domutation(LINK)
struct LOC_getlocations *LINK;
{
  struct LOC_domutation V;


  V.LINK = LINK;
  getmuthap(1, &V);
}  /* domutation */


static void setnumbers(LINK)
struct LOC_getlocations *LINK;
{
  int nhap1, nhap2, isys, FORLIM;

  LINK->ngene++;

  segstart[LINK->ngene - 1] = LINK->here + 1;
  probstart[LINK->ngene - 1] = LINK->there + 1;
  probend[LINK->ngene - 1] = LINK->there + LINK->nseg;

  LINK->there += LINK->nseg;

  nhap1 = gethapn(LINK->hap1, LINK);
  nhap2 = gethapn(LINK->hap2, LINK);
  FORLIM = nlocus;
  for (isys = 0; isys < FORLIM; isys++) {
    hapstore[isys][nhap1 - 1] = LINK->hap1[isys];
    hapstore[isys][nhap2 - 1] = LINK->hap2[isys];
  }
  genenumber[nhap1 - 1][nhap2 - 1] = LINK->ngene;
  genenumber[nhap2 - 1][nhap1 - 1] = LINK->ngene;

  if (minfreq != 0.0) {
    if (LINK->rarepresent)
      rare[LINK->ngene - 1] = true;
    else
      rare[LINK->ngene - 1] = false;
  } else
    rare[LINK->ngene - 1] = false;
  if (risk) {
    risk1[LINK->ngene - 1] = LINK->riskhet;
    risk2[LINK->ngene - 1] = LINK->riskhom;
  }

  LINK->thisseg++;
  invgenenum1[LINK->thisseg - 1] = nhap1;
  invgenenum2[LINK->thisseg - 1] = nhap2;
}  /* setnumbers */


static void hapscr(system, nscramble, LINK)
int system, nscramble;
struct LOC_getlocations *LINK;
{
  int i, j;

  if (LINK->thishet[system - 1])
    nscramble++;
  if (system != nlocus)
    hapscr(system + 1, nscramble, LINK);
  else
    setnumbers(LINK);
  if (nscramble <= 1)
    return;
  if (LINK->hap1[system - 1] == LINK->hap2[system - 1])
    return;
  i = LINK->hap1[system - 1];
  j = LINK->hap2[system - 1];
  LINK->hap1[system - 1] = j;
  LINK->hap2[system - 1] = i;
  if (system != nlocus)
    hapscr(system + 1, nscramble, LINK);
  else
    setnumbers(LINK);
  LINK->hap1[system - 1] = i;
  LINK->hap2[system - 1] = j;
}  /* hapscr */


static void sethap(system, LINK)
int system;
struct LOC_getlocations *LINK;
{
  int i, j;
  locusvalues *WITH;
  int FORLIM, FORLIM1;

  WITH = thislocus[system - 1];
  if (LINK->thishet[system - 1]) {
    FORLIM = WITH->nallele;
    for (i = 1; i < FORLIM; i++) {
      LINK->hap1[system - 1] = i;
      FORLIM1 = WITH->nallele;
      for (j = i + 1; j <= FORLIM1; j++) {
	LINK->hap2[system - 1] = j;
	if (system != nlocus)
	  sethap(system + 1, LINK);
	else {
	  LINK->rarepresent = checkrare(LINK);
	  if (risk)
	    checkrisk(&LINK->riskhet, &LINK->riskhom, LINK);
	  LINK->there = LINK->start;
	  LINK->thisseg = LINK->here;
	  hapscr(1, 0, LINK);
	  LINK->here += LINK->nseg;
	}
      }
    }
    return;
  }
  FORLIM = WITH->nallele;
  for (i = 1; i <= FORLIM; i++) {
    LINK->hap1[system - 1] = i;
    LINK->hap2[system - 1] = i;
    if (system != nlocus)
      sethap(system + 1, LINK);
    else {
      LINK->rarepresent = checkrare(LINK);
      if (risk)
	checkrisk(&LINK->riskhet, &LINK->riskhom, LINK);
      LINK->thisseg = LINK->here;
      LINK->there = LINK->start;
      hapscr(1, 0, LINK);
      LINK->here += LINK->nseg;
    }
  }
}  /* sethap */


static void starthap(LINK)
struct LOC_getlocations *LINK;
{
  int i, FORLIM;

  LINK->nseg = 1;
  FORLIM = LINK->nhet;
  for (i = 2; i <= FORLIM; i++)
    LINK->nseg *= 2;
  sethap(1, LINK);
  LINK->start = LINK->there;
}  /* starthap */


static void gethet1(system, LINK)
int system;
struct LOC_getlocations *LINK;
{
  LINK->thishet[system - 1] = false;
  if (system != nlocus)
    gethet1(system + 1, LINK);
  else
    starthap(LINK);
  LINK->thishet[system - 1] = true;
  LINK->nhet++;
  if (system != nlocus)
    gethet1(system + 1, LINK);
  else
    starthap(LINK);
  LINK->nhet--;
}  /* gethet1 */


void getlocations()
{
  struct LOC_getlocations V;


  V.nhet = 0;
  V.here = 0;
  V.there = 0;
  V.ngene = 0;
  V.start = 0;
  gethet1(1, &V);
  if (mutsys != 0)
    domutation(&V);
  if (sexlink && risk)
    domalerisk(&V);
}


/* Local variables for inputdata: */
struct LOC_inputdata {
  int nfactor[maxlocus];
} ;

static void inputerror(nerror, par1, par2)
int nerror, par1, par2;
{
  printf("Fatal error detected in procedure inputdata\n");
  switch (nerror) {

  case 0:
    printf("Number of loci %2d exceeds the constant maxlocus\n", par1);
    break;

  case 1:
    printf("Number of loci read %2d. Less than minimum of 1\n", par1);
    break;

  case 2:
    printf(
      "Error detected reading loci order. Locus number %2d in position %2d exceeds number of loci\n",
      par2, par1);
    break;

  case 3:
    printf(
      "Error detected reading loci order. Illegal locus number %2d in position %2d\n",
      par2, par1);
    break;

  case 4:
    printf(
      "Error detected reading loci order. Locus number repeated in positions %2d and %2d\n",
      par1, par2);
    break;

  case 5:
    printf(
      "Error detected reading locus description. Illegal locus type %2d for locus %2d\n",
      par2, par1);
    break;

  case 6:
    printf(
      "Error detected reading locus description for system %2d. Number of alleles  %2d exceeds maxall\n",
      par1, par2); /*bug fix, was par1 instead of par2*/
    break;

  case 7:
    printf(
      "Error detected reading locus description for system %2d. Illegal number of alleles  %2d\n",
      par1, par2);
    break;

  case 8:
    printf(
      "Error detected reading locus description for system %2d. Number of factors  %2d exceeds maxfact\n",
      par1, par2);
    break;

  case 9:
    printf(
      "Error detected reading locus description for system %2d. Illegal number of factors  %2d\n",
      par1, par2);
    break;

  case 10:
    printf(
      "Error detected reading locus description for system %2d. Alleles not codominant\n",
      par1);
    break;

  case 11:
    printf("Error detected reading pedigree record %2d. Illegal code for sex %2d\n",
	   par1, par2);
    break;

  case 12:
    printf(
      "Error detected reading pedigree record at pedigree%2d. Maximum number of pedigree records exceeded\n",
      par1);
    break;

  case 13:
    printf(
      "Error detected reading pedigree record %2d. Maximum number of individuals exceeded\n",
      par1);
    break;

  case 14:
    printf(
      "Error detected reading pedigree record %2d. Illegal binary factor code %2d\n",
      par1, par2);
    break;

  case 15:
    printf(
      "Error detected reading pedigree record %2d. No allelic pair for genotype\n",
      par1);
    break;

  case 16:
    printf(
      "Error detected reading pedigree record %2d. Allele number %2d exceeds maxall\n",
      par1, par2);
    break;

  case 17:
    printf(
      "Error detected reading pedigree record %2d. Illegal allele number %2d\n",
      par1, par2);
    break;

  case 18:
    printf("Number of systems after factorization (%3d) exceeds maxsystem\n",
	   par1);
    break;

  case 19:
    printf("Number of systems after factorization (%3d) less than minimum of 1\n",
	   par1);
    break;

  case 20:
    printf("Number of recombination types (%3d) exceeds maxrectype\n", par1);
    break;

  case 21:
    printf("Number of recombination types (%3d) less than minimum of 1\n",
	   par1);
    break;

  case 22:
    printf(
      "End of file detected in tempdat by procedure readthg before all data found\n");
    break;

  case 23:
    printf(
      "Error detected reading iterated locus in simdata. Value (%3d) greater than nlocus\n",
      par1);
    break;

  case 24:
    printf("Error detected reading iterated locus in simdata. Illegal value (%3d)\n",
	   par1);
    break;

  case 25:
    printf("Number of iterated parameters greater then maxn\n");
    break;

  case 26:
    printf(
      "Error detected reading pedigree record %2d. Liability class (%2d) exceeds nclass\n",
      par1, par2);
    break;

  case 27:
    printf(
      "Error detected reading pedigree record %2d. Illegal liability class (%2d)\n",
      par1, par2);
    break;

  case 28:
    printf(
      "Error detected reading locus description for system%2d. Liability classes (%3d) exceed maxliab\n",
      par1, par2);
    break;

  case 29:
    printf(
      "Error detected reading locus description for system%2d. Illegal number of liability classes (%3d)\n",
      par1, par2);
    break;

  case 30:
    printf(
      "Error detected reading locus description for system%2d. Penetrance out of range\n",
      par1);
    break;

  case 31:
    printf(
      "Error detected reading locus description for system%2d. Number of traits (%3d) can only be one.\n",
      par1, par2);
    break;

  case 32:
    printf(
      "Error detected reading locus description for system%2d. Number of traits out of range (%3d)\n",
      par1, par2);
    break;

  case 33:
    printf(
      "Error detected reading locus description for system%2d. Variance must be positive\n",
      par1);
    break;

  case 34:
    printf(
      "Error detected reading locus description for system%2d. Variance multiplier must be positive\n",
      par1);
    break;

  case 35:
    printf(
      "Error detected reading locus description for system%2d. Risk allele %3d) exceeds nallele\n",
      par1, par2);
    break;

  case 36:
    printf(
      "Error detected reading locus description for system%2d. Illegal risk allele (%3d)\n",
      par1, par2);
    break;

  case 37:
    printf("Error detected reading simdata. Risk locus %3d) exceeds nlocus\n",
	   par2);
    break;

  case 38:
    printf("Error detected reading simdata. Illegal value for risk locus %3d)\n",
	   par2);
    break;

  case 39:
    printf("Error detected reading simdata. Mutation locus %3d) exceeds nlocus\n",
	   par2);
    break;

  case 40:
    printf("Error detected reading simdata. Illegal value for mutation locus %3d)\n",
	   par2);
    break;
  }
  printf("Press return to halt...\007\n");
  scanf("%*[^\n]");
  getchar();
  exit(0);   /*SUN*/
}  /* inputerror */

/*inputerror*/

static void inputwarning(nwarning, par1, par2)
int nwarning, par1, par2;
{
  switch (nwarning) {   /*Version 5.1*/

  case 0:
    printf("Illegal sex difference parameter %2d Parameter should be 0, 1, or 2\n",
	   par1);
    break;

  case 1:
    printf("Illegal interference parameter %2d Lack of interference assumed\n",
	   par1);
    break;

  case 2:
    printf(
      "Illegal sex difference parameter %2d Parameter must be 0 with sex-linked data\n",
      par1);
    break;

  case 3:
    printf("Warning: No risk calculations are permitted during the simulation.\n");
    printf("However the risk locus is %3d\n", par2);
    break;

  case 4:
    printf(
      "Non-standard affection status%4d interpreted as normal in pedigree record%5d\n",
      par2, par1);
    break;
  }
  printf("Press return to continue...\007\n");
  scanf("%*[^\n]");
  getchar();
}  /* inputwarning */

/* Local variables for readped: */
struct LOC_readped {
  struct LOC_inputdata *LINK;
  int sequence;
} ;

/* Local variables for getphenotype: */
struct LOC_getphenotype {
  struct LOC_readped *LINK;
  thisperson **p;
  int system;
} ;


static void readbin(phen, thislocus, LINK)
phenotype **phen;
locusvalues *thislocus;
struct LOC_getphenotype *LINK;
{
  int i, j;
  phenotype *WITH;
  int FORLIM;

  WITH = *phen;
  WITH->which = binary_;
  WITH->phenf = 0;
  FORLIM = LINK->LINK->LINK->nfactor[LINK->system - 1];
  for (i = 1; i <= FORLIM; i++) {
    fscanf(simped, "%d", &j);
    if (j != 0 && j != 1)
      inputerror(14, (*LINK->p)->id, j);
    if (j == 1)
      WITH->phenf = ((int)WITH->phenf) | (1 << ((int)i));
  }
}  /* readbin */


static void readnumber(phen, thislocus, LINK)
phenotype **phen;
locusvalues *thislocus;
struct LOC_getphenotype *LINK;
{
  int i, j;
  phenotype *WITH;

  WITH = *phen;
  WITH->which = binary_;
  WITH->alleles[0] = 0; /*new representation*/
  WITH->alleles[1] = 0; /*new representation*/
  for (i = 1; i <= 2; i++) {
    fscanf(simped, "%d", &j);
    if (j > maxall)
      inputerror(16, (*LINK->p)->id, j);
    if (j < 0)
      inputerror(17, (*LINK->p)->id, j);
    WITH->alleles[i-1] = j;
  }
}  /* readnumber */


static void readaff(phen, thislocus, LINK)
phenotype **phen;
locusvalues *thislocus;
struct LOC_getphenotype *LINK;
{
  int thisval;   /*Version 5.1 has this as an integer*/
  phenotype *WITH;

  WITH = *phen;
  WITH->which = affection;
  fscanf(simped, "%d", &thisval);
  if (thisval == missaff)
    WITH->aff = 0;
  else if (thisval == affval)
    WITH->aff = 2;
  else {
    if (thisval != 1)
      inputwarning(4, (*LINK->p)->id, thisval);
    WITH->aff = 1;
  }
  if (thislocus->UU.U0.nclass == 1)
    WITH->liability = 1;
  else
    fscanf(simped, "%d", &WITH->liability);
  if (WITH->liability > thislocus->UU.U0.nclass)
    inputerror(26, (*LINK->p)->id, WITH->liability);
  if (WITH->liability <= 0)
    inputerror(27, (*LINK->p)->id, WITH->liability);

  /*Added warning => Version 5.1*/
}  /* readaff */


static void readquan(phen, thislocus, LINK)
phenotype **phen;
locusvalues *thislocus;
struct LOC_getphenotype *LINK;
{
  int i;
  double xval;
  phenotype *WITH;
  int FORLIM;

  WITH = *phen;
  if (!sexlink || !(*LINK->p)->male) {
    WITH->which = quantitative;
    FORLIM = thislocus->UU.U1.ntrait;
    for (i = 0; i < FORLIM; i++)
      fscanf(simped, "%lg", &WITH->x[i]);
    WITH->missing = true;
    FORLIM = thislocus->UU.U1.ntrait;
    for (i = 0; i < FORLIM; i++) {
      if (WITH->x[i] != missval)
	WITH->missing = false;
    }
    return;
  }
  WITH->which = affection;
  fscanf(simped, "%lg", &xval);
  if (xval == missval)
    WITH->aff = missaff;
  else if (xval == affall)
    WITH->aff = affall;
  else
    WITH->aff = unafqu;
  WITH->liability = 1;
  FORLIM = thislocus->UU.U1.ntrait;
  for (i = 2; i <= FORLIM; i++)
    fscanf(simped, "%lg", &xval);
}  /* readquan */

static void getphenotype(p_, LINK)
thisperson **p_;
struct LOC_readped *LINK;
{
  struct LOC_getphenotype V;
  int thisread;
  thisperson *WITH;
  int FORLIM;


  V.LINK = LINK;
  V.p = p_;
  WITH = *V.p;
  FORLIM = nlocus;
  for (thisread = 0; thisread < FORLIM; thisread++) {
    V.system = order[thisread];
    WITH->phen[V.system - 1] = NULL;
    if (thislocus[V.system - 1]->which != null_)
      WITH->phen[V.system - 1] = (phenotype *)Malloc(sizeof(phenotype));
    /*         phen[system]:=phenpoint(NewPtr(SizeOf(phenotype)));*/
    switch (thislocus[V.system - 1]->which) {

    case quantitative:
      readquan(&WITH->phen[V.system - 1], thislocus[V.system - 1], &V);
      break;

    case affection:
      readaff(&WITH->phen[V.system - 1], thislocus[V.system - 1], &V);
      break;

    case binary_:
      if (thislocus[V.system - 1]->format == 3)
	readnumber(&WITH->phen[V.system - 1], thislocus[V.system - 1], &V);
      else
	readbin(&WITH->phen[V.system - 1], thislocus[V.system - 1], &V);
      break;
    }
    if (lastpriv == V.system) {
      WITH->privphen = (phenotype *)Malloc(sizeof(phenotype));
      /*          privphen:=phenpoint(NewPtr(SizeOf(phenotype)));*/
      switch (thislocus[V.system - 1]->privlocus->which) {

      case quantitative:
	readquan(&WITH->privphen, thislocus[V.system - 1]->privlocus, &V);
	break;

      case affection:
	readaff(&WITH->privphen, thislocus[V.system - 1]->privlocus, &V);
	break;
      }
    }
  }
}  /* getphenotype */


static void getind(id, LINK)
int *id;
struct LOC_readped *LINK;
{
  int idorig;   /*getind*/

  fscanf(simped, "%d", id);
  if (*id == 0) {
    return;
  }  /*Store original ids*/
  idorig = *id;
  *id += LINK->sequence;
  if (*id > maxind)
    inputerror(13, *id, *id);
  if (person[*id] == NULL)
    person[*id] = (thisperson *)Malloc(sizeof(thisperson));
  /*        person[id]:=ind(NewPtr(SizeOf(thisperson))); ;*/
  person[*id]->origid = idorig;
}  /* getind */

/*getind*/

static void multimarriage(p, LINK)
thisperson **p;
struct LOC_readped *LINK;
{
  thisperson *q, *child;   /*multimarriage*/
  thisperson *WITH;

  if ((*p)->foff == NULL) {
    (*p)->multi = false;
    return;
  }
  WITH = *p;
  if (WITH->male)
    q = WITH->foff->ma;
  else
    q = WITH->foff->pa;
  child = WITH->foff;
  (*p)->multi = false;
  do {
    if (WITH->male) {
      WITH->multi = (q == child->ma);
      child = child->nextpa;
    } else {
      WITH->multi = (q == child->pa);
      child = child->nextma;
    }
  } while (!(child == NULL || WITH->multi));
}  /* multimarriage */


static void readped(LINK)
struct LOC_inputdata *LINK;
{
  struct LOC_readped V;
  int i, newid, sex_, profield, newped, nuperson, thisone, thisped;
  int startped[maxped], endped[maxped];
  thisperson *holdloop;
  int FORLIM;
  thisperson *WITH;

  /*multimarriage*/

  V.LINK = LINK;
  for (i = 0; i <= maxind; i++)
    person[i] = NULL;
  V.sequence = 0;
  nuperson = 0;
  nuped = 1;
  for (i = 0; i < maxloop; i++) {
    looppers[nuped - 1][i][0] = NULL;
    looppers[nuped - 1][i][1] = NULL;
  }
  proband[nuped - 1] = NULL;
  fscanf(simped, "%d", &newped);
  thisped = newped;
  startped[0] = 1;
  while (!P_eof(simped)) { 
    nuperson++;
    getind(&thisone, &V);
    if (proband[nuped - 1] == NULL)
      proband[nuped - 1] = person[thisone];
    WITH = person[thisone];
    /*               ped := thisped;*/
    WITH->ped = nuped;   /* L. Sandkuijl's suggestion */
    WITH->origped = thisped;   /* L. Sandkuijl's suggestion */
    WITH->id = thisone;
    getind(&newid, &V);
    WITH->pa = person[newid];
    getind(&newid, &V);
    WITH->ma = person[newid];
    getind(&newid, &V);
    WITH->foff = person[newid];
    getind(&newid, &V);
    WITH->nextpa = person[newid];
    getind(&newid, &V);
    WITH->nextma = person[newid];
    fscanf(simped, "%d", &sex_);
    if (sex_ != 1 && sex_ != 2)
      inputerror(11, WITH->id, sex_);
    if (sex_ == 1) {
      WITH->male = true;
      men++;
    } else {
      WITH->male = false;
      women++;
    }
    WITH->unknown = false;
    WITH->inloop = 0;
    fscanf(simped, "%d", &profield);
    person[thisone]->prbnum = profield;
    if (profield == 1)
      proband[nuped - 1] = person[thisone];
    else if (profield > 1 && profield - 1 <= maxloop) {
      if (looppers[nuped - 1][profield - 2][1] == NULL)
	looppers[nuped - 1][profield - 2][1] = person[thisone];
      else
	looppers[nuped - 1][profield - 2][0] = person[thisone];
    }
    getphenotype(&person[thisone], &V);
    /* We now require the last column to contain integers: a non-zero
    value in this column indicates that the person is available for
    marker-typing. Learned from Lodewijk Sandkuyl's code.*/
    fscanf(simped, "%d%*[^\n]", &person[thisone]->input_availnum);
    if (all_available) {
      person[thisone]->availnum = 2;
    }
    else {
      person[thisone]->availnum = person[thisone]->input_availnum;
    }
    getc(simped);
    WITH = person[thisone];   /*Count the number in each code*/
    if (WITH->input_availnum != 0 && WITH->input_availnum != 1 && WITH->input_availnum != 2 &&
	WITH->input_availnum != 3) {
      fprintf(simout,
	"ERROR: illegal availability code for person %6d in pedigree %6d\n",
	WITH->origid, WITH->origped);
      printf("ERROR: illegal availability code for person %6d in pedigree %6d\n",
	     WITH->origid, WITH->origped);
      printf("Press return to halt...\n");
      scanf("%*[^\n]");
      getchar();
      exit(0);   /*SUN*/
    }
    code[WITH->input_availnum]++;
    if (!P_eof(simped)) /*Testing EOF*/
      fscanf(simped, "%d", &newped);
    if (newped == thisped)
      continue;
    V.sequence += nuperson;
    endped[nuped - 1] = V.sequence;
    nuperson = 0;
    nuped++;
    if (nuped > maxped)
      inputerror(12, newped, nuped);
    startped[nuped - 1] = V.sequence + 1;
    for (i = 0; i < maxloop; i++) {
      looppers[nuped - 1][i][0] = NULL;
      looppers[nuped - 1][i][1] = NULL;
    }
    proband[nuped - 1] = NULL;
    thisped = newped;
  }
  totperson = V.sequence + nuperson;
  endped[nuped - 1] = totperson;
  FORLIM = nuped;
  for (newped = 0; newped < FORLIM; newped++) {
    if (looppers[newped][0][1] != NULL && looppers[newped][0][0] == NULL)
      looppers[newped][0][0] = proband[newped];
    for (i = 0; i < maxloop; i++) {
      if (looppers[newped][i][0] == NULL)
	looppers[newped][i][1] = NULL;
      else {
	looppers[newped][i][0]->inloop = i + 1;
	looppers[newped][i][1]->inloop = i + 1;
	if (looppers[newped][i][0]->pa == NULL &&
	    looppers[newped][i][1]->pa != NULL) {
	  holdloop = looppers[newped][i][0];
	  looppers[newped][i][0] = looppers[newped][i][1];
	  looppers[newped][i][1] = holdloop;
	}
      }
    }
  }
  FORLIM = totperson;
  for (thisone = 1; thisone <= FORLIM; thisone++)
    multimarriage(&person[thisone], &V);
}  /* readped */

/* Local variables for readloci: */
struct LOC_readloci {
  struct LOC_inputdata *LINK;
  int i, whichtype, nupriv;
} ;

/* Local variables for getlocus: */
struct LOC_getlocus {
  struct LOC_readloci *LINK;
  int system;
} ;

static void getbin(locus, system, LINK)
locusvalues **locus;
int *system;
struct LOC_getlocus *LINK;
{
  int i, j, k;
  locusvalues *WITH;
  int FORLIM, FORLIM1;

  fscanf(simdata, "%d%*[^\n]", &LINK->LINK->LINK->nfactor[*system - 1]);
  getc(simdata);
  if (LINK->LINK->LINK->nfactor[*system - 1] > maxfact)
    inputerror(8, *system, LINK->LINK->LINK->nfactor[*system - 1]);
  if (LINK->LINK->LINK->nfactor[*system - 1] <= 0)
    inputerror(9, *system, LINK->LINK->LINK->nfactor[*system - 1]);
  WITH = *locus;
  WITH->UU.U2.numfact = LINK->LINK->LINK->nfactor[*system - 1];
  FORLIM = WITH->nallele;
  for (i = 0; i < FORLIM; i++){
    WITH->UU.U2.bin_allele[i] = 0; /*new representation*/
  }
  FORLIM = WITH->nallele;
  for (i = 0; i < FORLIM; i++) {
    FORLIM1 = LINK->LINK->LINK->nfactor[*system - 1];
    for (j = 1; j <= FORLIM1; j++) {
      fscanf(simdata, "%d", &k);
      if (k == 1)
    	WITH->UU.U2.bin_allele[i] = ((int)WITH->UU.U2.bin_allele[i]) | (1 << ((int)j)); /*new representation*/
    }
  }
  fscanf(simdata, "%*[^\n]");
  getc(simdata);
}  /* getbin */


static void getnumber(locus, system, LINK)
locusvalues **locus;
int *system;
struct LOC_getlocus *LINK;
{
  int i;
  locusvalues *WITH;
  int FORLIM;

  WITH = *locus;
  FORLIM = WITH->nallele;
  for (i = 1; i <= FORLIM; i++) {
    if (i > maxall) {
      fprintf(simout, " *** ERROR: maxall not big enough ***\n");
      printf(" *** ERROR: maxall not big enough ***\n");
      printf(" Press return to halt...\n");
      scanf("%*[^\n]");
      getchar();
      exit(0);   /*SUN*/
    }
    WITH->UU.U2.num_allele[i - 1][0] = i; /*new representation*/
  }
}  /* getnumber */


static void getpen(locus, LINK)
locusvalues **locus;
struct LOC_getlocus *LINK;
{
  int i, j, k, l;
  locusvalues *WITH;
  int FORLIM, FORLIM1, FORLIM2;

  WITH = *locus;
  fscanf(simdata, "%d%*[^\n]", &WITH->UU.U0.nclass);
  getc(simdata);
  if (WITH->UU.U0.nclass > maxliab)
    inputerror(28, LINK->system, WITH->UU.U0.nclass);
  if (WITH->UU.U0.nclass <= 0)
    inputerror(29, LINK->system, WITH->UU.U0.nclass);
  FORLIM = WITH->UU.U0.nclass;
  for (l = 0; l < FORLIM; l++) {
    FORLIM1 = WITH->nallele;
    for (i = 1; i <= FORLIM1; i++) {
      FORLIM2 = WITH->nallele;
      for (j = i - 1; j < FORLIM2; j++) {
	fscanf(simdata, "%lg", &WITH->UU.U0.pen[i][j][2][l]);
	if ((unsigned)WITH->UU.U0.pen[i][j][2][l] > 1.0)
	  inputerror(30, LINK->system, LINK->system);
	WITH->UU.U0.pen[i][j][1][l] = 1 - WITH->UU.U0.pen[i][j][2][l];
	WITH->UU.U0.pen[i][j][0][l] = 1.0;
	for (k = 0; k <= 2; k++)
	  WITH->UU.U0.pen[j + 1][i - 1][k][l] = WITH->UU.U0.pen[i][j][k][l];
      }
    }
    fscanf(simdata, "%*[^\n]");
    getc(simdata);
    FORLIM1 = WITH->nallele;
    for (i = 0; i < FORLIM1; i++)
      WITH->UU.U0.pen[0][i][0][l] = 1.0;
    if (sexlink) {
      FORLIM1 = WITH->nallele;
      for (i = 0; i < FORLIM1; i++) {   /*Corrected to match Version 5.1*/
	fscanf(simdata, "%lg", &WITH->UU.U0.pen[0][i][2][l]);
	if ((unsigned)WITH->UU.U0.pen[0][i][2][l] > 1.0)
	  inputerror(30, LINK->system, LINK->system);
      }
      FORLIM1 = WITH->nallele;
      for (i = 0; i < FORLIM1; i++)
	WITH->UU.U0.pen[0][i][1][l] = 1.0 - WITH->UU.U0.pen[0][i][2][l];
      fscanf(simdata, "%*[^\n]");
      getc(simdata);
    }
  }
}  /* getpen */


static void getquan(locus, privelege, LINK)
locusvalues **locus;
boolean privelege;
struct LOC_getlocus *LINK;
{
  /*Get information of a quantitative locus. Privelege says whether it is
  a priveleged locus or not*/
  int i, j, k;
  locusvalues *WITH;
  int FORLIM, FORLIM1, FORLIM2;

  WITH = *locus;
  fscanf(simdata, "%d%*[^\n]", &WITH->UU.U1.ntrait);
  getc(simdata);
  if (WITH->UU.U1.ntrait > 1)
    inputerror(31, LINK->system, WITH->UU.U1.ntrait);
  if (WITH->UU.U1.ntrait <= 0)
    inputerror(32, LINK->system, WITH->UU.U0.nclass);
  FORLIM = WITH->UU.U1.ntrait;
  for (i = 0; i < FORLIM; i++) {
    FORLIM1 = WITH->nallele;
    for (j = 1; j <= FORLIM1; j++) {
      FORLIM2 = WITH->nallele;
      for (k = j; k <= FORLIM2; k++) {
	fscanf(simdata, "%lg", &WITH->UU.U1.pm[j][k - 1][i]);
	WITH->UU.U1.pm[k][j - 1][i] = WITH->UU.U1.pm[j][k - 1][i];
      }
    }
  }
  fscanf(simdata, "%*[^\n]");
  getc(simdata);
  if (privelege && LINK->LINK->nupriv != lastpriv)
    return;
  FORLIM = WITH->UU.U1.ntrait;
  for (i = 0; i < FORLIM; i++) {
    FORLIM1 = WITH->UU.U1.ntrait;
    for (j = i; j < FORLIM1; j++) {
      fscanf(simdata, "%lg", &WITH->UU.U1.vmat[i][j]);
      if (i + 1 == j + 1 && WITH->UU.U1.vmat[i][j] <= 0.0)
	inputerror(33, LINK->system, LINK->system);
      WITH->UU.U1.vmat[j][i] = WITH->UU.U1.vmat[i][j];
    }
  }
  fscanf(simdata, "%*[^\n]");
  getc(simdata);
  invertmat(WITH->UU.U1.vmat, WITH->UU.U1.ntrait, &WITH->UU.U1.det);
  WITH->UU.U1.det = 1 / sqrt(WITH->UU.U1.det);
  fscanf(simdata, "%lg%*[^\n]", &WITH->UU.U1.conmat);
  getc(simdata);
  if (WITH->UU.U1.conmat <= 0.0)
    inputerror(34, LINK->system, LINK->system);
  WITH->UU.U1.conmat = 1 / WITH->UU.U1.conmat;
  WITH->UU.U1.contrait = 1.0;
  FORLIM = WITH->UU.U1.ntrait;
  for (i = 1; i <= FORLIM; i++)
    WITH->UU.U1.contrait *= WITH->UU.U1.conmat;
  WITH->UU.U1.contrait = sqrt(WITH->UU.U1.contrait);
}  /* getquan */

static void getlocus(system_, LINK)
int system_;
struct LOC_readloci *LINK;
{
  struct LOC_getlocus V;
  int j;
  locusvalues *WITH, *WITH1;
  int FORLIM;
  int TEMP;


  V.LINK = LINK;
  V.system = system_;
  thislocus[V.system - 1] = (locusvalues *)Malloc(sizeof(locusvalues));
  WITH = thislocus[V.system - 1];
  /*     thislocus[system]:=locuspoint(NewPtr(SizeOf(locusvalues)));*/
  WITH->privlocus = NULL;
  fscanf(simdata, "%d%d", &LINK->whichtype, &WITH->nallele);
  if (LINK->whichtype < 0 && LINK->whichtype > 4)
    inputerror(5, V.system, LINK->whichtype);
  if (WITH->nallele > maxall)
    inputerror(6, V.system, WITH->nallele);
  if (WITH->nallele <= 0)
    inputerror(7, V.system, WITH->nallele);
  switch (LINK->whichtype) {

  case 0:
    WITH->which = quantitative;
    break;

  case 1:
    WITH->which = affection;
    break;

  case 2:
    WITH->which = binary_;
    break;

  case 3:
    WITH->which = binary_;
    break;
  }
  WITH->format = LINK->whichtype;
  if (lastpriv == 0) {
    fscanf(simdata, "%*[^\n]");
    getc(simdata);
  } else {
    fscanf(simdata, "%d%*[^\n]", &LINK->whichtype);
    getc(simdata);
    if (LINK->whichtype == 0 || LINK->whichtype == 1) {
      LINK->nupriv++;
      WITH->privlocus = (locusvalues *)Malloc(sizeof(locusvalues));
      /*           privlocus:=locuspoint(NewPtr(SizeOf(locusvalues)));*/
      WITH->privlocus->nallele = WITH->nallele;
      WITH1 = WITH->privlocus;
      switch (LINK->whichtype) {

      case 0:
	WITH1->which = quantitative;
	break;

      case 1:
	WITH1->which = affection;
	break;
      }
    }
  }
  if (!disequi) {
    FORLIM = WITH->nallele;
    for (j = 0; j < FORLIM; j++)
      fscanf(simdata, "%lg", &WITH->freq[j]);
    fscanf(simdata, "%*[^\n]");
    getc(simdata);
  }
  switch (WITH->which) {

  case binary_:
    if (WITH->format == 3)
      getnumber(&thislocus[V.system - 1], &V.system, &V);
    else
      getbin(&thislocus[V.system - 1], &V.system, &V);
    break;

  case affection:
    getpen(&thislocus[V.system - 1], &V);
    break;

  case quantitative:
    getquan(&thislocus[V.system - 1], false, &V);
    break;
  }
  if (WITH->privlocus != NULL) {
    switch (WITH->privlocus->which) {

    case affection:
      getpen(&WITH->privlocus, &V);
      break;

    case quantitative:
      getquan(&WITH->privlocus, true, &V);
      break;
    }
  }
  if (LINK->nupriv == lastpriv && lastpriv != 0)
    lastpriv = V.system;
  if (!risk)
    return;
  if (LINK->i == risksys) {
    fscanf(simdata, "%d%*[^\n]", &TEMP);
    getc(simdata);
    riskall = TEMP;
  }
  if (riskall > thislocus[LINK->i - 1]->nallele)
    inputerror(35, LINK->i, (int)riskall);
  if (riskall < 0)
    inputerror(36, LINK->i, (int)riskall);
}  /* getlocus */


static void gettheta(sex_, LINK)
thetavalues **sex_;
struct LOC_readloci *LINK;
{
  thetarray oldtheta;
  int i;
  thetavalues *WITH;
  int FORLIM;

  *sex_ = (thetavalues *)Malloc(sizeof(thetavalues));
  for (i = 0; i < maxlocus; i++)
    (*sex_)->theta[i] = 0.0;
  nuneed = 7;
  for(i = 2; i < mlocus; i++)
    nuneed = 5 * nuneed - 3;
  if (*sex_ == NULL)
    malloc_err("item of type thetavalues");
   /*Next line added by A. A. Schaffer*/
  (*sex_)->segprob = (double*) malloc(nuneed * sizeof(double));
  if ((*sex_)->segprob == NULL)
    malloc_err("a segprob array in procedure gettheta");
  WITH = *sex_;
  if (*sex_ == maletheta || readfemale) {
    FORLIM = mlocus - 2;
    for (i = 0; i <= FORLIM; i++)
      fscanf(simdata, "%lf", &(*sex_)->theta[i]);
    if (interfer && !mapping)
      fscanf(simdata, "%lf", &(*sex_)->theta[mlocus - 1]);
    fscanf(simdata, "%*[^\n]");
    getc(simdata);
  } else {
    fscanf(simdata, "%lf%*[^\n]", &distratio);
    getc(simdata);
  }
  /*         FOR j:=1 TO maxneed DO segprob[j]:=0.0;*/
  if (!interfer || mapping) {
    return;
  }  /*=ln(1/0.0001-1.0)*/
  memcpy(oldtheta, WITH->theta, sizeof(thetarray));
  WITH->theta[0] = (oldtheta[0] + oldtheta[nlocus - 1] - oldtheta[nlocus - 2]) / 2.0;
  WITH->theta[nlocus - 2] =
    (oldtheta[nlocus - 2] + oldtheta[nlocus - 1] - oldtheta[0]) / 2.0;
  WITH->theta[nlocus - 1] =
    (oldtheta[0] + oldtheta[nlocus - 2] - oldtheta[nlocus - 1]) / 2.0;
  FORLIM = nlocus;
  for (i = 0; i < FORLIM; i++) {
    if (WITH->theta[i] > 0.0) {
      if (WITH->theta[i] < 0.999)
	WITH->theta[i] = log(1.0 / WITH->theta[i] - 1.0);
      else
	WITH->theta[i] = -6.9;   /*=ln(1/0.999-1.0)*/
    } else
      WITH->theta[i] = 9.21;
  }
}  /* gettheta */

/*readped*/

static void readloci(LINK)
struct LOC_inputdata *LINK;
{
  struct LOC_readloci V;
  int j, coupling, autosomal, independent, difference, FORLIM, FORLIM1;
  locusvalues *WITH;


  V.LINK = LINK;
  lastpriv = 0;
  fscanf(simdata, "%d%d%d%*[^\n]", &nlocus, &risksys, &autosomal);
  getc(simdata);
  mlocus = nlocus; /*A. A. Schaffer*/
  /*Replace the above line with the next when using epistasis*/
  /*readln(simdata,nlocus,risksys,autosomal,lastpriv);*/
  if (nlocus > maxlocus)
    inputerror(0, nlocus, nlocus);
  if (nlocus <= 0)
    inputerror(1, nlocus, nlocus);
  if (risksys > maxlocus)
    inputerror(37, risksys, risksys);
  if (risksys < 0)
    inputerror(38, risksys, risksys);
  if (risksys != 0)
    inputwarning(3, risksys, risksys);
  risk = false;   /*risksys<>0*/
  sexlink = (autosomal == 1);
  printf("YOU ARE USING LINKAGE/SLINK (V%s) WITH%3d-POINT", version, nlocus);
  if (sexlink)
    printf(" SEXLINKED DATA\n");
  else
    printf(" AUTOSOMAL DATA\n");
  fscanf(simdata, "%d%lg%lg%d%*[^\n]", &mutsys, &mutmale, &mutfemale,
	 &coupling);
  getc(simdata);
  if (mutsys > maxlocus)
    inputerror(39, mutsys, mutsys);
  if (mutsys < 0)
    inputerror(40, mutsys, mutsys);
  if (coupling == 1)
    disequi = true;
  else
    disequi = false;
  if (disequi)
    hapfreq = (thisarray *)Malloc(sizeof(thisarray));
  FORLIM = nlocus;
  /*     hapfreq:=genpoint(NewPtr(SizeOf(thisarray)));*/
  for (V.i = 1; V.i <= FORLIM; V.i++) {
    fscanf(simdata, "%d", &j);
    if (j > nlocus)
      inputerror(2, V.i, j);
    if (j <= 0)
      inputerror(3, V.i, j);
    order[j - 1] = V.i;
  }
  FORLIM = nlocus;
  for (V.i = 1; V.i <= FORLIM; V.i++) {
    FORLIM1 = V.i;
    for (j = 1; j < FORLIM1; j++) {
      if (order[V.i - 1] == order[j - 1])
	inputerror(4, V.i, j);
    }
  }
  fscanf(simdata, "%*[^\n]");
  getc(simdata);
  if (mutsys != 0)
    mutsys = order[mutsys - 1];
  V.nupriv = 0;
  FORLIM = nlocus;
  for (V.i = 1; V.i <= FORLIM; V.i++)
    getlocus(order[V.i - 1], &V);
  if (risksys != 0)
    risksys = order[risksys - 1];
  increment[nlocus - 1] = 1;
  for (V.i = nlocus - 1; V.i >= 1; V.i--)
    increment[V.i - 1] = increment[V.i] * thislocus[V.i]->nallele;
  fgeno = 1;
  FORLIM = nlocus;
  for (j = 0; j < FORLIM; j++)
    fgeno *= thislocus[j]->nallele;
  mgeno = fgeno;
  nuhap = fgeno;
  FORLIM = nlocus;
  for (V.i = 1; V.i <= FORLIM; V.i++)
    nohom[V.i - 1] = false;
  if (disequi) {
    allocate_thisarray(hapfreq, mgeno);
    for (V.i = 1; V.i <= mgeno; V.i++)
      fscanf(simdata, "%lg", &hapfreq->genarray[V.i - 1]);
    fscanf(simdata, "%*[^\n]");
    getc(simdata);
  } else {
    FORLIM = nlocus;
    for (V.i = 1; V.i <= FORLIM; V.i++) {
      WITH = thislocus[V.i - 1];
      if (WITH->which == affection || WITH->which == quantitative) {
	if (WITH->freq[affall - 1] < minfreq)
	  nohom[V.i - 1] = true;
      }
    }
  }
  fgeno = fgeno * (fgeno + 1) / 2;
  if (!sexlink)
    mgeno = fgeno;
  fscanf(simdata, "%d", &difference);
  if ((unsigned int)difference > 2) {
    inputwarning(0, difference, difference);
    difference = 0;
  }
  sexdif = (difference != 0);
  readfemale = (difference == 2);
  fscanf(simdata, "%d%*[^\n]", &independent);
  getc(simdata);
  if ((unsigned int)independent > 2) {
    inputwarning(1, independent, independent);
    independent = 0;
  }
  interfer = (independent != 0);
  mapping = (independent == 2);
  gettheta(&maletheta, &V);
  if (sexdif)
    gettheta(&femaletheta, &V);
  else
    femaletheta = maletheta;
  if (sexlink && sexdif) {
    inputwarning(2, difference, difference);
    sexdif = false;
    readfemale = false;
  }
  fscanf(simdata, "%d%lg%lg%*[^\n]", &which, &inc, &finish);
  getc(simdata);
  finish += 0.0001;
  if (!sexlink) {
    if (mutsys == 0)
      thispath = auto_;
    else
      thispath = mauto;
  } else if (mutsys == 0)
    thispath = sex;
  else
    thispath = msex;
  segscale = scale;
  FORLIM = nlocus;
  for (V.i = 1; V.i <= FORLIM; V.i++) ;
}  /* readloci */


void inputdata()
{
  struct LOC_inputdata V;
  int i;
  int maxonelocusgeno;

  printf("Reading SIMDATA.DAT\n");
  if (simdata != NULL)
    simdata = freopen("simdata.dat", "r", simdata);
  else
    simdata = fopen("simdata.dat", "r");
  if (simdata == NULL)
    exit(FileNotFound);
  /*Changed for SUN Pascal*/
  hapstore = (int **) malloc(maxlocus * sizeof(int *));
  maxonelocusgeno = (maxall * (maxall +1)/2);
  for (i = 0 ; i < maxlocus; i++)
    hapstore[i] = (int *) malloc(maxonelocusgeno * sizeof(int));
  readloci(&V);
  printf("Reading SIMPED.DAT\n");
  if (simped != NULL)
    simped = freopen("simped.dat", "r", simped);
  else
    simped = fopen("simped.dat", "r");
  if (simped == NULL)
    exit(FileNotFound);
  /*SUN*/
  readped(&V);
}  /* inputdata */



