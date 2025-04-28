#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/time.h>

/* Output from p2c, the Pascal-to-C translator */
/* From input file "msim.p" */


/* p2c: msim.p, line 1: 
 * Note: Unexpected name "outfile" in program header [262] */
/* p2c: msim.p, line 1: 
 * Note: Unexpected name "ipedfile" in program header [262] */
/* p2c: msim.p, line 1: 
 * Note: Unexpected name "datafile" in program header [262] */
/* p2c: msim.p, line 1: 
 * Note: Unexpected name "speedfile" in program header [262] */
/* p2c: msim.p, line 1: 
 * Note: Unexpected name "msim" in program header [262] */
/* p2c: msim.p, line 1: 
 * Note: Unexpected name "lodfile" in program header [262] */
/* p2c: msim.p, line 1: 
 * Note: Unexpected name "simout" in program header [262] */


/*3 July 1992*/

/*MSIM is a modified version of MLINK, modified by Daniel E. Weeks with
  the help of Mark Lathrop during June 1989.  It is a companion program
  to the generalized simulation program SLINK.
  Partial Update History:
   6/21/90 - Fixed msegsex and msegsexf according to Lodewijk.
             Added logical variable lodprint.
   8/90    - Fixed to be consistent with Version 5.1.
             Added checking of size of maxfact.
  11/90    - Fixed error with clods (percentage).
   9/90    - Now reads number of pedigrees/replicate from the file 'simout.dat'
  11/91    - QUADMAX function replaced (J. Ott).

  INPUT:
  The three typical MLINK-type datafiles, where ipedfile and speedfile have
  been produced by UNKNOWN.
  datafile.dat
  ipedfile.dat
  speedfile.dat
  simout.dat

  OUTPUT:
  msim.dat      => Summary statistics
  outfile.dat
  lodfile.dat   => Listing of lods by male theta (only made if lodprint=TRUE)

  NOTES:

  1) The stream file outputs have been removed.

  2) This program has been modified to process a replicate of the
  original pedigrees at a time.  This should make it possible
  to analyze a large number of replicates even on a computer
  without much memory.

  3) The statistical summary is placed in the file named 'msim.dat',
  while the usual output of mlink is written to outfile.

  4) A list of lod scores by pedigree is kept in the fie 'lodfile.dat'
  However, right now it just writes out the male thetas.
  */

/*include <p2c/p2c.h>*/

#define true            1
#define false           0
#define FileNotFound   10  /*From p2c.h*/

#ifndef NULL
#define NULL            0
#endif

#ifndef EXIT_SUCCESS
#  ifdef vms
#define EXIT_SUCCESS 1
#define EXIT_FAILURE (02000000000L)
#else
#define EXIT_SUCCESS 0
#define EXIT_FAILURE 1
#  endif     /*From p2c.h*/
#endif
#ifndef Void
#define Void void
#endif

#ifndef Static
#define Static static
#endif

#ifndef Local
#define Local static
#endif

#ifndef Char
#define Char char
#endif

#define Free free
#define Malloc malloc

#define version         "2.51"
    /*VERSION OF MSIM, based on version 4.91 of LINKAGE*/

#define print           false
    /*print = FALSE turns off most output for faster execution*/
#define lodprint        true
    /*FALSE turns off output to the lodfile for faster execution*/

#define maxpnt          27
    /*Max. number of points at which the lik. can be evaluated*/
/* SOME USER DEFINED CONSTANTS */
/*THE PROGRAM WILL TELL YOU IF THE FOLLOWING TWO CONSTANTS CAN BE REDUCED*/
/*IF THE PROGRAM TERMINATES WITH AN ERROR IN RECOMBINATION INCREASE MAXNEED*/
#define maxneed         32
    /*MAXIMUM NUMBER OF RECOMBINATION PROBABILITIES*/
/*THE FOLLOWING SHOULD BE LARGER THAN MININT*/
#define maxcensor       10000   /*MAXIMUM FOR CENSORING ARRAY*/
#define maxlocus        4   /*MAXIMUM NUMBER OF LOCI */
#define maxseg          8   /*2 TO THE POWER maxlocus-1 */
#define maxall          8   /*MAX NUMBER OF ALLELES AT A SINGLE LOCUS*/
#define maxhap          16   /*MAXIMUM NUMBER OF HAPLOTYPES*/
/* = n1 x n2 x ... ETC. WHERE ni = NUMBER OF ALLELES LOCUS I*/

#define maxfem          (maxhap * (maxhap + 1) / 2)
    /*MAX. NO. OF JOINT GENOTYPES FOR A FEMALE*/
#define maxmal          maxfem
    /*MAXIMUM NUMBER OF JOINT GENOTYPES FOR A MALE*/
/* = maxfem (AUTOSOMAL) OR maxhap (SEXLINKED)*/

#define maxind          600   /*MAXIMUM NUMBER OF INDIVIDUALS*/
#define maxped          65   /*MAXIMUM NUMBER OF PEDIGREES*/
#define maxchild        50   /*MAXIMUM NUMBER OF FULLSIBS IN A SIBSHIP*/
#define maxloop         3   /*MAXIMUM NUMBER OF LOOPS PER PEDIGREE*/

#define minfreq         0.0

#define affall          2
    /*DISEASE ALLELE FOR QUANTITATIVE TRAITS OR AFFECTION STATUS*/
/* QUANTITATIVE TRAIT */
#define maxtrait        1
    /*MAXIMUM NUMBER OF QUANTITATIVE FACTORS AT A SINGLE LOCUS*/

#define missval         0.0   /*MISSING VALUES FOR QUANTITATIVE TRAITS */
/* AFFECTION STATUS */

#define missaff         0   /*MISSING VALUE FOR AFFECTION STATUS */
#define affval          2   /*CODE FOR AFFECTED INDIVIDUAL*/
#define maxliab         400   /*MAXIMUM NUMBER OF LIABILITY CLASSES */
/* BINARY (FACTOR UNION) SYSTEM */
#define maxfact         maxall
    /*MAXIMUM NUMBER OF BINARY CODES AT A SINGLE LOCUS*/
/* OTHERS */

#define scale           1.0   /*SCALE FACTOR*/
#define scalemult       1.0   /*SCALE WEIGHT FOR EACH LOCUS INCLUDED*/

#define fitmodel        false
    /*TRUE IF ESTIMATING PARAMETERS OTHER THAN REC*/
#define score           true   /*LEAVE score TRUE FOR msim TO WORK CORRECTLY*/
/*CALCULATE LOD SCORES*/
#define byfamily        true
    /*LEAVE byfamily TRUE FOR msim TO WORK CORRECTLY*/
/*GIVE LOD SCORES BY FAMILY*/

#define zerolike        (-1.0e20)
/*FOR INCONSISTENT DATA OR RECOMBINATION */
/*SUN*/
#define log10_          2.30259

#define minint          (-32767)   /*MINIMUM ALLOWED INTEGER*/
/*GRADIENT APPROXIMATIONS*/

#define approximate     false
    /*do not change or else change approxrec below*/

#define epsilon         0.0

/*next two typedefs come from p2c*/
typedef char boolean;
typedef short uchar;

typedef struct gennurec {
  long genenumber[maxhap][maxhap];
} gennurec;

typedef struct approxrec {
  boolean approxarray[1][1];
} approxrec;

typedef struct censorrec {
  boolean censor[maxcensor - minint + 1];
} censorrec;

typedef double genotype[maxfem];
typedef double covmatrix[maxtrait][maxtrait];
typedef char hapvector[maxlocus];
typedef long mutarray[3][2];
typedef double thesemeans[maxtrait];
typedef thesemeans means[maxall + 1][maxall];
typedef enum {
  auto_, mauto, sex, msex
} pathway;
typedef enum {
  peelup, peeldown
} direction;

typedef long binset;

typedef binset phenarray[maxall];

typedef enum {
  affection, quantitative, binary_, null_
} locustype;

typedef struct locusvalues {
  long nallele, format;
  double freq[maxall];
  struct locusvalues *privlocus;
  locustype which;
  union {
    struct {
      double pen[maxall + 1][maxall][3][maxliab];
      long nclass;
    } U0;
    struct {
      long ntrait;
      means pm;
      covmatrix vmat;
      double det, contrait, conmat;
    } U1;
    phenarray allele;
  } UU;
} locusvalues;

typedef struct phenotype {
  locustype which;
  union {
    binset phenf;
    struct {
      double x[maxtrait];
      boolean missing;
    } U1;
    struct {
      long aff, liability;
    } U0;
  } UU;
} phenotype;

typedef phenotype *indphen[maxlocus];

typedef struct thisarray {
  genotype genarray;
} thisarray;

typedef boolean possvect[maxall][maxall];
typedef possvect possarray[maxlocus];

typedef struct information {
  possarray possible;
} information;

typedef struct thisperson {
  long id, ped, inloop;
  struct thisperson *pa, *ma, *foff, *nextpa, *nextma;
  thisarray *gen;
  indphen phen;
  phenotype *privphen;
  boolean unknown, multi, done, up, male, firstpass;
  information *store;
} thisperson;

typedef long subhap[maxseg];
typedef double thetarray[maxlocus];
typedef double happrob[maxneed];

typedef struct thetavalues {
  thetarray theta;
  happrob segprob;
} thetavalues;


Static double zlod[maxped + 1][maxpnt];
Static double ztheta[maxpnt];
Static long jmax;
Static double rep, zmax;
Static long clods[maxped + 1][3];   /*counts of lods >= constants 1,2,3*/
Static double lodlimit[3];
Static long nrep, ip, npt, j;
    /*nrep is a counter of number of replications*/
Static double aver[maxped + 1][maxpnt], variance[maxped + 1][maxpnt],
	      big[maxped + 1][maxpnt], small[maxped + 1][maxpnt];
Static double avermaxlod[maxped + 1], variancemaxlod[maxped + 1],
	      bigmaxlod[maxped + 1], smallmaxlod[maxped + 1];
/*These are arrays where various statistics are stored*/
Static long nfactor[maxlocus];
Static long opeds;   /*Number of original pedigrees*/
Static boolean infile;
    /*TRUE if currently in the middle of the pedigree file*/
Static long oldped;   /*oldped stores number of next pedigree to be read.*/
Static double lods[maxped], stand[maxped];
Static boolean zeromale[maxlocus], zerofemale[maxlocus];
Static pathway thispath;
Static boolean informative[maxped];
Static boolean rare[maxfem], risk1[maxfem], risk2[maxfem];
Static boolean riskmale[maxhap];
Static approxrec *approxstruct;
Static censorrec *censorstruct;
Static long thisc;
Static boolean malechild[maxchild];
Static thisarray *thischild[maxchild];
Static long nchild;
Static char seghap1[maxfem], seghap2[maxfem];
Static short segstart[maxfem];
Static short probstart[maxfem], probend[maxfem];
Static boolean nohom[maxlocus];
Static locusvalues *thislocus[maxlocus];
Static long increment[maxlocus], order[maxlocus];
/*PEOPLE*/
Static thisperson *person[maxind + 1];
Static thisperson *proband[maxped];
Static thisperson *looppers[maxped][maxloop][2];
/*MUTATION */
Static char muthap[maxhap];
Static gennurec *gennustruct;
/*RECOMBINATION*/
Static thetavalues *maletheta, *femaletheta;
/*FREQUENCIES */
Static thisarray *hapfreq;
/*RISK*/
Static char riskall;
/*OTHERS*/
Static long i, risksys, mutsys, nlocus, which, lastpriv, thissystem,
	    lastspeed, lastseg, segperson;
/*variables used for reading speedfile: lastspeed,lastseg,segperson*/
Static char nuhap;
Static short fgeno, mgeno;
Static long nuped, totperson;
Static double segscale, mutmale, mutfemale, like, alike, tlike, finish, inc,
	      distratio, scorevalue, holdtheta;
Static boolean interfer, disequi, sexlink, risk, sexdif, readfemale, mapping,
	       dolod, firstapprox, lasttime, firsttime, firsteff;
Static FILE *outfile, *ipedfile, *datafile, *speedfile, *msim, *lodfile,
	    *simout;
Static char chtemp;
Static long k;

/* Two routines taken from */
/* "p2c"  Copyright (C) 1989, 1990, 1991 Free Software Foundation.
 * By Dave Gillespie, daveg@csvax.cs.caltech.edu.  Version --VERSION--.
 * This file may be copied, modified, etc. in any way.  It is not restricted
 * by the licence agreement accompanying p2c itself.
 */


/* Check if at end of file, using Pascal "eof" semantics.  End-of-file for
   stdin is broken; remove the special case for it to be broken in a
   different way. */
#include <stdio.h>

int P_eof(f) 
FILE *f;
{
    register int ch;

    if (feof(f))
        return 1;
    if (f == stdin)
        return 0;    /* not safe to look-ahead on the keyboard! */
    ch = getc(f);
    if (ch == EOF)
        return 1;
    ungetc(ch, f);
    return 0;
}


/* Check if at end of line (or end of entire file). */

int P_eoln(f)
FILE *f;
{
    register int ch;

    ch = getc(f);
    if (ch == EOF)
        return 1;
    ungetc(ch, f);
    return (ch == '\n');
}


Static double quadmax(x1, x2, x3, f1, f2, f3)
double x1, x2, x3, f1, f2, f3;
{
  /*Method based on formulas (4.9) - (4.14) in Ott (1985) Analysis of Human
    Genetic Linkage, first edition.  Johns Hopkins Univ. Press, Baltimore.
    Computes quadratic curve for function values fi over arguments xi.  If
    one of the flanking fi is largest that value is taken to be quadmax.
    If interval lengths differ by more than a factor of 4, quadmax is
    simply the max. of the 3 fi values*/
  double d, d1, d2, ex, TEMP;

  if (x1 == x2) {   /*worry about x values*/
    if (f2 > f3)
      return f2;
    else
      return f3;
  } else if (x2 == x3) {
    if (f1 > f2)
      return f1;
    else
      return f2;
  } else {
    d = (x2 - x1) / (x3 - x2);
    if (d > 4.0 || d < 0.25 || f1 < -100.0) {
	  /*then simply look at max(fi)*/
	    if (f1 >= f2 && f1 >= f3)
	return f1;
      else if (f3 >= f2 && f2 >= f1)
	return f3;
      else
	return f2;
    } else if (f1 == f2 && f2 == f3)
      return f1;
    else if (f2 >= f1 && f2 >= f3) {
      d = x3 - x1;
      d1 = (f3 - f1) / d;
      d2 = 0.5 * d / ((f2 - f1) / (x2 - x1) - (f3 - f2) / (x3 - x2));
      ex = 0.5 * (x1 + x3) + d1 * d2;
      TEMP = x2 - ex;
      return (f2 + 0.5 * (TEMP * TEMP) / d2);
    } else if (f1 > f3)
      return f1;
    else
      return f3;
  }

  /*verify f values*/
}


/* quadmax */
/*quadmax*/


Static double min(a, b)
double a, b;
{
  if (a < b)
    return a;
  else
    return b;
}  /* min */


Static double max(a, b)
double a, b;
{
  if (a > b)
    return a;
  else
    return b;
}  /* max */


Static double mapfunction(theta1, theta2)
double theta1, theta2;
{
  /*User defined function giving recombination between
flanking markers as a function of recombination
between adjacent markers*/
  return ((theta1 + theta2) / (1 + 4 * theta1 * theta2));
}  /* mapfunction */


Static double getdist(theta)
double *theta;
{
  if (*theta < 0.5)
    return (log(1.0 - 2.0 * *theta) / -2.0);
  else
    return 10.0;
}  /* getdist */


Static double invdist(dist)
double *dist;
{
  if (*dist < 10.0)
    return ((1.0 - exp(-2.0 * *dist)) / 2.0);
  else
    return 0.5;
}  /* invdist */


Static Void invert(m, n, det)
double (*m)[maxtrait];
long n;
double *det;
{
  covmatrix v;
  double val;
  long i, j, k;

  *det = 1.0;
  for (i = 1; i <= n; i++) {
    val = m[i - 1][i - 1];
    for (k = 0; k <= i - 2; k++)
      val -= v[k][i - 1] * v[k][i - 1];
    *det *= val;
    v[i - 1][i - 1] = sqrt(val);
    for (j = i; j < n; j++) {
      val = m[i - 1][j];
      for (k = 0; k <= i - 2; k++)
	val -= v[k][i - 1] * v[k][j];
      v[i - 1][j] = val / v[i - 1][i - 1];
      v[j][i - 1] = 0.0;
    }
  }
  for (i = 1; i <= n; i++) {
    m[i - 1][i - 1] = 1 / v[i - 1][i - 1];
    for (j = i + 1; j <= n; j++) {
      val = 0.0;
      for (k = i - 1; k <= j - 2; k++)
	val -= m[k][i - 1] * v[k][j - 1];
      m[j - 1][i - 1] = val / v[j - 1][j - 1];
      m[i - 1][j - 1] = 0.0;
    }
  }
  for (i = 0; i < n; i++) {
    for (j = 0; j < i + 1; j++) {
      val = 0.0;
      for (k = j; k < n; k++)
	val += m[k][i] * m[k][j];
      v[i][j] = val;
      v[j][i] = val;
    }
  }
  memcpy(m, v, sizeof(covmatrix));
}  /* invert */


/* Local variables for recombination: */
struct LOC_recombination {
  long here;
  double p1, p2, p3, p4;
} ;

/* Local variables for recombine: */
struct LOC_recombine {
  struct LOC_recombination *LINK;
  double *theta;
  double *segprob;
  long there, nhap;
  boolean thishet[maxlocus];
  hapvector hap1, hap2;
  hapvector thishap1[maxseg], thishap2[maxseg];
} ;

Local Void scramble(LINK)
struct LOC_recombine *LINK;
{
  long whichhap, start, length, i, j, k;
  double recval, val;
  long FORLIM2;

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


Local Void setrec(val, LINK)
double val;
struct LOC_recombine *LINK;
{
  LINK->nhap++;
  memcpy(LINK->thishap1[LINK->nhap - 1], LINK->hap1, sizeof(hapvector));
  memcpy(LINK->thishap2[LINK->nhap - 1], LINK->hap2, sizeof(hapvector));
  LINK->there++;
  LINK->segprob[LINK->there - 1] = val;
}  /* setrec */


Local Void dointer(LINK)
struct LOC_recombine *LINK;
{
  long i;
  boolean temphet[3];
  long FORLIM;

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


Local Void nexthet(i, val, inphase, LINK)
long i;
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


Local Void getrecprob(LINK)
struct LOC_recombine *LINK;
{
  long i;

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


Local Void gethet(system, LINK)
long *system;
struct LOC_recombine *LINK;
{
  long newsystem;

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

Local Void recombine(theta_, segprob_, LINK)
double *theta_;
double *segprob_;
struct LOC_recombination *LINK;
{
  struct LOC_recombine V;
  long system;


  V.LINK = LINK;
  V.theta = theta_;
  V.segprob = segprob_;
  LINK->here = 0;
  system = 1;
  gethet(&system, &V);
}  /* recombine */


Local Void getfemaletheta(LINK)
struct LOC_recombination *LINK;
{
  double dist;
  long ntheta, i;

  if (interfer)
    ntheta = nlocus;
  else
    ntheta = nlocus - 1;
  for (i = 0; i < ntheta; i++) {
    dist = getdist(&maletheta->theta[i]) * distratio;
    femaletheta->theta[i] = invdist(&dist);
  }
}  /* getfemaletheta */


Static Void recombination()
{
  struct LOC_recombination V;
  long i;
  thetarray oldtheta;
  thetavalues *WITH;
  long FORLIM;


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
  if (!print)
    return;
  if (firsteff) {
    if (V.here < maxneed)
      printf(" Maxneed may be reduced to %12ld\n", V.here);
  }
}  /* recombination */


/* Local variables for getlocations: */
struct LOC_getlocations {
  long ngene, nseg, here, there, start, nhet, thisseg;
  boolean rarepresent, riskhom, riskhet;
  hapvector hap1, hap2;
  boolean thishet[maxlocus];
} ;

Local boolean checkrare(LINK)
struct LOC_getlocations *LINK;
{
  long i;
  boolean check;
  long FORLIM;
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


Local Void checkrisk(riskhet, riskhom, LINK)
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


Local long gethapn(hap, LINK)
char *hap;
struct LOC_getlocations *LINK;
{
  long i, n, FORLIM;

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

Local Void setrisk(LINK)
struct LOC_domalerisk *LINK;
{
  long n;

  n = gethapn(LINK->LINK->hap1, LINK->LINK);
  if (LINK->LINK->hap1[risksys - 1] == riskall)
    riskmale[n - 1] = true;
  else
    riskmale[n - 1] = false;
}  /* setrisk */


Local Void getriskhap(system, LINK)
long system;
struct LOC_domalerisk *LINK;
{
  long i;
  locusvalues *WITH;
  long FORLIM;

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


Local Void domalerisk(LINK)
struct LOC_getlocations *LINK;
{
  struct LOC_domalerisk V;


  V.LINK = LINK;
  getriskhap(1L, &V);
}  /* domalerisk */

/* Local variables for domutation: */
struct LOC_domutation {
  struct LOC_getlocations *LINK;
} ;

Local Void setmutation(LINK)
struct LOC_domutation *LINK;
{
  long i, n;

  n = gethapn(LINK->LINK->hap1, LINK->LINK);
  if (LINK->LINK->hap1[mutsys - 1] == thislocus[mutsys - 1]->nallele) {
    muthap[n - 1] = n;
    return;
  }
  i = LINK->LINK->hap1[mutsys - 1];
  LINK->LINK->hap1[mutsys - 1] = thislocus[mutsys - 1]->nallele;
  muthap[n - 1] = gethapn(LINK->LINK->hap1, LINK->LINK);
  LINK->LINK->hap1[mutsys - 1] = i;
}  /* setmutation */


Local Void getmuthap(system, LINK)
long system;
struct LOC_domutation *LINK;
{
  long i;
  locusvalues *WITH;
  long FORLIM;

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


Local Void domutation(LINK)
struct LOC_getlocations *LINK;
{
  struct LOC_domutation V;


  V.LINK = LINK;
  getmuthap(1L, &V);
}  /* domutation */


Local Void setnumbers(LINK)
struct LOC_getlocations *LINK;
{
  long nhap1, nhap2;

  LINK->ngene++;

  segstart[LINK->ngene - 1] = LINK->here + 1;
  probstart[LINK->ngene - 1] = LINK->there + 1;
  probend[LINK->ngene - 1] = LINK->there + LINK->nseg;

  LINK->there += LINK->nseg;

  nhap1 = gethapn(LINK->hap1, LINK);
  nhap2 = gethapn(LINK->hap2, LINK);
  gennustruct->genenumber[nhap1 - 1][nhap2 - 1] = LINK->ngene;
  gennustruct->genenumber[nhap2 - 1][nhap1 - 1] = LINK->ngene;

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
  seghap1[LINK->thisseg - 1] = nhap1;
  seghap2[LINK->thisseg - 1] = nhap2;
}  /* setnumbers */


Local Void hapscr(system, nscramble, LINK)
long system, nscramble;
struct LOC_getlocations *LINK;
{
  long i, j;

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


Local Void sethap(system, LINK)
long system;
struct LOC_getlocations *LINK;
{
  long i, j;
  locusvalues *WITH;
  long FORLIM, FORLIM1;

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
	  hapscr(1L, 0L, LINK);
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
      hapscr(1L, 0L, LINK);
      LINK->here += LINK->nseg;
    }
  }
}  /* sethap */


Local Void starthap(LINK)
struct LOC_getlocations *LINK;
{
  long i, FORLIM;

  LINK->nseg = 1;
  FORLIM = LINK->nhet;
  for (i = 2; i <= FORLIM; i++)
    LINK->nseg *= 2;
  sethap(1L, LINK);
  LINK->start = LINK->there;
}  /* starthap */


Local Void gethet1(system, LINK)
long system;
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


Static Void getlocations()
{
  struct LOC_getlocations V;


  V.nhet = 0;
  V.here = 0;
  V.there = 0;
  V.ngene = 0;
  V.start = 0;
  gethet1(1L, &V);
  if (mutsys != 0)
    domutation(&V);
  if (sexlink && risk)
    domalerisk(&V);
}  /* getlocations */


Local Void inputerror(nerror, par1, par2)
long nerror, par1, par2;
{
  printf("Fatal error detected in procedure inputdata\n");
  switch (nerror) {

  case 0:
    printf("Number of loci %2ld exceeds the constant maxlocus\n", par1);
    break;

  case 1:
    printf("Number of loci read %2ld. Less than minimum of 1\n", par1);
    break;

  case 2:
    printf(
      "Error detected reading loci order. Locus number %2ld in position %2ld exceeds number of loci\n",
      par2, par1);
    break;

  case 3:
    printf(
      "Error detected reading loci order. Illegal locus number %2ld in position %2ld\n",
      par2, par1);
    break;

  case 4:
    printf(
      "Error detected reading loci order. Locus number repeated in positions %2ld and %2ld\n",
      par1, par2);
    break;

  case 5:
    printf(
      "Error detected reading locus description. Illegal locus type %2ld for locus %2ld\n",
      par2, par1);
    break;

  case 6:
    printf(
      "Error detected reading locus description for system %2ld. Number of alleles  %2ld exceeds maxall\n",
      par1, par1);
    break;

  case 7:
    printf(
      "Error detected reading locus description for system %2ld. Illegal number of alleles  %2ld\n",
      par1, par2);
    break;

  case 8:
    printf(
      "Error detected reading locus description for system %2ld. Number of factors  %2ld exceeds maxfact\n",
      par1, par2);
    break;

  case 9:
    printf(
      "Error detected reading locus description for system %2ld. Illegal number of factors  %2ld\n",
      par1, par2);
    break;

  case 10:
    printf(
      "Error detected reading locus description for system %2ld. Alleles not codominant\n",
      par1);
    break;

  case 11:
    printf("Error detected reading pedigree record %2ld. Illegal code for sex %2ld\n",
	   par1, par2);
    break;

  case 12:
    printf(
      "Error detected reading pedigree record at pedigree%2ld. Maximum number of pedigree records exceeded\n",
      par1);
    break;

  case 13:
    printf(
      "Error detected reading pedigree record %2ld. Maximum number of individuals exceeded\n",
      par1);
    break;

  case 14:
    printf(
      "Error detected reading pedigree record %2ld. Illegal binary factor code %2ld\n",
      par1, par2);
    break;

  case 15:
    printf(
      "Error detected reading pedigree record %2ld. No allelic pair for genotype\n",
      par1);
    break;

  case 16:
    printf(
      "Error detected reading pedigree record %2ld. Allele number %2ld exceeds maxall\n",
      par1, par2);
    break;

  case 17:
    printf(
      "Error detected reading pedigree record %2ld. Illegal allele number %2ld\n",
      par1, par2);
    break;

  case 18:
    printf("Number of systems after factorization (%3ld) exceeds maxsystem\n",
	   par1);
    break;

  case 19:
    printf("Number of systems after factorization (%3ld) less than minimum of 1\n",
	   par1);
    break;

  case 20:
    printf("Number of recombination types (%3ld) exceeds maxrectype\n", par1);
    break;

  case 21:
    printf("Number of recombination types (%3ld) less than minimum of 1\n",
	   par1);
    break;

  case 22:
    printf(
      "End of file detected in tempdat by procedure readthg before all data found\n");
    break;

  case 23:
    printf(
      "Error detected reading iterated locus in datafile. Value (%3ld) greater than nlocus\n",
      par1);
    break;

  case 24:
    printf(
      "Error detected reading iterated locus in datafile. Illegal value (%3ld)\n",
      par1);
    break;

  case 25:
    printf("Number of iterated parameters greater then maxn\n");
    break;

  case 26:
    printf(
      "Error detected reading pedigree record %2ld. Liability class (%2ld) exceeds nclass\n",
      par1, par2);
    break;

  case 27:
    printf(
      "Error detected reading pedigree record %2ld. Illegal liability class (%2ld)\n",
      par1, par2);
    break;

  case 28:
    printf(
      "Error detected reading locus description for system%2ld. Liability classes (%3ld) exceed maxliab\n",
      par1, par2);
    break;

  case 29:
    printf(
      "Error detected reading locus description for system%2ld. Illegal number of liability classes (%3ld)\n",
      par1, par2);
    break;

  case 30:
    printf(
      "Error detected reading locus description for system%2ld. Penetrance out of range\n",
      par1);
    break;

  case 31:
    printf(
      "Error detected reading locus description for system%2ld. Number of traits (%3ld) exceeds maxtrait\n",
      par1, par2);
    break;

  case 32:
    printf(
      "Error detected reading locus description for system%2ld. Number of traits out of range (%3ld)\n",
      par1, par2);
    break;

  case 33:
    printf(
      "Error detected reading locus description for system%2ld. Variance must be positive\n",
      par1);
    break;

  case 34:
    printf(
      "Error detected reading locus description for system%2ld. Variance multiplier must be positive\n",
      par1);
    break;

  case 35:
    printf(
      "Error detected reading locus description for system%2ld. Risk allele %3ld) exceeds nallele\n",
      par1, par2);
    break;

  case 36:
    printf(
      "Error detected reading locus description for system%2ld. Illegal risk allele (%3ld)\n",
      par1, par2);
    break;

  case 37:
    printf("Error detected reading datafile. Risk locus %3ld) exceeds nlocus\n",
	   par2);
    break;

  case 38:
    printf("Error detected reading datafile. Illegal value for risk locus %3ld)\n",
	   par2);
    break;

  case 39:
    printf("Error detected reading datafile. Mutation locus %3ld) exceeds nlocus\n",
	   par2);
    break;

  case 40:
    printf(
      "Error detected reading datafile. Illegal value for mutation locus %3ld)\n",
      par2);
    break;

  case 41:
    printf("Error in datafile, line 1, item 4: Program code not for MSIM/MLINK\n");
    break;
  }
  printf("Press return to halt...\007\n");
  scanf("%*[^\n]");
  getchar();
  exit(0);   /*SUN*/
}  /* inputerror */

/*inputerror*/

Local Void inputwarning(nwarning, par1, par2)
long nwarning, par1, par2;
{
  printf("Warning number from procedure inputdata\n");
  switch (nwarning) {

  case 0:
    printf("Illegal sex difference parameter %2ld Parameter should be 0, 1, or 2\n",
	   par1);
    break;

  case 1:
    printf("Illegal interference parameter %2ld Lack of interference assumed\n",
	   par1);
    break;

  case 2:
    printf(
      "Illegal sex difference parameter %2ld Parameter must be 0 with sex-linked data\n",
      par1);
    break;
  }
  printf("Press return to continue...\n");
  scanf("%*[^\n]");
  getchar();
}  /* inputwarning */

/* Local variables for readloci: */
struct LOC_readloci {
  long i, whichtype, nupriv;
} ;

/* Local variables for getlocus: */
struct LOC_getlocus {
  struct LOC_readloci *LINK;
  long system;
} ;

Local Void getbin(locus, system, LINK)
locusvalues **locus;
long *system;
struct LOC_getlocus *LINK;
{
  long i, j, k;
  locusvalues *WITH;
  long FORLIM, FORLIM1;

  fscanf(datafile, "%ld%*[^\n]", &nfactor[*system - 1]);
  getc(datafile);
  if (nfactor[*system - 1] > maxfact)
    inputerror(8L, *system, nfactor[*system - 1]);
  if (nfactor[*system - 1] <= 0)
    inputerror(9L, *system, nfactor[*system - 1]);
  WITH = *locus;
  FORLIM = WITH->nallele;
  for (i = 0; i < FORLIM; i++)
    WITH->UU.allele[i] = 0;
  FORLIM = WITH->nallele;
  for (i = 0; i < FORLIM; i++) {
    FORLIM1 = nfactor[*system - 1];
    for (j = 1; j <= FORLIM1; j++) {
      fscanf(datafile, "%ld", &k);
      if (k == 1)
	WITH->UU.allele[i] = ((long)WITH->UU.allele[i]) | (1L << ((int)j));
    }
  }
  fscanf(datafile, "%*[^\n]");
  getc(datafile);
}  /* getbin */


Local Void getnumber(locus, system, LINK)
locusvalues **locus;
long *system;
struct LOC_getlocus *LINK;
{
  long i;
  locusvalues *WITH;
  long FORLIM;

  WITH = *locus;
  FORLIM = WITH->nallele;
  for (i = 1; i <= FORLIM; i++) {
    if (i > maxfact) {
      fprintf(msim, " *** ERROR: maxfact not big enough ***\n");
      fprintf(outfile, " *** ERROR: maxfact not big enough ***\n");
      printf(" *** ERROR: maxfact not big enough ***\n");
      printf(" Press return to halt...\n");
      scanf("%*[^\n]");
      getchar();
      exit(0);   /*SUN*/
    }
    WITH->UU.allele[i - 1] = 1L << ((int)i);
  }
}  /* getnumber */


Local Void getpen(locus, LINK)
locusvalues **locus;
struct LOC_getlocus *LINK;
{
  long i, j, k, l;
  locusvalues *WITH;
  long FORLIM, FORLIM1, FORLIM2;

  WITH = *locus;
  fscanf(datafile, "%ld%*[^\n]", &WITH->UU.U0.nclass);
  getc(datafile);
  if (WITH->UU.U0.nclass > maxliab)
    inputerror(28L, LINK->system, WITH->UU.U0.nclass);
  if (WITH->UU.U0.nclass <= 0)
    inputerror(29L, LINK->system, WITH->UU.U0.nclass);
  FORLIM = WITH->UU.U0.nclass;
  for (l = 0; l < FORLIM; l++) {
    FORLIM1 = WITH->nallele;
    for (i = 1; i <= FORLIM1; i++) {
      FORLIM2 = WITH->nallele;
      for (j = i - 1; j < FORLIM2; j++) {
	fscanf(datafile, "%lg", &WITH->UU.U0.pen[i][j][2][l]);
	if ((unsigned)WITH->UU.U0.pen[i][j][2][l] > 1.0)
	  inputerror(30L, LINK->system, LINK->system);
	WITH->UU.U0.pen[i][j][1][l] = 1 - WITH->UU.U0.pen[i][j][2][l];
	WITH->UU.U0.pen[i][j][0][l] = 1.0;
	for (k = 0; k <= 2; k++)
	  WITH->UU.U0.pen[j + 1][i - 1][k][l] = WITH->UU.U0.pen[i][j][k][l];
      }
    }
    fscanf(datafile, "%*[^\n]");
    getc(datafile);
    FORLIM1 = WITH->nallele;
    for (i = 0; i < FORLIM1; i++)
      WITH->UU.U0.pen[0][i][0][l] = 1.0;
    if (sexlink) {
      FORLIM1 = WITH->nallele;
      for (i = 0; i < FORLIM1; i++) {
	fscanf(datafile, "%lg", &WITH->UU.U0.pen[0][i][2][l]);
	if ((unsigned)WITH->UU.U0.pen[0][i][2][l] > 1.0)
	  inputerror(30L, LINK->system, LINK->system);
      }
      FORLIM1 = WITH->nallele;
      for (i = 0; i < FORLIM1; i++)
	WITH->UU.U0.pen[0][i][1][l] = 1.0 - WITH->UU.U0.pen[0][i][2][l];
      fscanf(datafile, "%*[^\n]");
      getc(datafile);
    }
  }
}  /* getpen */


Local Void getquan(locus, privelege, LINK)
locusvalues **locus;
boolean privelege;
struct LOC_getlocus *LINK;
{
  /*Get information of a quantitative locus. Privelege says whether it is
  a priveleged locus or not*/
  long i, j, k;
  locusvalues *WITH;
  long FORLIM, FORLIM1, FORLIM2;

  WITH = *locus;
  fscanf(datafile, "%ld%*[^\n]", &WITH->UU.U1.ntrait);
  getc(datafile);
  if (WITH->UU.U1.ntrait > maxtrait)
    inputerror(31L, LINK->system, WITH->UU.U1.ntrait);
  if (WITH->UU.U1.ntrait <= 0)
    inputerror(32L, LINK->system, WITH->UU.U0.nclass);
  FORLIM = WITH->UU.U1.ntrait;
  for (i = 0; i < FORLIM; i++) {
    FORLIM1 = WITH->nallele;
    for (j = 1; j <= FORLIM1; j++) {
      FORLIM2 = WITH->nallele;
      for (k = j; k <= FORLIM2; k++) {
	fscanf(datafile, "%lg", &WITH->UU.U1.pm[j][k - 1][i]);
	WITH->UU.U1.pm[k][j - 1][i] = WITH->UU.U1.pm[j][k - 1][i];
      }
    }
  }
  fscanf(datafile, "%*[^\n]");
  getc(datafile);
  if (privelege && LINK->LINK->nupriv != lastpriv)
    return;
  FORLIM = WITH->UU.U1.ntrait;
  for (i = 0; i < FORLIM; i++) {
    FORLIM1 = WITH->UU.U1.ntrait;
    for (j = i; j < FORLIM1; j++) {
      fscanf(datafile, "%lg", &WITH->UU.U1.vmat[i][j]);
      if (i + 1 == j + 1 && WITH->UU.U1.vmat[i][j] <= 0.0)
	inputerror(33L, LINK->system, LINK->system);
      WITH->UU.U1.vmat[j][i] = WITH->UU.U1.vmat[i][j];
    }
  }
  fscanf(datafile, "%*[^\n]");
  getc(datafile);
  invert(WITH->UU.U1.vmat, WITH->UU.U1.ntrait, &WITH->UU.U1.det);
  WITH->UU.U1.det = 1 / sqrt(WITH->UU.U1.det);
  fscanf(datafile, "%lg%*[^\n]", &WITH->UU.U1.conmat);
  getc(datafile);
  if (WITH->UU.U1.conmat <= 0.0)
    inputerror(34L, LINK->system, LINK->system);
  WITH->UU.U1.conmat = 1 / WITH->UU.U1.conmat;
  WITH->UU.U1.contrait = 1.0;
  FORLIM = WITH->UU.U1.ntrait;
  for (i = 1; i <= FORLIM; i++)
    WITH->UU.U1.contrait *= WITH->UU.U1.conmat;
  WITH->UU.U1.contrait = sqrt(WITH->UU.U1.contrait);
}  /* getquan */

Local Void getlocus(system_, LINK)
long system_;
struct LOC_readloci *LINK;
{
  struct LOC_getlocus V;
  long j;
  locusvalues *WITH, *WITH1;
  long FORLIM;
  int TEMP;


  V.LINK = LINK;
  V.system = system_;
  thislocus[V.system - 1] = (locusvalues *)Malloc(sizeof(locusvalues));
  WITH = thislocus[V.system - 1];
  /*     thislocus[system]:=locuspoint(NewPtr(SizeOf(locusvalues)));*/
  WITH->privlocus = NULL;
  fscanf(datafile, "%ld%ld", &LINK->whichtype, &WITH->nallele);
  if (LINK->whichtype < 0 && LINK->whichtype > 4)
    inputerror(5L, V.system, LINK->whichtype);
  if (WITH->nallele > maxall)
    inputerror(6L, V.system, WITH->nallele);
  if (WITH->nallele <= 0)
    inputerror(7L, V.system, WITH->nallele);
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
    fscanf(datafile, "%*[^\n]");
    getc(datafile);
  } else {
    fscanf(datafile, "%ld%*[^\n]", &LINK->whichtype);
    getc(datafile);
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
      fscanf(datafile, "%lg", &WITH->freq[j]);
    fscanf(datafile, "%*[^\n]");
    getc(datafile);
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
    fscanf(datafile, "%d%*[^\n]", &TEMP);
    getc(datafile);
    riskall = TEMP;
  }
  if (riskall > thislocus[LINK->i - 1]->nallele)
    inputerror(35L, LINK->i, (long)riskall);
  if (riskall < 0)
    inputerror(36L, LINK->i, (long)riskall);
}  /* getlocus */


Local Void gettheta(sex_, LINK)
thetavalues **sex_;
struct LOC_readloci *LINK;
{
  thetarray oldtheta;
  long i;
  thetavalues *WITH;
  long FORLIM;

  *sex_ = (thetavalues *)Malloc(sizeof(thetavalues));
  WITH = *sex_;
  /*     sex:=thetapoint(NewPtr(SizeOf(thetavalues)));*/
  if (*sex_ == maletheta || readfemale) {
    FORLIM = nlocus - 2;
    for (i = 0; i <= FORLIM; i++)
      fscanf(datafile, "%lg", &WITH->theta[i]);
    if (interfer && !mapping)
      fscanf(datafile, "%lg", &WITH->theta[nlocus - 1]);
    fscanf(datafile, "%*[^\n]");
    getc(datafile);
  } else {
    fscanf(datafile, "%lg%*[^\n]", &distratio);
    getc(datafile);
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


Local Void readloci()
{
  struct LOC_readloci V;
  long j, coupling, autosomal, independent, difference, FORLIM, FORLIM1;
  locusvalues *WITH;


  lastpriv = 0;
  fscanf(datafile, "%ld%ld%ld%ld%*[^\n]", &nlocus, &risksys, &autosomal, &V.i);
  getc(datafile);   /*i=prog.code*/
  /*Replace the above line with the next when using epistasis*/
  /*readln(datafile,nlocus,risksys,autosomal,lastpriv);*/
  if (V.i != 5)
    inputerror(41L, nlocus, nlocus);
  if (nlocus > maxlocus)
    inputerror(0L, nlocus, nlocus);
  if (nlocus <= 0)
    inputerror(1L, nlocus, nlocus);
  if (risksys > maxlocus)
    inputerror(37L, risksys, risksys);
  if (risksys < 0)
    inputerror(38L, risksys, risksys);
  risk = (risksys != 0);
  sexlink = (autosomal == 1);
  printf("YOU ARE USING LINKAGE/MSIM (V%s) WITH%3ld-POINT", version, nlocus);
  if (sexlink)
    printf(" SEXLINKED DATA\n");
  else
    printf(" AUTOSOMAL DATA\n");
  fscanf(datafile, "%ld%lg%lg%ld%*[^\n]", &mutsys, &mutmale, &mutfemale,
	 &coupling);
  getc(datafile);
  if (mutsys > maxlocus)
    inputerror(39L, mutsys, mutsys);
  if (mutsys < 0)
    inputerror(40L, mutsys, mutsys);
  if (coupling == 1)
    disequi = true;
  else
    disequi = false;
  if (disequi)
    hapfreq = (thisarray *)Malloc(sizeof(thisarray));
  FORLIM = nlocus;
  /*     hapfreq:=genpoint(NewPtr(SizeOf(thisarray)));*/
  for (V.i = 1; V.i <= FORLIM; V.i++) {
    fscanf(datafile, "%ld", &j);
    if (j > nlocus)
      inputerror(2L, V.i, j);
    if (j <= 0)
      inputerror(3L, V.i, j);
    order[j - 1] = V.i;
  }
  FORLIM = nlocus;
  for (V.i = 1; V.i <= FORLIM; V.i++) {
    FORLIM1 = V.i;
    for (j = 1; j < FORLIM1; j++) {
      if (order[V.i - 1] == order[j - 1])
	inputerror(4L, V.i, j);
    }
  }
  fscanf(datafile, "%*[^\n]");
  getc(datafile);
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
    FORLIM = mgeno;
    for (V.i = 1; V.i <= FORLIM; V.i++)
      fscanf(datafile, "%lg", &hapfreq->genarray[V.i - 1]);
    fscanf(datafile, "%*[^\n]");
    getc(datafile);
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
  fscanf(datafile, "%ld", &difference);
  if ((unsigned long)difference > 2) {
    inputwarning(0L, difference, difference);
    difference = 0;
  }
  sexdif = (difference != 0);
  readfemale = (difference == 2);
  fscanf(datafile, "%ld%*[^\n]", &independent);
  getc(datafile);
  if ((unsigned long)independent > 2) {
    inputwarning(1L, independent, independent);
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
    inputwarning(2L, difference, difference);
    sexdif = false;
    readfemale = false;
  }
  fscanf(datafile, "%ld%lg%lg%*[^\n]", &which, &inc, &finish);
  getc(datafile);
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


Static Void inputdata()
{
  /*readloci*/

  printf("Reading IPEDFILE.DAT\n");
  if (ipedfile != NULL)
    ipedfile = freopen("ipedfile.dat", "r", ipedfile);
  else
    ipedfile = fopen("ipedfile.dat", "r");
  if (ipedfile == NULL)
    exit(FileNotFound);
  /*SUN*/
  printf("Reading SPEEDFILE.DAT\n");   /*SUN*/
  if (speedfile != NULL)
    speedfile = freopen("speedfile.dat", "r", speedfile);
  else
    speedfile = fopen("speedfile.dat", "r");
  if (speedfile == NULL)
    exit(FileNotFound);
  /*SUN*/
  printf("Reading DATAFILE.DAT\n");
  if (datafile != NULL)
    datafile = freopen("datafile.dat", "r", datafile);
  else
    datafile = fopen("datafile.dat", "r");
  if (datafile == NULL)
    exit(FileNotFound);
  if (outfile != NULL) {
    /*SUN*/
    outfile = freopen("outfile.dat", "w", outfile);
  } else
    outfile = fopen("outfile.dat", "w");
  if (outfile == NULL)
    exit(FileNotFound);
  if (msim != NULL) {
    /*SUN*/
    msim = freopen("msim.dat", "w", msim);
  } else
    msim = fopen("msim.dat", "w");
  if (msim == NULL)
    exit(FileNotFound);
  /*SUN*/
  if (lodprint) {   /*SUN*/
    if (lodfile != NULL)
      lodfile = freopen("lodfile.dat", "w", lodfile);
    else
      lodfile = fopen("lodfile.dat", "w");
    if (lodfile == NULL)
      exit(FileNotFound);
  }
  readloci();
}  /* inputdata */


Static Void readspseg()
{
  /*Reads from the speedfile in appropriate segments*/
  long i, j, a, b, sys;
  Char ch;
  information *WITH;
  long FORLIM;

  /*Note that each person[]^.unknown is set to FALSE as the person is read in*/
  i = lastspeed;
  j = lastspeed - lastseg;
  if (j > 0 && i <= segperson) {
    person[j]->unknown = true;
    person[j]->store = (information *)Malloc(sizeof(information));
    WITH = person[j]->store;
    FORLIM = nlocus;
    /*     person[j]^.store:=infoptr(NewPtr(SizeOf(information)));*/
    for (sys = 0; sys < FORLIM; sys++) {
      for (a = 0; a < maxall; a++) {
	for (b = 0; b < maxall; b++)
	  WITH->possible[sys][a][b] = false;
      }
    }
  }
  while (!P_eof(speedfile) && i <= segperson) {
    ch = getc(speedfile);
    if (ch == '\n')
      ch = ' ';
    if (ch == 'i' || ch == 'I') {
      fscanf(speedfile, "%c%ld%*[^\n]", &ch, &i);
      getc(speedfile);
      if (ch == '\n')
	ch = ' ';
      if (i > segperson)
	break;
      j = i - lastseg;
      person[j]->unknown = true;
      person[j]->store = (information *)Malloc(sizeof(information));
      WITH = person[j]->store;
      FORLIM = nlocus;
      /*         person[j]^.store:=infoptr(NewPtr(SizeOf(information)));*/
      for (sys = 0; sys < FORLIM; sys++) {
	for (a = 0; a < maxall; a++) {
	  for (b = 0; b < maxall; b++)
	    WITH->possible[sys][a][b] = false;
	}
      }
      continue;
    }
    WITH = person[j]->store;
    fscanf(speedfile, "%ld%ld%ld%*[^\n]", &sys, &a, &b);
    getc(speedfile);
    if (sys <= nlocus)
      WITH->possible[order[sys - 1] - 1][a - 1][b - 1] = true;
  }
  lastspeed = i;
  lastseg = segperson;
}


/* Local variables for readpedseg: */
struct LOC_readpedseg {
  long sequence;
  long startped[maxped], endped[maxped];
} ;

/*Exit from the procedure readpedseg*/

Local Void inputerror_(nerror, par1, par2, LINK)
long nerror, par1, par2;
struct LOC_readpedseg *LINK;
{
  printf("Fatal error detected in procedure readpedseg\n");
  switch (nerror) {

  case 11:
    printf("Error detected reading pedigree record %2ld. Illegal code for sex %2ld\n",
	   par1, par2);
    break;

  case 12:
    printf(
      "Error detected reading pedigree record at pedigree%2ld. Maximum number of pedigree records exceeded\n",
      par1);
    break;

  case 13:
    printf(
      "Error detected reading pedigree record %2ld. Maximum number of individuals exceeded\n",
      par1);
    break;

  case 14:
    printf(
      "Error detected reading pedigree record %2ld. Illegal binary factor code %2ld\n",
      par1, par2);
    break;

  case 15:
    printf(
      "Error detected reading pedigree record %2ld. No allelic pair for genotype\n",
      par1);
    break;

  case 16:
    printf(
      "Error detected reading pedigree record %2ld. Allele number %2ld exceeds maxall\n",
      par1, par2);
    break;

  case 17:
    printf(
      "Error detected reading pedigree record %2ld. Illegal allele number %2ld\n",
      par1, par2);
    break;

  case 20:
    printf("Number of recombination types (%3ld) exceeds maxrectype\n", par1);
    break;

  case 21:
    printf("Number of recombination types (%3ld) less than minimum of 1\n",
	   par1);
    break;

  case 22:
    printf(
      "End of file detected in tempdat by procedure readthg before all data found\n");
    break;

  case 23:
    printf(
      "Error detected reading iterated locus in datafile. Value (%3ld) greater than nlocus\n",
      par1);
    break;

  case 24:
    printf(
      "Error detected reading iterated locus in datafile. Illegal value (%3ld)\n",
      par1);
    break;

  case 25:
    printf("Number of iterated parameters greater then maxn\n");
    break;

  case 26:
    printf(
      "Error detected reading pedigree record %2ld. Liability class (%2ld) exceeds nclass\n",
      par1, par2);
    break;

  case 27:
    printf(
      "Error detected reading pedigree record %2ld. Illegal liability class (%2ld)\n",
      par1, par2);
    break;
  }
  printf("Press return to halt...\n");
  scanf("%*[^\n]");
  getchar();
  exit(0);   /*SUN*/
}  /* inputerror */


Local Void inputwarning_(nwarning, par1, par2, LINK)
long nwarning, par1, par2;
struct LOC_readpedseg *LINK;
{
  printf("Warning number from procedure readpedseg\n");
  switch (nwarning) {

  case 3:
    printf(
      "Non-standard affection status%4ld interpreted as normal in pedigree record%5ld\n",
      par2, par1);
    break;
  }
}  /* inputwarning */

/* Local variables for getphenotype: */
struct LOC_getphenotype {
  struct LOC_readpedseg *LINK;
  thisperson **p;
  long system;
} ;

Local Void readbin(phen, thislocus, LINK)
phenotype **phen;
locusvalues *thislocus;
struct LOC_getphenotype *LINK;
{
  long i, j;
  phenotype *WITH;
  long FORLIM;

  WITH = *phen;
  WITH->which = binary_;
  WITH->UU.phenf = 0;
  FORLIM = nfactor[LINK->system - 1];
  for (i = 1; i <= FORLIM; i++) {
    fscanf(ipedfile, "%ld", &j);
    if (j != 0 && j != 1)
      inputerror_(14L, (*LINK->p)->id, j, LINK->LINK);
    if (j == 1)
      WITH->UU.phenf = ((long)WITH->UU.phenf) | (1L << ((int)i));
  }
}  /* readbin */


Local Void readnumber(phen, thislocus, LINK)
phenotype **phen;
locusvalues *thislocus;
struct LOC_getphenotype *LINK;
{
  long i, j;
  phenotype *WITH;

  WITH = *phen;
  WITH->which = binary_;
  WITH->UU.phenf = 0;
  for (i = 1; i <= 2; i++) {
    fscanf(ipedfile, "%ld", &j);
    if (j > maxall)
      inputerror_(16L, (*LINK->p)->id, j, LINK->LINK);
    if (j < 0)
      inputerror_(17L, (*LINK->p)->id, j, LINK->LINK);
    if (j != 0)
      WITH->UU.phenf = ((long)WITH->UU.phenf) | (1L << ((int)j));
  }
}  /* readnumber */


Local Void readaff(phen, thislocus, LINK)
phenotype **phen;
locusvalues *thislocus;
struct LOC_getphenotype *LINK;
{
  long thisval;   /* Corrected to match Version 5.1 */
  phenotype *WITH;

  WITH = *phen;
  WITH->which = affection;
  fscanf(ipedfile, "%ld", &thisval);
  if (thisval == missaff)
    WITH->UU.U0.aff = 0;
  else if (thisval == affval)
    WITH->UU.U0.aff = 2;
  else {
    if (thisval != 1)
      inputwarning_(3L, (*LINK->p)->id, thisval, LINK->LINK);
    WITH->UU.U0.aff = 1;
  }
  if (thislocus->UU.U0.nclass == 1)
    WITH->UU.U0.liability = 1;
  else
    fscanf(ipedfile, "%ld", &WITH->UU.U0.liability);
  if (WITH->UU.U0.liability > thislocus->UU.U0.nclass)
    inputerror_(26L, (*LINK->p)->id, WITH->UU.U0.liability, LINK->LINK);
  if (WITH->UU.U0.liability <= 0)
    inputerror_(27L, (*LINK->p)->id, WITH->UU.U0.liability, LINK->LINK);
}  /* readaff */


Local Void readquan(phen, thislocus, LINK)
phenotype **phen;
locusvalues *thislocus;
struct LOC_getphenotype *LINK;
{
  long i;
  double xval;
  phenotype *WITH;
  long FORLIM;

  WITH = *phen;
  if (!sexlink || !(*LINK->p)->male) {
    WITH->which = quantitative;
    FORLIM = thislocus->UU.U1.ntrait;
    for (i = 0; i < FORLIM; i++)
      fscanf(ipedfile, "%lg", &WITH->UU.U1.x[i]);
    WITH->UU.U1.missing = true;
    FORLIM = thislocus->UU.U1.ntrait;
    for (i = 0; i < FORLIM; i++) {
      if (WITH->UU.U1.x[i] != missval)
	WITH->UU.U1.missing = false;
    }
    return;
  }
  WITH->which = affection;
  fscanf(ipedfile, "%lg", &xval);
  if (xval == missval)
    WITH->UU.U0.aff = missaff;
  else if (xval == affall)
    WITH->UU.U0.aff = affall;
  else
    WITH->UU.U0.aff = -11;
  WITH->UU.U0.liability = 1;
  FORLIM = thislocus->UU.U1.ntrait;
  for (i = 2; i <= FORLIM; i++)
    fscanf(ipedfile, "%lg", &xval);
}  /* readquan */


Local Void getphenotype(p_, LINK)
thisperson **p_;
struct LOC_readpedseg *LINK;
{
  struct LOC_getphenotype V;
  long thisread;
  thisperson *WITH;
  long FORLIM;


  V.LINK = LINK;
  V.p = p_;
  WITH = *V.p;
  FORLIM = nlocus;
  for (thisread = 0; thisread < FORLIM; thisread++) {
    V.system = order[thisread];
    WITH->phen[V.system - 1] = NULL;
    if (thislocus[V.system - 1]->which != null_)
      WITH->phen[V.system - 1] = (phenotype *)Malloc(sizeof(phenotype));
    /*        phen[system]:=phenpoint(NewPtr(SizeOf(phenotype)));*/
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
      /*         privphen:=phenpoint(NewPtr(SizeOf(phenotype)));*/
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


Local Void getinformative(LINK)
struct LOC_readpedseg *LINK;
{
  long i, j, k, l, m, count, nchild;
  thisperson *child;
  long FORLIM, FORLIM1, FORLIM2;
  phenotype *WITH;
  locusvalues *WITH1;
  long FORLIM3;

  if (risk) {
    FORLIM = nuped;
    for (i = 0; i < FORLIM; i++)
      informative[i] = true;
    return;
  }
  FORLIM = nuped;
  for (i = 0; i < FORLIM; i++) {
    informative[i] = false;
    FORLIM1 = LINK->endped[i];
    for (j = LINK->startped[i]; j <= FORLIM1; j++) {
      if (person[j]->foff != NULL) {
	nchild = 0;
	child = person[j]->foff;
	do {
	  nchild++;
	  if (person[j]->male)
	    child = child->nextpa;
	  else
	    child = child->nextma;
	} while (child != NULL);
	count = 0;   /* Moved to this position in Version 5.1 */
	if (nchild > 1 || nchild == 1 && person[j]->pa != NULL) {
	  FORLIM2 = nlocus;
	  for (k = 0; k < FORLIM2; k++) {
	    WITH = person[j]->phen[k];
	    WITH1 = thislocus[k];
	    if (WITH1->which != binary_)
	      count++;
	    else if (WITH->UU.phenf == 0)
	      count++;
	    else {
	      l = 0;
	      FORLIM3 = WITH1->nallele;
	      for (m = 1; m <= FORLIM3; m++) {
		if ((unsigned long)m < 32 && ((1L << m) & WITH->UU.phenf) != 0)
		  l++;
	      }
	      if (l > 1)
		count++;
	    }
	  }
	}
	if (count > 1)
	  informative[i] = true;
      }
    }
  }
}  /* getinformative */


Local Void getind(id, LINK)
long *id;
struct LOC_readpedseg *LINK;
{
  fscanf(ipedfile, "%ld", id);
  if (*id == 0)
    return;
  *id += LINK->sequence;
  if (*id > maxind)
    inputerror_(13L, *id, *id, LINK);
  if (person[*id] == NULL)
    person[*id] = (thisperson *)Malloc(sizeof(thisperson));
  /*       person[id]:=ind(NewPtr(SizeOf(thisperson))); ;*/

}  /* getind */

/*getind*/

Local Void multimarriage(p, LINK)
thisperson **p;
struct LOC_readpedseg *LINK;
{
  thisperson *q, *child, *WITH;

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


/* readspseg */
/*readspseg*/


Static Void readpedseg()
{
  struct LOC_readpedseg V;
  long i, newid, sex_, profield, newped, nuperson, thisone, thisped;
  thisperson *holdloop;
  long FORLIM;
  thisperson *WITH;

  /*multimarriage*/

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
  if (infile)
    newped = oldped;
  else
    fscanf(ipedfile, "%ld", &newped);
  thisped = newped;
  V.startped[0] = 1;
  while (!P_eof(ipedfile)) {
    nuperson++;
    getind(&thisone, &V);
    if (proband[nuped - 1] == NULL)
      proband[nuped - 1] = person[thisone];
    WITH = person[thisone];
    WITH->ped = thisped;
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
    fscanf(ipedfile, "%ld", &sex_);
    if (sex_ != 1 && sex_ != 2)
      inputerror_(11L, WITH->id, sex_, &V);
    if (sex_ == 1)
      WITH->male = true;
    else
      WITH->male = false;
    WITH->unknown = false;
    WITH->inloop = 0;
    fscanf(ipedfile, "%ld", &profield);
    if (profield == 1)
      proband[nuped - 1] = person[thisone];
    else if (profield > 1 && profield - 1 <= maxloop) {
      if (looppers[nuped - 1][profield - 2][1] == NULL)
	looppers[nuped - 1][profield - 2][1] = person[thisone];
      else
	looppers[nuped - 1][profield - 2][0] = person[thisone];
    }
    getphenotype(&person[thisone], &V);
    fscanf(ipedfile, "%*[^\n]");
    getc(ipedfile);
    if (!P_eof(ipedfile))
      fscanf(ipedfile, "%ld", &newped);
    /*We only know that we have reached the end of the current pedigree
    by the fact that the pedigree number changes i.e. newped<>thisped
    So we have to read one more record than I wish to read.
    Test for: newped<>thisped and nuped=opeds => save newped in oldped
    and don't read it on the next call to readpedseg*/
    if (newped != thisped && nuped == opeds)
    {   /*Exit the procedure readpedseg*/
      oldped = newped;
      infile = true;
      goto _L10;
    }
    if (newped == thisped)
      continue;
    V.sequence += nuperson;
    V.endped[nuped - 1] = V.sequence;
    nuperson = 0;
    nuped++;
    if (nuped > maxped)
      inputerror_(12L, newped, nuped, &V);
    V.startped[nuped - 1] = V.sequence + 1;
    for (i = 0; i < maxloop; i++) {
      looppers[nuped - 1][i][0] = NULL;
      looppers[nuped - 1][i][1] = NULL;
    }
    proband[nuped - 1] = NULL;
    thisped = newped;
  }  /*while not eof() loop*/
_L10:
  totperson = V.sequence + nuperson;
  V.endped[nuped - 1] = totperson;
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
  /*seg person is used in reading from the speedfile */
  segperson += totperson;
  FORLIM = totperson;
  for (thisone = 1; thisone <= FORLIM; thisone++)
    multimarriage(&person[thisone], &V);
  getinformative(&V);
}


/* readpedseg */
/*readpedseg*/


Static Void cleanup(p)
thisperson **p;
{
  thisperson *WITH;

  WITH = *p;
  Free(WITH->gen);
  /*     DisposPtr(Ptr(gen));*/
  WITH->gen = NULL;
}


/* Local variables for getvect: */
struct LOC_getvect {
  thisperson *p;
  hapvector hap1, hap2;
} ;

/*Local Void getgene PP((long syste, double val, struct LOC_getvect *LINK));
  Local Void ugetgene PP((long syste, double val, struct LOC_getvect *LINK));*/

static void getgene();
static void ugetgene();

Local double quanfun(phen, thislocus, i, j, mean, LINK)
phenotype *phen;
locusvalues *thislocus;
long i, j;
double *mean;
struct LOC_getvect *LINK;
{
  double val;
  long it, jt, FORLIM, FORLIM1;

  val = 1.0;
  if (phen->UU.U1.missing)
    return val;
  val = 0.0;
  FORLIM = thislocus->UU.U1.ntrait;
  for (it = 0; it < FORLIM; it++) {
    FORLIM1 = thislocus->UU.U1.ntrait;
    for (jt = 0; jt < FORLIM1; jt++) {
      if (i == j)
	val += thislocus->UU.U1.vmat[it]
	       [jt] * (phen->UU.U1.x[jt] - mean[jt]) *
	       (phen->UU.U1.x[it] - mean[it]);
      else
	val += thislocus->UU.U1.conmat * thislocus->UU.U1.vmat[it]
	       [jt] * (phen->UU.U1.x[jt] - mean[jt]) *
	       (phen->UU.U1.x[it] - mean[it]);
    }
  }
  val = thislocus->UU.U1.det * exp(-val * 0.5);
  if (i != j)
    val *= thislocus->UU.U1.contrait;
  return val;
}  /* quanfun */

/*quanfun*/

Local Void getval(syste, i, j, val, LINK)
long syste, i, j;
double *val;
struct LOC_getvect *LINK;
{
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
    *val *= WITH->UU.U0.pen[i][j - 1][WITH1->UU.U0.aff]
      [WITH1->UU.U0.liability - 1];
    break;
  }
}  /* getval */

/* Local variables for setval: */
struct LOC_setval {
  struct LOC_getvect *LINK;
  double val;
  long nhap1, nhap2;
} ;

Local Void prior(LINK)
struct LOC_setval *LINK;
{
  long i;   /*prior*/
  long FORLIM;
  locusvalues *WITH;
  thisarray *WITH1;

  if (!disequi) {
    if (sexlink && LINK->LINK->p->male) {
      FORLIM = nlocus;
      for (i = 0; i < FORLIM; i++) {
	WITH = thislocus[i];
	LINK->val *= WITH->freq[LINK->LINK->hap1[i] - 1];
      }
    } else {
      FORLIM = nlocus;
      for (i = 0; i < FORLIM; i++) {
	WITH = thislocus[i];
	LINK->val *= WITH->freq[LINK->LINK->hap1[i] - 1] *
		     WITH->freq[LINK->LINK->hap2[i] - 1];
      }
      if (LINK->nhap1 != LINK->nhap2)
	LINK->val = 2.0 * LINK->val;
    }
  } else {
    WITH1 = hapfreq;
    if (sexlink && LINK->LINK->p->male)
      LINK->val *= WITH1->genarray[LINK->nhap1 - 1];
    else {
      LINK->val *= WITH1->genarray[LINK->nhap1 - 1] *
		   WITH1->genarray[LINK->nhap2 - 1];
      if (LINK->nhap1 != LINK->nhap2)
	LINK->val = 2.0 * LINK->val;
    }
  }
  LINK->val *= segscale;
}  /* prior */

/*getval*/

Local Void setval(val_, LINK)
double val_;
struct LOC_getvect *LINK;
{
  struct LOC_setval V;
  long here, count, i, FORLIM;
  gennurec *WITH;
  thisarray *WITH1;

  /*prior*/

  V.LINK = LINK;
  V.val = val_;
  count = 1;
  V.nhap1 = 1;
  V.nhap2 = 1;
  if (sexlink && LINK->p->male) {
    FORLIM = nlocus;
    for (i = 0; i < FORLIM; i++)
      V.nhap1 += increment[i] * (LINK->hap1[i] - 1);
    here = V.nhap1;
  } else {
    WITH = gennustruct;
    FORLIM = nlocus;
    for (i = 0; i < FORLIM; i++) {
      V.nhap1 += increment[i] * (LINK->hap1[i] - 1);
      V.nhap2 += increment[i] * (LINK->hap2[i] - 1);
      if (LINK->hap1[i] != LINK->hap2[i])
	count *= 2;
      here = WITH->genenumber[V.nhap1 - 1][V.nhap2 - 1];
    }
  }
  if (LINK->p->pa == NULL)
    prior(&V);
  WITH1 = LINK->p->gen;
  if (disequi) {
    WITH1->genarray[here - 1] = V.val;
    return;
  }
  if (count != 1)
    count /= 2;
  for (i = 1; i <= count; i++) {
    WITH1->genarray[here - 1] = V.val;
    here++;
  }
}  /* setval */

/* Local variables for getgene: */
struct LOC_getgene {
  struct LOC_getvect *LINK;
  long syste;
  double val, newval;
} ;

Local Void facmale(LINK)
struct LOC_getgene *LINK;
{
  long j;   /*facmale*/
  thisperson *WITH;
  locusvalues *WITH1;
  long FORLIM;

  WITH = LINK->LINK->p;
  WITH1 = thislocus[LINK->syste - 1];
  FORLIM = WITH1->nallele;
  for (j = 1; j <= FORLIM; j++) {
    if (WITH->phen[LINK->syste - 1]->UU.phenf == WITH1->UU.allele[j - 1] ||
	WITH->phen[LINK->syste - 1]->UU.phenf == 0) {
      LINK->LINK->hap1[LINK->syste - 1] = j;
      if (LINK->syste != nlocus)
	getgene(LINK->syste + 1, LINK->val, LINK->LINK);
      else
	setval(LINK->val, LINK->LINK);
    }
  }
}  /* facmale */

/*facmale*/

Local Void affmale(LINK)
struct LOC_getgene *LINK;
{
  long j;   /*affmale*/
  locusvalues *WITH;
  long FORLIM;

  WITH = thislocus[LINK->syste - 1];
  FORLIM = WITH->nallele;
  for (j = 1; j <= FORLIM; j++) {
    LINK->newval = LINK->val;
    getval(LINK->syste, 0L, j, &LINK->newval, LINK->LINK);
    LINK->LINK->hap1[LINK->syste - 1] = j;
    if (LINK->newval != 0.0) {
      if (LINK->syste != nlocus)
	getgene(LINK->syste + 1, LINK->newval, LINK->LINK);
      else
	setval(LINK->newval, LINK->LINK);
    }
  }
}  /* affmale */

/*affmale*/

Local Void quanmale(LINK)
struct LOC_getgene *LINK;
{
  long j;   /*quanmale*/
  thisperson *WITH;
  locusvalues *WITH1;
  long FORLIM;

  WITH = LINK->LINK->p;
  WITH1 = thislocus[LINK->syste - 1];
  if (WITH->phen[LINK->syste - 1]->UU.U0.aff == affall ||
      WITH->phen[LINK->syste - 1]->UU.U0.aff == missaff) {
    LINK->newval = LINK->val;
    LINK->LINK->hap1[LINK->syste - 1] = affall;
    if (LINK->newval != 0.0) {
      if (LINK->syste != nlocus)
	getgene(LINK->syste + 1, LINK->newval, LINK->LINK);
      else
	setval(LINK->newval, LINK->LINK);
    }
  }
  if (WITH->phen[LINK->syste - 1]->UU.U0.aff == affall &&
      WITH->phen[LINK->syste - 1]->UU.U0.aff != missaff)
    return;
  FORLIM = WITH1->nallele;
  for (j = 1; j <= FORLIM; j++) {
    if (j != affall) {
      LINK->newval = LINK->val;
      LINK->LINK->hap1[LINK->syste - 1] = j;
      if (LINK->newval != 0.0) {
	if (LINK->syste != nlocus)
	  getgene(LINK->syste + 1, LINK->newval, LINK->LINK);
	else
	  setval(LINK->newval, LINK->LINK);
      }
    }
  }
}  /* quanmale */

/*quanmale*/

Local Void fac(LINK)
struct LOC_getgene *LINK;
{
  long i, j;   /*fac*/
  thisperson *WITH;
  locusvalues *WITH1;
  long FORLIM, FORLIM1;

  WITH = LINK->LINK->p;
  WITH1 = thislocus[LINK->syste - 1];
  FORLIM = WITH1->nallele;
  for (i = 1; i <= FORLIM; i++) {
    LINK->LINK->hap1[LINK->syste - 1] = i;
    FORLIM1 = WITH1->nallele;
    for (j = i; j <= FORLIM1; j++) {
      if (WITH->phen[LINK->syste - 1]->UU.phenf ==
	  (WITH1->UU.allele[i - 1] | WITH1->UU.allele[j - 1]) ||
	  WITH->phen[LINK->syste - 1]->UU.phenf == 0) {
	LINK->LINK->hap2[LINK->syste - 1] = j;
	if (LINK->syste != nlocus)
	  getgene(LINK->syste + 1, LINK->val, LINK->LINK);
	else
	  setval(LINK->val, LINK->LINK);
      }
    }
  }
  if (!disequi)
    return;
  WITH = LINK->LINK->p;
  WITH1 = thislocus[LINK->syste - 1];
  FORLIM = WITH1->nallele;
  for (i = 1; i <= FORLIM; i++) {
    LINK->LINK->hap2[LINK->syste - 1] = i;
    FORLIM1 = WITH1->nallele;
    for (j = i; j <= FORLIM1; j++) {
      if (WITH->phen[LINK->syste - 1]->UU.phenf ==
	  (WITH1->UU.allele[i - 1] | WITH1->UU.allele[j - 1]) ||
	  WITH->phen[LINK->syste - 1]->UU.phenf == 0) {
	LINK->LINK->hap1[LINK->syste - 1] = j;
	if (LINK->syste != nlocus)
	  getgene(LINK->syste + 1, LINK->val, LINK->LINK);
	else
	  setval(LINK->val, LINK->LINK);
      }
    }
  }
}  /* fac */

/*fac*/

Local Void aff(LINK)
struct LOC_getgene *LINK;
{
  /*Used with an affection status phenotype or when
  thislocus[syste]^which is null*/
  long i, j;
  locusvalues *WITH;
  long FORLIM, FORLIM1;

  WITH = thislocus[LINK->syste - 1];
  FORLIM = WITH->nallele;
  for (i = 1; i <= FORLIM; i++) {
    LINK->LINK->hap1[LINK->syste - 1] = i;
    FORLIM1 = WITH->nallele;
    for (j = i; j <= FORLIM1; j++) {
      if (!nohom[LINK->syste - 1] || i != affall || j != affall) {
	LINK->LINK->hap2[LINK->syste - 1] = j;
	LINK->newval = LINK->val;
	getval(LINK->syste, i, j, &LINK->newval, LINK->LINK);
	if (LINK->newval != 0.0) {
	  if (LINK->syste != nlocus)
	    getgene(LINK->syste + 1, LINK->newval, LINK->LINK);
	  else
	    setval(LINK->newval, LINK->LINK);
	}
      }
    }
  }
  if (!disequi)
    return;
  WITH = thislocus[LINK->syste - 1];
  FORLIM = WITH->nallele;
  for (i = 1; i <= FORLIM; i++) {
    LINK->LINK->hap2[LINK->syste - 1] = i;
    FORLIM1 = WITH->nallele;
    for (j = i; j <= FORLIM1; j++) {
      if (!nohom[LINK->syste - 1] || i != affall || j != affall) {
	LINK->LINK->hap1[LINK->syste - 1] = j;
	LINK->newval = LINK->val;
	getval(LINK->syste, i, j, &LINK->newval, LINK->LINK);
	if (LINK->newval != 0.0) {
	  if (LINK->syste != nlocus)
	    getgene(LINK->syste + 1, LINK->newval, LINK->LINK);
	  else
	    setval(LINK->newval, LINK->LINK);
	}
      }
    }
  }
}  /* aff */


Local Void quanval(LINK)
struct LOC_getgene *LINK;
{
  /*Uses this only when thislocus[syste]^.which is not null*/
  long i, j;
  locusvalues *WITH;
  long FORLIM, FORLIM1;

  WITH = thislocus[LINK->syste - 1];
  FORLIM = WITH->nallele;
  for (i = 1; i <= FORLIM; i++) {
    LINK->LINK->hap1[LINK->syste - 1] = i;
    FORLIM1 = WITH->nallele;
    for (j = i; j <= FORLIM1; j++) {
      if (!nohom[LINK->syste - 1] || i != affall || j != affall) {
	LINK->LINK->hap2[LINK->syste - 1] = j;
	LINK->newval = LINK->val;
	getval(LINK->syste, i, j, &LINK->newval, LINK->LINK);
	if (LINK->newval != 0.0) {
	  if (LINK->syste != nlocus)
	    getgene(LINK->syste + 1, LINK->newval, LINK->LINK);
	  else
	    setval(LINK->newval, LINK->LINK);
	}
      }
    }
  }
  if (!disequi)
    return;
  WITH = thislocus[LINK->syste - 1];
  FORLIM = WITH->nallele;
  for (i = 1; i <= FORLIM; i++) {
    LINK->LINK->hap2[LINK->syste - 1] = i;
    FORLIM1 = WITH->nallele;
    for (j = i; j <= FORLIM1; j++) {
      if (!nohom[LINK->syste - 1] || i != affall || j != affall) {
	LINK->LINK->hap1[LINK->syste - 1] = j;
	LINK->newval = LINK->val;
	getval(LINK->syste, i, j, &LINK->newval, LINK->LINK);
	if (LINK->newval != 0.0) {
	  if (LINK->syste != nlocus)
	    getgene(LINK->syste + 1, LINK->newval, LINK->LINK);
	  else
	    setval(LINK->newval, LINK->LINK);
	}
      }
    }
  }
}  /* quanval */

/*setval*/

Local Void getgene(syste_, val_, LINK)
long syste_;
double val_;
struct LOC_getvect *LINK;
{
  struct LOC_getgene V;
  locusvalues *WITH;


  V.LINK = LINK;
  V.syste = syste_;
  V.val = val_;
  WITH = thislocus[V.syste - 1];
  if (sexlink && LINK->p->male) {
    switch (WITH->which) {

    case binary_:
      facmale(&V);
      break;

    case affection:
      affmale(&V);
      break;

    case quantitative:
      quanmale(&V);
      break;

    case null_:
      if (WITH->privlocus->which == affection)
	affmale(&V);
      else
	quanmale(&V);
      break;
    }
    return;
  }
  switch (WITH->which) {

  case binary_:
    fac(&V);
    break;

  case affection:
    aff(&V);
    break;

  case quantitative:
    quanval(&V);
    break;

  case null_:
    aff(&V);
    break;
  }
}  /* getgene */

/* Local variables for ugetgene: */
struct LOC_ugetgene {
  struct LOC_getvect *LINK;
  long syste;
  double val, newval;
} ;

Local Void facmale_(LINK)
struct LOC_ugetgene *LINK;
{
  long j;   /*facmale*/
  thisperson *WITH;
  information *WITH1;
  locusvalues *WITH2;
  long FORLIM;

  WITH = LINK->LINK->p;
  WITH1 = LINK->LINK->p->store;
  WITH2 = thislocus[LINK->syste - 1];
  FORLIM = WITH2->nallele;
  for (j = 1; j <= FORLIM; j++) {
    if (WITH->phen[LINK->syste - 1]->UU.phenf == WITH2->UU.allele[j - 1] ||
	WITH->phen[LINK->syste - 1]->UU.phenf == 0) {
      if (WITH1->possible[LINK->syste - 1][0][j - 1]) {
	LINK->LINK->hap1[LINK->syste - 1] = j;
	if (LINK->syste != nlocus)
	  ugetgene(LINK->syste + 1, LINK->val, LINK->LINK);
	else
	  setval(LINK->val, LINK->LINK);
      }
    }
  }
}  /* facmale */

/*facmale*/

Local Void affmale_(LINK)
struct LOC_ugetgene *LINK;
{
  long j;
  information *WITH;
  locusvalues *WITH1;
  long FORLIM;

  WITH = LINK->LINK->p->store;
  WITH1 = thislocus[LINK->syste - 1];
  FORLIM = WITH1->nallele;
  for (j = 1; j <= FORLIM; j++) {
    if (WITH->possible[LINK->syste - 1][0][j - 1]) {
      LINK->newval = LINK->val;
      getval(LINK->syste, 0L, j, &LINK->newval, LINK->LINK);
      LINK->LINK->hap1[LINK->syste - 1] = j;
      if (LINK->newval != 0.0) {
	if (LINK->syste != nlocus)
	  ugetgene(LINK->syste + 1, LINK->newval, LINK->LINK);
	else
	  setval(LINK->newval, LINK->LINK);
      }
    }
  }
}  /* affmale */


Local Void quanmale_(LINK)
struct LOC_ugetgene *LINK;
{
  long j;
  thisperson *WITH;
  information *WITH1;
  locusvalues *WITH2;
  information *WITH3;
  long FORLIM;

  WITH = LINK->LINK->p;
  WITH1 = LINK->LINK->p->store;
  WITH2 = thislocus[LINK->syste - 1];
  if (WITH->phen[LINK->syste - 1]->UU.U0.aff == affall ||
      WITH->phen[LINK->syste - 1]->UU.U0.aff == missaff) {
    if (WITH1->possible[LINK->syste - 1][0][affall - 1]) {
      LINK->newval = LINK->val;
      LINK->LINK->hap1[LINK->syste - 1] = affall;
      if (LINK->newval != 0.0) {
	if (LINK->syste != nlocus)
	  ugetgene(LINK->syste + 1, LINK->newval, LINK->LINK);
	else
	  setval(LINK->newval, LINK->LINK);
      }
    }
  }
  WITH3 = LINK->LINK->p->store;
  if (WITH->phen[LINK->syste - 1]->UU.U0.aff == affall &&
      WITH->phen[LINK->syste - 1]->UU.U0.aff != missaff)
    return;
  FORLIM = WITH2->nallele;
  for (j = 1; j <= FORLIM; j++) {
    if (j != affall) {
      if (WITH3->possible[LINK->syste - 1][0][j - 1]) {
	LINK->newval = LINK->val;
	LINK->LINK->hap1[LINK->syste - 1] = j;
	if (LINK->newval != 0.0) {
	  if (LINK->syste != nlocus)
	    ugetgene(LINK->syste + 1, LINK->newval, LINK->LINK);
	  else
	    setval(LINK->newval, LINK->LINK);
	}
      }
    }
  }
}  /* quanmale */


Local Void fac_(LINK)
struct LOC_ugetgene *LINK;
{
  long i, j;
  thisperson *WITH;
  information *WITH1;
  locusvalues *WITH2;
  long FORLIM, FORLIM1;

  WITH = LINK->LINK->p;
  WITH1 = LINK->LINK->p->store;
  WITH2 = thislocus[LINK->syste - 1];
  FORLIM = WITH2->nallele;
  for (i = 1; i <= FORLIM; i++) {
    LINK->LINK->hap1[LINK->syste - 1] = i;
    FORLIM1 = WITH2->nallele;
    for (j = i; j <= FORLIM1; j++) {
      if (WITH->phen[LINK->syste - 1]->UU.phenf ==
	  (WITH2->UU.allele[i - 1] | WITH2->UU.allele[j - 1]) ||
	  WITH->phen[LINK->syste - 1]->UU.phenf == 0) {
	if (WITH1->possible[LINK->syste - 1][i - 1][j - 1]) {
	  LINK->LINK->hap2[LINK->syste - 1] = j;
	  if (LINK->syste != nlocus)
	    ugetgene(LINK->syste + 1, LINK->val, LINK->LINK);
	  else
	    setval(LINK->val, LINK->LINK);
	}
      }
    }
  }
  if (!disequi)
    return;
  WITH = LINK->LINK->p;
  WITH1 = LINK->LINK->p->store;
  WITH2 = thislocus[LINK->syste - 1];
  FORLIM = WITH2->nallele;
  for (i = 1; i <= FORLIM; i++) {
    LINK->LINK->hap2[LINK->syste - 1] = i;
    FORLIM1 = WITH2->nallele;
    for (j = i; j <= FORLIM1; j++) {
      if (WITH->phen[LINK->syste - 1]->UU.phenf ==
	  (WITH2->UU.allele[i - 1] | WITH2->UU.allele[j - 1]) ||
	  WITH->phen[LINK->syste - 1]->UU.phenf == 0) {
	if (WITH1->possible[LINK->syste - 1][i - 1][j - 1]) {
	  LINK->LINK->hap1[LINK->syste - 1] = j;
	  if (LINK->syste != nlocus)
	    ugetgene(LINK->syste + 1, LINK->val, LINK->LINK);
	  else
	    setval(LINK->val, LINK->LINK);
	}
      }
    }
  }
}  /* fac */


Local Void aff_(LINK)
struct LOC_ugetgene *LINK;
{
  /*Used with an affection status phenotype or when
  thislocus[syste]^which is null*/
  long i, j;
  thisperson *WITH;
  information *WITH1;
  locusvalues *WITH2;
  long FORLIM, FORLIM1;

  WITH = LINK->LINK->p;
  WITH1 = LINK->LINK->p->store;
  WITH2 = thislocus[LINK->syste - 1];
  FORLIM = WITH2->nallele;
  for (i = 1; i <= FORLIM; i++) {
    LINK->LINK->hap1[LINK->syste - 1] = i;
    FORLIM1 = WITH2->nallele;
    for (j = i; j <= FORLIM1; j++) {
      if (!nohom[LINK->syste - 1] || i != affall || j != affall) {
	if (WITH1->possible[LINK->syste - 1][i - 1][j - 1]) {
	  LINK->LINK->hap2[LINK->syste - 1] = j;
	  LINK->newval = LINK->val;
	  getval(LINK->syste, i, j, &LINK->newval, LINK->LINK);
	  if (LINK->newval != 0.0) {
	    if (LINK->syste != nlocus)
	      ugetgene(LINK->syste + 1, LINK->newval, LINK->LINK);
	    else
	      setval(LINK->newval, LINK->LINK);
	  }
	}
      }
    }
  }
  if (!disequi)
    return;
  WITH = LINK->LINK->p;
  WITH1 = LINK->LINK->p->store;
  WITH2 = thislocus[LINK->syste - 1];
  FORLIM = WITH2->nallele;
  for (i = 1; i <= FORLIM; i++) {
    LINK->LINK->hap2[LINK->syste - 1] = i;
    FORLIM1 = WITH2->nallele;
    for (j = i; j <= FORLIM1; j++) {
      if (!nohom[LINK->syste - 1] || i != affall || j != affall) {
	if (WITH1->possible[LINK->syste - 1][i - 1][j - 1]) {
	  LINK->LINK->hap1[LINK->syste - 1] = j;
	  LINK->newval = LINK->val;
	  getval(LINK->syste, i, j, &LINK->newval, LINK->LINK);
	  if (LINK->newval != 0.0) {
	    if (LINK->syste != nlocus)
	      ugetgene(LINK->syste + 1, LINK->newval, LINK->LINK);
	    else
	      setval(LINK->newval, LINK->LINK);
	  }
	}
      }
    }
  }
}  /* aff */


Local Void quanval_(LINK)
struct LOC_ugetgene *LINK;
{
  /*Uses this only when thislocus[syste]^.which is not null*/
  long i, j;
  thisperson *WITH;
  information *WITH1;
  locusvalues *WITH2;
  long FORLIM, FORLIM1;

  WITH = LINK->LINK->p;
  WITH1 = LINK->LINK->p->store;
  WITH2 = thislocus[LINK->syste - 1];
  FORLIM = WITH2->nallele;
  for (i = 1; i <= FORLIM; i++) {
    LINK->LINK->hap1[LINK->syste - 1] = i;
    FORLIM1 = WITH2->nallele;
    for (j = i; j <= FORLIM1; j++) {
      if (!nohom[LINK->syste - 1] || i != affall || j != affall) {
	if (WITH1->possible[LINK->syste - 1][i - 1][j - 1]) {
	  LINK->LINK->hap2[LINK->syste - 1] = j;
	  LINK->newval = LINK->val;
	  getval(LINK->syste, i, j, &LINK->newval, LINK->LINK);
	  if (LINK->newval != 0.0) {
	    if (LINK->syste != nlocus)
	      ugetgene(LINK->syste + 1, LINK->newval, LINK->LINK);
	    else
	      setval(LINK->newval, LINK->LINK);
	  }
	}
      }
    }
  }
  if (!disequi)
    return;
  WITH = LINK->LINK->p;
  WITH1 = LINK->LINK->p->store;
  WITH2 = thislocus[LINK->syste - 1];
  FORLIM = WITH2->nallele;
  for (i = 1; i <= FORLIM; i++) {
    LINK->LINK->hap2[LINK->syste - 1] = i;
    FORLIM1 = WITH2->nallele;
    for (j = i; j <= FORLIM1; j++) {
      if (!nohom[LINK->syste - 1] || i != affall || j != affall) {
	if (WITH1->possible[LINK->syste - 1][i - 1][j - 1]) {
	  LINK->LINK->hap1[LINK->syste - 1] = j;
	  LINK->newval = LINK->val;
	  getval(LINK->syste, i, j, &LINK->newval, LINK->LINK);
	  if (LINK->newval != 0.0) {
	    if (LINK->syste != nlocus)
	      ugetgene(LINK->syste + 1, LINK->newval, LINK->LINK);
	    else
	      setval(LINK->newval, LINK->LINK);
	  }
	}
      }
    }
  }
}  /* quanval */


Local Void ugetgene(syste_, val_, LINK)
long syste_;
double val_;
struct LOC_getvect *LINK;
{
  struct LOC_ugetgene V;
  locusvalues *WITH;


  V.LINK = LINK;
  V.syste = syste_;
  V.val = val_;
  WITH = thislocus[V.syste - 1];
  if (sexlink && LINK->p->male) {
    switch (WITH->which) {

    case binary_:
      facmale_(&V);
      break;

    case affection:
      affmale_(&V);
      break;

    case quantitative:
      quanmale_(&V);
      break;

    case null_:
      if (WITH->privlocus->which == affection)
	affmale_(&V);
      else
	quanmale_(&V);
      break;
    }
    return;
  }
  switch (WITH->which) {

  case binary_:
    fac_(&V);
    break;

  case affection:
    aff_(&V);
    break;

  case quantitative:
    quanval_(&V);
    break;

  case null_:
    aff_(&V);
    break;
  }
}  /* ugetgene */


/* cleanup */
/*cleanup*/

Static Void getvect(p_)
thisperson *p_;
{
  struct LOC_getvect V;


  V.p = p_;
  if (V.p->unknown)
    ugetgene(1L, 1.0, &V);
  else
    getgene(1L, 1.0, &V);
}  /* getvect */


/* Local variables for likelihood: */
struct LOC_likelihood {
  long thisped;
  thisperson *proband;
  long loopgen[maxloop];
  double homo, hetero;
  long nuscale;
  thisarray *holdpoint[maxloop];
} ;

/*Local Void collapsedown PP((thisperson *p, struct LOC_likelihood *LINK));*/

static void collapsedown();

Local Void prob(p, LINK)
thisperson **p;
struct LOC_likelihood *LINK;
{
  long i;
  thisperson *WITH;
  thisarray *WITH1;
  long FORLIM;

  WITH = *p;
  if (WITH->gen != NULL)
    return;
  WITH->gen = (thisarray *)Malloc(sizeof(thisarray));
  WITH1 = WITH->gen;
  FORLIM = fgeno;
  /*       gen:=genpoint(NewPtr(SizeOf(thisarray)));*/
  for (i = 0; i < FORLIM; i++)
    WITH1->genarray[i] = 0.0;
  if (WITH->inloop > 0) {
    if (looppers[LINK->thisped - 1][WITH->inloop - 1][0] == *p)
      WITH->gen->genarray[LINK->loopgen[WITH->inloop - 1] - 1] = LINK->
	    holdpoint[WITH->inloop - 1]->
	  genarray[LINK->loopgen[WITH->inloop - 1] - 1];
    else
      WITH->gen->genarray[LINK->loopgen[WITH->inloop - 1] - 1] = 1.0;
    return;
  }
  getvect(*p);
  if ((*p)->pa == NULL)
    LINK->nuscale++;
}  /* prob */

/* Local variables for seg: */
struct LOC_seg {
  struct LOC_likelihood *LINK;
  thisperson **p, **q, **r, *child, *father, *mother;
  long fseg, sseg, sstart, send, fstart, fend, nfirst, nsecond, firstseg,
       secondseg;
  double pf, ps;
  thetavalues *firstsex, *secondsex;
} ;

Local Void getapprox(LINK)
struct LOC_seg *LINK;
{
  long first;
  double maxval;   /*getapprox*/
  thisarray *WITH;
  long FORLIM;
  approxrec *WITH1;

  maxval = (*LINK->p)->gen->genarray[0];
  WITH = (*LINK->p)->gen;
  FORLIM = fgeno;
  for (first = 0; first < FORLIM; first++) {
    if (WITH->genarray[first] > maxval)
      maxval = WITH->genarray[first];
  }
  WITH1 = approxstruct;
  WITH = (*LINK->p)->gen;
  FORLIM = fgeno;
  for (first = 0; first < FORLIM; first++)
    WITH1->approxarray[LINK->LINK->thisped - 1]
      [first] = (WITH->genarray[first] > epsilon);
  if (lasttime)
    return;
  WITH1 = approxstruct;
  WITH = (*LINK->p)->gen;
  FORLIM = fgeno;
  for (first = 0; first < FORLIM; first++) {
    if (!WITH1->approxarray[LINK->LINK->thisped - 1][first])
      WITH->genarray[first] = 0.0;
  }
}  /* getapprox */

/*getapprox*/

Local double msegsex(LINK)
struct LOC_seg *LINK;
{
  long g1, g2, g3, g4, g5, g6, g7, g8, ms, mf, ms1, ms2, mf1, mf2, j, k, f1,
       f2, s1, s2;
  double val, temp2;
  double temp[maxchild];   /*msegsex*/
  long FORLIM;
  gennurec *WITH;
  thetavalues *WITH1;
  long FORLIM1;
  thisarray *WITH2;

  FORLIM = nchild;
  for (k = 0; k < FORLIM; k++)
    temp[k] = 0.0;
  WITH = gennustruct;
  if ((*LINK->p)->male) {
    mf = muthap[LINK->fseg - 1];
    LINK->secondseg = LINK->sseg;
    WITH1 = LINK->secondsex;
    FORLIM = LINK->send;
    for (j = LINK->sstart - 1; j < FORLIM; j++) {
      temp2 = WITH1->segprob[j];
      s1 = seghap1[LINK->secondseg - 1];
      s2 = seghap2[LINK->secondseg - 1];
      ms1 = muthap[s1 - 1];
      ms2 = muthap[s2 - 1];
      if (s1 != s2) {
	g1 = WITH->genenumber[LINK->fseg - 1][s1 - 1];
	g2 = WITH->genenumber[LINK->fseg - 1][s2 - 1];
	g3 = WITH->genenumber[LINK->fseg - 1][ms1 - 1];
	g4 = WITH->genenumber[LINK->fseg - 1][ms2 - 1];
	g5 = WITH->genenumber[mf - 1][s1 - 1];
	g6 = WITH->genenumber[mf - 1][s2 - 1];
	g7 = WITH->genenumber[mf - 1][ms1 - 1];
	g8 = WITH->genenumber[mf - 1][ms2 - 1];
	FORLIM1 = nchild;
	for (k = 0; k < FORLIM1; k++) {
	  WITH2 = thischild[k];
	  if (malechild[k])
	    temp[k] += temp2 *
		((1 - LINK->ps) * (WITH2->genarray[s1 - 1] + WITH2->genarray[s2 - 1]) +
		 LINK->ps * (WITH2->genarray[ms1 - 1] + WITH2->genarray[ms2 - 1]));
	  else
	    temp[k] += temp2 * ((1 - LINK->pf) * (1 - LINK->ps) *
		    (WITH2->genarray[g1 - 1] + WITH2->genarray[g2 - 1]) +
		  (1 - LINK->pf) * LINK->ps * (WITH2->genarray[g3 - 1] +
		      WITH2->genarray[g4 - 1]) + LINK->pf * (1 - LINK->ps) *
		    (WITH2->genarray[g5 - 1] + WITH2->genarray[g6 - 1]) +
		  LINK->pf * LINK->ps *
		    (WITH2->genarray[g7 - 1] + WITH2->genarray[g8 - 1]));
/* p2c: msim.p, line 4352: 
 * Note: Line breaker spent 0.0 seconds, 5000 tries on line 3556 [251] */
	}
      } else {
	g1 = WITH->genenumber[LINK->fseg - 1][s1 - 1];
	g3 = WITH->genenumber[LINK->fseg - 1][ms1 - 1];
	g5 = WITH->genenumber[mf - 1][s1 - 1];
	g7 = WITH->genenumber[mf - 1][ms1 - 1];
	FORLIM1 = nchild;
	for (k = 0; k < FORLIM1; k++) {
	  WITH2 = thischild[k];
	  if (malechild[k])
	    temp[k] += 2.0 * temp2 * ((1 - LINK->ps) * WITH2->genarray[s1 - 1] +
				      LINK->ps * WITH2->genarray[ms1 - 1]);
	  else
	    temp[k] += 2.0 * temp2 *
		((1 - LINK->pf) * (1 - LINK->ps) * WITH2->genarray[g1 - 1] +
		 LINK->ps * (1 - LINK->pf) * WITH2->genarray[g3 - 1] +
		 LINK->pf * (1 - LINK->ps) * WITH2->genarray[g5 - 1] +
		 LINK->pf * LINK->ps * WITH2->genarray[g7 - 1]);
		/*Ask Mark*/
	}
      }
      LINK->secondseg++;
    }
  } else {
    LINK->firstseg = LINK->fseg;
    ms = muthap[LINK->sseg - 1];
    WITH1 = LINK->firstsex;
    FORLIM = LINK->fend;
    for (j = LINK->fstart - 1; j < FORLIM; j++) {
      temp2 = WITH1->segprob[j];
      f1 = seghap1[LINK->firstseg - 1];
      f2 = seghap2[LINK->firstseg - 1];
      mf1 = muthap[f1 - 1];
      mf2 = muthap[f2 - 1];
      if (f1 != f2) {
	g1 = WITH->genenumber[LINK->sseg - 1][f1 - 1];
	g2 = WITH->genenumber[LINK->sseg - 1][f2 - 1];
	g3 = WITH->genenumber[LINK->sseg - 1][mf1 - 1];
	g4 = WITH->genenumber[LINK->sseg - 1][mf2 - 1];
	g5 = WITH->genenumber[ms - 1][f1 - 1];
	g6 = WITH->genenumber[ms - 1][f2 - 1];
	g7 = WITH->genenumber[ms - 1][mf1 - 1];
	g8 = WITH->genenumber[ms - 1][mf2 - 1];
	FORLIM1 = nchild;
	for (k = 0; k < FORLIM1; k++) {
	  WITH2 = thischild[k];
	  if (malechild[k])
	    temp[k] += temp2 *
		((1 - LINK->pf) * (WITH2->genarray[f1 - 1] + WITH2->genarray[f2 - 1]) +
		 LINK->pf * (WITH2->genarray[mf1 - 1] + WITH2->genarray[mf2 - 1]));
	  else
	    temp[k] += temp2 * ((1 - LINK->pf) * (1 - LINK->ps) *
		    (WITH2->genarray[g1 - 1] + WITH2->genarray[g2 - 1]) +
		  (1 - LINK->ps) * LINK->pf * (WITH2->genarray[g3 - 1] +
		      WITH2->genarray[g4 - 1]) + LINK->ps * (1 - LINK->pf) *
		    (WITH2->genarray[g5 - 1] + WITH2->genarray[g6 - 1]) +
		  LINK->pf * LINK->ps *
		    (WITH2->genarray[g7 - 1] + WITH2->genarray[g8 - 1]));
/* p2c: msim.p, line 4352: 
 * Note: Line breaker spent 0.0 seconds, 5000 tries on line 3616 [251] */
	}
      } else {
	g1 = WITH->genenumber[LINK->sseg - 1][f1 - 1];
	g3 = WITH->genenumber[LINK->sseg - 1][mf1 - 1];
	g5 = WITH->genenumber[ms - 1][f1 - 1];
	g7 = WITH->genenumber[ms - 1][mf1 - 1];
	FORLIM1 = nchild;
	for (k = 0; k < FORLIM1; k++) {
	  WITH2 = thischild[k];
	  if (malechild[k])
	    temp[k] += 2.0 * temp2 * ((1 - LINK->pf) * WITH2->genarray[f1 - 1] +
				      LINK->pf * WITH2->genarray[mf1 - 1]);
	  else
	    temp[k] += 2.0 * temp2 *
		((1 - LINK->pf) * (1 - LINK->ps) * WITH2->genarray[g1 - 1] +
		 LINK->pf * (1 - LINK->ps) * WITH2->genarray[g3 - 1] +
		 LINK->ps * (1 - LINK->pf) * WITH2->genarray[g5 - 1] +
		 LINK->pf * LINK->ps * WITH2->genarray[g7 - 1]);
	}
      }
      LINK->firstseg++;
    }
  }
  val = 1.0;
  FORLIM = nchild;
  for (k = 0; k < FORLIM; k++)
    val *= temp[k];
  return val;
}  /* msegsex */

/*msegsex*/

Local double msegsexf(LINK)
struct LOC_seg *LINK;
{
  long g1, g2, g3, g4, g5, g6, g7, g8, mf, ms1, ms2, j, k, l, s1, s2, slength;
  double val, temp1;
  double temp[maxchild][maxseg];
  double temp2[maxseg];   /*msegsexf*/
  long FORLIM, FORLIM1;
  gennurec *WITH;
  thetavalues *WITH1;
  thisarray *WITH2;

  mf = muthap[LINK->fseg - 1];
  slength = LINK->send - LINK->sstart + 1;
  FORLIM = nchild;
  for (k = 0; k < FORLIM; k++) {
    for (l = 0; l < slength; l++)
      temp[k][l] = 0.0;
  }
  LINK->secondseg = LINK->sseg;
  WITH = gennustruct;
  WITH1 = LINK->secondsex;
  FORLIM = LINK->send;
  for (j = LINK->sstart - 1; j < FORLIM; j++) {
    for (l = 0; l < slength; l++)
      temp2[l] = WITH1->segprob[j + l * slength];
    s1 = seghap1[LINK->secondseg - 1];
    s2 = seghap2[LINK->secondseg - 1];
    ms1 = muthap[s1 - 1];
    ms2 = muthap[s2 - 1];
    if (s1 != s2) {
      g1 = WITH->genenumber[LINK->fseg - 1][s1 - 1];
      g2 = WITH->genenumber[LINK->fseg - 1][s2 - 1];
      g3 = WITH->genenumber[LINK->fseg - 1][ms1 - 1];
      g4 = WITH->genenumber[LINK->fseg - 1][ms2 - 1];
      g5 = WITH->genenumber[mf - 1][s1 - 1];
      g6 = WITH->genenumber[mf - 1][s2 - 1];
      g7 = WITH->genenumber[mf - 1][ms1 - 1];
      g8 = WITH->genenumber[mf - 1][ms2 - 1];
      FORLIM1 = nchild;
      for (k = 0; k < FORLIM1; k++) {
	WITH2 = thischild[k];
	if (malechild[k])
	  val = (1 - LINK->ps) * (WITH2->genarray[s1 - 1] + WITH2->genarray[s2 - 1]) +
		LINK->ps * (WITH2->genarray[ms1 - 1] + WITH2->genarray[ms2 - 1]);
	else
	  val = (1 - LINK->pf) * (1 - LINK->ps) *
		(WITH2->genarray[g1 - 1] + WITH2->genarray[g2 - 1]) +
		(1 - LINK->pf) * LINK->ps *
		(WITH2->genarray[g3 - 1] + WITH2->genarray[g4 - 1]) +
		LINK->pf * (1 - LINK->ps) *
		(WITH2->genarray[g5 - 1] + WITH2->genarray[g6 - 1]) +
		LINK->pf * LINK->ps *
		(WITH2->genarray[g7 - 1] + WITH2->genarray[g8 - 1]);
/* p2c: msim.p, line 4352: 
 * Note: Line breaker spent 0.0 seconds, 5000 tries on line 3704 [251] */
	for (l = 0; l < slength; l++)
	  temp[k][l] += temp2[l] * val;
      }
    } else {
      g1 = WITH->genenumber[LINK->fseg - 1][s1 - 1];
      g3 = WITH->genenumber[LINK->fseg - 1][ms1 - 1];
      g5 = WITH->genenumber[mf - 1][s1 - 1];
      g7 = WITH->genenumber[mf - 1][ms1 - 1];
      FORLIM1 = nchild;
      for (k = 0; k < FORLIM1; k++) {
	WITH2 = thischild[k];
	if (malechild[k])
	  val = 2.0 * ((1 - LINK->ps) * WITH2->genarray[s1 - 1] +
		       LINK->ps * WITH2->genarray[ms1 - 1]);
	else
	  val = 2.0 * ((1 - LINK->pf) * (1 - LINK->ps) * WITH2->genarray[g1 - 1] +
		       LINK->ps * (1 - LINK->pf) * WITH2->genarray[g3 - 1] +
		       LINK->pf * (1 - LINK->ps) * WITH2->genarray[g5 - 1] +
		       LINK->pf * LINK->ps * WITH2->genarray[g7 - 1]);
	for (l = 0; l < slength; l++)
	  temp[k][l] += temp2[l] * val;
      }
    }
    LINK->secondseg++;
  }
  temp1 = 0.0;
  for (l = 0; l < slength; l++) {
    val = 1.0;
    FORLIM1 = nchild;
    for (k = 0; k < FORLIM1; k++)
      val *= temp[k][l];
    temp1 += val;
  }
  return temp1;
}  /* msegsexf */

/*msegsexf*/

Local double segsex(LINK)
struct LOC_seg *LINK;
{
  long g1, g2, j, k, f1, f2, s1, s2;
  double val, temp2;
  double temp[maxchild];   /*segsex*/
  long FORLIM;
  gennurec *WITH;
  thetavalues *WITH1;
  long FORLIM1;
  thisarray *WITH2;

  FORLIM = nchild;
  for (k = 0; k < FORLIM; k++)
    temp[k] = 0.0;
  WITH = gennustruct;
  if ((*LINK->p)->male) {
    LINK->secondseg = LINK->sseg;
    WITH1 = LINK->secondsex;
    FORLIM = LINK->send;
    for (j = LINK->sstart - 1; j < FORLIM; j++) {
      temp2 = WITH1->segprob[j];
      s1 = seghap1[LINK->secondseg - 1];
      s2 = seghap2[LINK->secondseg - 1];
      if (s1 != s2) {
	g1 = WITH->genenumber[LINK->fseg - 1][s1 - 1];
	g2 = WITH->genenumber[LINK->fseg - 1][s2 - 1];
	FORLIM1 = nchild;
	for (k = 0; k < FORLIM1; k++) {
	  WITH2 = thischild[k];
	  if (malechild[k])
	    temp[k] += temp2 * (WITH2->genarray[s1 - 1] + WITH2->genarray[s2 - 1]);
	  else
	    temp[k] += temp2 * (WITH2->genarray[g1 - 1] + WITH2->genarray[g2 - 1]);
	}
      } else {
	g1 = WITH->genenumber[LINK->fseg - 1][s1 - 1];
	FORLIM1 = nchild;
	for (k = 0; k < FORLIM1; k++) {
	  WITH2 = thischild[k];
	  if (malechild[k])
	    temp[k] += 2.0 * temp2 * WITH2->genarray[s1 - 1];
	  else
	    temp[k] += 2.0 * temp2 * WITH2->genarray[g1 - 1];
	}
      }
      LINK->secondseg++;
    }
  } else {
    LINK->firstseg = LINK->fseg;
    WITH1 = LINK->firstsex;
    FORLIM = LINK->fend;
    for (j = LINK->fstart - 1; j < FORLIM; j++) {
      temp2 = WITH1->segprob[j];
      f1 = seghap1[LINK->firstseg - 1];
      f2 = seghap2[LINK->firstseg - 1];
      if (f1 != f2) {
	g1 = WITH->genenumber[LINK->sseg - 1][f1 - 1];
	g2 = WITH->genenumber[LINK->sseg - 1][f2 - 1];
	FORLIM1 = nchild;
	for (k = 0; k < FORLIM1; k++) {
	  WITH2 = thischild[k];
	  if (malechild[k])
	    temp[k] += temp2 * (WITH2->genarray[f1 - 1] + WITH2->genarray[f2 - 1]);
	  else
	    temp[k] += temp2 * (WITH2->genarray[g1 - 1] + WITH2->genarray[g2 - 1]);
	}
      } else {
	g1 = WITH->genenumber[LINK->sseg - 1][f1 - 1];
	FORLIM1 = nchild;
	for (k = 0; k < FORLIM1; k++) {
	  WITH2 = thischild[k];
	  if (malechild[k])
	    temp[k] += 2.0 * temp2 * WITH2->genarray[f1 - 1];
	  else
	    temp[k] += 2.0 * temp2 * WITH2->genarray[g1 - 1];
	}
      }
      LINK->firstseg++;
    }
  }
  val = 1.0;
  FORLIM = nchild;
  for (k = 0; k < FORLIM; k++)
    val *= temp[k];
  return val;
}  /* segsex */

/*segsex*/

Local double segsexf(LINK)
struct LOC_seg *LINK;
{
  long g1, g2, j, k, l, s1, s2, slength;
  double val, temp1;
  double temp[maxchild][maxseg];
  double temp2[maxseg];
  long FORLIM, FORLIM1;
  gennurec *WITH;
  thetavalues *WITH1;
  thisarray *WITH2;

  slength = LINK->send - LINK->sstart + 1;
  FORLIM = nchild;
  for (k = 0; k < FORLIM; k++) {
    for (l = 0; l < slength; l++)
      temp[k][l] = 0.0;
  }
  LINK->secondseg = LINK->sseg;
  WITH = gennustruct;
  WITH1 = LINK->secondsex;
  FORLIM = LINK->send;
  for (j = LINK->sstart - 1; j < FORLIM; j++) {
    for (l = 0; l < slength; l++)
      temp2[l] = WITH1->segprob[j + l * slength];
    s1 = seghap1[LINK->secondseg - 1];
    s2 = seghap2[LINK->secondseg - 1];
    if (s1 != s2) {
      g1 = WITH->genenumber[LINK->fseg - 1][s1 - 1];
      g2 = WITH->genenumber[LINK->fseg - 1][s2 - 1];
      FORLIM1 = nchild;
      for (k = 0; k < FORLIM1; k++) {
	WITH2 = thischild[k];
	if (malechild[k])
	  val = WITH2->genarray[s1 - 1] + WITH2->genarray[s2 - 1];
	else
	  val = WITH2->genarray[g1 - 1] + WITH2->genarray[g2 - 1];
	for (l = 0; l < slength; l++)
	  temp[k][l] += temp2[l] * val;
      }
    } else {
      g1 = WITH->genenumber[LINK->fseg - 1][s1 - 1];
      FORLIM1 = nchild;
      for (k = 0; k < FORLIM1; k++) {
	WITH2 = thischild[k];
	if (malechild[k])
	  val = 2.0 * WITH2->genarray[s1 - 1];
	else
	  val = 2.0 * WITH2->genarray[g1 - 1];
	for (l = 0; l < slength; l++)
	  temp[k][l] += temp2[l] * val;
      }
    }
    LINK->secondseg++;
  }
  temp1 = 0.0;
  for (l = 0; l < slength; l++) {
    val = 1.0;
    FORLIM1 = nchild;
    for (k = 0; k < FORLIM1; k++)
      val *= temp[k][l];
    temp1 += val;
  }
  return temp1;
}  /* segsexf */


Local double segfun(LINK)
struct LOC_seg *LINK;
{
  long g1, g2, g3, g4, i, j, k, f1, f2, s1, s2;
  double val, temp1, temp2;
  double temp[maxchild];
  long FORLIM;
  gennurec *WITH;
  thetavalues *WITH1, *WITH2;
  long FORLIM1, FORLIM2;
  thisarray *WITH3;

  LINK->firstseg = LINK->fseg;
  FORLIM = nchild;
  for (k = 0; k < FORLIM; k++)
    temp[k] = 0.0;
  WITH = gennustruct;
  WITH1 = LINK->firstsex;
  FORLIM = LINK->fend;
  for (i = LINK->fstart - 1; i < FORLIM; i++) {
    temp1 = WITH1->segprob[i];
    f1 = seghap1[LINK->firstseg - 1];
    f2 = seghap2[LINK->firstseg - 1];
    LINK->secondseg = LINK->sseg;
    WITH2 = LINK->secondsex;
    if (f1 != f2) {
      FORLIM1 = LINK->send;
      for (j = LINK->sstart - 1; j < FORLIM1; j++) {
	temp2 = temp1 * WITH2->segprob[j];
	s1 = seghap1[LINK->secondseg - 1];
	s2 = seghap2[LINK->secondseg - 1];
	if (s1 != s2) {
	  g1 = WITH->genenumber[f1 - 1][s1 - 1];
	  g2 = WITH->genenumber[f1 - 1][s2 - 1];
	  g3 = WITH->genenumber[f2 - 1][s1 - 1];
	  g4 = WITH->genenumber[f2 - 1][s2 - 1];
	  FORLIM2 = nchild;
	  for (k = 0; k < FORLIM2; k++) {
	    WITH3 = thischild[k];
	    temp[k] += temp2 *
		       (WITH3->genarray[g1 - 1] + WITH3->genarray[g2 - 1] +
			WITH3->genarray[g3 - 1] + WITH3->genarray[g4 - 1]);
	  }
	} else {
	  g1 = WITH->genenumber[f1 - 1][s1 - 1];
	  g3 = WITH->genenumber[f2 - 1][s1 - 1];
	  FORLIM2 = nchild;
	  for (k = 0; k < FORLIM2; k++) {
	    WITH3 = thischild[k];
	    temp[k] += 2 * temp2 *
		       (WITH3->genarray[g1 - 1] + WITH3->genarray[g3 - 1]);
	  }
	}
	LINK->secondseg++;
      }
    } else {
      FORLIM1 = LINK->send;
      for (j = LINK->sstart - 1; j < FORLIM1; j++) {
	temp2 = temp1 * WITH2->segprob[j];
	s1 = seghap1[LINK->secondseg - 1];
	s2 = seghap2[LINK->secondseg - 1];
	if (s1 != s2) {
	  g1 = WITH->genenumber[f1 - 1][s1 - 1];
	  g2 = WITH->genenumber[f1 - 1][s2 - 1];
	  FORLIM2 = nchild;
	  for (k = 0; k < FORLIM2; k++) {
	    WITH3 = thischild[k];
	    temp[k] += 2 * temp2 *
		       (WITH3->genarray[g1 - 1] + WITH3->genarray[g2 - 1]);
	  }
	} else {
	  g1 = WITH->genenumber[f1 - 1][s1 - 1];
	  FORLIM2 = nchild;
	  for (k = 0; k < FORLIM2; k++) {
	    WITH3 = thischild[k];
	    temp[k] += 4 * temp2 * WITH3->genarray[g1 - 1];
	  }
	}
	LINK->secondseg++;
      }
    }
    LINK->firstseg++;
  }
  val = 1.0;
  FORLIM = nchild;
  for (k = 0; k < FORLIM; k++)
    val *= temp[k];
  return val;
}  /* segfun */


Local double msegfast(LINK)
struct LOC_seg *LINK;
{
  long g1, g2, g3, g4, g5, g6, g7, g8, g9, g10, g11, g12, g13, g14, g15, g16,
       i, j, k, l, f1, f2, s1, s2, ms1, ms2, mf1, mf2, slength;
  double val, temp1;
  double temp[maxchild][maxseg];
  double temp2[maxseg];   /*msegfast*/
  long FORLIM, FORLIM1;
  gennurec *WITH;
  thetavalues *WITH1, *WITH2;
  long FORLIM2;
  thisarray *WITH3;

  LINK->firstseg = LINK->fseg;
  slength = LINK->send - LINK->sstart + 1;
  FORLIM = nchild;
  for (k = 0; k < FORLIM; k++) {
    for (l = 0; l < slength; l++)
      temp[k][l] = 0.0;
  }
  WITH = gennustruct;
  WITH1 = LINK->firstsex;
  FORLIM = LINK->fend;
  for (i = LINK->fstart - 1; i < FORLIM; i++) {
    temp1 = WITH1->segprob[i];
    f1 = seghap1[LINK->firstseg - 1];
    f2 = seghap2[LINK->firstseg - 1];
    mf1 = muthap[f1 - 1];
    mf2 = muthap[f2 - 1];
    LINK->secondseg = LINK->sseg;
    WITH2 = LINK->secondsex;
    if (f1 != f2) {
      FORLIM1 = LINK->send;
      for (j = LINK->sstart - 1; j < FORLIM1; j++) {
	for (l = 0; l < slength; l++)
	  temp2[l] = temp1 * WITH2->segprob[j + l * slength];
	s1 = seghap1[LINK->secondseg - 1];
	s2 = seghap2[LINK->secondseg - 1];
	ms1 = muthap[s1 - 1];
	ms2 = muthap[s2 - 1];
	if (s1 != s2) {
	  g1 = WITH->genenumber[f1 - 1][s1 - 1];
	  g2 = WITH->genenumber[f1 - 1][s2 - 1];
	  g3 = WITH->genenumber[f2 - 1][s1 - 1];
	  g4 = WITH->genenumber[f2 - 1][s2 - 1];
	  g5 = WITH->genenumber[f1 - 1][ms1 - 1];
	  g6 = WITH->genenumber[f1 - 1][ms2 - 1];
	  g7 = WITH->genenumber[f2 - 1][ms1 - 1];
	  g8 = WITH->genenumber[f2 - 1][ms2 - 1];
	  g9 = WITH->genenumber[mf1 - 1][s1 - 1];
	  g10 = WITH->genenumber[mf1 - 1][s2 - 1];
	  g11 = WITH->genenumber[mf2 - 1][s1 - 1];
	  g12 = WITH->genenumber[mf2 - 1][s2 - 1];
	  g13 = WITH->genenumber[mf1 - 1][ms1 - 1];
	  g14 = WITH->genenumber[mf1 - 1][ms2 - 1];
	  g15 = WITH->genenumber[mf2 - 1][ms1 - 1];
	  g16 = WITH->genenumber[mf2 - 1][ms2 - 1];
	  FORLIM2 = nchild;
	  for (k = 0; k < FORLIM2; k++) {
	    WITH3 = thischild[k];
	    val = (1 - LINK->ps) * (1 - LINK->pf) * (WITH3->genarray[g1 - 1] +
		    WITH3->genarray[g2 - 1] + WITH3->genarray[g3 - 1] +
		    WITH3->genarray[g4 - 1]) + LINK->ps * (1 - LINK->pf) *
		  (WITH3->genarray[g5 - 1] + WITH3->genarray[g6 - 1] + WITH3->
		      genarray[g7 - 1] + WITH3->genarray[g8 - 1]) +
		LINK->pf * (1 - LINK->ps) * (WITH3->genarray[g9 - 1] +
		    WITH3->genarray[g10 - 1] + WITH3->genarray[g11 - 1] +
		    WITH3->genarray[g12 - 1]) + LINK->pf * LINK->ps *
		  (WITH3->genarray[g13 - 1] + WITH3->genarray[g14 - 1] +
		    WITH3->genarray[g15 - 1] + WITH3->genarray[g16 - 1]);
/* p2c: msim.p, line 4352: 
 * Note: Line breaker spent 0.0 seconds, 5000 tries on line 4063 [251] */
	    for (l = 0; l < slength; l++)
	      temp[k][l] += temp2[l] * val;
	  }
	} else {
	  g1 = WITH->genenumber[f1 - 1][s1 - 1];
	  g3 = WITH->genenumber[f2 - 1][s1 - 1];
	  g5 = WITH->genenumber[f1 - 1][ms1 - 1];
	  g7 = WITH->genenumber[f2 - 1][ms1 - 1];
	  g9 = WITH->genenumber[mf1 - 1][s1 - 1];
	  g11 = WITH->genenumber[mf2 - 1][s1 - 1];
	  g13 = WITH->genenumber[mf1 - 1][ms1 - 1];
	  g15 = WITH->genenumber[mf2 - 1][ms1 - 1];
	  FORLIM2 = nchild;
	  for (k = 0; k < FORLIM2; k++) {
	    WITH3 = thischild[k];
	    val = 2 * ((1 - LINK->ps) * (1 - LINK->pf) *
		    (WITH3->genarray[g1 - 1] + WITH3->genarray[g3 - 1]) +
		  LINK->ps * (1 - LINK->pf) * (WITH3->genarray[g5 - 1] +
		      WITH3->genarray[g7 - 1]) + LINK->pf * (1 - LINK->ps) *
		    (WITH3->genarray[g9 - 1] + WITH3->genarray[g11 - 1]) +
		  LINK->pf * LINK->ps *
		    (WITH3->genarray[g13 - 1] + WITH3->genarray[g15 - 1]));
/* p2c: msim.p, line 4352: 
 * Note: Line breaker spent 1.0 seconds, 5000 tries on line 4087 [251] */
	    for (l = 0; l < slength; l++)
	      temp[k][l] += temp2[l] * val;
	  }
	}
	LINK->secondseg++;
      }
    } else {
      FORLIM1 = LINK->send;
      for (j = LINK->sstart - 1; j < FORLIM1; j++) {
	for (l = 0; l < slength; l++)
	  temp2[l] = temp1 * WITH2->segprob[j + l * slength];
	s1 = seghap1[LINK->secondseg - 1];
	s2 = seghap2[LINK->secondseg - 1];
	ms1 = muthap[s1 - 1];
	ms2 = muthap[s2 - 1];
	if (s1 != s2) {
	  g1 = WITH->genenumber[f1 - 1][s1 - 1];
	  g2 = WITH->genenumber[f1 - 1][s2 - 1];
	  g5 = WITH->genenumber[f1 - 1][ms1 - 1];
	  g6 = WITH->genenumber[f1 - 1][ms2 - 1];
	  g9 = WITH->genenumber[mf1 - 1][s1 - 1];
	  g10 = WITH->genenumber[mf1 - 1][s2 - 1];
	  g13 = WITH->genenumber[mf1 - 1][ms1 - 1];
	  g14 = WITH->genenumber[mf1 - 1][ms2 - 1];
	  FORLIM2 = nchild;
	  for (k = 0; k < FORLIM2; k++) {
	    WITH3 = thischild[k];
	    val = 2 * ((1 - LINK->ps) * (1 - LINK->pf) *
		    (WITH3->genarray[g1 - 1] + WITH3->genarray[g2 - 1]) +
		  LINK->ps * (1 - LINK->pf) * (WITH3->genarray[g5 - 1] +
		      WITH3->genarray[g6 - 1]) + LINK->pf * (1 - LINK->ps) *
		    (WITH3->genarray[g9 - 1] + WITH3->genarray[g10 - 1]) +
		  LINK->pf * LINK->ps *
		    (WITH3->genarray[g13 - 1] + WITH3->genarray[g14 - 1]));
/* p2c: msim.p, line 4352: 
 * Note: Line breaker spent 0.0 seconds, 5000 tries on line 4123 [251] */
	    for (l = 0; l < slength; l++)
	      temp[k][l] += temp2[l] * val;
	  }
	} else {
	  g1 = WITH->genenumber[f1 - 1][s1 - 1];
	  g5 = WITH->genenumber[f1 - 1][ms1 - 1];
	  g9 = WITH->genenumber[mf1 - 1][s1 - 1];
	  g13 = WITH->genenumber[mf1 - 1][ms1 - 1];
	  FORLIM2 = nchild;
	  for (k = 0; k < FORLIM2; k++) {
	    WITH3 = thischild[k];
	    val = 4 * ((1 - LINK->ps) * (1 - LINK->pf) * WITH3->genarray[g1 - 1] +
		       LINK->ps * (1 - LINK->pf) * WITH3->genarray[g5 - 1] +
		       LINK->pf * (1 - LINK->ps) * WITH3->genarray[g9 - 1] +
		       LINK->pf * LINK->ps * WITH3->genarray[g13 - 1]);
	    for (l = 0; l < slength; l++)
	      temp[k][l] += temp2[l] * val;
	  }
	}
	LINK->secondseg++;
      }
    }
    LINK->firstseg++;
  }
  temp1 = 0.0;
  for (l = 0; l < slength; l++) {
    val = 1.0;
    FORLIM1 = nchild;
    for (k = 0; k < FORLIM1; k++)
      val *= temp[k][l];
    temp1 += val;
  }
  return temp1;
}  /* msegfast */

/*msegfast*/

Local double msegfun(LINK)
struct LOC_seg *LINK;
{
  long g1, g2, g3, g4, g5, g6, g7, g8, g9, g10, g11, g12, g13, g14, g15, g16,
       i, j, k, f1, f2, s1, s2, ms1, ms2, mf1, mf2;
  double val, temp1, temp2;
  double temp[maxchild];   /*msegfun*/
  long FORLIM;
  gennurec *WITH;
  thetavalues *WITH1, *WITH2;
  long FORLIM1, FORLIM2;
  thisarray *WITH3;

  LINK->firstseg = LINK->fseg;
  FORLIM = nchild;
  for (k = 0; k < FORLIM; k++)
    temp[k] = 0.0;
  WITH = gennustruct;
  WITH1 = LINK->firstsex;
  FORLIM = LINK->fend;
  for (i = LINK->fstart - 1; i < FORLIM; i++) {
    temp1 = WITH1->segprob[i];
    f1 = seghap1[LINK->firstseg - 1];
    f2 = seghap2[LINK->firstseg - 1];
    mf1 = muthap[f1 - 1];
    mf2 = muthap[f2 - 1];
    LINK->secondseg = LINK->sseg;
    WITH2 = LINK->secondsex;
    if (f1 != f2) {
      FORLIM1 = LINK->send;
      for (j = LINK->sstart - 1; j < FORLIM1; j++) {
	temp2 = temp1 * WITH2->segprob[j];
	s1 = seghap1[LINK->secondseg - 1];
	s2 = seghap2[LINK->secondseg - 1];
	ms1 = muthap[s1 - 1];
	ms2 = muthap[s2 - 1];
	if (s1 != s2) {
	  g1 = WITH->genenumber[f1 - 1][s1 - 1];
	  g2 = WITH->genenumber[f1 - 1][s2 - 1];
	  g3 = WITH->genenumber[f2 - 1][s1 - 1];
	  g4 = WITH->genenumber[f2 - 1][s2 - 1];
	  g5 = WITH->genenumber[f1 - 1][ms1 - 1];
	  g6 = WITH->genenumber[f1 - 1][ms2 - 1];
	  g7 = WITH->genenumber[f2 - 1][ms1 - 1];
	  g8 = WITH->genenumber[f2 - 1][ms2 - 1];
	  g9 = WITH->genenumber[mf1 - 1][s1 - 1];
	  g10 = WITH->genenumber[mf1 - 1][s2 - 1];
	  g11 = WITH->genenumber[mf2 - 1][s1 - 1];
	  g12 = WITH->genenumber[mf2 - 1][s2 - 1];
	  g13 = WITH->genenumber[mf1 - 1][ms1 - 1];
	  g14 = WITH->genenumber[mf1 - 1][ms2 - 1];
	  g15 = WITH->genenumber[mf2 - 1][ms1 - 1];
	  g16 = WITH->genenumber[mf2 - 1][ms2 - 1];
	  FORLIM2 = nchild;
	  for (k = 0; k < FORLIM2; k++) {
	    WITH3 = thischild[k];
	    temp[k] += temp2 * ((1 - LINK->ps) * (1 - LINK->pf) * (WITH3->
			genarray[g1 - 1] + WITH3->genarray[g2 - 1] + WITH3->
			genarray[g3 - 1] + WITH3->genarray[g4 - 1]) +
		  LINK->ps * (1 - LINK->pf) * (WITH3->genarray[g5 - 1] +
		      WITH3->genarray[g6 - 1] + WITH3->genarray[g7 - 1] +
		      WITH3->genarray[g8 - 1]) + LINK->pf * (1 - LINK->ps) *
		    (WITH3->genarray[g9 - 1] + WITH3->genarray[g10 - 1] +
		      WITH3->genarray[g11 - 1] + WITH3->genarray[g12 - 1]) +
		  LINK->pf * LINK->ps * (WITH3->genarray[g13 - 1] + WITH3->
			genarray[g14 - 1] + WITH3->genarray[g15 - 1] + WITH3->
			genarray[g16 - 1]));
/* p2c: msim.p, line 4352: 
 * Note: Line breaker spent 0.0 seconds, 5000 tries on line 4229 [251] */
	  }
	} else {
	  g1 = WITH->genenumber[f1 - 1][s1 - 1];
	  g3 = WITH->genenumber[f2 - 1][s1 - 1];
	  g5 = WITH->genenumber[f1 - 1][ms1 - 1];
	  g7 = WITH->genenumber[f2 - 1][ms1 - 1];
	  g9 = WITH->genenumber[mf1 - 1][s1 - 1];
	  g11 = WITH->genenumber[mf2 - 1][s1 - 1];
	  g13 = WITH->genenumber[mf1 - 1][ms1 - 1];
	  g15 = WITH->genenumber[mf2 - 1][ms1 - 1];
	  FORLIM2 = nchild;
	  for (k = 0; k < FORLIM2; k++) {
	    WITH3 = thischild[k];
	    temp[k] += 2 * temp2 * ((1 - LINK->ps) * (1 - LINK->pf) *
		    (WITH3->genarray[g1 - 1] + WITH3->genarray[g3 - 1]) +
		  LINK->ps * (1 - LINK->pf) * (WITH3->genarray[g5 - 1] +
		      WITH3->genarray[g7 - 1]) + LINK->pf * (1 - LINK->ps) *
		    (WITH3->genarray[g9 - 1] + WITH3->genarray[g11 - 1]) +
		  LINK->pf * LINK->ps *
		    (WITH3->genarray[g13 - 1] + WITH3->genarray[g15 - 1]));
/* p2c: msim.p, line 4352: 
 * Note: Line breaker spent 0.0 seconds, 5000 tries on line 4251 [251] */
	  }
	}
	LINK->secondseg++;
      }
    } else {
      FORLIM1 = LINK->send;
      for (j = LINK->sstart - 1; j < FORLIM1; j++) {
	temp2 = temp1 * WITH2->segprob[j];
	s1 = seghap1[LINK->secondseg - 1];
	s2 = seghap2[LINK->secondseg - 1];
	ms1 = muthap[s1 - 1];
	ms2 = muthap[s2 - 1];
	if (s1 != s2) {
	  g1 = WITH->genenumber[f1 - 1][s1 - 1];
	  g2 = WITH->genenumber[f1 - 1][s2 - 1];
	  g5 = WITH->genenumber[f1 - 1][ms1 - 1];
	  g6 = WITH->genenumber[f1 - 1][ms2 - 1];
	  g9 = WITH->genenumber[mf1 - 1][s1 - 1];
	  g10 = WITH->genenumber[mf1 - 1][s2 - 1];
	  g13 = WITH->genenumber[mf1 - 1][ms1 - 1];
	  g14 = WITH->genenumber[mf1 - 1][ms2 - 1];
	  FORLIM2 = nchild;
	  for (k = 0; k < FORLIM2; k++) {
	    WITH3 = thischild[k];
	    temp[k] += 2 * temp2 * ((1 - LINK->ps) * (1 - LINK->pf) *
		    (WITH3->genarray[g1 - 1] + WITH3->genarray[g2 - 1]) +
		  LINK->ps * (1 - LINK->pf) * (WITH3->genarray[g5 - 1] +
		      WITH3->genarray[g6 - 1]) + LINK->pf * (1 - LINK->ps) *
		    (WITH3->genarray[g9 - 1] + WITH3->genarray[g10 - 1]) +
		  LINK->pf * LINK->ps *
		    (WITH3->genarray[g13 - 1] + WITH3->genarray[g14 - 1]));
/* p2c: msim.p, line 4352: 
 * Note: Line breaker spent 0.0 seconds, 5000 tries on line 4284 [251] */
	  }
	} else {
	  g1 = WITH->genenumber[f1 - 1][s1 - 1];
	  g5 = WITH->genenumber[f1 - 1][ms1 - 1];
	  g9 = WITH->genenumber[mf1 - 1][s1 - 1];
	  g13 = WITH->genenumber[mf1 - 1][ms1 - 1];
	  FORLIM2 = nchild;
	  for (k = 0; k < FORLIM2; k++) {
	    WITH3 = thischild[k];
	    temp[k] += 4 * temp2 *
		((1 - LINK->ps) * (1 - LINK->pf) * WITH3->genarray[g1 - 1] +
		 LINK->ps * (1 - LINK->pf) * WITH3->genarray[g5 - 1] +
		 LINK->pf * (1 - LINK->ps) * WITH3->genarray[g9 - 1] +
		 LINK->pf * LINK->ps * WITH3->genarray[g13 - 1]);
	  }
	}
	LINK->secondseg++;
      }
    }
    LINK->firstseg++;
  }
  val = 1.0;
  FORLIM = nchild;
  for (k = 0; k < FORLIM; k++)
    val *= temp[k];
  return val;
}  /* msegfun */

/*msegfun*/

Local double segfast(LINK)
struct LOC_seg *LINK;
{
  long g1, g2, g3, g4, i, j, k, l, f1, f2, s1, s2, slength;
  double val, temp1;
  double temp[maxchild][maxseg];
  double temp2[maxseg];   /*segfast*/
  long FORLIM, FORLIM1;
  gennurec *WITH;
  thetavalues *WITH1, *WITH2;
  long FORLIM2;
  thisarray *WITH3;

  LINK->firstseg = LINK->fseg;
  slength = LINK->send - LINK->sstart + 1;
  FORLIM = nchild;
  for (k = 0; k < FORLIM; k++) {
    for (l = 0; l < slength; l++)
      temp[k][l] = 0.0;
  }
  WITH = gennustruct;
  WITH1 = LINK->firstsex;
  FORLIM = LINK->fend;
  for (i = LINK->fstart - 1; i < FORLIM; i++) {
    temp1 = WITH1->segprob[i];
    f1 = seghap1[LINK->firstseg - 1];
    f2 = seghap2[LINK->firstseg - 1];
    LINK->secondseg = LINK->sseg;
    WITH2 = LINK->secondsex;
    if (f1 != f2) {
      FORLIM1 = LINK->send;
      for (j = LINK->sstart - 1; j < FORLIM1; j++) {
	for (l = 0; l < slength; l++)
	  temp2[l] = temp1 * WITH2->segprob[j + l * slength];
	s1 = seghap1[LINK->secondseg - 1];
	s2 = seghap2[LINK->secondseg - 1];
	if (s1 != s2) {
	  g1 = WITH->genenumber[f1 - 1][s1 - 1];
	  g2 = WITH->genenumber[f1 - 1][s2 - 1];
	  g3 = WITH->genenumber[f2 - 1][s1 - 1];
	  g4 = WITH->genenumber[f2 - 1][s2 - 1];
	  FORLIM2 = nchild;
	  for (k = 0; k < FORLIM2; k++) {
	    WITH3 = thischild[k];
	    val = WITH3->genarray[g1 - 1] + WITH3->genarray[g2 - 1] +
		  WITH3->genarray[g3 - 1] + WITH3->genarray[g4 - 1];
	    for (l = 0; l < slength; l++)
	      temp[k][l] += temp2[l] * val;
	  }
	} else {
	  g1 = WITH->genenumber[f1 - 1][s1 - 1];
	  g3 = WITH->genenumber[f2 - 1][s1 - 1];
	  FORLIM2 = nchild;
	  for (k = 0; k < FORLIM2; k++) {
	    WITH3 = thischild[k];
	    val = 2.0 * (WITH3->genarray[g1 - 1] + WITH3->genarray[g3 - 1]);
	    for (l = 0; l < slength; l++)
	      temp[k][l] += temp2[l] * val;
	  }
	}
	LINK->secondseg++;
      }
    } else {
      FORLIM1 = LINK->send;
      for (j = LINK->sstart - 1; j < FORLIM1; j++) {
	for (l = 0; l < slength; l++)
	  temp2[l] = temp1 * WITH2->segprob[j + l * slength];
	s1 = seghap1[LINK->secondseg - 1];
	s2 = seghap2[LINK->secondseg - 1];
	if (s1 != s2) {
	  g1 = WITH->genenumber[f1 - 1][s1 - 1];
	  g2 = WITH->genenumber[f1 - 1][s2 - 1];
	  FORLIM2 = nchild;
	  for (k = 0; k < FORLIM2; k++) {
	    WITH3 = thischild[k];
	    val = 2.0 * (WITH3->genarray[g1 - 1] + WITH3->genarray[g2 - 1]);
	    for (l = 0; l < slength; l++)
	      temp[k][l] += temp2[l] * val;
	  }
	} else {
	  g1 = WITH->genenumber[f1 - 1][s1 - 1];
	  FORLIM2 = nchild;
	  for (k = 0; k < FORLIM2; k++) {
	    WITH3 = thischild[k];
	    val = 4.0 * WITH3->genarray[g1 - 1];
	    for (l = 0; l < slength; l++)
	      temp[k][l] += temp2[l] * val;
	  }
	}
	LINK->secondseg++;
      }
    }
    LINK->firstseg++;
  }
  temp1 = 0.0;
  for (l = 0; l < slength; l++) {
    val = 1.0;
    FORLIM1 = nchild;
    for (k = 0; k < FORLIM1; k++)
      val *= temp[k][l];
    temp1 += val;
  }
  return temp1;
}  /* segfast */

/*segfast*/

Local Void initseg(LINK)
struct LOC_seg *LINK;
{
  /*initseg*/
  if ((*LINK->p)->male) {
    LINK->nfirst = mgeno;
    LINK->nsecond = fgeno;
    LINK->firstsex = maletheta;
    LINK->secondsex = femaletheta;
    LINK->pf = mutmale;
    LINK->ps = mutfemale;
  } else {
    LINK->nfirst = fgeno;
    LINK->nsecond = mgeno;
    LINK->firstsex = femaletheta;
    LINK->secondsex = maletheta;
    LINK->pf = mutfemale;
    LINK->ps = mutmale;
  }
  prob(&LINK->father, LINK->LINK);
  prob(&LINK->mother, LINK->LINK);
  LINK->child = LINK->father->foff;
  nchild = 0;
  do {
    prob(&LINK->child, LINK->LINK);
    if (LINK->child->ma == LINK->mother && !LINK->child->up) {
      nchild++;
      thischild[nchild - 1] = LINK->child->gen;
      malechild[nchild - 1] = LINK->child->male;
    }
    LINK->child = LINK->child->nextpa;
  } while (LINK->child != NULL);   /* initseg */
}

/*initseg*/

Local Void exitseg(LINK)
struct LOC_seg *LINK;
{
  /*exitseg*/
  LINK->LINK->nuscale++;
  LINK->child = LINK->father->foff;
  do {
    if (LINK->child->ma == LINK->mother && !LINK->child->up)
      cleanup(&LINK->child);
    LINK->child = LINK->child->nextpa;
  } while (LINK->child != NULL);   /* exitseg */
}

/*exitseg*/

Local Void segsexctop(LINK)
struct LOC_seg *LINK;
{
  double segval;
  long first, second, sstop;
  boolean thiscensor, thisrare;   /*segsexctop*/
  censorrec *WITH;
  thisarray *WITH1;
  long FORLIM;
  thisarray *WITH2;

  initseg(LINK);
  WITH = censorstruct;
  if ((*LINK->p)->male) {
    WITH1 = (*LINK->p)->gen;
    FORLIM = LINK->nfirst;
    for (first = 0; first < FORLIM; first++) {
      if (WITH1->genarray[first] != 0.0) {
	segval = 0.0;
	LINK->fseg = first + 1;
	thisrare = rare[first];
	second = 1;
	WITH2 = (*LINK->q)->gen;
	do {
	  LINK->sstart = probstart[second - 1];
	  LINK->send = probend[second - 1];
	  LINK->sseg = segstart[second - 1];
	  sstop = second + LINK->send - LINK->sstart + 1;
	  if (thisc < maxcensor) {
	    thisc++;
	    thiscensor = WITH->censor[thisc - minint];
	  } else
	    thiscensor = (WITH2->genarray[second - 1] == 0.0 ||
			  thisrare && rare[second - 1]);
	  if (!thiscensor) {
	    if (mutsys != 0)
	      segval += WITH2->genarray[second - 1] * msegsexf(LINK);
	    else
	      segval += WITH2->genarray[second - 1] * segsexf(LINK);
	  }
	  second = sstop;
	} while (second <= LINK->nsecond);
	WITH1->genarray[first] *= segval * segscale;
      }
    }
  } else {
    WITH1 = (*LINK->p)->gen;
    FORLIM = LINK->nfirst;
    for (first = 0; first < FORLIM; first++) {
      if (WITH1->genarray[first] != 0.0) {
	segval = 0.0;
	LINK->fstart = probstart[first];
	LINK->fend = probend[first];
	LINK->fseg = segstart[first];
	thisrare = rare[first];
	second = 1;
	WITH2 = (*LINK->q)->gen;
	do {
	  LINK->sseg = second;
	  if (thisc < maxcensor) {
	    thisc++;
	    thiscensor = WITH->censor[thisc - minint];
	  } else
	    thiscensor = (WITH2->genarray[second - 1] == 0.0 ||
			  thisrare && rare[second - 1]);
	  if (!thiscensor) {
	    if (mutsys != 0)
	      segval += WITH2->genarray[second - 1] * msegsex(LINK);
	    else
	      segval += WITH2->genarray[second - 1] * segsex(LINK);
	  }
	  second++;
	} while (second <= LINK->nsecond);
	WITH1->genarray[first] *= segval * segscale;
      }
    }
  }
  cleanup(LINK->q);
  exitseg(LINK);
}  /* segsexctop */

/*segsexctop*/

Local Void segsextop(LINK)
struct LOC_seg *LINK;
{
  double segval, val;
  long first, second, sstop;   /*segsextop*/
  censorrec *WITH;
  thisarray *WITH1;
  long FORLIM;
  thisarray *WITH2;

  initseg(LINK);
  WITH = censorstruct;
  if ((*LINK->p)->male) {
    WITH1 = (*LINK->p)->gen;
    FORLIM = LINK->nfirst;
    for (first = 1; first <= FORLIM; first++) {
      if (WITH1->genarray[first - 1] != 0.0) {
	segval = 0.0;
	LINK->fseg = first;
	second = 1;
	WITH2 = (*LINK->q)->gen;
	do {
	  LINK->sstart = probstart[second - 1];
	  LINK->send = probend[second - 1];
	  LINK->sseg = segstart[second - 1];
	  sstop = second + LINK->send - LINK->sstart + 1;
	  if (thisc <= maxcensor)
	    thisc++;
	  if (WITH2->genarray[second - 1] == 0.0) {
	    second = sstop;
	    if (thisc <= maxcensor)
	      WITH->censor[thisc - minint] = true;
	  } else {
	    if (mutsys != 0)
	      val = msegsexf(LINK);
	    else
	      val = segsexf(LINK);
	    if (thisc <= maxcensor)
	      WITH->censor[thisc - minint] = (val == 0.0);
	    segval += WITH2->genarray[second - 1] * val;
	    second = sstop;
	  }
	} while (second <= LINK->nsecond);
	WITH1->genarray[first - 1] *= segval * segscale;
      }
    }
  } else {
    WITH1 = (*LINK->p)->gen;
    FORLIM = LINK->nfirst;
    for (first = 0; first < FORLIM; first++) {
      if (WITH1->genarray[first] != 0.0) {
	segval = 0.0;
	LINK->fstart = probstart[first];
	LINK->fend = probend[first];
	LINK->fseg = segstart[first];
	second = 1;
	WITH2 = (*LINK->q)->gen;
	do {
	  LINK->sseg = second;
	  if (thisc <= maxcensor)
	    thisc++;
	  if (WITH2->genarray[second - 1] == 0.0) {
	    second++;
	    if (thisc <= maxcensor)
	      WITH->censor[thisc - minint] = true;
	  } else {
	    if (mutsys != 0)
	      val = msegsex(LINK);
	    else
	      val = segsex(LINK);
	    if (thisc <= maxcensor)
	      WITH->censor[thisc - minint] = (val == 0.0);
	    segval += WITH2->genarray[second - 1] * val;
	    second++;
	  }
	} while (second <= LINK->nsecond);
	WITH1->genarray[first] *= segval * segscale;
      }
    }
  }
  cleanup(LINK->q);
  exitseg(LINK);
}  /* segsextop */

/*segsextop*/

Local Void segsexup(LINK)
struct LOC_seg *LINK;
{
  double segval;
  long first, second;   /*segsexup*/
  censorrec *WITH;
  thisarray *WITH1;
  long FORLIM;
  thisarray *WITH2;
  long FORLIM1;

  initseg(LINK);
  WITH = censorstruct;
  if ((*LINK->p)->male) {
    WITH1 = (*LINK->p)->gen;
    FORLIM = LINK->nfirst;
    for (first = 1; first <= FORLIM; first++) {
      if (WITH1->genarray[first - 1] != 0.0) {
	segval = 0.0;
	LINK->fseg = first;
	second = 1;
	WITH2 = (*LINK->q)->gen;
	FORLIM1 = LINK->nsecond;
	for (second = 0; second < FORLIM1; second++) {
	  if (WITH2->genarray[second] != 0.0) {
	    LINK->sstart = probstart[second];
	    LINK->send = probend[second];
	    LINK->sseg = segstart[second];
	    if (mutsys != 0)
	      segval += WITH2->genarray[second] * msegsex(LINK);
	    else
	      segval += WITH2->genarray[second] * segsex(LINK);
	  }
	}
	WITH1->genarray[first - 1] *= segval * segscale;
      }
    }
  } else {
    WITH1 = (*LINK->p)->gen;
    FORLIM = LINK->nfirst;
    for (first = 0; first < FORLIM; first++) {
      if (WITH1->genarray[first] != 0.0) {
	segval = 0.0;
	LINK->fstart = probstart[first];
	LINK->fend = probend[first];
	LINK->fseg = segstart[first];
	WITH2 = (*LINK->q)->gen;
	FORLIM1 = LINK->nsecond;
	for (second = 1; second <= FORLIM1; second++) {
	  if (WITH2->genarray[second - 1] != 0.0) {
	    LINK->sseg = second;
	    if (mutsys != 0)
	      segval += WITH2->genarray[second - 1] * msegsex(LINK);
	    else
	      segval += WITH2->genarray[second - 1] * segsex(LINK);
	  }
	}
	WITH1->genarray[first] *= segval * segscale;
      }
    }
  }
  cleanup(LINK->q);
  exitseg(LINK);
}  /* segsexup */

/*segsexup*/

Local Void segsexdown(LINK)
struct LOC_seg *LINK;
{
  short here;
  genotype gene;
  double val, temp2;
  long j, first, second;   /*segsexdown*/
  short FORLIM;
  gennurec *WITH;
  censorrec *WITH1;
  thisarray *WITH2;
  long FORLIM1;
  thisarray *WITH3;
  long FORLIM2;
  thetavalues *WITH4;
  long FORLIM3;

  initseg(LINK);
  FORLIM = fgeno;
  for (here = 0; here < FORLIM; here++)
    gene[here] = 0.0;
  WITH = gennustruct;
  WITH1 = censorstruct;
  WITH2 = (*LINK->p)->gen;
  FORLIM1 = LINK->nfirst;
  for (first = 0; first < FORLIM1; first++) {
    if (WITH2->genarray[first] != 0.0) {
      LINK->fseg = first + 1;
      second = 1;
      WITH3 = (*LINK->q)->gen;
      FORLIM2 = LINK->nsecond;
      for (second = 0; second < FORLIM2; second++) {
	if (WITH3->genarray[second] != 0.0) {
	  val = WITH3->genarray[second] * (*LINK->p)->gen->genarray[first] *
		segscale;
	  LINK->sstart = probstart[second];
	  LINK->send = probend[second];
	  LINK->sseg = segstart[second];
	  if (nchild != 0)
	    val *= segsex(LINK);
	  if (val != 0.0) {
	    LINK->secondseg = LINK->sseg;
	    WITH4 = femaletheta;
	    FORLIM3 = LINK->send;
	    for (j = LINK->sstart - 1; j < FORLIM3; j++) {
	      temp2 = WITH4->segprob[j];
	      if ((*LINK->r)->male) {
		here = seghap1[LINK->secondseg - 1];
		gene[here - 1] += temp2 * val;
		here = seghap2[LINK->secondseg - 1];
		gene[here - 1] += temp2 * val;
	      } else {
		here = WITH->genenumber[seghap1[LINK->secondseg - 1] - 1]
		  [first];
		gene[here - 1] += temp2 * val;
		here = WITH->genenumber[seghap2[LINK->secondseg - 1] - 1]
		  [first];
		gene[here - 1] += temp2 * val;
	      }
	      LINK->secondseg++;
	    }
	  }  /*second*/
	}
      }
    }  /*first*/
  }
  memcpy((*LINK->p)->gen->genarray, (*LINK->r)->gen->genarray,
	 sizeof(genotype));
  memcpy((*LINK->r)->gen->genarray, gene, sizeof(genotype));
  WITH2 = (*LINK->r)->gen;
  FORLIM1 = fgeno;
  for (first = 0; first < FORLIM1; first++)
    WITH2->genarray[first] *= (*LINK->p)->gen->genarray[first];
  cleanup(LINK->p);
  cleanup(LINK->q);
  exitseg(LINK);
}  /* segsexdown */

/*segsexdown*/

Local Void segctop(LINK)
struct LOC_seg *LINK;
{
  double segval;
  long first, second, sstop;
  boolean thiscensor, thisrare;   /*segctop*/
  censorrec *WITH;
  thisarray *WITH1;
  long FORLIM;
  thisarray *WITH2;

  initseg(LINK);
  WITH = censorstruct;
  WITH1 = (*LINK->p)->gen;
  FORLIM = fgeno;
  for (first = 0; first < FORLIM; first++) {
    if (WITH1->genarray[first] != 0.0) {
      segval = 0.0;
      LINK->fstart = probstart[first];
      LINK->fend = probend[first];
      LINK->fseg = segstart[first];
      thisrare = rare[first];
      second = 1;
      WITH2 = (*LINK->q)->gen;
      do {
	LINK->sstart = probstart[second - 1];
	LINK->send = probend[second - 1];
	LINK->sseg = segstart[second - 1];
	sstop = second + LINK->send - LINK->sstart + 1;
	if (thisc < maxcensor) {
	  thisc++;
	  thiscensor = WITH->censor[thisc - minint];
	} else
	  thiscensor = (WITH2->genarray[second - 1] == 0.0 ||
			thisrare && rare[second - 1]);
	if (thiscensor)
	  second = sstop;
	else {
	  if (mutsys != 0)
	    segval += WITH2->genarray[second - 1] * msegfast(LINK);
	  else
	    segval += WITH2->genarray[second - 1] * segfast(LINK);
	  second = sstop;
	}
      } while (second <= fgeno);
      WITH1->genarray[first] *= segval * segscale;
    }
  }
  exitseg(LINK);
  if (approximate)
    getapprox(LINK);
}  /* segctop */

/*segctop*/

Local Void segtop(LINK)
struct LOC_seg *LINK;
{
  double segval, val;
  long first, second, sstop;
  boolean thatrare, thisrare;   /*segtop*/
  censorrec *WITH;
  thisarray *WITH1;
  long FORLIM;
  thisarray *WITH2;

  initseg(LINK);
  WITH = censorstruct;
  WITH1 = (*LINK->p)->gen;
  FORLIM = fgeno;
  for (first = 0; first < FORLIM; first++) {
    if (WITH1->genarray[first] != 0.0) {
      thisrare = rare[first];
      segval = 0.0;
      LINK->fstart = probstart[first];
      LINK->fend = probend[first];
      LINK->fseg = segstart[first];
      second = 1;
      WITH2 = (*LINK->q)->gen;
      do {
	LINK->sstart = probstart[second - 1];
	LINK->send = probend[second - 1];
	LINK->sseg = segstart[second - 1];
	sstop = second + LINK->send - LINK->sstart + 1;
	if (thisc <= maxcensor)
	  thisc++;
	thatrare = (thisrare && rare[second - 1]);
	if (WITH2->genarray[second - 1] == 0.0 || thatrare) {
	  second = sstop;
	  if (thisc <= maxcensor)
	    WITH->censor[thisc - minint] = true;
	} else {
	  if (mutsys != 0)
	    val = msegfast(LINK);
	  else
	    val = segfast(LINK);
	  if (thisc <= maxcensor)
	    WITH->censor[thisc - minint] = (val == 0.0);
	  segval += WITH2->genarray[second - 1] * val;
	  second = sstop;
	}
      } while (second <= fgeno);
      WITH1->genarray[first] *= segval * segscale;
    }
  }
  cleanup(LINK->q);
  exitseg(LINK);
  if (approximate)
    getapprox(LINK);
}  /* segtop */

/*segtop*/

Local Void segcapprox(LINK)
struct LOC_seg *LINK;
{
  double segval;
  long first, second, sstop;
  boolean thiscensor;   /*segcapprox*/
  censorrec *WITH;
  approxrec *WITH1;
  thisarray *WITH2;
  long FORLIM;
  thisarray *WITH3;

  initseg(LINK);
  WITH = censorstruct;
  WITH1 = approxstruct;
  WITH2 = (*LINK->p)->gen;
  FORLIM = fgeno;
  for (first = 0; first < FORLIM; first++) {
    if (WITH2->genarray[first] != 0.0) {
      segval = 0.0;
      LINK->fstart = probstart[first];
      LINK->fend = probend[first];
      LINK->fseg = segstart[first];
      second = 1;
      WITH3 = (*LINK->q)->gen;
      do {
	LINK->sstart = probstart[second - 1];
	LINK->send = probend[second - 1];
	LINK->sseg = segstart[second - 1];
	sstop = second + LINK->send - LINK->sstart + 1;
	if (thisc < maxcensor) {
	  thisc++;
	  thiscensor = WITH->censor[thisc - minint];
	} else
	  thiscensor = (WITH3->genarray[second - 1] == 0.0);
	if (thiscensor || !WITH1->approxarray[LINK->LINK->thisped - 1][first])
	  second = sstop;
	else {
	  if (mutsys != 0)
	    segval += WITH3->genarray[second - 1] * msegfast(LINK);
	  else
	    segval += WITH3->genarray[second - 1] * segfast(LINK);
	  second = sstop;
	}
      } while (second <= fgeno);
      WITH2->genarray[first] *= segval * segscale;
    }
  }
  cleanup(LINK->q);
  exitseg(LINK);
}  /* segcapprox */

/*segcapprox*/

Local Void segup(LINK)
struct LOC_seg *LINK;
{
  double segval;
  long first, second;   /*segup*/
  censorrec *WITH;
  thisarray *WITH1;
  long FORLIM;
  thisarray *WITH2;
  long FORLIM1;

  initseg(LINK);
  WITH = censorstruct;
  WITH1 = (*LINK->p)->gen;
  FORLIM = fgeno;
  for (first = 0; first < FORLIM; first++) {
    if (WITH1->genarray[first] != 0.0) {
      segval = 0.0;
      LINK->fstart = probstart[first];
      LINK->fend = probend[first];
      LINK->fseg = segstart[first];
      WITH2 = (*LINK->q)->gen;
      FORLIM1 = fgeno;
      for (second = 0; second < FORLIM1; second++) {
	if (WITH2->genarray[second] != 0.0) {
	  LINK->sstart = probstart[second];
	  LINK->send = probend[second];
	  LINK->sseg = segstart[second];
	  if (mutsys != 0)
	    segval += WITH2->genarray[second] * msegfun(LINK);
	  else
	    segval += WITH2->genarray[second] * segfun(LINK);
	}
      }
      WITH1->genarray[first] *= segval * segscale;
    }
  }
  cleanup(LINK->q);
  exitseg(LINK);
}  /* segup */

/*segup*/

Local Void segdown(LINK)
struct LOC_seg *LINK;
{
  short here;
  genotype gene;
  double val, temp1, temp2;
  long f1, f2, s1, s2, i, j, first, second;
  short FORLIM;
  gennurec *WITH;
  censorrec *WITH1;
  thisarray *WITH2;
  long FORLIM1;
  thisarray *WITH3;
  long FORLIM2;
  thetavalues *WITH4;
  long FORLIM3;
  thetavalues *WITH5;
  long FORLIM4;

  initseg(LINK);
  FORLIM = fgeno;
  for (here = 0; here < FORLIM; here++)
    gene[here] = 0.0;
  WITH = gennustruct;
  WITH1 = censorstruct;
  WITH2 = (*LINK->p)->gen;
  FORLIM1 = mgeno;
  for (first = 0; first < FORLIM1; first++) {
    if (WITH2->genarray[first] != 0.0) {
      LINK->fstart = probstart[first];
      LINK->fend = probend[first];
      LINK->fseg = segstart[first];
      WITH3 = (*LINK->q)->gen;
      FORLIM2 = fgeno;
      for (second = 0; second < FORLIM2; second++) {
	if (WITH3->genarray[second] != 0.0) {
	  LINK->sstart = probstart[second];
	  LINK->send = probend[second];
	  LINK->sseg = segstart[second];
	  val = segscale * WITH3->genarray[second] *
		(*LINK->p)->gen->genarray[first];
	  if (nchild != 0)
	    val = segfun(LINK) * val;
	  if (val != 0.0) {
	    LINK->firstseg = LINK->fseg;
	    WITH4 = maletheta;
	    FORLIM3 = LINK->fend;
	    for (i = LINK->fstart - 1; i < FORLIM3; i++) {
	      temp1 = WITH4->segprob[i];
	      LINK->secondseg = LINK->sseg;
	      WITH5 = femaletheta;
	      FORLIM4 = LINK->send;
	      for (j = LINK->sstart - 1; j < FORLIM4; j++) {
		temp2 = WITH5->segprob[j] * temp1 * val;
		f1 = seghap1[LINK->firstseg - 1];
		f2 = seghap2[LINK->firstseg - 1];
		s1 = seghap1[LINK->secondseg - 1];
		s2 = seghap2[LINK->secondseg - 1];
		here = WITH->genenumber[s1 - 1][f1 - 1];
		gene[here - 1] += temp2;
		here = WITH->genenumber[s1 - 1][f2 - 1];
		gene[here - 1] += temp2;
		here = WITH->genenumber[s2 - 1][f1 - 1];
		gene[here - 1] += temp2;
		here = WITH->genenumber[s2 - 1][f2 - 1];
		gene[here - 1] += temp2;
		LINK->secondseg++;
	      }
	      LINK->firstseg++;
	    }
	  }  /*second*/
	}
      }
    }  /*first*/
  }
  memcpy((*LINK->p)->gen->genarray, (*LINK->r)->gen->genarray,
	 sizeof(genotype));
  memcpy((*LINK->r)->gen->genarray, gene, sizeof(genotype));
  WITH2 = (*LINK->r)->gen;
  FORLIM1 = fgeno;
  for (first = 0; first < FORLIM1; first++)
    WITH2->genarray[first] *= (*LINK->p)->gen->genarray[first];
  cleanup(LINK->p);
  cleanup(LINK->q);
  exitseg(LINK);
}  /* segdown */

/*segdown*/

Local Void msegsexdown(LINK)
struct LOC_seg *LINK;
{
  short here;
  genotype gene;
  double val, temp2;
  long ms2, ms1, mf, j, first, second;   /*msegsexdown*/
  short FORLIM;
  gennurec *WITH;
  censorrec *WITH1;
  thisarray *WITH2;
  long FORLIM1;
  thisarray *WITH3;
  long FORLIM2;
  thetavalues *WITH4;
  long FORLIM3;

  initseg(LINK);
  FORLIM = fgeno;
  for (here = 0; here < FORLIM; here++)
    gene[here] = 0.0;
  WITH = gennustruct;
  WITH1 = censorstruct;
  WITH2 = (*LINK->p)->gen;
  FORLIM1 = LINK->nfirst;
  for (first = 0; first < FORLIM1; first++) {
    if (WITH2->genarray[first] != 0.0) {
      LINK->fseg = first + 1;
      second = 1;
      WITH3 = (*LINK->q)->gen;
      FORLIM2 = LINK->nsecond;
      for (second = 0; second < FORLIM2; second++) {
	if (WITH3->genarray[second] != 0.0) {
	  val = WITH3->genarray[second] * (*LINK->p)->gen->genarray[first] *
		segscale;
	  LINK->sstart = probstart[second];
	  LINK->send = probend[second];
	  LINK->sseg = segstart[second];
	  if (nchild != 0)
	    val *= msegsex(LINK);
	  if (val != 0.0) {
	    mf = muthap[seghap1[LINK->fseg - 1] - 1];
	    LINK->secondseg = LINK->sseg;
	    WITH4 = femaletheta;
	    FORLIM3 = LINK->send;
	    for (j = LINK->sstart - 1; j < FORLIM3; j++) {
	      ms1 = muthap[seghap1[LINK->secondseg - 1] - 1];
	      ms2 = muthap[seghap2[LINK->secondseg - 1] - 1];
	      temp2 = WITH4->segprob[j];
	      if ((*LINK->r)->male) {
		here = seghap1[LINK->secondseg - 1];
		gene[here - 1] += (1 - LINK->ps) * temp2 * val;
		here = seghap2[LINK->secondseg - 1];
		gene[here - 1] += (1 - LINK->ps) * temp2 * val;
		here = ms1;
		gene[here - 1] += LINK->ps * temp2 * val;
		here = ms2;
		gene[here - 1] += LINK->ps * temp2 * val;
	      } else {
		here = WITH->genenumber[seghap1[LINK->secondseg - 1] - 1]
		  [first];
		gene[here - 1] += (1 - LINK->pf) * (1 - LINK->ps) * temp2 * val;
		here = WITH->genenumber[seghap2[LINK->secondseg - 1] - 1]
		  [first];
		gene[here - 1] += (1 - LINK->pf) * (1 - LINK->ps) * temp2 * val;
		here = WITH->genenumber[seghap1[LINK->secondseg - 1] - 1]
		  [mf - 1];
		gene[here - 1] += LINK->pf * (1 - LINK->ps) * temp2 * val;
		here = WITH->genenumber[seghap2[LINK->secondseg - 1] - 1]
		  [mf - 1];
		gene[here - 1] += LINK->pf * (1 - LINK->ps) * temp2 * val;
		here = WITH->genenumber[ms1 - 1][first];
		gene[here - 1] += (1 - LINK->pf) * LINK->ps * temp2 * val;
		here = WITH->genenumber[ms2 - 1][first];
		gene[here - 1] += (1 - LINK->pf) * LINK->ps * temp2 * val;
		here = WITH->genenumber[ms1 - 1][mf - 1];
		gene[here - 1] += LINK->pf * LINK->ps * temp2 * val;
		here = WITH->genenumber[ms2 - 1][mf - 1];
		gene[here - 1] += LINK->pf * LINK->ps * temp2 * val;
	      }
	      LINK->secondseg++;
	    }
	  }  /*second*/
	}
      }
    }  /*first*/
  }
  memcpy((*LINK->p)->gen->genarray, (*LINK->r)->gen->genarray,
	 sizeof(genotype));
  memcpy((*LINK->r)->gen->genarray, gene, sizeof(genotype));
  WITH2 = (*LINK->r)->gen;
  FORLIM1 = fgeno;
  for (first = 0; first < FORLIM1; first++)
    WITH2->genarray[first] *= (*LINK->p)->gen->genarray[first];
  cleanup(LINK->p);
  cleanup(LINK->q);
  exitseg(LINK);
}  /* msegsexdown */

/*msegsexdown*/

Local Void msegdown(LINK)
struct LOC_seg *LINK;
{
  short here;
  genotype gene;
  double val, temp, temp1, temp2;
  long i, j, first, second, f1, f2, s1, s2, mf1, mf2, ms1, ms2;   /*msegdown*/
  short FORLIM;
  gennurec *WITH;
  censorrec *WITH1;
  thisarray *WITH2;
  long FORLIM1;
  thisarray *WITH3;
  long FORLIM2;
  thetavalues *WITH4;
  long FORLIM3;
  thetavalues *WITH5;
  long FORLIM4;

  initseg(LINK);
  FORLIM = fgeno;
  for (here = 0; here < FORLIM; here++)
    gene[here] = 0.0;
  WITH = gennustruct;
  WITH1 = censorstruct;
  WITH2 = (*LINK->p)->gen;
  FORLIM1 = fgeno;
  for (first = 0; first < FORLIM1; first++) {
    if (WITH2->genarray[first] != 0.0) {
      LINK->fstart = probstart[first];
      LINK->fend = probend[first];
      LINK->fseg = segstart[first];
      WITH3 = (*LINK->q)->gen;
      FORLIM2 = fgeno;
      for (second = 0; second < FORLIM2; second++) {
	if (WITH3->genarray[second] != 0.0) {
	  LINK->sstart = probstart[second];
	  LINK->send = probend[second];
	  LINK->sseg = segstart[second];
	  val = WITH3->genarray[second] * (*LINK->p)->gen->genarray[first] *
		segscale;
	  if (nchild != 0)
	    val = msegfun(LINK) * val;
	  if (val != 0.0) {
	    LINK->firstseg = LINK->fseg;
	    WITH4 = maletheta;
	    FORLIM3 = LINK->fend;
	    for (i = LINK->fstart - 1; i < FORLIM3; i++) {
	      temp1 = WITH4->segprob[i];
	      LINK->secondseg = LINK->sseg;
	      WITH5 = femaletheta;
	      FORLIM4 = LINK->send;
	      for (j = LINK->sstart - 1; j < FORLIM4; j++) {
		temp2 = WITH5->segprob[j];
		f1 = seghap1[LINK->firstseg - 1];
		f2 = seghap2[LINK->firstseg - 1];
		s1 = seghap1[LINK->secondseg - 1];
		s2 = seghap2[LINK->secondseg - 1];
		temp = (1 - LINK->pf) * (1 - LINK->ps) * temp1 * temp2 * val;
		here = WITH->genenumber[s1 - 1][f1 - 1];
		gene[here - 1] += temp;
		here = WITH->genenumber[s1 - 1][f2 - 1];
		gene[here - 1] += temp;
		here = WITH->genenumber[s2 - 1][f1 - 1];
		gene[here - 1] += temp;
		here = WITH->genenumber[s2 - 1][f2 - 1];
		gene[here - 1] += temp;
		ms1 = muthap[s1 - 1];
		ms2 = muthap[s2 - 1];
		temp = (1 - LINK->pf) * LINK->ps * temp1 * temp2 * val;
		here = WITH->genenumber[ms1 - 1][f1 - 1];
		gene[here - 1] += temp;
		here = WITH->genenumber[ms1 - 1][f2 - 1];
		gene[here - 1] += temp;
		here = WITH->genenumber[ms2 - 1][f1 - 1];
		gene[here - 1] += temp;
		here = WITH->genenumber[ms2 - 1][f2 - 1];
		gene[here - 1] += temp;
		mf1 = muthap[f1 - 1];
		mf2 = muthap[f2 - 1];
		temp = LINK->pf * (1 - LINK->ps) * temp1 * temp2 * val;
		here = WITH->genenumber[mf1 - 1][s1 - 1];
		gene[here - 1] += temp;
		here = WITH->genenumber[mf1 - 1][s2 - 1];
		gene[here - 1] += temp;
		here = WITH->genenumber[mf2 - 1][s1 - 1];
		gene[here - 1] += temp;
		here = WITH->genenumber[mf2 - 1][s2 - 1];
		gene[here - 1] += temp;
		temp = LINK->pf * LINK->ps * temp1 * temp2 * val;
		here = WITH->genenumber[mf1 - 1][ms1 - 1];
		gene[here - 1] += temp;
		here = WITH->genenumber[mf1 - 1][ms2 - 1];
		gene[here - 1] += temp;
		here = WITH->genenumber[mf2 - 1][ms1 - 1];
		gene[here - 1] += temp;
		here = WITH->genenumber[mf2 - 1][ms2 - 1];
		gene[here - 1] += temp;
		LINK->secondseg++;
	      }
	      LINK->firstseg++;
	    }
	  }  /*second*/
	}
      }
    }  /*first*/
  }
  memcpy((*LINK->p)->gen->genarray, (*LINK->r)->gen->genarray,
	 sizeof(genotype));
  memcpy((*LINK->r)->gen->genarray, gene, sizeof(genotype));
  WITH2 = (*LINK->r)->gen;
  FORLIM1 = fgeno;
  for (first = 0; first < FORLIM1; first++)
    WITH2->genarray[first] *= (*LINK->p)->gen->genarray[first];
  cleanup(LINK->p);
  cleanup(LINK->q);
  exitseg(LINK);
}  /* msegdown */

/*prob*/

Local Void seg(p_, q_, r_, peel, LINK)
thisperson **p_, **q_, **r_;
direction peel;
struct LOC_likelihood *LINK;
{
  struct LOC_seg V;
  boolean phaseunkn;

  /*msegdown*/

  V.LINK = LINK;
  V.p = p_;
  V.q = q_;
  V.r = r_;
  phaseunkn = ((*V.q)->pa == NULL && (*V.q)->firstpass &&
	       (*V.q)->inloop == 0 && !disequi);
  if ((*V.p)->male) {
    V.father = *V.p;
    V.mother = *V.q;
  } else {
    V.father = *V.q;
    V.mother = *V.p;
  }
  if (peel == peelup) {
    if (approximate) {
      if (firstapprox && firsttime) {
	switch (thispath) {

	case auto_:
	  segtop(&V);
	  break;

	case mauto:
	  segtop(&V);
	  break;

	case sex:
	  segsextop(&V);
	  break;

	case msex:
	  segsextop(&V);
	  break;
	}
      } else if (firstapprox) {
	switch (thispath) {

	case auto_:
	  segctop(&V);
	  break;

	case mauto:
	  segctop(&V);
	  break;

	case sex:
	  segsexctop(&V);
	  break;

	case msex:
	  segsexctop(&V);
	  break;
	}
      } else {
	switch (thispath) {

	case auto_:
	  segcapprox(&V);
	  break;

	case mauto:
	  segcapprox(&V);
	  break;

	case sex:
	  segsexctop(&V);
	  break;

	case msex:
	  segsexctop(&V);
	  break;
	}
      }
    } else if (phaseunkn) {
      if (firsttime) {
	switch (thispath) {

	case auto_:
	  segtop(&V);
	  break;

	case mauto:
	  segtop(&V);
	  break;

	case sex:
	  segsextop(&V);
	  break;

	case msex:
	  segsextop(&V);
	  break;
	}
      } else {  /*not firsttime*/
	switch (thispath) {

	case auto_:
	  segctop(&V);
	  break;

	case mauto:
	  segctop(&V);
	  break;

	case sex:
	  segsexctop(&V);
	  break;

	case msex:
	  segsexctop(&V);
	  break;
	}
      }
    } else {
      switch (thispath) {

      case auto_:
	segup(&V);
	break;

      case mauto:
	segup(&V);
	break;

      case sex:
	segsexup(&V);
	break;

      case msex:
	segsexup(&V);
	break;
      }
    }
  } else {  /*not peelup*/
    switch (thispath) {

    case auto_:
      segdown(&V);
      break;

    case mauto:
      msegdown(&V);
      break;

    case sex:
      segsexdown(&V);
      break;

    case msex:
      msegsexdown(&V);
      break;
    }
  }
  (*V.q)->firstpass = false;
  (*V.p)->firstpass = false;

  /*first approximate not first time*/
  /*approximate*/
  /*do not approximate*/
  /*phaseinfo*/
}  /* seg */

Local Void collapseup(p, LINK)
thisperson *p;
struct LOC_likelihood *LINK;
{
  thisperson *q, *child, *nextchild;
  boolean down;

  p->done = true;
  if (p->foff == NULL)
    return;
  down = false;
  child = p->foff;
  while (child != NULL) {
    down = false;
    if (p->male)
      q = child->ma;
    else
      q = child->pa;
    if (!q->done) {
      collapsedown(q, LINK);
      nextchild = child;
      while (nextchild != NULL) {
	if (nextchild->pa == q || nextchild->ma == q) {
	  if (!nextchild->up)
	    collapseup(nextchild, LINK);
	  else
	    down = true;
	}
	if (p->male)
	  nextchild = nextchild->nextpa;
	else
	  nextchild = nextchild->nextma;
      }
      if (q->multi)
	collapseup(q, LINK);
      if (!down)
	seg(&p, &q, &child, peelup, LINK);
      else
	collapsedown(p, LINK);
    }
    if (p->male)
      child = child->nextpa;
    else
      child = child->nextma;
  }
}  /* collapseup */

/*collapseup*/

Local Void collapsedown(p, LINK)
thisperson *p;
struct LOC_likelihood *LINK;
{
  if (p->pa == NULL)
    return;
  p->up = true;
  collapseup(p->pa, LINK);
  seg(&p->pa, &p->ma, &p, peeldown, LINK);
}  /* collapsedown */

/*collapsedown*/

Local Void riskcumul(LINK)
struct LOC_likelihood *LINK;
{
  long i;
  thisarray *WITH;
  long FORLIM;

  WITH = LINK->proband->gen;
  if (sexlink && LINK->proband->male) {
    FORLIM = mgeno;
    for (i = 0; i < FORLIM; i++) {
      if (riskmale[i])
	LINK->hetero += WITH->genarray[i];
    }
    return;
  }
  FORLIM = fgeno;
  for (i = 0; i < FORLIM; i++) {
    if (risk2[i])
      LINK->homo += WITH->genarray[i];
    else if (risk1[i])
      LINK->hetero += WITH->genarray[i];
  }
}  /* riskcumul */

/*riskcumul*/

Local Void riskcalc(LINK)
struct LOC_likelihood *LINK;
{
  double normal;

  LINK->homo /= like;
  LINK->hetero /= like;
  normal = 1 - LINK->homo - LINK->hetero;
  if (!print)
    return;
  fprintf(outfile, "RISK FOR PERSON %6ld IN PEDIGREE %7ld\n",
	  LINK->proband->id, LINK->proband->ped);
  if (!LINK->proband->male || !sexlink)
    fprintf(outfile, "HOMOZYGOTE CARRIER   : %8.5f\n", LINK->homo);
  if (!LINK->proband->male || !sexlink)
    fprintf(outfile, "HETEROZYGOTE CARRIER : %8.5f\n", LINK->hetero);
  else
    fprintf(outfile, "MALE CARRIER         : %8.5f\n", LINK->hetero);
  fprintf(outfile, "NORMAL               : %8.5f\n", normal);
  printf("RISK FOR PERSON %6ld IN PEDIGREE %7ld\n",
	 LINK->proband->id, LINK->proband->ped);
  if (!LINK->proband->male || !sexlink)
    printf("HOMOZYGOTE CARRIER   : %8.5f\n", LINK->homo);
  if (!LINK->proband->male || !sexlink)
    printf("HETEROZYGOTE CARRIER : %8.5f\n", LINK->hetero);
  else
    printf("MALE CARRIER         : %8.5f\n", LINK->hetero);
  printf("NORMAL               : %8.5f\n", normal);
}  /* riskcalc */


/*seg*/

Static Void likelihood(thisped_, proband_)
long thisped_;
thisperson *proband_;
{
  struct LOC_likelihood V;
  long loopmax[maxloop];
  double tmplike;
  long i, j;
  boolean gocalc, alldone;
  thisperson *WITH;
  thisarray *WITH1;
  long FORLIM1;

  /*riskcalc*/

  V.thisped = thisped_;
  V.proband = proband_;
  if (!informative[V.thisped - 1]) {
    like = 0.0;
    return;
  }
  V.homo = 0.0;
  V.hetero = 0.0;
  tmplike = 0.0;
  alldone = false;
  V.nuscale = 0;
  for (i = 0; i < maxloop; i++) {
    V.loopgen[i] = 1;
    loopmax[i] = 1;
    V.holdpoint[i] = NULL;
    if (looppers[V.thisped - 1][i][0] != NULL) {
      WITH = looppers[V.thisped - 1][i][0];
      WITH->gen = (thisarray *)Malloc(sizeof(thisarray));
      WITH1 = WITH->gen;
      FORLIM1 = fgeno;
      /*          gen:=genpoint(NewPtr(SizeOf(thisarray)));*/
      for (j = 0; j < FORLIM1; j++)
	WITH1->genarray[j] = 0.0;
      getvect(looppers[V.thisped - 1][i][0]);
      if (looppers[V.thisped - 1][i][0]->pa == NULL)
	V.nuscale++;
      V.holdpoint[i] = WITH->gen;
      WITH->gen = NULL;
      if (WITH->male)
	loopmax[i] = mgeno;
      else
	loopmax[i] = fgeno;
    }
  }
  V.loopgen[0] = 0;
  do {
    i = 1;
    do {
      V.loopgen[i - 1]++;
      if (V.loopgen[i - 1] > loopmax[i - 1])
	V.loopgen[i - 1] = 1;
      else
	i = maxloop;
      i++;
    } while (i <= maxloop);
    gocalc = true;
    for (i = 0; i < maxloop; i++) {
      /*ML change*/
      if (V.holdpoint[i] != NULL) {
	if (V.holdpoint[i]->genarray[V.loopgen[i] - 1] == 0.0)
	  gocalc = false;
      }
    }
    if (gocalc) {
      FORLIM1 = totperson;
      for (i = 1; i <= FORLIM1; i++) {
	WITH = person[i];
	WITH->gen = NULL;
	WITH->done = false;
	WITH->up = false;
      }
      collapseup(V.proband, &V);
      collapsedown(V.proband, &V);
      like = 0.0;
      WITH1 = V.proband->gen;
      if (V.proband->male) {
	FORLIM1 = mgeno;
	for (i = 0; i < FORLIM1; i++)
	  like += WITH1->genarray[i];
      } else {
	FORLIM1 = fgeno;
	for (i = 0; i < FORLIM1; i++)
	  like += WITH1->genarray[i];
      }
      tmplike += like;
      if (risk && like != 0.0)
	riskcumul(&V);
      FORLIM1 = totperson;
      for (i = 1; i <= FORLIM1; i++) {
	if (person[i]->gen != NULL)
	  cleanup(&person[i]);
      }
    }
    alldone = true;
    for (i = 0; i < maxloop; i++)
      alldone = (alldone && V.loopgen[i] == loopmax[i]);
  } while (!alldone);
  like = tmplike;
  if (risk && like != 0.0)
    riskcalc(&V);
  if (like == 0.0)
    like = zerolike;
  else
    like = log(like) - V.nuscale * log(segscale);
  for (i = 0; i < maxloop; i++) {
    if (V.holdpoint[i] != NULL) {
      Free(V.holdpoint[i]);
      /*        DisposPtr(Ptr(holdpoint[i]));*/
      V.holdpoint[i] = NULL;
    }
  }
}


/* likelihood */
/*likelihood*/


Static Void checkzero()
{
  long i;
  thetavalues *WITH;
  long FORLIM;

  if (!firsttime) {
    WITH = maletheta;
    FORLIM = nlocus - 2;
    for (i = 0; i <= FORLIM; i++) {
      if (WITH->theta[i] != 0.0 && zeromale[i])
	firsttime = true;
    }
    WITH = femaletheta;
    FORLIM = nlocus - 2;
    for (i = 0; i <= FORLIM; i++) {
      if (WITH->theta[i] != 0.0 && zerofemale[i])
	firsttime = true;
    }
  }
  if (maletheta->theta[which - 1] == 0.0)
    firsttime = true;
  if (femaletheta->theta[which - 1] == 0.0)
    firsttime = true;
  if (!firsttime)
    return;
  WITH = maletheta;
  FORLIM = nlocus - 2;
  for (i = 0; i <= FORLIM; i++)
    zeromale[i] = (WITH->theta[i] == 0.0);
  WITH = femaletheta;
  FORLIM = nlocus - 2;
  for (i = 0; i <= FORLIM; i++)
    zerofemale[i] = (WITH->theta[i] == 0.0);
}  /* checkzero */


Static Void iterpeds()
{
  long i, locip;
  long thisped;
  double rep;
  long FORLIM;
  char FORLIM1;
  double TEMP;

  tlike = 0.0;
  alike = 0.0;
  FORLIM = totperson;
  for (i = 1; i <= FORLIM; i++)
    person[i]->done = false;
  FORLIM = totperson;
  for (i = 1; i <= FORLIM; i++)
    person[i]->firstpass = true;
  thisc = minint;
  recombination();
  checkzero();
  if (print) {
    for (i = 1; i <= 48; i++)
      putc('-', outfile);
    putc('\n', outfile);
    for (i = 1; i <= 48; i++)
      putc('-', outfile);
    putc('\n', outfile);
    if (sexdif)
      fprintf(outfile, "MALE THETAS   ");
    else
      fprintf(outfile, "THETAS ");
    FORLIM = nlocus - 2;
    for (i = 0; i <= FORLIM; i++)
      fprintf(outfile, "%6.3f", maletheta->theta[i]);
    if (interfer)
      fprintf(outfile, "%6.3f", maletheta->theta[nlocus - 1]);
    putc('\n', outfile);
    if (sexdif) {
      fprintf(outfile, "FEMALE THETAS ");
      FORLIM = nlocus - 2;
      for (i = 0; i <= FORLIM; i++)
	fprintf(outfile, "%6.3f", femaletheta->theta[i]);
      if (interfer)
	fprintf(outfile, "%6.3f", femaletheta->theta[nlocus - 1]);
      putc('\n', outfile);
    }
    for (i = 1; i <= 48; i++)
      putc('-', outfile);
    fprintf(outfile, "\nPEDIGREE  |  LN LIKE  | LOG 10 LIKE|");
    if (nlocus == 2)
      fprintf(outfile, " LOD SCORE\n");
    else
      fprintf(outfile, " MULTIPOINT LOD SCORE\n");
    for (i = 1; i <= 48; i++)
      putc('-', outfile);
    putc('\n', outfile);
    for (i = 1; i <= 48; i++)
      putchar('-');
    putchar('\n');
    for (i = 1; i <= 48; i++)
      putchar('-');
    putchar('\n');
    if (sexdif)
      printf("MALE THETAS   ");
    else
      printf("THETAS ");
    FORLIM = nlocus - 2;
    for (i = 0; i <= FORLIM; i++)
      printf("%6.3f", maletheta->theta[i]);
    if (interfer)
      printf("%6.3f", maletheta->theta[nlocus - 1]);
    putchar('\n');
    if (sexdif) {
      printf("FEMALE THETAS ");
      FORLIM = nlocus - 2;
      for (i = 0; i <= FORLIM; i++)
	printf("%6.3f", femaletheta->theta[i]);
      if (interfer)
	printf("%6.3f", femaletheta->theta[nlocus - 1]);
      putchar('\n');
    }
    for (i = 1; i <= 48; i++)
      putchar('-');
    printf("\nPEDIGREE  |  LN LIKE  | LOG 10 LIKE|");
    if (nlocus == 2)
      printf(" LOD SCORE\n");
    else
      printf(" MULTIPOINT LOD SCORE\n");
    for (i = 1; i <= 48; i++)
      putchar('-');
    putchar('\n');
  }  /* print out loop */
  FORLIM1 = nuped;
  for (thisped = 0; thisped < FORLIM1; thisped++) {
    likelihood((long)(thisped + 1), proband[thisped]);
    if (!risk) {
      if (maletheta->theta[which - 1] == 0.5)
	stand[thisped] = like / log10_;
    }
    if (print)
      fprintf(outfile, "%9ld %12.6f ", proband[thisped]->ped, like);
    if (print)
      printf("%9ld %12.6f ", proband[thisped]->ped, like);
    alike += like;
    like /= log10_;
    if (byfamily) {
      if (!risk) {
	if (nlocus == 2)   /*Use lod scores*/
	  lods[thisped] = like - stand[thisped];
	      /*Use multipoint lod scores*/
	else
	  lods[thisped] = like - stand[thisped];
      }
      /*  lods[thisped]:=-2.0*(stand[thisped]-like)*log10;  location scores*/
      if (print)
	fprintf(outfile, "%12.6f", like);
      if (print)
	fprintf(outfile, "%12.6f \n", lods[thisped]);
    }
    if (print)
      printf("%12.6f", like);
    if (print)
      printf("%12.6f \n", lods[thisped]);
    tlike += like;
  }
  for (i = 1; i <= 48; i++) {
    if (print)
      putc('-', outfile);
  }
  if (print)
    putc('\n', outfile);
  /*nrep is the current replication being analyzed; npt indicates
  the current theta being analyzed; opeds is the number of
  original pedigrees.*/
  /*Write out male theta(s) to lod score file; Don't write if
  only two systems and theta = 0.50 */
  if (lodprint) {
    if (nlocus != 2 || maletheta->theta[0] != 0.50) {
      FORLIM = nlocus - 2;
      for (i = 0; i <= FORLIM; i++)
	fprintf(lodfile, "%6.3f", maletheta->theta[i]);
    }
  }
  FORLIM = opeds;
  for (locip = 0; locip < FORLIM; locip++) {
    /*Write out lodscores by pedigree on the same line as the theta(s) */
    if (lodprint) {
      if (nlocus != 2 || maletheta->theta[0] != 0.50)
	fprintf(lodfile, "%12.6f", lods[locip]);
    }
    /* Store the lod scores for the point npt and the pedigree locip
    in the array zlod */
    zlod[locip][npt - 1] = lods[locip];
    if (nrep == 1) {
      aver[locip][npt - 1] = lods[locip];
      variance[locip][npt - 1] = 0.0;
      small[locip][npt - 1] = lods[locip];
      big[locip][npt - 1] = lods[locip];
    } else {
      rep = nrep;
      /* rep := double(nrep);*/
      aver[locip]
	[npt - 1] = (rep - 1.0) * aver[locip][npt - 1] / rep + lods[locip] / rep;
      TEMP = lods[locip] - aver[locip][npt - 1];
      variance[locip]
	[npt - 1] = (rep - 1.0) * variance[locip][npt - 1] / rep +
		    TEMP * TEMP / (rep - 1.0);
      small[locip][npt - 1] = min(small[locip][npt - 1], lods[locip]);
      big[locip][npt - 1] = max(big[locip][npt - 1], lods[locip]);
    }
  }  /*loop on locip*/
  /*Put a line feed into the lodfile */
  if (lodprint) {
    if (nlocus != 2 || maletheta->theta[0] != 0.50)
      putc('\n', lodfile);
  }
  FORLIM = opeds;
  /* Put the statistics for each study into the opeds+1 position of the arrays*/
  for (locip = 1; locip < FORLIM; locip++)
    lods[0] += lods[locip];
  locip = opeds + 1;
  zlod[locip - 1][npt - 1] = lods[0];
  if (nrep == 1) {
    aver[locip - 1][npt - 1] = lods[0];
    variance[locip - 1][npt - 1] = 0.0;
    small[locip - 1][npt - 1] = lods[0];
    big[locip - 1][npt - 1] = lods[0];
  } else {
    rep = nrep;
    /*rep := double(nrep);*/
    aver[locip - 1]
      [npt - 1] = (rep - 1.0) * aver[locip - 1][npt - 1] / rep + lods[0] / rep;
    TEMP = lods[0] - aver[locip - 1][npt - 1];
    variance[locip - 1][npt - 1] = (rep - 1.0) * variance[locip - 1]
				   [npt - 1] / rep + TEMP * TEMP / (rep - 1.0);
    small[locip - 1][npt - 1] = min(small[locip - 1][npt - 1], lods[0]);
    big[locip - 1][npt - 1] = max(big[locip - 1][npt - 1], lods[0]);
  }

  for (i = 1; i <= 48; i++) {
    if (print)
      putc('-', outfile);
  }
  if (print)
    putc('\n', outfile);
  if (print)
    fprintf(outfile, "TOTALS    %12.6f %12.6f\n", alike, tlike);
  for (i = 1; i <= 48; i++) {
    if (print)
      putchar('-');
  }
  if (print)
    putchar('\n');
  if (print)
    printf("TOTALS    %12.6f %12.6f\n", alike, tlike);
  alike = -2 * alike;
  if (print)
    fprintf(outfile, "-2 LN(LIKE) = % .5E", alike);
  if (print)
    printf("-2 LN(LIKE) = % .5E", alike);
  if (!risk) {
    if (nlocus == 2) {
      if (maletheta->theta[which - 1] == 0.5)
	scorevalue = tlike;
      alike = tlike - scorevalue;
      if (print)
	fprintf(outfile, " LOD SCORE = %12.6f", alike);
      if (print)
	printf(" LOD SCORE = %12.6f", alike);
    } else {
      if (maletheta->theta[which - 1] == 0.5)
	scorevalue = alike;
      alike = scorevalue - alike;
      if (print)
	fprintf(outfile, " MULTIPOINT LOD SCORE = %12.6f", alike);
      if (print)
	printf(" MULTIPOINT LOD SCORE = %12.6f", alike);
    }
  }
  if (print)
    putc('\n', outfile);
  if (print)
    putchar('\n');
  if (print) {
    if (firsteff) {
      if (thisc < maxcensor)
	printf("Maxcensor can be reduced to %12ld\n", thisc);
      else if (thisc > maxcensor)
	printf("you may gain efficiency by increasing maxcensor\n");
    }
  }
  firsttime = false;
  firsteff = false;
}  /* iterpeds */


Static Void initialize()
{
  long i;

  printf("\nProgram MSIM version %s\n\n", version);
  printf("The program constants are set to the following maxima:\n");
  printf("%6ld loci in mapping problem (maxlocus)\n", (long)maxlocus);
  printf("%6ld alleles at a single locus (maxall)\n", (long)maxall);
  printf("%6ld recombination probabilities (maxneed)\n", (long)maxneed);
  printf("%6ld maximum of censoring array (maxcensor)\n", (long)maxcensor);
  printf("%6ld haplotypes = n1 x n2 x ... where ni = current # alleles at locus i\n",
	 (long)maxhap);
  printf("%6ld joint genotypes for a female\n", (long)maxfem);
  printf("%6ld joint genotypes for a male\n", (long)maxmal);
  printf("%6ld individuals in all pedigrees combined (maxind)\n",
	 (long)maxind);
  printf("%6ld pedigrees (maxped)\n", (long)maxped);
  printf("%6ld binary codes at a single locus (maxfact)\n", (long)maxfact);
  printf("%6ld quantitative factor(s) at a single locus\n", (long)maxtrait);
  printf("%6ld liability classes\n", (long)maxliab);
  printf("%6ld binary codes at a single locus\n", (long)maxfact);
  printf("%8.2f base scaling factor for likelihood (scale)\n", scale);
  printf("%8.2f scale multiplier for each locus (scalemult)\n", scalemult);
  printf("%8.5f frequency for elimination of heterozygotes (minfreq)\n",
	 minfreq);
  if (minfreq != 0.0) {
    printf("IMPORTANT : RECOMPILE THIS PROGRAM WITH MINFREQ=0.0\n");
    printf("FOR THE ANALYSIS OF RECESSIVE TRAITS\n");
  }
  putchar('\n');
  /* assign(speedfile,'speedfil.dat');
     assign(datafile,'datafile.dat');
     assign(outfile,'outfile.dat');
     assign(lodfile,'lodfile.dat');
     assign(ipedfile,'ipedfile.dat');
     assign(msim,'msim.dat');
     assign(simout,'simout.dat'); */
  /*SUN*/
  firsteff = true;
  for (i = 1; i <= 66; i++)
    putchar('*');
  printf("\n This is MSIM, which is a modified version of MLINK\n");
  printf(" designed to be used on simulated data produced by\n");
  printf(" the companion program SLINK.\n\n");
  printf(" The pedigree data are read from ipedfile.dat.\n");
  printf(" The statistical summary is written to msim.dat\n");
  for (i = 1; i <= 66; i++)
    putchar('*');
  putchar('\n');

  infile = false;   /*Indicates at the beginning of the pedigree*/
  /*Read the number of pedigrees per replicate from the most recent simout.dat*/
  printf("Reading SIMOUT.DAT\n");
  if (simout != NULL)
    simout = freopen("simout.dat", "r", simout);
  else
    simout = fopen("simout.dat", "r");
  if (simout == NULL)
    exit(FileNotFound);
  /*SUN*/
  for (i = 1; i <= 5; i++) {   /*The no. of pedigrees is on the sixth line*/
    fscanf(simout, "%*[^\n]");
    getc(simout);
  }
  do {
    chtemp = getc(simout);
    if (chtemp == '\n')
      chtemp = ' ';
  } while (chtemp != 's');   /*Read to the 's' in 'pedigrees' */
  fscanf(simout, "%ld%*[^\n]", &opeds);
  getc(simout);   /*Read the number of pedigrees */
  if (opeds > 1)
    printf(" There are %5ld pedigrees per replicate.\n", opeds);
  else
    printf(" There is %5ld pedigree per replicate.\n", opeds);
  if (opeds <= maxped)
    return;
  printf("ERROR: maxped is too small\n");
  fprintf(outfile, "ERROR: maxped is too small\n");
  printf("Press return to halt...\n");
  scanf("%*[^\n]");
  getchar();
  exit(0);
}


/* initialize */
/*initialize*/


Static Void quadout(fout)
FILE **fout;
{
  long i, ip, FORLIM, FORLIM1;

  fprintf(*fout,
	  " Average Maximum Lod Scores based on quadratic interpolation \n");
  for (i = 1; i <= 66; i++)
    putc('-', *fout);
  fprintf(*fout,
	  "\nPedigree  Average        StdDev         Min            Max\n");
  FORLIM = opeds + 1;
  for (ip = 1; ip <= FORLIM; ip++) {
    if (ip <= opeds)
      fprintf(*fout, "%6ld%15.6f%15.6f%15.6f%15.6f\n",
	      ip, avermaxlod[ip - 1], sqrt(variancemaxlod[ip - 1]),
	      smallmaxlod[ip - 1], bigmaxlod[ip - 1]);
    else {
      fprintf(*fout, "\nStudy %15.6f%15.6f%15.6f%15.6f\n",
	      avermaxlod[ip - 1], sqrt(variancemaxlod[ip - 1]),
	      smallmaxlod[ip - 1], bigmaxlod[ip - 1]);
      for (i = 1; i <= 66; i++)
	putc('-', *fout);
      putc('\n', *fout);
    }
  }
  putc('\n', *fout);
  for (i = 1; i <= 66; i++)
    putc('-', *fout);
  fprintf(*fout,
	  "\n   Number of (interpolated) maximum lod scores greater than\n");
  fprintf(*fout, "   a given constant\n");
  for (i = 1; i <= 66; i++)
    putc('-', *fout);
  fprintf(*fout, "\n Constant  Pedigree  Number  Percent\n");
  for (i = 0; i <= 2; i++) {
    FORLIM1 = opeds + 1;
    for (ip = 1; ip <= FORLIM1; ip++) {
      if (ip <= opeds)
	fprintf(*fout, "%8.3f%10ld%8ld%9.3f\n",
		lodlimit[i], ip, clods[ip - 1]
		[i], 100.0 * clods[ip - 1][i] / nrep);
      else
	fprintf(*fout, "%8.3f%10s%8ld%9.3f\n",
		lodlimit[i], "Study", clods[ip - 1]
		[i], 100.0 * clods[ip - 1][i] / nrep);
    }
  }
  for (i = 1; i <= 66; i++)
    putc('-', *fout);
  putc('\n', *fout);
}


/* quadout */
/*quadout*/


Static Void acknowl(ff)
FILE **ff;
{
  fprintf(*ff,
    "Please use these two references when reporting results based on SLINK:\n\n");
  fprintf(*ff,
    "Ott J (1989) Computer-simulation methods in human linkage analysis.\n");
  fprintf(*ff, "Proc Natl Acad Sci USA 86:4175-4178\n\n");
  fprintf(*ff,
	  "Weeks DE, Ott J, Lathrop GM (1990) SLINK: a general simulation\n");
  fprintf(*ff,
    "program for linkage analysis.  Am J Hum Genet 47(3):A204 (Supplement)\n");
}  /* acknowl */


Static Void getparam()
{
  FILE *limit;
  long i;

  /* assign(limit,'LIMIT.DAT'); */
  /*SUN*/
  limit = NULL;
  printf("Reading LIMIT.DAT\n");
  if (limit != NULL)
    limit = freopen("limit.dat", "r", limit);
  else
    limit = fopen("limit.dat", "r");
  if (limit == NULL)
    exit(FileNotFound);
  /*SUN*/
  for (i = 1; i <= 3; i++)   /*Jurg 10/19/91*/
    lodlimit[i - 1] = i;
  if (!(P_eoln(limit) | P_eof(limit)))
    fscanf(limit, "%lg%lg%lg", lodlimit, &lodlimit[1], &lodlimit[2]);
  if (limit != NULL)
    fclose(limit);
}


/* getparam */
/*getparam*/


int main(argc, argv)
int argc;
Char *argv[];
{
  long FORLIM, FORLIM1;
  thisperson *WITH;
  long FORLIM2;
  double TEMP;
  FILE *TEMP1;

  /*  PASCAL_MAIN(argc, argv);*/
  simout = NULL;
  lodfile = NULL;
  msim = NULL;
  speedfile = NULL;
  datafile = NULL;
  ipedfile = NULL;
  outfile = NULL;
  initialize();
  inputdata();   /* Took out readped */
  getparam();
  /*Initialize lod counters */
  for (i = 1; i <= 3; i++) {
    for (ip = 1; ip <= maxped; ip++)
      clods[ip - 1][i - 1] = 0;
  }
  /*Initialize variables for reading from speedfile */
  lastseg = 0;
  segperson = 0;
  lastspeed = 0;
  firsttime = true;
  lasttime = false;
  dolod = false;
  /* Note: MLINK Version 5.1 contains no new statement for the approxstruct */
  /*       while Version 5.03 contains none of the approx stuff. */
  /*       This new statement was added as a safety measure
  => the Mac crashes without it. */
  if (approximate)
    approxstruct = (approxrec *)Malloc(sizeof(approxrec));
  /*   if approximate then approxstruct:=approxpnt(NewPtr(SizeOf(approxrec))); */
  censorstruct = (censorrec *)Malloc(sizeof(censorrec));
  /*  censorstruct:=censorpnt(NewPtr(SizeOf(censorrec)));*/
  for (k = minint; k <= maxcensor; k++)
    censorstruct->censor[k - minint] = false;
  gennustruct = (gennurec *)Malloc(sizeof(gennurec));
  /*  gennustruct:=gennuptr(NewPtr(SizeOf(gennurec)));*/
  firstapprox = true;
  getlocations();

  fprintf(outfile, "LINKAGE/MSIM (V%s) WITH%3ld-POINT", version, nlocus);
  if (sexlink)
    fprintf(outfile, " SEXLINKED DATA\n");
  else
    fprintf(outfile, " AUTOSOMAL DATA\n");
  fprintf(outfile, " ORDER OF LOCI: ");
  FORLIM = nlocus;
  for (i = 1; i <= FORLIM; i++) {
    thissystem = 1;
    while (order[thissystem - 1] != i)
      thissystem++;
    fprintf(outfile, "%3ld", thissystem);
  }
  putc('\n', outfile);
  holdtheta = maletheta->theta[which - 1];
  nrep = 0;
  for (i = 1; i <= maxind; i++)
    person[i] = NULL;
  while (!P_eof(ipedfile)) {   /*for i  begin*/
    nrep++;
    if (nrep % 10 == 0)
      putchar('*');
    else
      putchar('.');
/* p2c: msim.p, line 4809:
 * Note: Using % for possibly-negative arguments [317] */
    fflush(stdout);
    /*  P_ioresult = 0;   *SUN*/
    if (nrep % 50 == 0)
      printf(" %ld replicates done\n", nrep);
/* p2c: msim.p, line 4814:
 * Note: Using % for possibly-negative arguments [317] */
    firsttime = true;
    readpedseg();
    FORLIM1 = totperson;
    for (i = 1; i <= FORLIM1; i++)
      person[i]->store = NULL;
    /* Set all store pointers to NIL in order to be able to use Dispose later */
    readspseg();
    npt = 1;
    if (!risk) {
      maletheta->theta[which - 1] = 0.5;
      /* Store thetas in ztheta array for quadratic interpolation */
      ztheta[npt - 1] = maletheta->theta[which - 1];
      iterpeds();
      maletheta->theta[which - 1] = holdtheta;
      npt++;
    }
    /* repeat
    */
    while (maletheta->theta[which - 1] <= finish) {
      iterpeds();
      /* Store thetas in ztheta array for quadratic interpolation */
      ztheta[npt - 1] = maletheta->theta[which - 1];
      maletheta->theta[which - 1] += inc;
      npt++;
      if (npt <= maxpnt)
	continue;
      printf("\nERROR: maxpnt is too small (%ld)\n", (long)maxpnt);
      fprintf(outfile, "ERROR: maxpnt is too small (%ld)\n", (long)maxpnt);
      printf(" Either increase the constant maxpnt and re-compile\n");
      printf(" or reduce the number of thetas by changing the datafile.dat\n");
      printf("Press return to halt...\n");
      scanf("%*[^\n]");
      getchar();
      exit(0);   /*SUN*/
    }  /* increment on theta[which] */
    /* This resets which to 0 => causes problems. Also don't want to
    read the datafile too many times....
    IF not eof(datafile)
    THEN
    readln(datafile,maletheta^.theta[which],inc,finish)
    ELSE which:=0;
    finish:=finish+0.0001;
    until which=0; */

    FORLIM1 = totperson;
    /* Need to reuse memory once the pedigrees have been processed */
    for (i = 1; i <= FORLIM1; i++) {   /*all persons in this replicate*/
      if (person[i] != NULL) {
	WITH = person[i];
	FORLIM2 = nlocus;
	for (j = 1; j <= FORLIM2; j++) {
	  if (WITH->phen[j - 1] != NULL)
	    Free(WITH->phen[j - 1]);
	}
	if (person[i]->store != NULL)
	  Free(person[i]->store);
	Free(person[i]);
	person[i] = NULL;
      }  /*IF person[i]<>NIL*/
    }

    /* Now update the maxlods for pedigree i based on the zlod[i,ip] for
    pedigree i at points ip => Need not store all lod scores across
    replicates. */
    /* Note that npt is one too big, so */
    npt--;

    /* Now do quadratic interpolation if appropriate */
    if (sexdif || nlocus != 2 || npt <= 2)
      continue;
    FORLIM1 = opeds + 1;
    for (i = 1; i <= FORLIM1; i++) {
      zmax = -100000.0;
      FORLIM2 = npt;
      for (ip = 1; ip <= FORLIM2; ip++) {
	if (zmax <= zlod[i - 1][ip - 1]) {
	  jmax = ip;
	  zmax = zlod[i - 1][ip - 1];
	}
      }
      if (jmax == 1)   /*Maximum at ztheta=0.50*/
	zmax = 0.0;
      else {
	if (jmax == 2)   /*jmax=2 represents smallest theta*/
	  jmax = 3;
	if (jmax == npt)   /*jmax=npt represent largest theta*/
	  jmax = npt - 1;
	zmax = quadmax(ztheta[jmax - 2], ztheta[jmax - 1], ztheta[jmax],
		       zlod[i - 1][jmax - 2], zlod[i - 1]
		       [jmax - 1], zlod[i - 1][jmax]);
      }

      /* Update counts of maximum lod scores greater than 1, 2 or 3 */
      if (zmax >= lodlimit[0])
	clods[i - 1][0]++;
      if (zmax >= lodlimit[1])
	clods[i - 1][1]++;
      if (zmax >= lodlimit[2])
	clods[i - 1][2]++;

      /* Update the average of the maxlods over all replications,
      so averMaxLod[i] gives the average maximum lod score for
      pedigree i so far */
      if (nrep == 1) {
	avermaxlod[i - 1] = zmax;
	variancemaxlod[i - 1] = 0.0;
	smallmaxlod[i - 1] = zmax;
	bigmaxlod[i - 1] = zmax;
      } else {
	rep = nrep;
	/* rep := double(nrep);*/
	avermaxlod[i - 1] = (rep - 1.0) * avermaxlod[i - 1] / rep + zmax / rep;
	TEMP = zmax - avermaxlod[i - 1];
	variancemaxlod[i - 1] = (rep - 1.0) * variancemaxlod[i - 1] / rep +
				TEMP * TEMP / (rep - 1.0);
	smallmaxlod[i - 1] = min(smallmaxlod[i - 1], zmax);
	bigmaxlod[i - 1] = max(bigmaxlod[i - 1], zmax);
      }
    }
  }
  /* of while not eof(ipedfile) */
  if (ipedfile != NULL)
    fclose(ipedfile);
  ipedfile = NULL;   /*SUN*/
  if (speedfile != NULL)
    fclose(speedfile);
  speedfile = NULL;   /*SUN*/
  if (datafile != NULL)
    fclose(datafile);
  datafile = NULL;   /*SUN*/
  if (outfile != NULL)
    fclose(outfile);
  outfile = NULL;   /*SUN*/
  /* Copy most recent simout.dat to msim.dat */
  fprintf(msim, " ********* Data from most recent SIMOUT.DAT *********\n");
  /*assign(simout,'simout.dat');*/
  if (simout != NULL)
    simout = freopen("simout.dat", "r", simout);
  else
    simout = fopen("simout.dat", "r");
  if (simout == NULL)
    exit(FileNotFound);
  /*SUN*/
  while (!P_eof(simout)) {
    chtemp = getc(simout);
    if (chtemp == '\n')
      chtemp = ' ';
    putc(chtemp, msim);
    if (P_eoln(simout)) {
      fscanf(simout, "%*[^\n]");
      getc(simout);
      putc('\n', msim);
    }
  }
  fprintf(msim, " ********* End of most recent SIMOUT.DAT *********\n\n");
  if (simout != NULL)
    fclose(simout);
  simout = NULL;   /*SUN*/
  /* Output actual order to 'msim.dat' */
  fprintf(msim, " ORDER OF LOCI: ");
  FORLIM1 = nlocus;
  for (i = 1; i <= FORLIM1; i++) {
    thissystem = 1;
    while (order[thissystem - 1] != i)
      thissystem++;
    fprintf(msim, "%3ld", thissystem);
  }
  putc('\n', msim);
  putchar('\n');
  if (nlocus == 2)
    fprintf(msim, " Average Lod Scores at Given Thetas\n");
  else
    fprintf(msim, " Average Multipoint Lod Scores at Given Thetas\n");
  putc('\n', msim);
  if (nlocus == 2)
    printf(" Average Lod Scores at Given Thetas\n");
  else
    printf(" Average Multipoint Lod Scores at Given Thetas\n");
  printf("\nNumber of replicates = %6ld\n", nrep);
  fprintf(msim, "Number of replicates = %6ld\n", nrep);
  putchar('\n');
  putc('\n', msim);
  for (i = 1; i <= 66; i++)
    putc('-', msim);
  putc('\n', msim);
  for (i = 1; i <= 66; i++)
    putchar('-');
  putchar('\n');
  maletheta->theta[which - 1] = holdtheta;
  /*Don't need to write out the lods at
  theta = 0.5*/
  npt = 2;
  while (maletheta->theta[which - 1] <= finish) {
    if (sexdif)
      printf("MALE THETAS   ");
    else
      printf("THETAS ");
    FORLIM1 = nlocus;
    for (i = 1; i < FORLIM1; i++)
      printf("%6.3f", maletheta->theta[i - 1]);
    if (interfer)
      printf("%6.3f", maletheta->theta[nlocus - 1]);
    putchar('\n');
    if (sexdif) {
      printf("FEMALE THETAS ");
      FORLIM1 = nlocus;
      for (i = 1; i < FORLIM1; i++)
	printf("%6.3f", femaletheta->theta[i - 1]);
      if (interfer)
	printf("%6.3f", femaletheta->theta[nlocus - 1]);
      putchar('\n');
    }
    if (sexdif)
      fprintf(msim, "MALE THETAS   ");
    else
      fprintf(msim, "THETAS ");
    FORLIM1 = nlocus;
    for (i = 1; i < FORLIM1; i++)
      fprintf(msim, "%6.3f", maletheta->theta[i - 1]);
    if (!sexdif && nlocus == 2) {
      if (interfer)
	fprintf(msim, "%6.3f", maletheta->theta[nlocus - 1]);
    }
    putc('\n', msim);
    if (sexdif) {
      fprintf(msim, "FEMALE THETAS ");
      FORLIM1 = nlocus;
      for (i = 1; i < FORLIM1; i++)
	fprintf(msim, "%6.3f", femaletheta->theta[i - 1]);
      if (interfer)
	fprintf(msim, "%6.3f", femaletheta->theta[nlocus - 1]);
      putc('\n', msim);
    }
    for (i = 1; i <= 66; i++)
      putc('-', msim);
    putc('\n', msim);
    for (i = 1; i <= 66; i++)
      putchar('-');
    printf("\nPedigree  Average        StdDev         Min            Max\n");
    fprintf(msim,
	    "Pedigree  Average        StdDev         Min            Max\n");
    FORLIM1 = opeds + 1;
    for (ip = 1; ip <= FORLIM1; ip++) {
      if (ip <= opeds) {
	printf("%6ld%15.6f%15.6f%15.6f%15.6f\n",
	       ip, aver[ip - 1][npt - 1], sqrt(variance[ip - 1][npt - 1]),
	       small[ip - 1][npt - 1], big[ip - 1][npt - 1]);
	fprintf(msim, "%6ld%15.6f%15.6f%15.6f%15.6f\n",
		ip, aver[ip - 1][npt - 1], sqrt(variance[ip - 1][npt - 1]),
		small[ip - 1][npt - 1], big[ip - 1][npt - 1]);
      } else {
	putc('\n', msim);
	printf("\nStudy %15.6f%15.6f%15.6f%15.6f\n",
	       aver[ip - 1][npt - 1], sqrt(variance[ip - 1][npt - 1]),
	       small[ip - 1][npt - 1], big[ip - 1][npt - 1]);
	fprintf(msim, "Study %15.6f%15.6f%15.6f%15.6f\n",
		aver[ip - 1][npt - 1], sqrt(variance[ip - 1][npt - 1]),
		small[ip - 1][npt - 1], big[ip - 1][npt - 1]);
	for (i = 1; i <= 66; i++)
	  putc('-', msim);
	putc('\n', msim);
	for (i = 1; i <= 66; i++)
	  putchar('-');
	putchar('\n');
      }
    }
    maletheta->theta[which - 1] += inc;
    npt++;
  }
  /* npt is one too big, so...*/
  npt--;
  for (i = 1; i <= 66; i++)
    putc('-', msim);
  putc('\n', msim);
  for (i = 1; i <= 66; i++)
    putchar('-');
  putchar('\n');
  /* Now output the quadratic interpolation results if appropriate */
  if (!sexdif && nlocus == 2 && npt > 2) {
    quadout(&msim);
    TEMP1 = stdout;
/* p2c: msim.p, line 5051:
 * Note: Taking address of stdout; consider setting VarFiles = 0 [144] */
    quadout(&TEMP1);
  }
  acknowl(&msim);
  if (outfile != NULL)
    fclose(outfile);
  if (ipedfile != NULL)
    fclose(ipedfile);
  if (datafile != NULL)
    fclose(datafile);
  if (speedfile != NULL)
    fclose(speedfile);
  if (msim != NULL)
    fclose(msim);
  if (lodfile != NULL)
    fclose(lodfile);
  if (simout != NULL)
    fclose(simout);
  exit(EXIT_SUCCESS);
}  /* msim */



/* End. */
