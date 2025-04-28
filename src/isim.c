#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/time.h>

/* Output from p2c, the Pascal-to-C translator */
/* From input file "isim.p" */


/*#include <p2c/p2c.h>*/
/* p2c: isim.p, line 1: 
 * Note: Unexpected name "isim" in program header [262] */
/* p2c: isim.p, line 1: 
 * Note: Unexpected name "outfile" in program header [262] */
/* p2c: isim.p, line 1: 
 * Note: Unexpected name "datafile" in program header [262] */
/* p2c: isim.p, line 1: 
 * Note: Unexpected name "speedfile" in program header [262] */
/* p2c: isim.p, line 1: 
 * Note: Unexpected name "ipedfile" in program header [262] */
/* p2c: isim.p, line 1: 
 * Note: Unexpected name "simout" in program header [262] */
/* p2c: isim.p, line 1: 
 * Note: Unexpected name "final" in program header [262] */


/*3 July 1992*/

/*ISIM is a modified version of ILINK, modified by Daniel E. Weeks with
  the help of Mark Lathrop during June 1989.  It is a companion program
  to the generalized simulation program SLINK.

  Update History:
   8/90 - Fixed msegsex and msegsexf. Made compatible with version 5.1.

  INPUT:
  The three typical ILINK-type datafiles, where ipedfile and speedfile have
  been produced by UNKNOWN.
   datafile.dat
   ipedfile.dat
   speedfile.dat
   simout.dat

  OUTPUT:
   isim.dat      => Summary statistics
   outfile.dat
   final.dat

  NOTES:

  1) The stream file outputs have been removed.

  2) This program has been modified to process a replicate of the
  original pedigrees at a time.  This should make it possible
  to analyze a large number of replicates even on a computer
  without much memory. In addition, this allows calculation
  of the maximum lod score on each pedigree.

  3) The statistical summary is placed in the file named 'isim.dat',
  while the usual output of ILINK is written to outfile.
  */


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


typedef double real;


#define version         "2.51"
    /*VERSION OF ISIM, based on version 4.91 of LINKAGE*/

#define prn             false
    /* Turns off most outfile output for more efficient running */
#define finalprn        false
    /*Turns off output to final.dat for more efficient running */
/* SOME USER DEFINED CONSTANTS */
/*THE PROGRAM WILL TELL YOU IF THE FOLLOWING TWO CONSTANTS CAN BE REDUCED*/
/*IF THE PROGRAM TERMINATES WITH AN ERROR IN RECOMBINATION INCREASE MAXNEED*/

#define maxneed         800
    /*MAXIMUM NUMBER OF RECOMBINATION PROBABILITIES*/
/*THE FOLLOWING SHOULD BE LARGER THAN MININT*/
#define maxcensor       10000   /*MAXIMUM FOR CENSORING ARRAY*/
#define maxlocus        4   /*MAXIMUM NUMBER OF LOCI */
#define maxseg          8   /*2 TO THE POWER maxlocus-1 */
#define maxall          6   /*MAX NUMBER OF ALLELES AT A SINGLE LOCUS*/
#define maxhap          16   /*MAXIMUM NUMBER OF HAPLOTYPES*/
/* = n1 x n2 x ... ETC. WHERE ni = NUMBER OF ALLELES LOCUS I*/

#define maxfem          (maxhap * (maxhap + 1) / 2)
    /*MAX. NO. OF JOINT GENOTYPES FOR A FEMALE*/
#define maxmal          maxfem
    /*MAXIMUM NUMBER OF JOINT GENOTYPES FOR A MALE*/
/* = maxfem (AUTOSOMAL) OR maxhap (SEXLINKED)*/

#define maxind          600   /*MAXIMUM NUMBER OF INDIVIDUALS*/
#define maxped          65   /*MAXIMUM NUMBER OF PEDIGREES*/
#define maxchild        20   /*MAXIMUM NUMBER OF FULLSIBS IN A SIBSHIP*/
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
#define maxliab         20   /*MAXIMUM NUMBER OF LIABILITY CLASSES */
/* BINARY (FACTOR UNION) SYSTEM */
#define maxfact         10
    /*MAXIMUM NUMBER OF BINARY CODES AT A SINGLE LOCUS*/
/* OTHERS */

#define scale           3.0   /*SCALE FACTOR*/
#define scalemult       2.0   /*SCALE WEIGHT FOR EACH LOCUS INCLUDED*/

#define fitmodel        false
    /*TRUE IF ESTIMATING PARAMETERS OTHER THAN REC*/
#define byfamily        true   /*GIVE LOD SCORES BY FAMILY*/

#define zerolike        (-1.0e20)
/*FOR INCONSISTENT DATA OR RECOMBINATION */
/*SUN*/
#define log10_          2.30259

#define minint          (-32767)   /*MINIMUM ALLOWED INTEGER*/
/*GRADIENT APPROXIMATIONS*/

#define approximate     false
    /*do not change or else change approxrec below*/

#define epsilon         1.0e-4   /*SUN*/
/*GEMINI*/

#define maxn            20   /* MAXIMUM NUMBER OF ITERATED PARAMETERS */
#define nbit            50   /*NUMBER OF BITS OF MACHINE PRECISION*/

#define tolconst        1.0e-3
/*TOLERANCE ON LIKELIHOOD FOR STOP*/
/*SUN*/
/*DO NOT CHANGE THE FOLLOWING CONSTANTS*/
#define tolg            1.0e-35   /*SUN*/
#define xpmcon          1.0e30   /*SUN*/
#define clbcon          1.0e10   /*SUN*/
#define tbndcon         1.0e5   /*SUN*/
#define curv            0.75
#define small           1.0e-10   /*SUN*/

/*next two typedefs come from p2c*/
typedef char boolean;
typedef short uchar;


typedef struct gennurec {
  long genenumber[maxhap][maxhap];
} gennurec;

typedef struct approxrec {
  boolean approxarray[1][1];   /*Was packed*/
} approxrec;

typedef struct censorrec {
  boolean censor[maxcensor - minint + 1];   /*Was packed*/
} censorrec;

typedef real genotype[maxfem];
typedef real covmatrix[maxtrait][maxtrait];
typedef char hapvector[maxlocus];   /*Was packed*/
typedef long mutarray[3][2];
typedef real thesemeans[maxtrait];
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
  real freq[maxall];
  struct locusvalues *privlocus;
  locustype which;
  union {
    struct {
      real pen[maxall + 1][maxall][3][maxliab];
      long nclass;
    } U0;
    struct {
      long ntrait;
      means pm;
      covmatrix vmat;
      real det, contrait, conmat;
    } U1;
    phenarray allele;
  } UU;
} locusvalues;

typedef struct phenotype {
  locustype which;
  union {
    binset phenf;
    struct {
      real x[maxtrait];
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
typedef real thetarray[maxlocus];
typedef real happrob[maxneed];

typedef struct thetavalues {
  thetarray theta;
  happrob segprob;
} thetavalues;

/*GEMINI*/
typedef real vector[maxn];
typedef vector matrix[maxn];


Static boolean itmax;   /*Flag indicating iterated parameters*/
Static long nrep, j;   /*nrep is a counter of number of replications*/
Static real rep, ltemp, aver, variance, big, smll, ddist, segscale, mutmale,
	    mutfemale, like, alike, distratio, penlike, tol, tmin, gh, trupb,
	    ptg, tbnd, f, fsmf, fsav2, t, hx, xsave, fxph, fxmh, xp, xpm, ytp;
Static long clods[3];
Static real lodlimit[3];
/*These are where various statistics are stored*/
Static long nfactor[maxlocus];
Static long opeds;   /*Number of original pedigrees*/
Static boolean infile;
    /*TRUE if currently in the middle of the pedigree file*/
Static long oldped;   /*oldped stores number of next pedigree to be read.*/
Static pathway thispath;
Static boolean informative[maxped];
Static boolean rare[maxfem], risk1[maxfem], risk2[maxfem];   /*Was packed*/
Static boolean riskmale[maxhap];   /*Was packed*/
Static approxrec *approxstruct;
Static censorrec *censorstruct;
Static long thisc;
Static boolean malechild[maxchild];
Static thisarray *thischild[maxchild];
Static long nchild;
Static char seghap1[maxfem], seghap2[maxfem];
Static uchar segstart[maxfem];
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
Static long risksys, mutsys, nlocus, lastpriv, i, lastspeed, lastseg,
	    segperson;
/*speedfile modification above; added variable lastspeed,lastseg,segperson*/
Static char nuhap;
Static uchar fgeno, mgeno;
Static long nuped, totperson;
Static boolean interfer, disequi, sexlink, risk, sexdif, readfemale,
	       inconsistent, dolod, firstapprox, firsttime, lasttime;
Static FILE *isim, *outfile, *datafile, *speedfile, *ipedfile, *simout;
Static Char chtemp;
/*isim*/
Static FILE *final;
Static boolean mapping, penalty;
/*GEMINI*/
Static boolean active;
Static enum {
  go, quit, restart
} continue_;

/*VARIABLES FROM COMMON MIN1*/
Static matrix tl;
Static vector d, g, gsave, y, p;

/*VARIABLE FROM COMMON MIN2*/
Static long nit, nfe, idg, idif, isw, iret, ibnd, ihess, ivar, ihx, maxit;

/*DIMENSIONED IN GEMINI*/
Static vector holdxall, xall, x, gv, se;
Static long holditp[maxn], itp[maxn];
Static real bnd[maxn][2];
Static long nall, n, icall;

/*ARRAYS IN UPDATE*/
Static vector wtil, ztil, w, z, wtjp1, ztjp1, s, dp;

Static real nu, muj, lambj, lambj2, sb, sc, sd, fbcd, alpha, sa, thet1, thet2,
	    aa, bb, cc, del2, alph1, alph2, rdiv, eta, aj, thj, bj, gamlj,
	    betlj, del;

Static long iswup;

/*FROM STEP*/
Static vector xt;

Static real fsave, sumt, twot, ft, f2t, ft2, scal;

/*FROM CHKBND */
Static real clb, check, eh, teq;

/*OTHERS*/
Static long itsys;
/*USUALLY IN INIB*/
Static matrix bmat;

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


Static Void outcontrol(z)
double *z;
{
  long i, j, FORLIM;

  j = 0;
  FORLIM = n;
  for (i = 0; i < FORLIM; i++) {
    if (prn)
      fprintf(outfile, "% .5E", z[i]);
    j++;
    if (j == 4) {
      if (prn)
	putc('\n', outfile);
      j = 0;
    }
  }
  if (prn)
    putc('\n', outfile);
}  /* outcontrol */


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
  if (*dist != 10.0)
    return ((1.0 - exp(-2 * *dist)) / 2.0);
  else
    return 0.5;
}  /* invdist */


Static Void invert(m, n, det)
real (*m)[maxtrait];
long n;
double *det;
{
  covmatrix v;
  real val;
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
  real p1, p2, p3, p4;
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
  long whichhap, start, length, i, j, k, FORLIM2;
  real recval, val;

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
  real newval, recval;

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
  real dist;
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
  if (!prn)
    return;
  if (firsttime) {
    if (V.here < maxneed)
      printf("Maxneed can be reduced to %12ld\n", V.here);
  }
}  /* recombination */


Static Void setparam()
{
  long i, j, k, l, m;
  thetavalues *WITH;
  long FORLIM;
  locusvalues *WITH1;
  thisarray *WITH2;
  long FORLIM1, FORLIM2;

  penalty = false;
  k = 0;
  WITH = maletheta;
  FORLIM = nlocus - 2;
  for (i = 0; i <= FORLIM; i++) {
    k++;
    WITH->theta[i] = xall[k - 1];
  }
  if (interfer && !mapping) {
    k++;
    maletheta->theta[nlocus - 1] = xall[k - 1];
  }
  if (sexdif) {
    if (readfemale) {
      WITH = femaletheta;
      FORLIM = nlocus - 2;
      for (i = 0; i <= FORLIM; i++) {
	k++;
	WITH->theta[i] = xall[k - 1];
      }
      if (interfer && !mapping) {
	k++;
	femaletheta->theta[nlocus - 1] = xall[k - 1];
      }
    } else {
      k++;
      distratio = xall[k - 1];
    }
  }
  if (itsys == 0)
    return;
  WITH1 = thislocus[itsys - 1];
  if (disequi) {
    WITH2 = hapfreq;
    FORLIM = nuhap;
    for (i = 0; i < FORLIM; i++)
      WITH2->genarray[i] = xall[i + k];
    WITH2->genarray[nuhap - 1] = 1.0;
    FORLIM = nuhap - 2;
    for (i = 0; i <= FORLIM; i++)
      WITH2->genarray[nuhap - 1] -= WITH2->genarray[i];
    penalty = (WITH2->genarray[nuhap - 1] < 0.0);
    k += nuhap - 1;
  } else {
    FORLIM = WITH1->nallele - 2;
    for (i = 0; i <= FORLIM; i++)
      WITH1->freq[i] = xall[i + k];
    WITH1->freq[WITH1->nallele - 1] = 1.0;
    FORLIM = WITH1->nallele - 2;
    for (i = 0; i <= FORLIM; i++)
      WITH1->freq[WITH1->nallele - 1] -= WITH1->freq[i];
    penalty = (WITH1->freq[WITH1->nallele - 1] < 0.0);
    k += WITH1->nallele - 1;
  }
  if (WITH1->which == affection) {
    FORLIM = WITH1->UU.U0.nclass;
    for (l = 0; l < FORLIM; l++) {
      FORLIM1 = WITH1->nallele;
      for (i = 1; i <= FORLIM1; i++) {
	FORLIM2 = WITH1->nallele;
	for (j = i; j <= FORLIM2; j++) {
	  k++;
	  WITH1->UU.U0.pen[i][j - 1][2][l] = xall[k - 1];
	  WITH1->UU.U0.pen[i][j - 1][1][l] = 1.0 - xall[k - 1];
	  for (m = 0; m <= 2; m++)
	    WITH1->UU.U0.pen[j][i - 1][m][l] = WITH1->UU.U0.pen[i][j - 1][m]
	      [l];
	}
      }
      if (sexlink) {
	FORLIM1 = WITH1->nallele;
	for (i = 0; i < FORLIM1; i++) {
	  k++;
	  WITH1->UU.U0.pen[0][i][2][l] = xall[k - 1];
	  WITH1->UU.U0.pen[0][i][1][l] = 1.0 - xall[k - 1];
	}
      }
    }
    return;
  }
  if (WITH1->which != quantitative)
    return;
  FORLIM = WITH1->UU.U1.ntrait;
  /*USER MUST CHANGE FOR SYSTEMS WITH MORE THAN 2 ALLELES*/
  for (i = 0; i < FORLIM; i++) {
    WITH1->UU.U1.pm[1][0][i] = xall[k];
    WITH1->UU.U1.pm[2][1][i] = xall[k + 1] + xall[k];
    WITH1->UU.U1.pm[1][1][i] = xall[k + 2] * xall[k + 1] + xall[k];
    WITH1->UU.U1.pm[0][0][i] = xall[k];
    WITH1->UU.U1.pm[0][1][i] = WITH1->UU.U1.pm[1][1][i];
    WITH1->UU.U1.pm[2][0][i] = WITH1->UU.U1.pm[1][1][i];
    k += 3;
  }
  FORLIM = WITH1->UU.U1.ntrait;
  for (i = 0; i < FORLIM; i++) {
    FORLIM1 = WITH1->UU.U1.ntrait;
    for (j = i; j < FORLIM1; j++) {
      k++;
      WITH1->UU.U1.vmat[i][j] = xall[k - 1];
      WITH1->UU.U1.vmat[j][i] = xall[k - 1];
    }
  }
  invert(WITH1->UU.U1.vmat, WITH1->UU.U1.ntrait, &WITH1->UU.U1.det);
  WITH1->UU.U1.det = 1 / sqrt(WITH1->UU.U1.det);
}  /* setparam */


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
{  /*getlocations??*/
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
}


typedef enum {
  nobnd, probnd, zerobnd, logbnd
} boundary;


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
    printf("Error in datafile, line 1, item 4: Program code not for ISIM/ILINK\n");
    break;
  }
  printf("Press return to halt...\007\n");
  scanf("%*[^\n]");
  getchar();
  exit(0);   /*SUN*/
}  /* inputerror */


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

  case 3:
    printf(
      "Non-standard affection status%4ld interpreted as normal in pedigree record%5ld\n",
      par2, par1);
    break;
  }
  printf("Press return to continue...\007\n");
  scanf("%*[^\n]");
  getchar();
}  /* inputwarning */


Local Void setiterate()
{
  boundary bndtype[maxn];
  long i, j, k, l;
  thetavalues *WITH;
  long FORLIM;
  locusvalues *WITH1;
  long FORLIM1, FORLIM2;

  fscanf(datafile, "%ld%*[^\n]", &i);
  getc(datafile);
  if (i > nlocus)
    inputerror(23L, nlocus, i);
  if (i < 0)
    inputerror(24L, i, i);
  if (i == 0)
    itsys = 0;
  else {
    itsys = 1;
    while (order[itsys - 1] != i)
      itsys++;
  }
  k = 0;
  WITH = maletheta;
  FORLIM = nlocus - 2;
  for (i = 0; i <= FORLIM; i++) {
    k++;
    if (k > maxn)
      inputerror(25L, (long)maxn, k);
    xall[k - 1] = WITH->theta[i];
    if (interfer && !mapping)
      bndtype[k - 1] = logbnd;
    else
      bndtype[k - 1] = probnd;
  }
  if (interfer && !mapping) {
    k++;
    if (k > maxn)
      inputerror(25L, (long)maxn, k);
    xall[k - 1] = maletheta->theta[nlocus - 1];
    bndtype[k - 1] = logbnd;
  }
  if (sexdif) {
    if (readfemale) {
      WITH = femaletheta;
      FORLIM = nlocus - 2;
      for (i = 0; i <= FORLIM; i++) {
	k++;
	if (k > maxn)
	  inputerror(25L, (long)maxn, k);
	xall[k - 1] = WITH->theta[i];
	if (interfer && !mapping)
	  bndtype[k - 1] = logbnd;
	else
	  bndtype[k - 1] = probnd;
      }
      if (interfer && !mapping) {
	k++;
	if (k > maxn)
	  inputerror(25L, (long)maxn, k);
	xall[k - 1] = femaletheta->theta[nlocus - 1];
	bndtype[k - 1] = logbnd;
      }
    } else {
      k++;
      if (k > maxn)
	inputerror(25L, (long)maxn, k);
      xall[k - 1] = distratio;
      bndtype[k - 1] = zerobnd;
    }
  }
  nall = k;
  if (itsys != 0) {
    WITH1 = thislocus[itsys - 1];
    if (disequi) {
      FORLIM = nuhap - 2;
      for (i = 0; i <= FORLIM; i++)
	xall[i + k] = hapfreq->genarray[i];
      nall = nuhap + k - 1;
    } else {
      FORLIM = WITH1->nallele - 2;
      for (i = 0; i <= FORLIM; i++)
	xall[i + k] = WITH1->freq[i];
      nall = WITH1->nallele + k - 1;
    }
    if (nall > maxn)
      inputerror(25L, (long)maxn, nall);
    FORLIM = nall;
    for (i = k; i < FORLIM; i++)
      bndtype[i] = probnd;
    k = nall;

    if (WITH1->which == quantitative) {
      nall += WITH1->UU.U1.ntrait * (WITH1->UU.U1.ntrait + 1) / 2 +
	      WITH1->UU.U1.ntrait * 3;

    } else if (WITH1->which == affection) {
      nall += WITH1->UU.U0.nclass * WITH1->nallele * (WITH1->nallele + 1) / 2;
      if (sexlink)
	nall += WITH1->nallele;
    }
    if (nall > maxn)
      inputerror(25L, (long)maxn, nall);

    if (WITH1->which == affection) {
      FORLIM = WITH1->UU.U0.nclass;
      for (l = 0; l < FORLIM; l++) {
	FORLIM1 = WITH1->nallele;
	for (i = 1; i <= FORLIM1; i++) {
	  FORLIM2 = WITH1->nallele;
	  for (j = i - 1; j < FORLIM2; j++) {
	    k++;
	    if (k > maxn)
	      inputerror(25L, (long)maxn, k);
	    xall[k - 1] = WITH1->UU.U0.pen[i][j][2][l];
	    bndtype[k - 1] = probnd;
	  }
	}
	if (sexlink) {
	  FORLIM1 = WITH1->nallele;
	  for (i = 0; i < FORLIM1; i++) {
	    k++;
	    if (k > maxn)
	      inputerror(25L, (long)maxn, k);
	    xall[k - 1] = WITH1->UU.U0.pen[0][i][2][l];
	    bndtype[k - 1] = probnd;
	  }
	}
      }
    } else if (WITH1->which == quantitative) {
      invert(WITH1->UU.U1.vmat, WITH1->UU.U1.ntrait, &WITH1->UU.U1.det);
      FORLIM = WITH1->UU.U1.ntrait;
      /*USER MUST MODIFY FOR A SYSTEM WILL MORE THAN 2 ALLELES */
      for (l = 0; l < FORLIM; l++) {
	FORLIM1 = WITH1->nallele;
	for (i = 1; i <= FORLIM1; i++) {
	  FORLIM2 = WITH1->nallele;
	  for (j = i; j <= FORLIM2; j++) {
	    k++;
	    if (i == 1 && j == 1)
	      xall[k - 1] = WITH1->UU.U1.pm[i][j - 1][l];
	    else if (i == 1)
	      xall[k - 1] = WITH1->UU.U1.pm[j][j - 1][l] - xall[k - 2];
	    else
	      xall[k - 1] = (WITH1->UU.U1.pm[i - 1][j - 1][l] - xall[k - 3]) /
			    xall[k - 2];
	    bndtype[k - 1] = nobnd;
	  }
	}
      }
      FORLIM = WITH1->UU.U1.ntrait;
      for (i = 1; i <= FORLIM; i++) {
	FORLIM1 = WITH1->UU.U1.ntrait;
	for (j = i; j <= FORLIM1; j++) {
	  k++;
	  if (k > maxn)
	    inputerror(25L, (long)maxn, k);
	  xall[k - 1] = WITH1->UU.U1.vmat[i - 1][j - 1];
	  if (i == j)
	    bndtype[k - 1] = zerobnd;
	  else
	    bndtype[k - 1] = nobnd;
	}
      }
    }

  }
  FORLIM = nall;
  for (i = 0; i < FORLIM; i++)
    itp[i] = 0;
  i = 0;
  while ((i <= nall) & (!P_eoln(datafile))) {
    i++;
    fscanf(datafile, "%ld", &itp[i - 1]);
  }
  n = 0;
  FORLIM = nall;
  for (i = 0; i < FORLIM; i++) {
    if (itp[i] == 1) {
      n++;
      x[n - 1] = xall[i];
    }
  }
  FORLIM = nall;
  for (i = 0; i < FORLIM; i++) {
    bnd[i][0] = -1000.0;
    bnd[i][1] = 1000.0;
    if (bndtype[i] == probnd) {
      bnd[i][0] = 0.0;
      bnd[i][1] = 1.0;
    } else if (bndtype[i] == logbnd) {
      bnd[i][0] = -20.0;
      bnd[i][1] = 20.0;
    } else if (bndtype[i] == zerobnd)
      bnd[i][0] = 0.0;
  }


}  /* setiterate */

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
  for (i = 1; i <= FORLIM; i++)
    WITH->UU.allele[i - 1] = 1L << ((int)i);
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

/*setiterate*/

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
  if (V.i != 3)
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
  printf("YOU ARE USING LINKAGE/ISIM (V%s) WITH%3ld-POINT", version, nlocus);
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
  for (V.i = 1; V.i <= FORLIM; V.i++)
    segscale *= scalemult;
}  /* readloci */


/* getlocations */
/*getlocations??*/


Static Void inputdata()
{
  /*readloci*/

  printf("Reading DATAFILE.DAT\n");
  if (datafile != NULL)
    datafile = freopen("datafile.dat", "r", datafile);
  else
    datafile = fopen("datafile.dat", "r");
  if (datafile == NULL)
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
  printf("Reading IPEDFILE.DAT\n");
  if (ipedfile != NULL)
    ipedfile = freopen("ipedfile.dat", "r", ipedfile);
  else
    ipedfile = fopen("ipedfile.dat", "r");
  if (ipedfile == NULL)
    exit(FileNotFound);
  if (outfile != NULL) {
    /*SUN*/
    outfile = freopen("outfile.dat", "w", outfile);
  } else
    outfile = fopen("outfile.dat", "w");
  if (outfile == NULL)
    exit(FileNotFound);
  if (final != NULL) {
    /*SUN*/
    final = freopen("final.dat", "w", final);
  } else
    final = fopen("final.dat", "w");
  if (final == NULL)
    exit(FileNotFound);
  /*SUN*/
  readloci();
  setiterate();
}  /* inputdata */


Static Void readspseg()
{  /*readspeed*/
  /*speedfile modification*/
  long i, j, a, b, sys;
  Char ch;
  information *WITH;
  long FORLIM;

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
  long i, sequence;
  long startped[maxped], endped[maxped];
} ;

/*Exit from the procedure readpedseg*/

Local Void inputerror_(nerror, par1, par2, LINK)
long nerror, par1, par2;
struct LOC_readpedseg *LINK;
{
  printf("Fatal error detected in procedure inputdata\n");
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
      LINK->i);
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
  long thisval;
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
  real xval;
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
	count = 0;
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
  /*       person[id]:=ind(NewPtr(SizeOf(thisperson)));*/
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
/*readspeed*/


Static Void readpedseg()
{
  struct LOC_readpedseg V;
  long newid, sex_, profield, newped, nuperson, thisone, thisped;
  thisperson *holdloop;
  long FORLIM;
  thisperson *WITH;

  /*multimarriage*/

  for (V.i = 0; V.i <= maxind; V.i++)
    person[V.i] = NULL;
  V.sequence = 0;
  nuperson = 0;
  nuped = 1;
  for (V.i = 1; V.i <= maxloop; V.i++) {
    looppers[nuped - 1][V.i - 1][0] = NULL;
    looppers[nuped - 1][V.i - 1][1] = NULL;
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
    Test for: newped<>thisped and nuped=oped => save newped in oldped
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
    for (V.i = 1; V.i <= maxloop; V.i++) {
      looppers[nuped - 1][V.i - 1][0] = NULL;
      looppers[nuped - 1][V.i - 1][1] = NULL;
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
    for (V.i = 1; V.i <= maxloop; V.i++) {
      if (looppers[newped][V.i - 1][0] == NULL)
	looppers[newped][V.i - 1][1] = NULL;
      else {
	looppers[newped][V.i - 1][0]->inloop = V.i;
	looppers[newped][V.i - 1][1]->inloop = V.i;
	if (looppers[newped][V.i - 1][0]->pa == NULL &&
	    looppers[newped][V.i - 1][1]->pa != NULL) {
	  holdloop = looppers[newped][V.i - 1][0];
	  looppers[newped][V.i - 1][0] = looppers[newped][V.i - 1][1];
	  looppers[newped][V.i - 1][1] = holdloop;
	}
      }
    }
  }
  /*speedfile modification*/
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
  real val;
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
  double val;
  real newval;
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
  double val;
  real newval;
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
  real homo, hetero;
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
  real pf, ps;
  thetavalues *firstsex, *secondsex;
} ;

Local Void getapprox(LINK)
struct LOC_seg *LINK;
{
  long first;
  real maxval;   /*getapprox*/
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
      [first] = (WITH->genarray[first] > maxval);
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
  real val, temp2;
  real temp[maxchild];   /*msegsex*/
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
/* p2c: isim.p, line 4619: Note:
 * Line breaker spent 0.0+1.00 seconds, 5000 tries on line 3873 [251] */
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
/* p2c: isim.p, line 4619: 
 * Note: Line breaker spent 0.0 seconds, 5000 tries on line 3932 [251] */
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
  real val, temp1;
  real temp[maxchild][maxseg];
  real temp2[maxseg];   /*msegsexf*/
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
	  val = (1 - LINK->pf) * (1 - LINK->ps) * (WITH2->genarray[g1 - 1] +
		  WITH2->genarray[g2 - 1]) + (1 - LINK->pf) * LINK->ps *
		(WITH2->genarray[g3 - 1] + WITH2->genarray[g4 - 1]) +
	      LINK->pf * (1 - LINK->ps) * (WITH2->genarray[g5 - 1] + WITH2->
		    genarray[g6 - 1]) + LINK->pf * LINK->ps *
		(WITH2->genarray[g7 - 1] + WITH2->genarray[g8 - 1]);
/* p2c: isim.p, line 4619: 
 * Note: Line breaker spent 0.0 seconds, 5000 tries on line 4018 [251] */
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
  real val, temp2;
  real temp[maxchild];   /*segsex*/
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
  real val, temp1;
  real temp[maxchild][maxseg];
  real temp2[maxseg];
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
  real val, temp1, temp2;
  real temp[maxchild];
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
  real val, temp1;
  real temp[maxchild][maxseg];
  real temp2[maxseg];   /*msegfast*/
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
/* p2c: isim.p, line 4619: 
 * Note: Line breaker spent 0.0 seconds, 5000 tries on line 4377 [251] */
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
/* p2c: isim.p, line 4619: 
 * Note: Line breaker spent 0.0 seconds, 5000 tries on line 4401 [251] */
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
/* p2c: isim.p, line 4619: Note:
 * Line breaker spent 0.0+1.00 seconds, 5000 tries on line 4437 [251] */
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
  real val, temp1, temp2;
  real temp[maxchild];   /*msegfun*/
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
/* p2c: isim.p, line 4619: 
 * Note: Line breaker spent 0.0 seconds, 5000 tries on line 4543 [251] */
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
/* p2c: isim.p, line 4619: 
 * Note: Line breaker spent 0.0 seconds, 5000 tries on line 4565 [251] */
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
/* p2c: isim.p, line 4619: 
 * Note: Line breaker spent 0.0 seconds, 5000 tries on line 4598 [251] */
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
  real val, temp1;
  real temp[maxchild][maxseg];
  real temp2[maxseg];   /*segfast*/
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
  real segval;
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
  real segval, val;
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
  real segval;
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
  uchar here;
  genotype gene;
  real val, temp2;
  long j, first, second;   /*segsexdown*/
  uchar FORLIM;
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
  real segval;
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
  real segval, val;
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
  real segval;
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
  real segval;
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
  uchar here;
  genotype gene;
  real val, temp1, temp2;
  long f1, f2, s1, s2, i, j, first, second;
  uchar FORLIM;
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
  uchar here;
  genotype gene;
  real val, temp2;
  long ms2, ms1, mf, j, first, second;   /*msegsexdown*/
  uchar FORLIM;
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
  uchar here;
  genotype gene;
  real val, temp, temp1, temp2;
  long i, j, first, second, f1, f2, s1, s2, mf1, mf2, ms1, ms2;   /*msegdown*/
  uchar FORLIM;
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
  real normal;

  LINK->homo /= like;
  LINK->hetero /= like;
  normal = 1 - LINK->homo - LINK->hetero;
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
  real tmplike;
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


Static Void fun(f, x)
double *f;
double *x;
{
  long i, k, FORLIM;

  k = 0;
  FORLIM = nall;
  for (i = 0; i < FORLIM; i++) {
    if (itp[i] == 1) {
      k++;
      xall[i] = x[k - 1];
    }
  }

  setparam();
  if (penalty) {
    *f = 1.5 * penlike;
    return;
  }
  FORLIM = totperson;
  for (i = 1; i <= FORLIM; i++)
    person[i]->done = false;
  FORLIM = totperson;
  for (i = 1; i <= FORLIM; i++)
    person[i]->firstpass = true;
  recombination();
  alike = 0.0;
  thisc = minint;
  FORLIM = nuped;
  for (i = 1; i <= FORLIM; i++) {
    likelihood(i, proband[i - 1]);
    alike += like;
  }
  if (prn) {
    if (firsttime) {
      if (thisc < maxcensor)
	printf("Maxcensor can be reduced to %12ld\n", thisc);
      else if (thisc > maxcensor)
	printf("You may gain efficiency by increasing maxcensor\n");
    }
  }
  firsttime = false;
  *f = -2 * alike;
  penlike = *f;
  firsttime = false;
}  /* fun */


/* Local variables for out: */
struct LOC_out {
  real thisval[maxped];
} ;

Local Void getlods(f, x, LINK)
double *f;
double *x;
struct LOC_out *LINK;
{
  long i, k, FORLIM;

  firstapprox = true;
  k = 0;
  FORLIM = nall;
  for (i = 0; i < FORLIM; i++) {
    if (itp[i] == 1) {
      k++;
      xall[i] = x[k - 1];
    }
  }
  FORLIM = nuped;
  for (i = 0; i < FORLIM; i++)
    LINK->thisval[i] = 0.0;
  setparam();
  if (!penalty) {
    if (byfamily) {
      for (i = 1; i <= 35; i++) {
	if (finalprn)
	  putc('-', final);
      }
      if (finalprn)
	putc('\n', final);
      if (finalprn)
	fprintf(final, "PEDIGREE |  LN LIKE  | LOG 10 LIKE\n");
      for (i = 1; i <= 35; i++) {
	if (finalprn)
	  putc('-', final);
      }
      if (finalprn)
	putc('\n', final);
    }
    FORLIM = totperson;
    for (i = 1; i <= FORLIM; i++)
      person[i]->done = false;
    FORLIM = totperson;
    for (i = 1; i <= FORLIM; i++)
      person[i]->firstpass = true;
    recombination();
    thisc = minint;
    FORLIM = nuped;
    for (i = 1; i <= FORLIM; i++) {
      likelihood(i, proband[i - 1]);
      LINK->thisval[i - 1] = like;
    }
  }
  firsttime = true;
  FORLIM = nlocus;
  for (i = 0; i < FORLIM; i++)
    maletheta->theta[i] = 0.5;
  FORLIM = nlocus;
  for (i = 0; i < FORLIM; i++)
    femaletheta->theta[i] = 0.5;
  if (penalty) {
    *f = 1.5 * penlike;
    return;
  }
  dolod = true;
  FORLIM = totperson;
  for (i = 1; i <= FORLIM; i++)
    person[i]->done = false;
  FORLIM = totperson;
  for (i = 1; i <= FORLIM; i++)
    person[i]->firstpass = true;
  recombination();
  alike = 0.0;
  thisc = minint;
  FORLIM = nuped;
  for (i = 0; i < FORLIM; i++) {
    likelihood(i + 1, proband[i]);
    alike += like;
    if (byfamily) {
      if (finalprn)
	fprintf(final, "%9ld %12.6f",
		proband[i]->ped, LINK->thisval[i] - like);
      if (finalprn)
	fprintf(final, "%12.6f", (LINK->thisval[i] - like) / log10_);
      if (like == zerolike) {
	if (finalprn)
	  fprintf(final, " inconsistent data\n");
	else if (finalprn)
	  putc('\n', final);
      }
      if (like != zerolike)
	LINK->thisval[i] -= like;
    }
  }
  *f = -2 * alike;
  penlike = *f;
  dolod = false;
}  /* getlods */


Static Void out()
{
  struct LOC_out V;
  long i, j, k, l;
  real deriv[3];
  thetarray oldtheta;
  real lods, dist;
  long FORLIM;
  locusvalues *WITH;
  long FORLIM1, FORLIM2;
  thetavalues *WITH1;
  double TEMP;


  firstapprox = true;
  lasttime = true;
  setparam();
  recombination();
  fun(&f, x);
  if (finalprn)
    fprintf(final, "CHROMOSOME ORDER OF LOCI : \n");
  FORLIM = nlocus;
  for (i = 1; i <= FORLIM; i++) {
    j = 0;
    do {
      j++;
    } while (order[j - 1] != i);
    if (finalprn)
      fprintf(final, "%3ld", j);
  }
  if (finalprn)
    putc('\n', final);
  if (itsys != 0) {
    if (finalprn)
      fprintf(final,
	      "****************** FINAL VALUES **********************\n");
    if (finalprn)
      fprintf(final, "PROVIDED FOR LOCUS %3ld (CHROMOSOME ORDER)\n",
	      order[itsys - 1]);
    if (finalprn)
      fprintf(final,
	      "******************************************************\n");
    WITH = thislocus[itsys - 1];
    if (disequi) {
      if (finalprn)
	fprintf(final, "HAPLOTYPE FREQUENCIES:\n");
      if (finalprn) {
	FORLIM = nuhap;
	for (i = 0; i < FORLIM; i++)
	  fprintf(final, "%9.6f", hapfreq->genarray[i]);
      }
    } else {
      if (finalprn)
	fprintf(final, "GENE FREQUENCIES :\n");
      if (finalprn) {
	FORLIM = WITH->nallele;
	for (i = 0; i < FORLIM; i++)
	  fprintf(final, "%9.6f", WITH->freq[i]);
      }
    }
    if (finalprn)
      putc('\n', final);
    if (WITH->which == quantitative) {
      if (finalprn)
	fprintf(final, "VALUES FOR GENOTYPE MEANS:\n");
      invert(WITH->UU.U1.vmat, WITH->UU.U1.ntrait, &WITH->UU.U1.det);
      if (finalprn) {
	FORLIM = WITH->UU.U1.ntrait;
	for (i = 0; i < FORLIM; i++) {
	  FORLIM1 = WITH->nallele;
	  for (j = 1; j <= FORLIM1; j++) {
	    FORLIM2 = WITH->nallele;
	    for (k = j - 1; k < FORLIM2; k++)
	      fprintf(final, "%6.3f ", WITH->UU.U1.pm[j][k][i]);
	  }
	  if (finalprn)
	    putc('\n', final);
	}
      }
      if (finalprn)
	fprintf(final, "COVARIANCE MATRIX:\n");
      if (finalprn) {
	FORLIM = WITH->UU.U1.ntrait;
	for (i = 1; i <= FORLIM; i++) {
	  for (j = 0; j < i; j++)
	    fprintf(final, "%6.3f ", WITH->UU.U1.vmat[i - 1][j]);
	  putc('\n', final);
	}
      }
    } else if (WITH->which == affection) {
      if (finalprn)
	fprintf(final, "PENETRANCES:\n");
      FORLIM = WITH->UU.U0.nclass;
      for (l = 0; l < FORLIM; l++) {
	if (finalprn) {
	  FORLIM1 = WITH->nallele;
	  for (i = 1; i <= FORLIM1; i++) {
	    FORLIM2 = WITH->nallele;
	    for (j = i - 1; j < FORLIM2; j++)
	      fprintf(final, "%5.3f ", WITH->UU.U0.pen[i][j][2][l]);
	  }
	}
	if (finalprn)
	  putc('\n', final);
	if (finalprn) {
	  if (sexlink) {
	    FORLIM1 = WITH->nallele;
	    for (i = 0; i < FORLIM1; i++)
	      fprintf(final, "%5.3f ", WITH->UU.U0.pen[0][i][2][l]);
	  }
	}
	if (finalprn)
	  putc('\n', final);
      }
    }
  }
  if (finalprn)
    fprintf(final, "******************************************************\n");
  if (interfer && !mapping) {
    if (finalprn)
      fprintf(final, "P VALUES:\n");
    WITH1 = maletheta;
    memcpy(oldtheta, WITH1->theta, sizeof(thetarray));
    WITH1->theta[0] =
      (oldtheta[0] + oldtheta[nlocus - 1] - oldtheta[nlocus - 2]) / 2.0;
    WITH1->theta[nlocus - 2] =
      (oldtheta[nlocus - 2] + oldtheta[nlocus - 1] - oldtheta[0]) / 2.0;
    WITH1->theta[nlocus - 1] =
      (oldtheta[0] + oldtheta[nlocus - 2] - oldtheta[nlocus - 1]) / 2.0;
    if (finalprn) {
      FORLIM = nlocus;
      for (i = 0; i < FORLIM; i++)
	fprintf(final, " %5.3f", WITH1->theta[i]);
    }
    if (finalprn)
      putc('\n', final);
    memcpy(WITH1->theta, oldtheta, sizeof(thetarray));
    if (sexdif) {
      WITH1 = femaletheta;
      memcpy(oldtheta, WITH1->theta, sizeof(thetarray));
      WITH1->theta[0] =
	(oldtheta[0] + oldtheta[nlocus - 1] - oldtheta[nlocus - 2]) / 2.0;
      WITH1->theta[nlocus - 2] =
	(oldtheta[nlocus - 2] + oldtheta[nlocus - 1] - oldtheta[0]) / 2.0;
      WITH1->theta[nlocus - 1] =
	(oldtheta[0] + oldtheta[nlocus - 2] - oldtheta[nlocus - 1]) / 2.0;
      if (finalprn)
	fprintf(final, "FEMALE:  \n");
      if (finalprn) {
	FORLIM = nlocus;
	for (i = 0; i < FORLIM; i++)
	  fprintf(final, " %5.3f", WITH1->theta[i]);
      }
      if (finalprn)
	putc('\n', final);
      memcpy(WITH1->theta, oldtheta, sizeof(thetarray));
    }
    if (ivar == 1 && icall == 0) {
      for (i = 0; i <= 2; i++) {
	TEMP = 1 + exp(xall[nall - i - 1]);
	deriv[i] = exp(xall[nall - i - 1]) / (TEMP * TEMP);
      }
      k = -1;
      for (i = 0; i <= 2; i++) {
	if (itp[nall - i - 1] == 1) {
	  k++;
	  l = -1;
	  for (j = 0; j <= 2; j++) {
	    if (itp[nall - j - 1] == 1) {
	      l++;
	      bmat[n - k - 1][n - l - 1] *= deriv[j] * deriv[i];
	    }
	  }
	  if (l == -1)
	    l = 0;
	  for (j = n - l - 1; j >= 0; j--) {
	    bmat[j][n - k - 1] *= deriv[i];
	    bmat[n - k - 1][j] = bmat[j][n - k - 1];
	  }
	}
      }
    }
  }
  if (finalprn)
    fprintf(final, "THETAS:\n");
  if (finalprn) {
    WITH1 = maletheta;
    FORLIM = nlocus - 2;
    for (i = 0; i <= FORLIM; i++)
      fprintf(final, " %5.3f", WITH1->theta[i]);
  }
  if (interfer) {
    if (finalprn)
      fprintf(final, " %5.3f\n", maletheta->theta[nlocus - 1]);
    else if (finalprn)
      putc('\n', final);
  }
  if (sexdif) {
    if (finalprn)
      fprintf(final, "FEMALE:\n");
    if (finalprn) {
      WITH1 = femaletheta;
      FORLIM = nlocus - 2;
      for (i = 0; i <= FORLIM; i++)
	fprintf(final, " %5.3f", WITH1->theta[i]);
    }
    if (interfer) {
      if (finalprn)
	fprintf(final, " %5.3f\n", femaletheta->theta[nlocus - 1]);
      else if (finalprn)
	putc('\n', final);
    }
    if (readfemale) {
      if (finalprn)
	fprintf(final, "FEMALE/MALE DIST RATIO :\n");
      FORLIM = nlocus - 2;
      for (i = 0; i <= FORLIM; i++) {
	dist = getdist(&maletheta->theta[i]);
	if (dist != 0.0)
	  dist = getdist(&femaletheta->theta[i]) / dist;
	else
	  dist = 0.0;
	if (finalprn)
	  fprintf(final, " %5.3f", dist);
      }
      if (interfer) {
	dist = getdist(&maletheta->theta[nlocus - 1]);
	if (dist != 0.0)
	  dist = getdist(&femaletheta->theta[nlocus - 1]) / dist;
	else
	  dist = 0.0;
	if (finalprn)
	  fprintf(final, " %5.3f", dist);
      }
      if (finalprn)
	putc('\n', final);
    } else {
      if (finalprn)
	fprintf(final, "CONSTANT FEMALE/MALE DIST RATIO :\n");
      if (finalprn)
	fprintf(final, "%5.3f\n", distratio);
    }
  }
  if (finalprn)
    fprintf(final, "******************************************************\n");
  if (finalprn)
    fprintf(final, "-2 LN(LIKE) = % .5E\n", f);
  if (!itmax && nrep == 1) {
    fprintf(isim, "THETAS:\n");
    WITH1 = maletheta;
    FORLIM = nlocus - 2;
    for (i = 0; i <= FORLIM; i++)
      fprintf(isim, " %5.3f", WITH1->theta[i]);
    if (interfer)
      fprintf(isim, " %5.3f\n", maletheta->theta[nlocus - 1]);
    else
      putc('\n', isim);
    if (sexdif) {
      fprintf(isim, "FEMALE:\n");
      WITH1 = femaletheta;
      FORLIM = nlocus - 2;
      for (i = 0; i <= FORLIM; i++)
	fprintf(isim, " %5.3f", WITH1->theta[i]);
      if (interfer)
	fprintf(isim, " %5.3f\n", femaletheta->theta[nlocus - 1]);
      else
	putc('\n', isim);
      if (readfemale) {
	fprintf(isim, "FEMALE/MALE DIST RATIO :\n");
	FORLIM = nlocus - 2;
	for (i = 0; i <= FORLIM; i++) {
	  ddist = getdist(&maletheta->theta[i]);
	  if (ddist != 0.0)
	    ddist = getdist(&femaletheta->theta[i]) / ddist;
	  else
	    ddist = 0.0;
	  fprintf(isim, " %5.3f", ddist);
	}
	if (interfer) {
	  ddist = getdist(&maletheta->theta[nlocus - 1]);
	  if (ddist != 0.0)
	    ddist = getdist(&femaletheta->theta[nlocus - 1]) / ddist;
	  else
	    ddist = 0.0;
	  fprintf(isim, " %5.3f", ddist);
	}
	putc('\n', isim);
      } else {
	fprintf(isim, "CONSTANT FEMALE/MALE DIST RATIO :\n");
	fprintf(isim, "%5.3f\n", distratio);
      }
    }
  }
  getlods(&lods, x, &V);
  lods = (lods - f) / (2 * log10_);
  if (nlocus != 2) {
    if (finalprn)
      fprintf(final, "OTTS GENERALIZED LOD SCORE  =% .5E\n", lods);
    else if (finalprn)
      fprintf(final, "LOD SCORE  =% .5E\n", lods);
  }
  if (lods >= lodlimit[0])
    clods[0]++;
  if (lods >= lodlimit[1])
    clods[1]++;
  if (lods >= lodlimit[2])
    clods[2]++;
  if (nrep == 1) {
    aver = lods;
    variance = 0.0;
    smll = aver;
    big = aver;
  } else {
    ltemp = lods;
    rep = nrep;
    /*  rep := real(nrep);*/
    aver = (rep - 1.0) * aver / rep + ltemp / rep;
    TEMP = ltemp - aver;
    variance = (rep - 1.0) * variance / rep + TEMP * TEMP / (rep - 1.0);
    smll = min(smll, ltemp);
    big = max(big, ltemp);
  }
  if (finalprn)
    fprintf(final, "NUMBER OF ITERATIONS = %5ld\n", nit);
  if (finalprn)
    fprintf(final, "NUMBER OF FUNCTION EVALUATIONS = %5ld\n", nfe);
  if (finalprn)
    fprintf(final, "PTG = % .5E\n", ptg);
  idg++;
  if (prn)
    fprintf(outfile, "EXIT CONDITION%12ld\n", idg);
  if (prn) {
    switch (idg) {

    case 1:
      fprintf(outfile, "Maximum possible accuracy reached\n");
      break;

    case 2:
      fprintf(outfile, "Search direction no longer downhill\n");
      break;

    case 3:
      fprintf(outfile,
	      "Accumulation of rounding error prevents further progress\n");
      break;

    case 4:
      fprintf(outfile,
	"All significant differences lost through cancellation in conditioning\n");
      break;

    case 5:
      fprintf(outfile, "Specified tolerance on normalized gradient met\n");
      break;

    case 6:
      fprintf(outfile, "Specified tolerance on gradient met\n");
      break;

    case 7:
      fprintf(outfile, "Maximum number of iterations reached\n");
      break;

    case 8:
      fprintf(outfile, "Excessive cancellation in gradient\n");
      break;
    }
  }
  if (ivar == 1 && icall == 0) {
    if (finalprn)
      fprintf(final, "VARIANCE-COVARIANCE OF THE ESTIMATES\n");
    if (interfer && !mapping) {
      if (finalprn)
	fprintf(final, "VALUES GIVEN FOR P VALUES\n");
    }
    FORLIM = n;
    for (i = 0; i < FORLIM; i++) {
      FORLIM1 = n;
      for (j = 0; j < FORLIM1; j++) {
	if (finalprn)
	  fprintf(final, "% .5E ", 2.0 * bmat[i][j]);
      }
      if (finalprn)
	putc('\n', final);
    }
  }
  if (finalprn)
    fprintf(final, "******************************************************\n");
  if (finalprn)
    fprintf(final, "******************************************************\n");
}  /* out */


/* Local variables for step: */
struct LOC_step {
  enum {
    go_, exit1, exit2, exit3
  } cont;
} ;

Local Void firststep(LINK)
struct LOC_step *LINK;
{
  long i, FORLIM;

  sumt = 0.0;
  idg = 0;
  FORLIM = n;
  for (i = 0; i < FORLIM; i++)
    xt[i] = x[i] + t * p[i];
  fun(&ft, xt);
  nfe++;
}  /* firststep */


Local Void decreaset(LINK)
struct LOC_step *LINK;
{
  long i, FORLIM;

  /*Modified by M. Lathrop 29/04/86 to trap tmin=0 problem*/
  if (tmin < small)
    tmin = small;
  LINK->cont = go_;
  while (LINK->cont == go_) {
    if (f - fsave == 0.0 && idif == 2)
      idg = 2;
    if (t < tmin) {
      LINK->cont = exit3;
      break;
    }
    t = 0.5 * t;
    FORLIM = n;
    for (i = 0; i < FORLIM; i++)
      xt[i] = x[i] + t * p[i];
    fun(&ft2, xt);
    nfe++;
    if (ft2 < f) {
      sumt += t;
      f = ft2;
      idg = 0;
      memcpy(x, xt, sizeof(vector));
      if (sumt < tmin)
	LINK->cont = exit3;
      else
	LINK->cont = exit2;
      continue;
    }
    if (ft + f - ft2 - ft2 <= 0.0)
      scal = 0.1;
    else {
      scal = 1.0 + 0.5 * (f - ft) / (f + ft - ft2 - ft2);
      if (scal < 0.1)
	scal = 0.1;
    }
    t = scal * t;
    FORLIM = n;
    for (i = 0; i < FORLIM; i++)
      xt[i] = x[i] + t * p[i];
    fun(&ft, xt);
    nfe++;
    if (f <= ft)
      continue;
    sumt += t;
    idg = 0;
    f = ft;
    memcpy(x, xt, sizeof(vector));
    if (t < tmin)
      LINK->cont = exit3;
    else
      LINK->cont = exit2;
  }
}  /* decreaset */


Local Void increaset(LINK)
struct LOC_step *LINK;
{
  long i, FORLIM;

  sumt += t;
  LINK->cont = go_;
  while (LINK->cont == go_) {
    twot = t + t;
    if (ibnd > 0 && tbnd >= 0.0 && twot > tbnd) {
      f = ft;
      memcpy(x, xt, sizeof(vector));
      if (prn)
	fprintf(outfile, "****** ACTIVE BOUNDARY CONSTRAINT *****\n");
      LINK->cont = exit1;
      continue;
    }
    memcpy(x, xt, sizeof(vector));
    FORLIM = n;
    for (i = 0; i < FORLIM; i++)
      xt[i] = x[i] + t * p[i];
    fun(&f2t, xt);
    nfe++;
    if (f2t > ft) {
      f = ft;
      LINK->cont = exit2;
      break;
    }
    if (f2t + f - ft - ft < 0.0 || ft - f2t + curv * (ft - f) >= 0.0) {
      sumt += t;
      t += t;
      ft = f2t;
    } else {
      sumt += t;
      f = f2t;
      memcpy(x, xt, sizeof(vector));
      LINK->cont = exit2;
    }
  }
}  /* increaset */


Static Void step()
{
  struct LOC_step V;


  firstapprox = true;
  firststep(&V);
  if (f > ft)
    increaset(&V);
  else
    decreaset(&V);
  if (V.cont == exit2)
    t = sumt;
  if (V.cont == exit3) {
    iret = 1;
    return;
  }
  memcpy(gsave, g, sizeof(vector));
  isw = 0;
  iret = 2;
}  /* step */


/* Local variables for update: */
struct LOC_update {
  long j, jm1;
  boolean cont, dfp;
} ;

Local Void prep(LINK)
struct LOC_update *LINK;
{
  long i, k, FORLIM;

  iswup = 1;
  wtil[0] = -gsave[0];
  ztil[0] = y[0];
  FORLIM = n;
  for (k = 2; k <= FORLIM; k++) {
    wtil[k - 1] = -gsave[k - 1];
    ztil[k - 1] = y[k - 1];
    for (i = 0; i <= k - 2; i++) {
      tl[i][k - 1] = tl[k - 1][i];
      wtil[k - 1] -= tl[k - 1][i] * wtil[i];
      ztil[k - 1] -= tl[k - 1][i] * ztil[i];
    }
  }
  sb = 0.0;
  sc = 0.0;
  sd = 0.0;
  FORLIM = n;
  for (i = 0; i < FORLIM; i++) {
    sb += t * ztil[i] * wtil[i] / d[i];
    sc += t * t * wtil[i] * wtil[i] / d[i];
    sd += ztil[i] * ztil[i] / d[i];
  }
}  /* prep */


Local Void sr1update(LINK)
struct LOC_update *LINK;
{
  alpha = -1.0;
  if (sc - sb == 0.0 || sb - sd == 0.0 || sd - 2.0 * sb + sc == 0.0) {
    idg = 3;
    return;
  }
  sa = 1.0 / (sqrt(fabs(sc - sb)) * sqrt(fabs(sb - sd)));
  thet1 = ((sd - sb) * sa + 1.0) / (2.0 * sb - sd - sc);
  thet2 = sa + (sa * (sb - sc) + 1.0) / (sd - 2.0 * sb + sc);
}  /* sr1update */


Local Void sr2update(LINK)
struct LOC_update *LINK;
{
  aa = sb / sc - 2.0 * (sd / sb) + sd / sc;
  bb = sb / sc - 1.0;
  cc = 1.0 - sb / sd;
  del2 = bb * bb - aa * cc;
  LINK->dfp = true;
  if (del2 > 0.00000001) {
    LINK->dfp = false;
    del = sqrt(del2);
    alph1 = (del - bb) / aa;
    alph2 = (-bb - del) / aa;
    if (fabs(alph1) < fabs(alph2))
      alpha = alph1;
    else
      alpha = alph2;
    if (fabs(alpha) < 0.00001)
      LINK->dfp = true;
    else {
      sa = (alpha + 1.0) * (alpha + 1.0) + sc / sb -
	   alpha * alpha * (sc / sb) * (sd / sb) - 1.0 + alpha * alpha * sd / sb;
      if (sa <= 0.0)
	sa = 0.0;
      else {
	sa = sqrt(sa);
	sa = 1.0 / (sa * sb);
      }
      rdiv = 1.0 / (alpha * alpha * sd + 2.0 * alpha * sb + sc);
      thet1 = -(sa * alpha * (alpha * sd + sb) + 1.0) * rdiv;
      thet2 = sa + (alpha * sa * (sc + alpha * sb) - alpha) * rdiv;
    }
  }
  if (!LINK->dfp)
    return;
  alpha = 0.0;
  sa = 1.0 / (sqrt(sb) * sqrt(sc));
  thet1 = -1.0 / sc;
  thet2 = sa;
}  /* sr2update */


Local Void getwzs(LINK)
struct LOC_update *LINK;
{
  long i, k, FORLIM;

  FORLIM = n;
  for (i = 0; i < FORLIM; i++) {
    w[i] = t * wtil[i] + alpha * ztil[i];
    z[i] = t * thet1 * wtil[i] + thet2 * ztil[i];
  }
  s[n - 1] = 0.0;
  FORLIM = n;
  for (k = 1; k < FORLIM; k++) {
    LINK->j = n - k + 1;
    LINK->jm1 = LINK->j - 1;
    s[LINK->jm1 - 1] = s[LINK->j - 1] +
		       w[LINK->j - 1] * w[LINK->j - 1] / d[LINK->j - 1];
  }
  nu = 1.0;
  eta = 0.0;
}  /* getwzs */


Local Void recur(LINK)
struct LOC_update *LINK;
{
  long i, k, FORLIM;

  LINK->cont = true;
  while (LINK->cont) {
    LINK->cont = false;
    if (iswup < 2) {
      FORLIM = n;
      for (i = 0; i < FORLIM; i++)
	wtjp1[i] = -gsave[i];
      memcpy(ztjp1, y, sizeof(vector));
    } else {
      FORLIM = n;
      for (i = 0; i < FORLIM; i++) {
	wtil[i] = alpha * y[i] - t * g[i];
	ztil[i] = thet2 * y[i] - t * thet1 * g[i];
      }
    }
    LINK->j = 0;
    lambj2 = 0.0;
    while (LINK->j < n - 1) {
      LINK->j++;
      if (iswup < 2) {
	FORLIM = n;
	for (k = LINK->j; k < FORLIM; k++) {
	  wtjp1[k] -= wtil[LINK->j - 1] * tl[k][LINK->j - 1];
	  ztjp1[k] -= ztil[LINK->j - 1] * tl[k][LINK->j - 1];
	}
      } else {
	FORLIM = n;
	for (k = LINK->j; k < FORLIM; k++) {
	  wtjp1[k] = wtil[k] - w[LINK->j - 1] * tl[k][LINK->j - 1];
	  ztjp1[k] = ztil[k] - z[LINK->j - 1] * tl[k][LINK->j - 1];
	}
      }
      aj = nu * z[LINK->j - 1] - eta * w[LINK->j - 1];
      thj = 1.0 + aj * w[LINK->j - 1] / d[LINK->j - 1];
      lambj2 = thj * thj + aj * aj * s[LINK->j - 1] / d[LINK->j - 1];
      if (iswup < 2) {
	if (lambj2 > 10.0) {
	  LINK->cont = true;
	  iswup = 2;
	  FORLIM = n;
	  for (k = 2; k <= FORLIM; k++) {
	    for (i = 0; i <= k - 2; i++)
	      tl[k - 1][i] = tl[i][k - 1];
	  }
	  LINK->j = n;
	}
      }
      if (LINK->cont)
	continue;
      dp[LINK->j - 1] = d[LINK->j - 1] * lambj2;
      lambj = sqrt(lambj2);
      if (thj > 0.0)
	lambj = -lambj;
      muj = thj - lambj;
      bj = thj * w[LINK->j - 1] + aj * s[LINK->j - 1];
      gamlj = bj * nu / (lambj2 * d[LINK->j - 1]);
      betlj = (aj - bj * eta) / (lambj2 * d[LINK->j - 1]);
      nu = -(nu / lambj);
      eta = -((eta + aj * aj / (muj * d[LINK->j - 1])) / lambj);
      if (iswup < 2) {
	FORLIM = n;
	for (k = LINK->j; k < FORLIM; k++)
	  tl[k][LINK->j - 1] += t * (betlj + thet1 * gamlj) * wtjp1[k] +
				(alpha * betlj + thet2 * gamlj) * ztjp1[k];
      } else {
	FORLIM = n;
	for (k = LINK->j; k < FORLIM; k++) {
	  tl[k][LINK->j - 1] = tl[k][LINK->j - 1] / lambj2 + betlj * wtil[k] +
			       gamlj * ztil[k];
	  wtil[k] = wtjp1[k];
	  ztil[k] = ztjp1[k];
	}
      }
    }
  }
  aj = nu * z[n - 1] - eta * w[n - 1];
  lambj = 1.0 + aj * w[n - 1] / d[n - 1];
  dp[n - 1] = d[n - 1] * lambj * lambj;
  memcpy(d, dp, sizeof(vector));
}  /* recur */


/*UPDATE*/

Static Void update()
{
  struct LOC_update V;


  prep(&V);
  if (sb < sc)
    fbcd = sb;
  else
    fbcd = sc;
  if (sd < fbcd)
    fbcd = sd;
  if (fbcd > small) {
    fbcd = 2.0 * sc * (sd / sb) / (sc + sd);
    if (fbcd < 1.0)
      sr1update(&V);
    else
      sr2update(&V);
  } else
    sr1update(&V);
  if (idg != 3) {
    getwzs(&V);
    recur(&V);
  }
}  /* update */


Local Void bldlt(b)
vector *b;
{
  long i, j, ic, k;
  real su, tlic, ff, hh, temp, temp1;
  matrix s;
  vector tu, xvec;
  long FORLIM, FORLIM1, FORLIM2;

  firstapprox = true;
  ff = f;
  if (icall > 0 && ihess <= 0) {
    memcpy(xvec, xall, sizeof(vector));
    hh = 10 * gh;
    FORLIM = n;
    for (i = 0; i < FORLIM; i++) {
      temp = x[i];
      x[i] += hh;
      fun(&f, x);
      tu[i] = f;
      for (j = 0; j < i + 1; j++) {
	temp1 = x[j];
	x[j] += hh;
	fun(&f, x);
	b[i][j] = f;
	b[j][i] = b[i][j];
	x[j] = temp1;
      }
      x[i] = temp;
    }
    memcpy(xall, xvec, sizeof(vector));
    FORLIM = n;
    for (i = 0; i < FORLIM; i++) {
      for (j = 0; j < i + 1; j++) {
	b[i][j] = (ff + b[i][j] - tu[i] - tu[j]) / (hh * hh);
	b[j][i] = b[i][j];
      }
    }
  }
  ic = 1;
  while (b[ic - 1][ic - 1] > 0.0 && ic <= n) {
    temp = b[ic - 1][ic - 1];
    d[ic - 1] = b[ic - 1][ic - 1];
    tl[ic - 1][ic - 1] = 1.0;
    if (ic != n) {
      FORLIM = n;
      for (k = ic; k < FORLIM; k++)
	tl[k][ic - 1] = b[k][ic - 1] / temp;
      FORLIM = n;
      for (i = ic; i < FORLIM; i++) {
	tlic = tl[i][ic - 1];
	FORLIM1 = n;
	for (k = i; k < FORLIM1; k++)
	  b[k][i] -= tlic * tl[k][ic - 1] * temp;
      }
    }
    ic++;
  }
  if (ic > n) {
    icall--;
    if (prn)
      fprintf(outfile, "FACTORIZATION SUCCEEDED\n");
  } else {
    icall++;
    if (prn)
      fprintf(outfile, "FACTORIZATION FAILED\n");
  }
  if (icall != 0)
    return;
  s[0][0] = 1.0;
  FORLIM = n;
  for (i = 2; i <= FORLIM; i++) {
    for (k = 0; k <= i - 2; k++) {
      su = 0.0;
      for (j = k; j <= i - 2; j++)
	su += tl[i - 1][j] * s[j][k];
      s[i - 1][k] = -su;
    }
    s[i - 1][i - 1] = 1.0;
  }
  FORLIM = n;
  for (i = 0; i < FORLIM; i++) {
    FORLIM1 = n;
    for (j = i; j < FORLIM1; j++) {
      su = 0.0;
      FORLIM2 = n;
      for (k = j; k < FORLIM2; k++)
	su += s[k][i] * s[k][j] / d[k];
      b[i][j] = su;
      b[j][i] = su;
    }
  }
  if (prn)
    fprintf(outfile, "B-MATRIX\n");
  if (prn) {
    FORLIM = n;
    for (i = 1; i <= FORLIM; i++) {
      for (j = 0; j < i; j++)
	fprintf(outfile, "% .5E", b[i - 1][j]);
      putc('\n', outfile);
    }
  }
  if (ivar != 1)
    return;
  k = 0;
  FORLIM = nall;
  for (i = 0; i < FORLIM; i++) {
    se[i] = 0.0;
    if (itp[i] == 1) {
      k++;
      se[i] = sqrt(b[k - 1][k - 1]);
    }
  }
}  /* bldlt */


Local Void inib()
{
  long i, j;
  FILE *in1;
  long FORLIM;

  in1 = NULL;
  if (icall != 1 && ihess > 1) {
    printf("Reading IN1.DAT (covariance matrix)\n");
    if (in1 != NULL)
      in1 = freopen("in1.dat", "r", in1);
    else
      in1 = fopen("in1.dat", "r");
    if (in1 == NULL)
      exit(FileNotFound);
    FORLIM = n;
    /*SUN*/
    for (i = 1; i <= FORLIM; i++) {
      for (j = 0; j < i; j++)
	fscanf(in1, "%lg", &bmat[i - 1][j]);
    }
    FORLIM = n;
    for (i = 0; i < FORLIM; i++) {
      for (j = 0; j < i + 1; j++)
	bmat[j][i] = bmat[i][j];
    }
  }
  /*    assign(in1,'IN1.DAT'); */
  /*SUN*/
  bldlt(bmat);
  if (icall == 1) {
    FORLIM = n;
    for (i = 0; i < FORLIM; i++) {
      for (j = 0; j < i + 1; j++)
	tl[i][j] = 0.0;
      tl[i][i] = 1.0;
      d[i] = 1.0;
    }
  }
  if (in1 != NULL)
    fclose(in1);
}  /* inib */


Local Void initialize()
{
  long i, j, FORLIM;

  firstapprox = true;
  if (prn)
    fprintf(outfile, "DIFFER INTER = % .5E TRUNC UPPER = % .5E\n", gh, trupb);
  nit = 0;
  idg = 0;
  idif = 1;
  t = 0.1;
  tmin = 0.0;
  fsmf = 0.0;
  fun(&f, x);
  nfe = 1;
  FORLIM = n;
  for (i = 0; i < FORLIM; i++) {
    for (j = 0; j < i + 1; j++)
      tl[i][j] = 0.0;
    tl[i][i] = 1.0;
    d[i] = 1.0;
  }
  if (ihess <= 0)
    return;
  icall = 0;
  inib();
  icall = 1;
  t = 1.0;
}  /* initialize */

/*initialize*/

Local Void gforward()
{
  long i, FORLIM;

  firstapprox = true;
  FORLIM = n;
  for (i = 0; i < FORLIM; i++) {
    hx = gh;
    if (ihx == 1) {
      if (fabs(x[i]) < 0.1)
	hx = gh * 0.1;
      else
	hx = gh * fabs(x[i]);
    }
    xsave = x[i];
    x[i] += hx;
    fun(&fxph, x);
    firstapprox = false;
    g[i] = (fxph - f) / hx;
    x[i] = xsave;
  }
  nfe += n;
}  /* gforward */


Local Void gcentral()
{
  long i, FORLIM;

  firstapprox = true;
  FORLIM = n;
  for (i = 0; i < FORLIM; i++) {
    hx = gh;
    if (ihx == 1) {
      if (fabs(x[i]) < 0.1)
	hx = gh * 0.1;
      else
	hx = gh * fabs(x[i]);
    }
    xsave = x[i];
    x[i] += hx;
    fun(&fxph, x);
    firstapprox = false;
    x[i] = xsave - hx;
    fun(&fxmh, x);
    g[i] = (fxph - fxmh) / (hx + hx);
    x[i] = xsave;
  }
  nfe += n + n;
}  /* gcentral */


Local Void getp()
{
  long i, nmj, j, FORLIM, FORLIM1;

  nit++;
  fsav2 = fsave;
  fsave = f;
  if (prn)
    putc('\n', outfile);
  if (prn)
    fprintf(outfile, "ITERATION %5ld T = %10.3f NFE = %5ld F = % .5E\n",
	    nit, t, nfe, f);
  if (prn)
    printf("ITERATION %5ld T = %10.3f NFE = %5ld F = % .5E\n", nit, t, nfe, f);
  if (nit > maxit) {
    idg = 6;
    return;
  }
  if (n < 20) {
    if (prn)
      fprintf(outfile, "X= ");
    outcontrol(x);
    if (prn)
      fprintf(outfile, "G= ");
    outcontrol(g);
  }
  if (n == 1)
    p[0] = -(g[0] / d[0]);
  else {
    gv[0] = -g[0];
    FORLIM = n;
    for (i = 2; i <= FORLIM; i++) {
      gv[i - 1] = -g[i - 1];
      for (j = 0; j <= i - 2; j++)
	gv[i - 1] -= tl[i - 1][j] * gv[j];
    }
    p[n - 1] = gv[n - 1] / d[n - 1];
    FORLIM = n;
    for (j = 1; j < FORLIM; j++) {
      nmj = n - j;
      p[nmj - 1] = gv[nmj - 1] / d[nmj - 1];
      FORLIM1 = n;
      for (i = nmj; i < FORLIM1; i++)
	p[nmj - 1] -= tl[i][nmj - 1] * p[i];
    }
  }
  if (prn)
    fprintf(outfile, "P= ");
  outcontrol(p);
  if (ibnd == 1) {
    FORLIM = n;
    for (i = 0; i < FORLIM; i++) {
      if (p[i] == 0.0)
	idg = 7;
    }
  }
  if (idg == 7)
    return;
  ptg = 0.0;
  FORLIM = n;
  for (i = 0; i < FORLIM; i++)
    ptg += p[i] * g[i];
  if (ptg >= 0.0) {
    idg = 1;
    return;
  }
  if (fabs(ptg) < tol) {
    idg = 4;
    fsmf = fsav2 - f;
    if (prn)
      fprintf(outfile, "FSMF = % .5E PTG = % .5E TMIN = % .5E\n",
	      fsmf, ptg, tmin);
    return;
  }
  xpm = xpmcon;
  if (nit == 1)
    return;
  FORLIM = n;
  for (i = 0; i < FORLIM; i++) {
    if (ihx == 1) {
      xp = fabs(x[i]);
      if (xp < 0.1)
	xp = 0.1;
      xp /= fabs(p[i]);
      if (xp < xpm)
	xpm = xp;
    } else if (xpm > fabs(1.0 / p[i]))
      xpm = fabs(1.0 / p[i]);
  }
  tmin = 0.5 * xpm * gh / trupb;
  if (idif == 2)
    t = 1.0;
  else {
    t = -2.0 * fsmf / ptg;
    if (t <= 0.0)
      t = 1.0;
    else if (1.0 < t)
      t = 1.0;
  }
  if (prn)
    fprintf(outfile, "FSMF = % .5E PTG =% .5E TMIN=% .5E\n", fsmf, ptg, tmin);
  if (prn)
    fprintf(outfile, "INITIAL T = % .5E\n", t);
}  /* getp */


Local Void getytp()
{
  long i, FORLIM;

  ytp = 0.0;
  fsmf = fsave - f;
  FORLIM = n;
  for (i = 0; i < FORLIM; i++) {
    y[i] = g[i] - gsave[i];
    ytp += y[i] * p[i];
  }
}  /* getytp */


Local Void chkbnd()
{
  /*This procedure modified by M. Lathrop 29/04/86
  with the introduction of continue and repeat*/
  long i, j, k, ik, ii, jk;
  boolean continue_;
  long FORLIM, FORLIM1;

  do {
    clb = clbcon;
    for (j = 1; j <= 2; j++) {
      k = 0;
      FORLIM1 = nall;
      for (i = 1; i <= FORLIM1; i++) {
	if (itp[i - 1] == 1) {
	  k++;
	  check = fabs(x[k - 1] - bnd[i - 1][j - 1]);
	  if (check <= clb) {
	    clb = check;
	    ik = k;
	    ii = i;
	    jk = j;
	  }
	}
      }
    }
    if (clb < 0.1 * gh) {
      ihess = 0;
      if (jk == 2)
	eh = -gh - gh;
      else
	eh = gh + gh;
      x[ik - 1] = bnd[ii - 1][jk - 1] + eh;
      xall[ii - 1] = x[ik - 1];
      itp[ii - 1] = 0;
      k = 0;
      FORLIM1 = nall;
      for (i = 0; i < FORLIM1; i++) {
	if (itp[i] == 1) {
	  k++;
	  x[k - 1] = xall[i];
	}
      }
      n = k;
      tol = tolconst;
      /*            tol:=tolconst*sqrt(n);*/
      continue_ = true;
      if (prn)
	fprintf(outfile,
		"******* A VARIABLE WAS SET TO A BOUND ***********\n");
    } else {
      tbnd = tbndcon;
      for (j = 0; j <= 1; j++) {
	k = 0;
	FORLIM = nall;
	for (i = 0; i < FORLIM; i++) {
	  if (itp[i] == 1) {
	    k++;
	    teq = (bnd[i][j] - x[k - 1]) / p[k - 1];
	    if (teq < tbnd && teq > 0.0)
	      tbnd = teq;
	  }
	}
      }
      continue_ = false;
      if (t * (2.0 + gh) >= tbnd)
	t = tbnd * (0.5 - gh);
      if (prn)
	fprintf(outfile, "TBND = % .5E RESET T = % .5E\n", tbnd, t);
    }
  } while (continue_);   /* chkbnd */
}


Local Void iterate()
{
  while (continue_ == go) {
    active = false;
    getp();
    if (idg != 0)
      continue_ = quit;
    if (idg == 0) {
      iret = 2;
      if (ibnd == 1)
	chkbnd();
      if (n == 0)
	continue_ = quit;
      if (!active)
	step();
    }
    /* THE NEXT IS NOT TRUE IF ACTIVE */
    if (iret != 2) {
      if (idif != 1) {
	continue_ = quit;
	break;
      }
      idif = 2;
      isw = 1;
      gcentral();
      if (prn)
	fprintf(outfile, "***** SWITCHING TO CENTRAL DIFFERERENCE*****\n");
      continue;
    }
    if (active) {
      continue_ = restart;
      break;
      continue;
    }
    if (idif == 1)
      gforward();
    else
      gcentral();
    if (isw != 0)
      continue;
    getytp();
    if (ytp <= 0.0)
      continue;
    if (n == 1)
      d[0] = -(y[0] * d[0] / (t * gsave[0]));
    else
      update();
    if (idg == 3)
      continue_ = quit;
  }
}  /* iterate */


/*START GEMINI*/

Static Void gemini()
{
  continue_ = go;
  while (continue_ != quit) {
    initialize();
    isw = 0;
    gforward();
    continue_ = go;
    iterate();
  }
  if (ivar == 1)
    inib();
}


/* gemini */
/*gemini*/


Static Void initisim()
{
  long i;

  printf("Program ISIM version %s\n\n", version);
  printf("The program constants are set to the following maxima:\n");
  printf("%6ld loci in mapping problem\n", (long)maxlocus);
  printf("%6ld alleles at a single locus\n", (long)maxall);
  printf("%6ld recombination probabilities (maxneed)\n", (long)maxneed);
  printf("%6ld maximum of censoring array (maxcensor)\n", (long)maxcensor);
  printf("%6ld haplotypes = n1 x n2 x ... where ni = current # alleles at locus i\n",
	 (long)maxhap);
  printf("%6ld joint genotypes for a female\n", (long)maxfem);
  printf("%6ld joint genotypes for a male\n", (long)maxmal);
  printf("%6ld individuals in all pedigrees combined\n", (long)maxind);
  printf("%6ld pedigrees (maxped\n", (long)maxped);
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
  /* assign(datafile,'datafile.dat');
     assign(speedfile,'speedfil.dat');
     assign(ipedfile,'ipedfile.dat');
     assign(outfile,'outfile.dat');
     assign(final,'final.dat');*/
  /*SUN*/
  for (i = 1; i <= 48; i++)
    putchar('*');
  printf("\n This is ISIM, which is a modified version of ILINK\n");
  printf(" designed to be used on datafiles generated by the \n");
  printf(" companion program SLINK.\n\n");
  printf(" Pedigree data are read from the file ipedfile.dat\n");
  printf(" The statistical summary is written to isim.dat \n");
  for (i = 1; i <= 48; i++)
    putchar('*');
  printf("\nReading SIMOUT.DAT\n");

  /*Read the number of pedigrees per replicate from the most recent simout.dat*/
  /* assign(simout,'simout.dat');*/
  /*SUN*/
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
    printf(" There are %2ld pedigrees per replicate.\n", opeds);
  else
    printf(" There is %2ld pedigree per replicate.\n", opeds);
}  /* initisim */


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

  /* assign(limit,'LIMIT.DAT');*/
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
  long FORLIM;
  thisperson *WITH;
  long FORLIM1;

  /*PASCAL_MAIN(argc, argv);*/
  final = NULL;
  simout = NULL;
  ipedfile = NULL;
  speedfile = NULL;
  datafile = NULL;
  outfile = NULL;
  isim = NULL;
  initisim();
  inputdata();
  getparam();
  /*speedfile modification*/
  lastseg = 0;
  segperson = 0;
  lastspeed = 0;
  if (datafile != NULL)
    fclose(datafile);
  datafile = NULL;   /*SUN*/
  gh = sqrt(exp(-nbit * log(2.0)));
  tol = tolconst;
  /* tol:=tolconst*sqrt(n);*/
  trupb = sqrt(gh);
  censorstruct = (censorrec *)Malloc(sizeof(censorrec));
  /*  censorstruct:=censorpnt(NewPtr(SizeOf(censorrec)));*/
  gennustruct = (gennurec *)Malloc(sizeof(gennurec));
  /*  gennustruct:=gennuptr(NewPtr(SizeOf(gennurec)));*/
  if (approximate)
    approxstruct = (approxrec *)Malloc(sizeof(approxrec));
  /*   approxstruct:=approxpnt(NewPtr(SizeOf(approxrec)));*/
  getlocations();
  memcpy(holditp, itp, maxn * sizeof(long));
  itmax = false;   /*Flag indicating iterated parameters*/
  FORLIM = nall;
  for (i = 1; i <= FORLIM; i++) {
    if (itp[i - 1] == 1)
      itmax = true;
  }
  memcpy(holdxall, xall, sizeof(vector));
  for (i = 1; i <= 3; i++)
    clods[i - 1] = 0;
  nrep = 1;
  /*assign(isim,'isim.dat');*/
  /*SUN*/
  if (isim != NULL)
    isim = freopen("isim.dat", "w", isim);
  else
    isim = fopen("isim.dat", "w");
  if (isim == NULL)
    exit(FileNotFound);
  /*SUN*/
  /* Copy most recent simout.dat to isim.dat */
  fprintf(isim, " ********* Data from most recent SIMOUT.DAT *********\n");
  /*assign(simout,'simout.dat');*/
  /*SUN*/
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
    putc(chtemp, isim);
    if (P_eoln(simout)) {
      fscanf(simout, "%*[^\n]");
      getc(simout);
      putc('\n', isim);
    }
  }
  fprintf(isim, " ********* End of most recent SIMOUT.DAT *********\n\n");
  if (simout != NULL)
    fclose(simout);
  simout = NULL;   /*SUN*/
  fprintf(isim, " ORDER OF LOCI : \n");
  FORLIM = nlocus;
  for (i = 1; i <= FORLIM; i++) {
    j = 0;
    do {
      j++;
    } while (order[j - 1] != i);
    fprintf(isim, "%3ld", j);
  }
  putc('\n', isim);
  if (nlocus != 2) {
    if (itmax)
      fprintf(isim, " Average Maximum of Ott's Generalized Lod Score\n");
    else
      fprintf(isim, " Average of Ott's Generalized Lod Score\n");
  } else if (itmax)
    fprintf(isim, " Average Maximum Lod Score\n");
  else
    fprintf(isim, " Average Lod Score\n");
  putc('\n', isim);
  while (!P_eof(ipedfile)) {   /*IF person[i]<>NIL*/
    memcpy(itp, holditp, maxn * sizeof(long));
    memcpy(xall, holdxall, sizeof(vector));
    n = 0;
    FORLIM = nall;
    for (i = 1; i <= FORLIM; i++) {
      if (itp[i - 1] == 1) {
	n++;
	x[n - 1] = xall[i - 1];
      }
    }
    ihess = 0;
    ibnd = 1;
    icall = 1;
    ivar = 0;
    ihx = 1;
    f = 0.0;
    ptg = 0.0;
    inconsistent = false;
    firsttime = true;
    lasttime = false;
    dolod = false;
    firstapprox = true;
    maxit = n * 50;
    readpedseg();
    FORLIM = totperson;
    for (i = 1; i <= FORLIM; i++)
      person[i]->store = NULL;
    /* Set all store pointers to NIL in order to be able to use Dispose later */
    readspseg();
    gemini();
    out();
    nrep++;
    if (nrep % 10 == 0)
      putchar('*');
    else
      putchar('.');
/* p2c: isim.p, line 6067:
 * Note: Using % for possibly-negative arguments [317] */
    fflush(stdout);
    /*P_ioresult = 0;   *SUN*/
    if (nrep % 50 == 0)
      printf(" %ld replicates analyzed\n", nrep);
    FORLIM = totperson;
/* p2c: isim.p, line 6072:
 * Note: Using % for possibly-negative arguments [317] */
    /* Need to reuse memory once the pedigrees have been processed */
    for (i = 1; i <= FORLIM; i++) {
      if (person[i] != NULL) {
	WITH = person[i];
	FORLIM1 = nlocus;
	for (j = 1; j <= FORLIM1; j++) {
	  if (WITH->phen[j - 1] != NULL)
	    Free(WITH->phen[j - 1]);
	}
	/*           DisposPtr(Ptr(phen[j]));*/
	if (person[i]->store != NULL)
	  Free(person[i]->store);
	/*         DisposPtr(Ptr(person[i]^.store));*/
	Free(person[i]);
	/*        DisposPtr(Ptr(person[i]));*/
	person[i] = NULL;
      }
    }
  }
  if (ipedfile != NULL)
    fclose(ipedfile);
  ipedfile = NULL;   /*SUN*/
  if (speedfile != NULL)
    fclose(speedfile);
  speedfile = NULL;   /*SUN*/
  if (outfile != NULL)
    fclose(outfile);
  outfile = NULL;   /*SUN*/
  if (final != NULL)
    fclose(final);
  final = NULL;   /*SUN*/
  fprintf(isim, "The number of replicates is %6ld\n", nrep - 1);
  for (i = 1; i <= 66; i++)
    putc('-', isim);
  putc('\n', isim);
  if (nlocus != 2) {
    if (itmax)
      printf(" Average Maximum of Ott's Generalized Lod Score\n");
    else
      printf(" Average of Ott's Generalized Lod Score\n");
  } else if (itmax)
    printf(" Average Maximum Lod Score\n");
  else
    printf(" Average Lod Score\n");
  printf("\nThe number of replicates is %6ld\n", nrep - 1);
  for (i = 1; i <= 66; i++)
    putchar('-');
  putchar('\n');
  if (itmax)
    printf(" Average Maximum | StdDev       | Min          | Max\n");
  else
    printf("      Average    | StdDev       | Min          | Max\n");
  if (itmax)
    fprintf(isim, " Average Maximum | StdDev       | Min          | Max\n");
  else
    fprintf(isim, "      Average    | StdDev       | Min          | Max\n");
  printf("%15.6f%15.6f%15.6f%15.6f\n", aver, sqrt(variance), smll, big);
  fprintf(isim, "%15.6f%15.6f%15.6f%15.6f\n", aver, sqrt(variance), smll, big);
  for (i = 1; i <= 66; i++)
    putc('-', isim);
  putc('\n', isim);
  for (i = 1; i <= 66; i++)
    putchar('-');
  putchar('\n');
  for (i = 1; i <= 66; i++)
    putc('-', isim);
  putc('\n', isim);
  for (i = 1; i <= 66; i++)
    putchar('-');
  putchar('\n');
  if (itmax)
    printf("   Number of maximum lod scores greater than \n");
  else
    printf("   Number of lod scores greater than \n");
  printf("   a given constant \n");
  for (i = 1; i <= 66; i++)
    putchar('-');
  printf("\n Constant  Number  Percent \n");
  for (i = 1; i <= 3; i++)
    printf("%8.3f%8ld%9.3f\n",
	   lodlimit[i - 1], clods[i - 1], 100.0 * clods[i - 1] / (nrep - 1));
  if (itmax)
    fprintf(isim, "   Number of maximum lod scores greater than \n");
  else
    fprintf(isim, "   Number of lod scores greater than \n");
  fprintf(isim, "   a given constant \n");
  for (i = 1; i <= 66; i++)
    putc('-', isim);
  fprintf(isim, "\n Constant  Number  Percent \n");
  for (i = 1; i <= 3; i++)
    fprintf(isim, "%8.3f%8ld%9.3f\n",
	    lodlimit[i - 1], clods[i - 1], 100.0 * clods[i - 1] / (nrep - 1));
  for (i = 1; i <= 66; i++)
    putc('-', isim);
  putc('\n', isim);
  for (i = 1; i <= 66; i++)
    putchar('-');
  putchar('\n');
  acknowl(&isim);
  if (isim != NULL)
    fclose(isim);
  if (outfile != NULL)
    fclose(outfile);
  if (datafile != NULL)
    fclose(datafile);
  if (speedfile != NULL)
    fclose(speedfile);
  if (ipedfile != NULL)
    fclose(ipedfile);
  if (simout != NULL)
    fclose(simout);
  if (final != NULL)
    fclose(final);
  exit(EXIT_SUCCESS);
}  /* isim */



/* End. */
