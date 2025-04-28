/*This file is part of a C version of SLINK a simulation program based
on the LINKAGE package. SLINK was originally described in: D. E.
Weeks, J. Ott, G. M. Lathrop (1990), SLINK: a general simulation
program for linkage analysis, American Journal of Human Genetics
47(1990), A204 (abstract). This version of SLINK is adapted to use the
faster pedigree traversal algorithms described in: R. W. Cottingham
Jr., R. M. Idury, A. A. Schaffer (1993), Faster sequential genetic
linkage computations, American Journal of Human Genetics 53(1993), pp.
252-263.*/
/* This file contains definitions sharec with the othe LINKAGE programs*/
/*March 14, 1994, A.A. Schaffer changed to make declarations more portable*/
/* May 2008 Edited to allow > 32 alleles at a numbered allele locus and
   to remove dependence on maxhap*/

#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <stdio.h>
#include <math.h>
#include <sys/time.h>

#define true            1
#define false           0
#if !defined(NULL)
#define NULL            0
#endif
#define FileNotFound   10  /*From p2c.h*/

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

#define Free free
#define Malloc malloc

#define version         "3.02"   /*version of SLINK*/
/* SOME USER DEFINED CONSTANTS */
/*THE PROGRAM WILL TELL YOU IF THE FOLLOWING TWO CONSTANTS CAN
  BE REDUCED*/
/*IF THE PROGRAM TERMINATES WITH AN ERROR IN RECOMBINATION INCREASE MAXNEED*/


    /*MAXIMUM NUMBER OF RECOMBINATION PROBABILITIES*/
/*THE FOLLOWING SHOULD BE LARGER THAN MININT*/

#define maxcensor       50000   /*MAXIMUM FOR CENSORING ARRAY*/
#define maxsys          40   /*MAXIMUM NUMBER OF SYSTEMS*/
#define maxlocus        8   /*MAXIMUM NUMBER OF LOCI */
#define maxrec          maxlocus   /*MAXIMUM POSSIBLE NUMBER OF RECOMB/MEI*/
#define maxseg          64   /*    (maxlocus-1)        */
/* = 2              I.E. 2 TO THE POWER maxlocus-1 */
#define maxall          129   /*MAX NUMBER OF ALLELES AT A SINGLE LOCUS*/

#define maxind          6900   /*MAXIMUM NUMBER OF INDIVIDUALS*/
#define maxped          1600   /*MAXIMUM NUMBER OF PEDIGREES*/
#define maxchild        16   /*MAXIMUM NUMBER OF FULLSIBS IN A SIBSHIP*/
#define maxloop         2   /*MAXIMUM NUMBER OF LOOPS PER PEDIGREE*/

#define minfreq         0.0

#define affall          2
/*DISEASE ALLELE FOR QUANTITATIVE TRAITS
                               OR AFFECTION STATUS*/
/* QUANTITATIVE TRAIT */
#define maxtrait        3
    /*MAXIMUM NUMBER OF QUANTITATIVE FACTORS AT A SINGLE LOCUS*/

#define missval         0.0   /*MISSING VALUES FOR QUANTITATIVE TRAITS */
/* AFFECTION STATUS */

#define missaff         0   /*MISSING VALUE FOR AFFECTION STATUS */
#define affval          2   /*CODE FOR AFFECTED INDIVIDUAL*/
#define maxliab        120   /*MAXIMUM NUMBER OF LIABILITY CLASSES */
/* BINARY (FACTOR UNION) SYSTEM */
#define maxfact         11
    /*MAXIMUM NUMBER OF BINARY CODES AT A SINGLE LOCUS*/
/* OTHERS */

#define scalemult       2.0   /*SCALE WEIGHT FOR EACH LOCUS INCLUDED*/

#define fitmodel        false
    /*TRUE IF ESTIMATING PARAMETERS OTHER THAN REC*/
#define dostream        true   /*STREAM FILE OUTPUT*/
#define byfamily        false   /*GIVE LOD SCORES BY FAMILY*/

#define zerolike        (-1.0e20)
    /*FOR INCONSISTENT DATA OR RECOMBINATION */
#define log10_          2.30259

#define minint          (-32767)   /*MINIMUM ALLOWED INTEGER*/
/*GRADIENT APPROXIMATIONS*/

#define approximate     false

/*next two typedefs come from p2c*/
typedef char boolean;
typedef short unschar;

unschar **approxarray;

typedef struct censorrec {
  boolean censor[maxcensor - minint];
} censorrec;

typedef double *genotype;
typedef double covmatrix[maxtrait][maxtrait];
typedef int mutarray[3][2];
typedef double thesemeans[maxtrait];
typedef thesemeans means[maxall + 1][maxall];
typedef enum {
  auto_, mauto, sex, msex
} pathway;
typedef enum {
  peelup, peeldown
} direction;

typedef int allelePair[2]; 
typedef int binset;

typedef allelePair num_phenarray[maxall]; /*new representation of a numbered alleles genotype*/
typedef binset bin_phenarray[maxall]; /*new representation of a binary factors genotype*/



typedef enum {
  affection, quantitative, binary_, null_
} locustype;

/* Commented out here because the definition is in sldefs.h 
  typedef struct locusvalues {
  int nallele, format;
  double freq[maxall];
  struct locusvalues *privlocus;
  locustype which;
  union {
    struct {
      double pen[maxall + 1][maxall][3][maxliab];
      int nclass;
    } U0;
    struct {
      int ntrait;
      means pm;
      covmatrix vmat;
      double det, contrait, conmat;
    } U1;
    num_phenarray num_allele;
    bin_phenarray bin_allele;
  } UU;
} locusvalues; */

typedef struct phenotype {
  locustype which;
  binset phenf;
  allelePair alleles;/*new representation of numbered allele genotype*/
  double x[maxtrait];
  boolean missing;
  int aff;
  int liability;
} phenotype;

typedef phenotype *hindphen[maxsys];
typedef phenotype *indphen[maxlocus];


/* a structure of type 'thisarray' stores the conditional
genotype probabilities for an individual. The
probabilities are stored in the array genarray.
We expect this array to be sparse. The field sparseflag
is a boolean array such that sparseflag[i] is nonzero
if and only if genarray[i] is nonzero. */

typedef struct thisarray {
  unsigned char *sparseflag;
  genotype genarray;
} thisarray;

typedef boolean possvect[maxall][maxall];
typedef possvect possarray[maxlocus];

typedef struct information {
  possarray possible;
} information;

/* a record of type this person stores information about one person
   in one pedigree */

/*typedef struct thisperson {
  int id, ped, inloop;
  struct thisperson *pa, *ma, *foff, *nextpa, *nextma;
  thisarray *gen;
  hindphen holdphen;
  indphen phen;
  phenotype *privphen;
  boolean unknown, multi, done, up, male, firstpass;
  information *store;
} thisperson; */



typedef int subhap[maxseg];
typedef double thetarray[maxlocus];
/* typedef double happrob[maxneed]; */
unsigned nuneed; /* Introduced by R. M. Idury, actual size of segprob
                    arrays */
unsigned nuprobclass; /* Introduced by A. A. Schaffer, actual number
			 of probclasses */

typedef struct thetavalues {
  thetarray theta,theta1, theta2;
  double *segprob;
} thetavalues;

 thetavalues *maletheta, *femaletheta;


int maxhaplo, maxfemgen; /*A. A. Schaffer*/
int maxclasssize, maxisozygclass;

/*The following declarations help store a partial correspondence
between genotypes and haplotypes,
haps1 between base[i] and fence[i] stores the left haplotypes
that genotype i can pass on to a child; haps2 stores the right
haplotypes; hind stores the index into an array of recombination probabilities
indicating the probability of this haplotype getting passed on from
genotype i. currentfence is used as a counter to fill
haps1, haps2, hind, base, and fence. */
unsigned currentfence;           /* R. M. Idury */
unsigned int *base, *fence; /*R. M. Idury*/
unsigned short *haps1, *haps2;
unsigned int *hind;
/* The arrays invgenenum1 and invgenenum2 store an effective inverse
to the genenumber mapping above. They convert an index into
genenumber into two haplotypes, replacing seghap1 an seghap2 */
unsigned *invgenenum1, *invgenenum2; /* R. M. Idury */


typedef struct cache{
 unsigned first;
 unsigned last;
} cache;


 pathway thispath;
 boolean informative[maxped];
 boolean *rare, *risk1, *risk2;
 boolean *riskmale;
 censorrec *censorstruct;
 int thisc;
 boolean malechild[maxchild];
 thisarray *thischild[maxchild];
 int nchild;
 int *segstart;
 unsigned int *probstart, *probend;
 boolean nohom[maxlocus];
/* locusvalues *thislocus[maxlocus]; */
 int increment[maxlocus], order[maxlocus];
 unsigned int *nonzgens; /*used to hold nonzero genotypes in pedigree
                            traversal, A.A. Schaffer */
 boolean *flag;  /*R. M. Idury*/  
 double *gene; /*used for local computations in pedigree traversal routines*/
/* The arrays psumcache and qsumcache store conditional probabilities of
   different haplotypes being passed on from p and q respectively.
   They are used in some of the pedigree traversal routines. */
 double *psumcache, *qsumcache;

/*Used in segup to keep track of entries in indpool, invpool,
  and nextpool that correspond to different haplotypes of a child*/
 cache *phapcache1;


/*PEOPLE*/
/* thisperson *person[maxind + 1];
 thisperson *proband[maxped];
 thisperson *looppers[maxped][maxloop][2]; */
/*MUTATION */
 unschar *muthap;
int **genenumber;
/*RECOMBINATION*/
 unsigned *segindex;
 double *tempseg, *tempseg2;
 double *segval;
/*FREQUENCIES */
 thisarray *hapfreq;
/*RISK*/
 int riskall;
/*OTHERS*/
 int risksys, mutsys, mlocus, lastpriv;
 int nuhap;
 int fgeno, mgeno;
 int nuped, totperson;
 double segscale, mutmale, mutfemale, like, alike, distratio;
 boolean interfer, disequi, sexlink, risk, sexdif, readfemale,
	        dolod, firstapprox, firsttime, lasttime;
 FILE *outfile, *ipedfile, *datafile, *stream, *speedfile;
/*ILINK*/
 FILE *final;
 boolean mapping;

/* Local variables for seg: */
/*

struct LOC_seg {
  struct LOC_likelihood *LINK;
  thisperson **p, **q, **r, *child, *father, *mother;
  int fseg, sseg, sstart, send, fstart, fend, nfirst, nsecond, firstseg,
       secondseg;
  double pf, ps;
  thetavalues *firstsex, *secondsex;
} ; */

/*Local variables for likelihood*/
/*
struct LOC_likelihood {
  int thisped;
  thisperson *proband;
  int loopgen[maxloop];
  double homo, hetero;
  int nuscale;
  thisarray *holdpoint[maxloop];
} ;
*/

void malloc_err(char*);

extern void getprobtable();
extern void segup();
extern void segdown();
extern void segsexup();
extern void segsexdown();
extern void   performCheckpoint();
extern void  recoverCheckpoint ();
extern void getvect();
extern void recombination();
extern void getlocations();
extern void readspeed();
extern void invert();
extern void cleanup();
extern void initseg();
extern void exitseg();
extern void segsexctop();
extern void segsextop();
extern void segctop();
extern void segtop();
extern void segcapprox();
extern void msegsexdown();
extern void msegdown();
extern double mapfunction();
extern double getdist();
extern double invdist();
extern void setparam();
extern void inputdata();
extern void oldsegsexup();
extern void oldsegup();
extern void malloc_err();
extern void invertmat();
extern void allocate_thisarray();


