/* Output from p2c, the Pascal-to-C translator */
/* From input file "unknown.p" */
#include<stdio.h>
#include<stdlib.h>
#include<string.h>

/*Changes at Columbia U. indicated by "change"*/
/*10 July 1993*/

#define aVersion         "5.20"   /*PRESENT VERSION OF LINKAGE*/

#define maxlocus        20   /*MAXIMUM NUMBER OF LOCI*/
#define maxall          20   /*MAX NUMBER OF ALLELES AT A SINGLE LOCUS*/

#define maxgeno         (maxall * (maxall + 1) / 2)
    /*MAX NUMBER OF SINGLE LOCUS GENOTYPES*/
/*Different definition than in analysis programs!*/

#define maxind          700   /*MAXIMUM NUMBER OF INDIVIDUALS IN A PEDIGREE*/
#define maxmarriage     3   /*MAXIMUM NUMBER OF MARRIAGES FOR ONE MALE*/

#define maxfact         (maxall + 2)
    /*MAXIMUM NUMBER OF BINARY CODES AT A SINGLE LOCUS*/

#define maxtrait        3
    /*MAXIMUM NUMBER OF QUANTITATIVE FACTORS AT A SINGLE LOCUS*/

#define missval         0.0   /*MISSING VALUES FOR QUANTITATIVE TRAITS*/

#define affall          2
    /*DISEASE ALLELE FOR QUANT. TRAITS OR AFFECTION STATUS*/
#define missaff         0   /*MISSING VALUE FOR AFFECTION STATUS*/
#define affval          2   /*CODE FOR AFFECTED INDIVIDUAL*/
#define maxliab        120   /*MAXIMUM NUMBER OF LIABILITY CLASSES*/

#ifndef NULL
#define NULL             0
#endif
#define FileNotFound   10
#define true            1
#define false           0
#ifdef vms
#define EXIT_SUCCESS    1
#define EXIT_FAILURE    (02000000000L)
#endif
#define Malloc          malloc
#define Free            free

#if !defined(ONE_ERROR_ONLY)
#define ONE_ERROR_ONLY   0   /*print one sure error, or all possible errors*/
#endif

#if !defined(DEPTH_MULTIPLE)
#define DEPTH_MULTIPLE   3
#endif

/*from p2c*/

typedef char boolean;

typedef boolean genotype[maxgeno];
typedef enum {
  peelup, peeldown
} direction;
typedef long binset;

typedef binset phenarray[maxall];
typedef boolean possvect[maxall][maxall];
typedef possvect possarray[maxlocus];
typedef enum {
  affection, quantitative, binary_
} locustype;

typedef struct locusvalues {
  long nallele;
  locustype which;
  union {
    long ntrait;
    struct {
      double pen[maxall + 1][maxall][3][maxliab];
      long nclass;
    } U0;
    struct {
      phenarray allele;
      long nfactor, format;
    } U2;
  } UU;
} locusvalues;

typedef struct phenotype {
  locustype which;
  union {
    struct {
      double x[maxtrait];
      boolean missing;
    } U1;
    struct {
      long aff, liability;
    } U0;
    binset phenf;
  } UU;
} phenotype;

typedef phenotype *indphen[maxlocus];

typedef struct thisarray {
  genotype genarray;
} thisarray;

typedef struct information {
  possarray possible;
} information;

typedef struct thisperson {
  long id, paid, maid, offid, npaid, nmaid, sex, profield, oldped, nseq;
  struct thisperson *pa, *ma, *foff, *nextpa, *nextma;
  thisarray *gen;
  indphen phen;
  information *store;
  boolean thisunknown[maxlocus];
  boolean unknown, multi, done, up, male;
} thisperson;

typedef enum {
  a, b
} haplotype;
typedef long subhap[(long)b - (long)a + 1];


static subhap seghap[maxgeno];
static locusvalues *thislocus[maxlocus];
static thisperson *person[maxind + 1];
static thisperson *proband, *loop1, *loop2;
static long genenumber[maxgeno][maxgeno];
static long risksys, mutsys, nsystem;
static short fgeno, mgeno;
static long nsequence, newped, whichsys, totperson;
static boolean sexlink, risk, disequi;
static FILE *speedfile, *datafile, *pedfile, *ipedfile;
static genotype gene;
static double one;   /*changed*/
static boolean makehomozygous;   /* Change - Added 7/8/93 */
static int numind;  /*number of individuals in pedigree*/
static int depth;  /*depth of recursion*/

/*The following two variables were added by A. A. Schaffer to
 help pinpoint Mendelian incompatibilities*/
static boolean incompat_traversal; /*found incompatibility on this traversal*/
static boolean first_incompat; /*Is this the first incomptibility */
static boolean detected_incompat[maxind+1]; /* array to record whether
                                               an incompatibility with this
                                               person has already been
                                               reported*/

/* Two routines taken from */
/* "p2c"  Copyright (C) 1989, 1990, 1991 Free Software Foundation.
 * By Dave Gillespie, daveg@csvax.cs.caltech.edu.  Version --VERSION--.
 * This file may be copied, modified, etc. in any way.  It is not restricted
 * by the licence agreement accompanying p2c itself.
 */


/* Check if at end of file, using Pascal "eof" semantics.  End-of-file for
   stdin is broken; remove the special case for it to be broken in a
   different way. */

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


static void respond()
{
  /*Change - new*/
  printf("*** Press <Enter> to continue\n");
  scanf("%*[^\n]");
  getchar();
}


static void inputerror(nerror, par1, par2)
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
    printf(
      "Error detected reading datafile. Linkage disequilibrium is not allowed with this program\n");
    break;

  case 42:
    printf("Locus %5ld in lod score list exceeds nlocus %5ld\n", par1, par2);
    break;

  case 43:
    printf("Illegal locus number %5ld in lod score list\n", par1);
    break;

  case 44:
    printf("Error detected reading pedigree record %2ld. One 0 allele\n",
	   par1);
    break;
  }
  respond();   /*changed*/
}


static void inputwarning(nwarning, par1, par2)
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
  respond();   /*changed*/
}


static void writespeed()
{
  long i, j, a_, b_, FORLIM;
  thisperson *WITH;
  information *WITH1;
  long FORLIM1, FORLIM2, FORLIM3;

  FORLIM = totperson;
  for (i = 1; i <= FORLIM; i++) {
    WITH = person[i];
    if (WITH->unknown && WITH->foff != NULL) {
      WITH1 = WITH->store;
      fprintf(speedfile, "id%7ld\n", WITH->nseq);
      FORLIM1 = nsystem;
      for (j = 1; j <= FORLIM1; j++) {
	FORLIM2 = thislocus[j - 1]->nallele;
	for (a_ = 1; a_ <= FORLIM2; a_++) {
	  FORLIM3 = thislocus[j - 1]->nallele;
	  for (b_ = 1; b_ <= FORLIM3; b_++) {
	    if (WITH1->possible[j - 1][a_ - 1][b_ - 1])
	      fprintf(speedfile, "%3ld%3ld%3ld\n", j, a_, b_);
	  }
	}
      }
    }
  }
}  /*writespeed*/


static void writeped()
{
  long i, j, k, a_, b_, FORLIM;
  thisperson *WITH;
  long FORLIM1;
  phenotype *WITH1;
  locusvalues *WITH2;
  long FORLIM2;

  FORLIM = totperson;
  for (i = 1; i <= FORLIM; i++) {
    WITH = person[i];
    fprintf(ipedfile, "%7ld%5ld%5ld%5ld%5ld%5ld",
	    WITH->oldped, WITH->id, WITH->paid, WITH->maid, WITH->offid,
	    WITH->npaid);
    fprintf(ipedfile, "%5ld%2ld%2ld ", WITH->nmaid, WITH->sex, WITH->profield);
    FORLIM1 = nsystem;
    for (j = 1; j <= FORLIM1; j++) {
      WITH1 = WITH->phen[j - 1];
      WITH2 = thislocus[j - 1];
      if (WITH2->which == binary_) {
	if (WITH2->UU.U2.format == 2) {
	  FORLIM2 = WITH2->UU.U2.nfactor;
	  for (k = 1; k <= FORLIM2; k++) {
	    if ((unsigned long)k < 32 && ((1L << k) & WITH1->UU.phenf) != 0)
	      fprintf(ipedfile, " 1");
	    else
	      fprintf(ipedfile, " 0");
	  }
	} else {
	  a_ = 0;
	  b_ = 0;
	  FORLIM2 = WITH2->nallele;
	  for (k = 1; k <= FORLIM2; k++) {
	    if ((unsigned long)k < 32 && ((1L << k) & WITH1->UU.phenf) != 0) {
	      if (a_ == 0)
		a_ = k;
	      else
		b_ = k;
	    }
	  }
	  if (b_ == 0)
	    b_ = a_;
	  fprintf(ipedfile, "%3ld%3ld", a_, b_);
	}
      } else if (WITH2->which == quantitative) {
	if (!sexlink || !WITH->male) {
	  FORLIM2 = WITH2->UU.ntrait;
	  for (k = 0; k < FORLIM2; k++)
	    fprintf(ipedfile, " %9.4f", WITH1->UU.U1.x[k]);
	} else {
	  FORLIM2 = WITH2->UU.ntrait;
	  for (k = 1; k <= FORLIM2; k++)
	    fprintf(ipedfile, " %9ld", WITH1->UU.U0.aff);
	}
      } else {
	fprintf(ipedfile, "%2ld", WITH1->UU.U0.aff);
	if (WITH2->UU.U0.nclass != 1)
	  fprintf(ipedfile, "%4ld", WITH1->UU.U0.liability);
      }
      if (j != nsystem)
	putc(' ', ipedfile);
    }
    putc('\n', ipedfile);
  }
}  /*writeped*/


static void infer()
{
  long i, j, k, l, kposs, lposs, count, pacount, macount;
  boolean someknown;
  long FORLIM;
  thisperson *WITH;
  information *WITH1;
  long FORLIM1;
  locusvalues *WITH2;
  long FORLIM2, FORLIM3;

  FORLIM = totperson;
  for (i = 1; i <= FORLIM; i++) {
    if (person[i]->unknown) {
      WITH = person[i];
      WITH1 = WITH->store;
      FORLIM1 = nsystem;
      for (j = 0; j < FORLIM1; j++) {
	if (thislocus[j]->which == binary_) {
	  if (WITH->phen[j]->UU.phenf == 0) {
	    WITH2 = thislocus[j];
	    count = 0;
	    FORLIM2 = WITH2->nallele;
	    for (k = 1; k <= FORLIM2; k++) {
	      FORLIM3 = WITH2->nallele;
	      for (l = k; l <= FORLIM3; l++) {
		if (WITH1->possible[j][k - 1][l - 1]) {
		  kposs = k;
		  lposs = l;
		  count++;
		}
	      }
	    }
	    if (count == 1) {
	      if (sexlink && WITH->male)
		WITH->phen[j]->UU.phenf = WITH2->UU.U2.allele[lposs - 1];
	      else
		WITH->phen[j]->UU.phenf = WITH2->UU.U2.allele[kposs - 1] |
					  WITH2->UU.U2.allele[lposs - 1];
	    }
	  }
	}
      }
      count = 0;
      FORLIM1 = nsystem;
      for (j = 0; j < FORLIM1; j++) {
	if (thislocus[j]->which != binary_)
	  count++;
	else if (WITH->phen[j]->UU.phenf == 0)
	  count++;
      }
      WITH->unknown = (count != 0);
    }
  }
  FORLIM = totperson;
  /*Infer children when parents are homozygotes*/
  for (i = 1; i <= FORLIM; i++) {
    if (person[i]->foff == NULL) {
      WITH = person[i];
      FORLIM1 = nsystem;
      for (j = 0; j < FORLIM1; j++) {
	WITH2 = thislocus[j];
	if (WITH->phen[j]->which == binary_) {
	  if (WITH->phen[j]->UU.phenf == 0) {
	    if (WITH->pa != NULL) {
	      pacount = 0;
	      macount = 0;
	      FORLIM2 = thislocus[j]->nallele;
	      for (k = 1; k <= FORLIM2; k++) {
		if ((WITH2->UU.U2.allele[k - 1] &
		     (~WITH->pa->phen[j]->UU.phenf)) == 0) {
		  kposs = k;
		  pacount++;
		}
	      }
	      FORLIM2 = thislocus[j]->nallele;
	      for (l = 1; l <= FORLIM2; l++) {
		if ((WITH2->UU.U2.allele[l - 1] &
		     (~WITH->ma->phen[j]->UU.phenf)) == 0) {
		  lposs = l;
		  macount++;
		}
	      }
	      if (macount == 1 && pacount == 1 && !(WITH->male && sexlink))
		WITH->phen[j]->UU.phenf = WITH2->UU.U2.allele[kposs - 1] |
					  WITH2->UU.U2.allele[lposs - 1];
	      else if (macount == 1 && WITH->male && sexlink)
		WITH->phen[j]->UU.phenf = WITH2->UU.U2.allele[lposs - 1];
	    }
	  }
	}
      }
    }
  }
  /*Replace by homozygotes if all unknown in a pedigree*/
  if (!makehomozygous)   /*change - added*/
    return;
  FORLIM = nsystem;
  for (j = 0; j < FORLIM; j++) {
    WITH2 = thislocus[j];
    if (WITH2->which == binary_ && WITH2->UU.U2.format == 3)
    {   /*change - 'format=3' added*/
      someknown = false;
      FORLIM1 = totperson;
      for (i = 1; i <= FORLIM1; i++) {
	if (person[i]->phen[j]->UU.phenf != 0)
	  someknown = true;
      }
      if (!someknown) {
	FORLIM1 = totperson;
	for (i = 1; i <= FORLIM1; i++)
	  person[i]->phen[j]->UU.phenf = WITH2->UU.U2.allele[0];
      }
    }
  }
}  /*infer*/


static void getunknown()
{
  long i, j, n, ahap, bhap, FORLIM, FORLIM1;
  thisperson *WITH;
  information *WITH1;
  locusvalues *WITH2;
  long FORLIM2, FORLIM3;

  FORLIM = totperson;
  for (i = 1; i <= FORLIM; i++)
    person[i]->unknown = false;
  FORLIM = totperson;
  for (i = 1; i <= FORLIM; i++) {
    FORLIM1 = nsystem;
    for (j = 0; j < FORLIM1; j++)
      person[i]->thisunknown[j] = false;
  }
  FORLIM = totperson;
  for (i = 1; i <= FORLIM; i++) {
    WITH = person[i];
    FORLIM1 = nsystem;
    for (j = 0; j < FORLIM1; j++) {
      if (thislocus[j]->which == binary_) {
	if (WITH->phen[j]->UU.phenf == 0)
	  WITH->thisunknown[j] = true;
	else if (thislocus[j]->which == quantitative) {
	  if (WITH->phen[j]->UU.U1.x[0] == missval)
	    WITH->thisunknown[j] = true;
	  else if (WITH->phen[j]->UU.U0.aff == missaff)
	    WITH->thisunknown[j] = true;
	}
      }
      if (WITH->thisunknown[j])
	WITH->unknown = true;
    }
  }
  FORLIM = totperson;
  for (i = 1; i <= FORLIM; i++) {
    WITH = person[i];
    if (WITH->unknown) {
      WITH->store = (information *)Malloc(sizeof(information));
      WITH1 = WITH->store;
      FORLIM1 = nsystem;
      for (n = 0; n < FORLIM1; n++) {
	WITH2 = thislocus[n];
	FORLIM2 = WITH2->nallele;
	for (ahap = 0; ahap < FORLIM2; ahap++) {
	  FORLIM3 = WITH2->nallele;
	  for (bhap = 0; bhap < FORLIM3; bhap++)
	    WITH1->possible[n][ahap][bhap] = true;
	}
      }
    }
  }
}  /*getunknown*/


static void getlocation(thislocus)
locusvalues *thislocus;
{
  long ahap, bhap, here, FORLIM, FORLIM1;

  here = 0;
  FORLIM = thislocus->nallele;
  for (ahap = 1; ahap <= FORLIM; ahap++) {
    FORLIM1 = thislocus->nallele;
    for (bhap = ahap; bhap <= FORLIM1; bhap++) {
      here++;
      genenumber[ahap - 1][bhap - 1] = here;
      genenumber[bhap - 1][ahap - 1] = here;
      seghap[here - 1][0] = ahap;
      seghap[here - 1][(long)b - (long)a] = bhap;
    }
  }
}  /*getlocation*/


/*  Local variables for getphenotype: */
struct LOC_getphenotype {
  thisperson **p;
} ;

 void readbin(phen, thislocus, LINK)
phenotype **phen;
locusvalues *thislocus;
struct LOC_getphenotype *LINK;
{
  long i, j;
  phenotype *WITH1;
  long FORLIM;

  WITH1 = *phen;
  WITH1->which = binary_;
  WITH1->UU.phenf = 0;
  FORLIM = thislocus->UU.U2.nfactor;
  for (i = 1; i <= FORLIM; i++) {
    fscanf(pedfile, "%ld", &j);
    if (j != 0 && j != 1)
      inputerror(14L, (*LINK->p)->id, j);
    if (j == 1)
      WITH1->UU.phenf = ((long)WITH1->UU.phenf) | (1L << ((int)i));
  }
}


 void readnumber(phen, LINK)
phenotype **phen;
struct LOC_getphenotype *LINK;
{
  long j, k;
  phenotype *WITH;

  WITH = *phen;
  WITH->which = binary_;
  WITH->UU.phenf = 0;
  fscanf(pedfile, "%ld%ld", &j, &k);
  if (j > maxall)
    inputerror(16L, (*LINK->p)->id, j);
  if (j < 0)
    inputerror(17L, (*LINK->p)->id, j);
  if (k > maxall)
    inputerror(16L, (*LINK->p)->id, k);
  if (k < 0)
    inputerror(17L, (*LINK->p)->id, k);
  if ((j == 0 || k == 0) && j != k) {
    inputerror(44L, (*LINK->p)->id, j);
    return;
  }
  if (j != 0)
    WITH->UU.phenf = ((long)WITH->UU.phenf) | (1L << ((int)j));
  if (k != 0)
    WITH->UU.phenf = ((long)WITH->UU.phenf) | (1L << ((int)k));
}

 void readaff(phen, thislocus, LINK)
phenotype **phen;
locusvalues *thislocus;
struct LOC_getphenotype *LINK;
{
  long thisval;
  phenotype *WITH;

  WITH = *phen;
  WITH->which = affection;
  fscanf(pedfile, "%ld", &thisval);
  if (thisval == missaff)
    WITH->UU.U0.aff = 0;
  else {
    if (thisval == affval)
      WITH->UU.U0.aff = 2;
    else {
      if (thisval != 1)
	inputwarning(3L, (*LINK->p)->id, thisval);
      WITH->UU.U0.aff = 1;
    }
  }
  if (thislocus->UU.U0.nclass == 1)
    WITH->UU.U0.liability = 1;
  else
    fscanf(pedfile, "%ld", &WITH->UU.U0.liability);
  if (WITH->UU.U0.liability > thislocus->UU.U0.nclass)
    inputerror(26L, (*LINK->p)->id, WITH->UU.U0.liability);
  if (WITH->UU.U0.liability <= 0)
    inputerror(27L, (*LINK->p)->id, WITH->UU.U0.liability);
}

 void readquan(phen, thislocus, LINK)
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
    FORLIM = thislocus->UU.ntrait;
    for (i = 0; i < FORLIM; i++)
      fscanf(pedfile, "%lg", &WITH->UU.U1.x[i]);
    WITH->UU.U1.missing = true;
    FORLIM = thislocus->UU.ntrait;
    for (i = 0; i < FORLIM; i++) {
      if (WITH->UU.U1.x[i] != missval)
	WITH->UU.U1.missing = false;
    }
    return;
  }
  WITH->which = affection;
  fscanf(pedfile, "%lg", &xval);
  if (xval == missval)
    WITH->UU.U0.aff = missaff;
  else {
    if (xval == affall)
      WITH->UU.U0.aff = affall;
    else
      WITH->UU.U0.aff = -11;
  }
  WITH->UU.U0.liability = 1;
  FORLIM = thislocus->UU.ntrait;
  for (i = 2; i <= FORLIM; i++)
    fscanf(pedfile, "%lg", &xval);
}


 void getphenotype(p_)
thisperson **p_;
{
  struct LOC_getphenotype V;
  long thisread, system;
  thisperson *WITH;
  long FORLIM;

  V.p = p_;
  WITH = *V.p;
  FORLIM = nsystem;
  for (thisread = 1; thisread <= FORLIM; thisread++) {
    system = thisread;
    WITH->phen[system - 1] = NULL;
    WITH->phen[system - 1] = (phenotype *)Malloc(sizeof(phenotype));
    switch (thislocus[system - 1]->which) {

    case quantitative:
      readquan(&WITH->phen[system - 1], thislocus[system - 1], &V);
      break;

    case affection:
      readaff(&WITH->phen[system - 1], thislocus[system - 1], &V);
      break;

    case binary_:
      if (thislocus[system - 1]->UU.U2.format == 3)
	readnumber(&WITH->phen[system - 1], &V);
      else
	readbin(&WITH->phen[system - 1], thislocus[system - 1], &V);
      break;
    }
  }
}  /*getphenotype*/


void getind(id, seq)
long *id, *seq;
{
  thisperson *WITH;

  *id = 0;
  fscanf(pedfile, "%ld", seq);
  if (*seq == 0)
    return;
  *id = *seq;
  if (*id > maxind)
    inputerror(13L, *id, *id);
  if (person[*id] != NULL)
    return;
  numind++;
  person[*id] = (thisperson *)Malloc(sizeof(thisperson));
  WITH = person[*id];
  WITH->gen = (thisarray *)Malloc(sizeof(thisarray));
  WITH->nseq = *seq + nsequence;
}  /*getind*/


 void multimarriage(p)
thisperson **p;
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
}  /*multimarriage*/


 void getrest()
{
  char whichchr;

  whichchr = ' ';
  while (!(P_eoln(pedfile) || whichchr == '\f')) {
    whichchr = getc(pedfile);
    if (whichchr == '\n')
      whichchr = ' ';
  }
  fscanf(pedfile, "%*[^\n]");
  getc(pedfile);
}  /*getrest*/


static void readped()
{
  long i, newid, thisone, thisped, tempid, FORLIM;
  thisperson *WITH;


  for (i = 0; i <= maxind; i++)
    person[i] = NULL;
  totperson = 0;
  loop1 = NULL;
  loop2 = NULL;
  proband = NULL;
  thisped = newped;
  printf("Ped. %2ld\n", thisped);   /*change - added*/
  while (!P_eof(pedfile) && thisped == newped) {
    totperson++;
    getind(&thisone, &tempid);
    if (proband == NULL)
      proband = person[thisone];
    WITH = person[thisone];
    WITH->id = tempid;
    WITH->oldped = newped;
    getind(&newid, &WITH->paid);
    WITH->pa = person[newid];
    getind(&newid, &WITH->maid);
    WITH->ma = person[newid];
    getind(&newid, &WITH->offid);
    WITH->foff = person[newid];
    getind(&newid, &WITH->npaid);
    WITH->nextpa = person[newid];
    getind(&newid, &WITH->nmaid);
    WITH->nextma = person[newid];
    WITH->store = NULL;
    WITH->up = false;
    fscanf(pedfile, "%ld", &WITH->sex);
    if (WITH->sex != 1 && WITH->sex != 2)
      inputerror(11L, WITH->id, WITH->sex);
    WITH->male = (WITH->sex == 1);
    fscanf(pedfile, "%ld", &WITH->profield);
    if (WITH->profield == 1)
      proband = person[thisone];
    else if (WITH->profield == 2) {
      if (loop2 == NULL)
	loop2 = person[thisone];
      else
	loop1 = person[thisone];
    }
    getphenotype(&person[thisone]);
    getrest();
    if (!P_eof(pedfile))
      fscanf(pedfile, "%ld", &newped);
  }
  nsequence += totperson;
  if (loop2 != NULL && loop1 == NULL)
    loop1 = proband;
  FORLIM = totperson;
  for (thisone = 1; thisone <= FORLIM; thisone++)
    multimarriage(&person[thisone]);
}  /*readped*/


/* Local variables for readloci: */
struct LOC_readloci {
  long whichtype;
} ;

/* Local variables for getlocus: */
struct LOC_getlocus {
  struct LOC_readloci *LINK;
  long system;
} ;

 void getquan(locus, LINK)
locusvalues **locus;
struct LOC_getlocus *LINK;
{
  long i;
  locusvalues *WITH;

  WITH = *locus;
  fscanf(datafile, "%ld%*[^\n]", &WITH->UU.ntrait);
  getc(datafile);
  if (WITH->UU.ntrait > maxtrait)
    inputerror(31L, LINK->system, WITH->UU.ntrait);
  if (WITH->UU.ntrait <= 0)
    inputerror(32L, LINK->system, WITH->UU.U0.nclass);
  for (i = 1; i <= 3; i++) {
    fscanf(datafile, "%*[^\n]");
    getc(datafile);
  }
}  /*getquan*/


 void getpen(locus, LINK)
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
	if ((unsigned)WITH->UU.U0.pen[i][j][2][l] > one)
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
      for (i = 0; i < FORLIM1; i++)
	fscanf(datafile, "%lg", &WITH->UU.U0.pen[0][i][2][l]);
      if ((unsigned)WITH->UU.U0.pen[0][j - 1][2][l] > one)
	inputerror(30L, LINK->system, LINK->system);
      FORLIM1 = WITH->nallele;
      for (i = 0; i < FORLIM1; i++)
	WITH->UU.U0.pen[0][i][1][l] = 1.0 - WITH->UU.U0.pen[0][i][2][l];
      fscanf(datafile, "%*[^\n]");
      getc(datafile);
    }
  }
}  /*getpen*/


 void getbin(locus, LINK)
locusvalues **locus;
struct LOC_getlocus *LINK;
{
  long i, j, k;
  locusvalues *WITH;
  long FORLIM, FORLIM1;

  WITH = *locus;
  fscanf(datafile, "%ld%*[^\n]", &WITH->UU.U2.nfactor);
  getc(datafile);
  if (WITH->UU.U2.nfactor > maxfact)
    inputerror(8L, LINK->system, WITH->UU.U2.nfactor);
  if (WITH->UU.U2.nfactor <= 0)
    inputerror(9L, LINK->system, WITH->UU.U2.nfactor);
  FORLIM = WITH->nallele;
  for (i = 0; i < FORLIM; i++)
    WITH->UU.U2.allele[i] = 0;
  FORLIM = WITH->nallele;
  for (i = 0; i < FORLIM; i++) {
    FORLIM1 = WITH->UU.U2.nfactor;
    for (j = 1; j <= FORLIM1; j++) {
      fscanf(datafile, "%ld", &k);
      if (k == 1)
	WITH->UU.U2.allele[i] = ((long)WITH->UU.U2.allele[i]) | (1L << ((int)j));
    }
  }
  fscanf(datafile, "%*[^\n]");
  getc(datafile);
}  /*getbin*/


 void getnumber(locus, LINK)
locusvalues **locus;
struct LOC_getlocus *LINK;
{
  long i;
  locusvalues *WITH;
  long FORLIM;

  WITH = *locus;
  FORLIM = WITH->nallele;
  for (i = 1; i <= FORLIM; i++)
    WITH->UU.U2.allele[i - 1] = 1L << ((int)i);
}  /*getnumber*/


 void getlocus(system_, LINK)
long system_;
struct LOC_readloci *LINK;
{
  struct LOC_getlocus V;
  locusvalues *WITH;


  V.LINK = LINK;
  V.system = system_;
  thislocus[V.system - 1] = (locusvalues *)Malloc(sizeof(locusvalues));
  WITH = thislocus[V.system - 1];
  fscanf(datafile, "%ld%ld", &LINK->whichtype, &WITH->nallele);
  switch (LINK->whichtype) {

  case 0:
    WITH->which = quantitative;
    break;

  case 1:
    WITH->which = affection;
    break;

  case 2:
    WITH->which = binary_;
    WITH->UU.U2.format = 2;
    break;

  case 3:
    WITH->which = binary_;
    WITH->UU.U2.format = 3;
    break;
  }
  fscanf(datafile, "%*[^\n]");
  getc(datafile);
  if (!disequi) {
    fscanf(datafile, "%*[^\n]");
    getc(datafile);
  }
  switch (WITH->which) {

  case quantitative:
    getquan(&thislocus[V.system - 1], &V);
    break;

  case affection:
    getpen(&thislocus[V.system - 1], &V);
    break;

  case binary_:
    if (WITH->UU.U2.format == 2)
      getbin(&thislocus[V.system - 1], &V);
    else
      getnumber(&thislocus[V.system - 1], &V);
    break;
  }
  if (risk && V.system == risksys) {
    fscanf(datafile, "%*[^\n]");
    getc(datafile);
  }
}  /*getlocus*/


static void readloci()
{
  struct LOC_readloci V;
  long i, coupling, autosomal;
  double mu;
  long FORLIM;


  fscanf(datafile, "%ld%ld%ld%*[^\n]", &nsystem, &risksys, &autosomal);
  getc(datafile);
  if (nsystem > maxlocus)
    inputerror(0L, nsystem, nsystem);
  if (nsystem <= 0)
    inputerror(1L, nsystem, nsystem);
  risk = (risksys != 0);
  sexlink = (autosomal == 1);
  printf("YOU ARE USING LINKAGE (V%s) WITH%3ld-POINT", aVersion, nsystem);
  if (sexlink)
    printf(" SEXLINKED DATA\n");
  else
    printf(" AUTOSOMAL DATA\n");
  fscanf(datafile, "%ld%lg%lg%ld%*[^\n]", &mutsys, &mu, &mu, &coupling);
  getc(datafile);
  disequi = (coupling == 1);
  fscanf(datafile, "%*[^\n]");
  getc(datafile);
  FORLIM = nsystem;
  for (i = 1; i <= FORLIM; i++)
    getlocus(i, &V);
}  /*readloci*/


static void cleanup(p)
thisperson **p;
{
  long i, j;
  thisperson *WITH;
  information *WITH1;
  thisarray *WITH2;
  long FORLIM, FORLIM1;

  WITH = *p;
  if (!WITH->unknown)
    return;
  WITH1 = WITH->store;
  WITH2 = WITH->gen;
  if (sexlink && WITH->male) {
    FORLIM = mgeno;
    for (i = 0; i < FORLIM; i++) {
      if (!WITH2->genarray[i])
	WITH1->possible[whichsys - 1][0][i] = false;
    }
    FORLIM = mgeno;
    for (i = 1; i < FORLIM; i++) {
      FORLIM1 = mgeno;
      for (j = 0; j < FORLIM1; j++)
	WITH1->possible[whichsys - 1][i][j] = false;
    }
    return;
  }
  FORLIM = fgeno;
  for (i = 0; i < FORLIM; i++) {
    if (!WITH2->genarray[i])
      WITH1->possible[whichsys - 1][seghap[i][0] - 1]
	[seghap[i][(long)b - (long)a] - 1] = false;
    WITH1->possible[whichsys - 1][seghap[i][(long)b - (long)a] - 1]
      [seghap[i][0] - 1] = WITH1->possible[whichsys - 1][seghap[i][0] - 1]
      [seghap[i][(long)b - (long)a] - 1];
  }
}  /*cleanup*/


static void getgene(system, p, phen)
long system;
thisperson *p;
phenotype **phen;
{
  long here, i, j;
  locusvalues *WITH;
  long FORLIM;
  phenotype *WITH1;
  long FORLIM1;

  here = 0;
  WITH = thislocus[system - 1];
  if (sexlink && p->male) {
    FORLIM = WITH->nallele;
    for (i = 1; i <= FORLIM; i++) {
      here++;
      switch (WITH->which) {

      case quantitative:
	WITH1 = phen[system - 1];
	if (i == affall)
	  p->gen->genarray[here - 1] = (WITH1->UU.U0.aff == affall ||
					WITH1->UU.U0.aff == missaff);
	else
	  p->gen->genarray[here - 1] = (WITH1->UU.U0.aff != affall ||
					WITH1->UU.U0.aff == missaff);
	break;

      case affection:
	WITH1 = phen[system - 1];
	p->gen->genarray[here - 1] = (WITH->UU.U0.pen[0][i - 1]
				      [WITH1->UU.U0.aff]
				      [WITH1->UU.U0.liability - 1] > 0.0);
	break;

      case binary_:
	WITH1 = phen[system - 1];
	p->gen->genarray[here - 1] =
	  (WITH1->UU.phenf == WITH->UU.U2.allele[i - 1] ||
	   WITH1->UU.phenf == 0);
	break;
      }
    }
    return;
  }
  FORLIM = WITH->nallele;
  for (i = 1; i <= FORLIM; i++) {
    FORLIM1 = WITH->nallele;
    for (j = i - 1; j < FORLIM1; j++) {
      here++;
      switch (WITH->which) {

      case quantitative:
	p->gen->genarray[here - 1] = true;
	break;

      case affection:
	WITH1 = phen[system - 1];
	p->gen->genarray[here - 1] = (WITH->UU.U0.pen[i][j][WITH1->UU.U0.aff]
				      [WITH1->UU.U0.liability - 1] > 0.0);
	break;

      case binary_:
	WITH1 = phen[system - 1];
	p->gen->genarray[here - 1] = (WITH1->UU.phenf ==
	      (WITH->UU.U2.allele[i - 1] | WITH->UU.U2.allele[j]) ||
	    WITH1->UU.phenf == 0);
	break;
      }
    }
  }
}  /*getgene*/


static void collapsedown ();

/* Local variables for seg: */
struct LOC_seg {
  thisperson **p, **q, **r, *child, *father, *mother;
  subhap firsthap, secondhap;
  long nfirst, nsecond;
} ;


 boolean segfun(child, first, second, LINK)
thisperson **child;
long first, second;
struct LOC_seg *LINK;
{
  boolean temp;
  haplotype thishap, thathap;
  thisarray *WITH;

  temp = false;
  WITH = (*child)->gen;
  if (!sexlink) {
    for (thishap = a;
	 (long)thishap <= (long)b;
	 thishap = (haplotype)((long)thishap + 1)) {
      for (thathap = a;
	   (long)thathap <= (long)b;
	   thathap = (haplotype)((long)thathap + 1))
	temp = (temp ||
	    WITH->genarray[genenumber[LINK->secondhap[(long)thishap - (long)a] - 1]
			   [LINK->firsthap[(long)thathap - (long)a] - 1] - 1]);
    }
    return temp;
  }
  if ((*child)->male) {
    if ((*LINK->p)->male) {
      for (thathap = a;
	   (long)thathap <= (long)b;
	   thathap = (haplotype)((long)thathap + 1))
	temp = (temp ||
		WITH->genarray[LINK->secondhap[(long)thathap - (long)a] - 1]);
    } else {
      for (thathap = a;
	   (long)thathap <= (long)b;
	   thathap = (haplotype)((long)thathap + 1))
	temp = (temp ||
		WITH->genarray[LINK->firsthap[(long)thathap - (long)a] - 1]);
    }
    return temp;
  }
  if ((*LINK->p)->male) {
    for (thathap = a;
	 (long)thathap <= (long)b;
	 thathap = (haplotype)((long)thathap + 1))
      temp = (temp ||
	  WITH->genarray[genenumber[LINK->secondhap[(long)thathap - (long)a] - 1]
			 [first - 1] - 1]);
  } else {
    for (thathap = a;
	 (long)thathap <= (long)b;
	 thathap = (haplotype)((long)thathap + 1))
      temp = (temp ||
	  WITH->genarray[genenumber[LINK->firsthap[(long)thathap - (long)a] - 1]
			 [second - 1] - 1]);
  }
  return temp;
}  /*segfun*/


 void segup(LINK)
struct LOC_seg *LINK;
{
  long first, second;
  boolean segval, val;
  thisarray *WITH;
  long FORLIM;
  thisarray *WITH1;
  long FORLIM1;
  boolean compat;
  int j;

  if ((*LINK->p)->male) {
    LINK->nfirst = mgeno;
    LINK->nsecond = fgeno;
  } else {
    LINK->nfirst = fgeno;
    LINK->nsecond = mgeno;
  }
  WITH = (*LINK->p)->gen;
  FORLIM = LINK->nfirst;
  for (first = 1; first <= FORLIM; first++) {
    if (WITH->genarray[first - 1]) {
      segval = false;
      memcpy(LINK->firsthap, seghap[first - 1], sizeof(subhap));
      WITH1 = (*LINK->q)->gen;
      FORLIM1 = LINK->nsecond;
      for (second = 1; second <= FORLIM1; second++) {
	if (WITH1->genarray[second - 1]) {
	  memcpy(LINK->secondhap, seghap[second - 1], sizeof(subhap));
	  val = true;
	  LINK->child = LINK->father->foff;
	  do {
	    if (LINK->mother == LINK->child->ma)
	      val = segfun(&LINK->child, first, second, LINK);
	    LINK->child = LINK->child->nextpa;
	  } while (val && LINK->child != NULL);
	  segval = (val || segval);
	}
      }
      WITH->genarray[first - 1] = segval;
    }
  }
  if (!incompat_traversal) {
    compat = false;
    FORLIM1 = fgeno;
    for (j = 0; j < FORLIM1; j++)
       compat = (compat || WITH->genarray[j]);
    if (!compat && (!detected_incompat[(*LINK->p)->id]) &&
	(!detected_incompat[(*LINK->q)->id])) {
      printf("\n One incompatibility involves the family in which person %ld is a parent\n",(*LINK->p)->id);
      incompat_traversal = true;
      detected_incompat[(*LINK->p)->id] = true;
      detected_incompat[(*LINK->q)->id] = true;
      if (ONE_ERROR_ONLY)
	exit(EXIT_FAILURE);
      if (first_incompat)
	respond();
    }
  }
  cleanup(LINK->q);
  LINK->child = LINK->father->foff;
  do {
    if (LINK->child->ma == LINK->mother)
      cleanup(&LINK->child);
    LINK->child = LINK->child->nextpa;
  } while (LINK->child != NULL);
}  /*segup*/


 void segdown(LINK)
struct LOC_seg *LINK;
{
  long first, second, here;
  boolean val;
  haplotype thishap, thathap;
  long FORLIM;
  thisarray *WITH, *WITH1;
  long FORLIM1;
  boolean compat;
  int j;

  FORLIM = fgeno;
  for (first = 0; first < FORLIM; first++)
    gene[first] = false;
  WITH = (*LINK->p)->gen;
  FORLIM = mgeno;
  for (first = 1; first <= FORLIM; first++) {
    if (WITH->genarray[first - 1]) {
      memcpy(LINK->firsthap, seghap[first - 1], sizeof(subhap));
      WITH1 = (*LINK->q)->gen;
      FORLIM1 = fgeno;
      for (second = 1; second <= FORLIM1; second++) {
	if (WITH1->genarray[second - 1]) {
	  memcpy(LINK->secondhap, seghap[second - 1], sizeof(subhap));
	  val = (WITH1->genarray[second - 1] &&
		 (*LINK->p)->gen->genarray[first - 1]);
	  LINK->child = LINK->father->foff;
	  do {
	    if (LINK->child->ma == LINK->mother) {
	      if (!LINK->child->up)
		val = segfun(&LINK->child, first, second, LINK);
	    }
	    LINK->child = LINK->child->nextpa;
	  } while (val && LINK->child != NULL);
	  if (val) {
	    if (!sexlink) {
	      for (thishap = a;
		   (long)thishap <= (long)b;
		   thishap = (haplotype)((long)thishap + 1)) {
		for (thathap = a;
		     (long)thathap <= (long)b;
		     thathap = (haplotype)((long)thathap + 1)) {
		  here = genenumber[LINK->secondhap[(long)thishap - (long)a] - 1]
		    [LINK->firsthap[(long)thathap - (long)a] - 1];
		  gene[here - 1] = (gene[here - 1] || val);
		}
	      }
	    } else if ((*LINK->r)->male) {
	      for (thathap = a;
		   (long)thathap <= (long)b;
		   thathap = (haplotype)((long)thathap + 1)) {
		here = LINK->secondhap[(long)thathap - (long)a];
		gene[here - 1] = (gene[here - 1] || val);
	      }
	    } else {
	      for (thathap = a;
		   (long)thathap <= (long)b;
		   thathap = (haplotype)((long)thathap + 1)) {
		here = genenumber[LINK->secondhap[(long)thathap - (long)a] - 1]
		  [first - 1];
		gene[here - 1] = (gene[here - 1] || val);
	      }
	    }
	  }
	}
      }
    }
  }
  WITH = (*LINK->r)->gen;
  FORLIM = fgeno;
  for (first = 0; first < FORLIM; first++)
    WITH->genarray[first] = (WITH->genarray[first] && gene[first]);

/*compatibility test*/
  if (!incompat_traversal) {
    compat = false;
    FORLIM1 = fgeno;
    for (j = 0; j < FORLIM1; j++)
       compat = (compat || WITH->genarray[j]);
    if (!compat && (!detected_incompat[(*LINK->p)->id]) &&
	(!detected_incompat[(*LINK->q)->id])) {
      printf("\n One incompatibility involves the family in which person %ld is a child\n",(*LINK->r)->id);
      incompat_traversal = true;
      detected_incompat[(*LINK->p)->id] = true;
      detected_incompat[(*LINK->q)->id] = true;
      if (ONE_ERROR_ONLY)
	exit(EXIT_FAILURE);
      if (first_incompat)
	respond();
    }
  }
  cleanup(LINK->p);
  cleanup(LINK->q);
  LINK->child = LINK->father->foff;
  do {
    if (LINK->child->ma == LINK->mother) {
      if (!LINK->child->up)
	cleanup(&LINK->child);
    }
    LINK->child = LINK->child->nextpa;
  } while (LINK->child != NULL);   /*segdown*/
}


 void seg(p_, q_, r_, peel)
thisperson **p_, **q_, **r_;
direction peel;
{
  struct LOC_seg V;


  V.p = p_;
  V.q = q_;
  V.r = r_;
  if ((*V.p)->male) {
    V.father = *V.p;
    V.mother = *V.q;
  } else {
    V.father = *V.q;
    V.mother = *V.p;
  }
  if (peel == peelup)
    segup(&V);
  else
    segdown(&V);
}  /*seg*/


 void collapseup(p)
thisperson *p;
{
  thisperson *q, *child, *nextchild;
  boolean down;

  depth++;
  if (depth > (DEPTH_MULTIPLE * numind)) {
    printf("The next pedigree appears to have an unbroken loop\n");
    printf("If you do not think so, increase DEPTH_MULTIPLE in unknown.c\n");
    exit(EXIT_FAILURE);
  }
  p->done = true;
  if (p->foff == NULL) {
    depth--;
    return;
  }
  down = false;
  child = p->foff;
  while (child != NULL) {
    down = false;
    if (p->male)
      q = child->ma;
    else
      q = child->pa;
    if (!q->done) {
      collapsedown(q);
      nextchild = child;
      while (nextchild != NULL) {
	if (nextchild->pa == q || nextchild->ma == q) {
	  if (!nextchild->up)
	    collapseup(nextchild);
	  else
	    down = true;
	}
	if (p->male)
	  nextchild = nextchild->nextpa;
	else
	  nextchild = nextchild->nextma;
      }
      if (q->multi)
	collapseup(q);
      if (!down)
	seg(&p, &q, &child, peelup);
      else
	collapsedown(p);
    }
    if (p->male)
      child = child->nextpa;
    else
      child = child->nextma;
  }
  depth--;
}  /*collapseup*/


static void collapsedown(p)
thisperson *p;
{
  depth++;
  if (depth > (DEPTH_MULTIPLE * numind)) {
    printf("This pedigree appears to have an unbroken loop\n");
    printf("If you do not think so, increase DEPTH_MULTIPLE in unknown.c\n");
    exit(EXIT_FAILURE);
  }
  if (p->pa == NULL) {
    depth--;
    return;
  }
  p->up = true;
  collapseup(p->pa);
  seg(&p->pa, &p->ma, &p, peeldown);
  depth--;
}  /*collapsedown*/


static void likelihood(proband)
thisperson **proband;
{
  long i, j;
  thisperson *WITH;
  information *WITH1;
  thisarray *WITH2;
  locusvalues *WITH3;
  long FORLIM, FORLIM1;


  collapsedown(*proband);
  collapseup(*proband);
  if (!(*proband)->thisunknown[whichsys - 1])
    return;
  WITH = *proband;
  WITH1 = WITH->store;
  WITH2 = WITH->gen;
  WITH3 = thislocus[whichsys - 1];
  if (sexlink && WITH->male) {
    FORLIM = WITH3->nallele;
    for (j = 0; j < FORLIM; j++)
      WITH1->possible[whichsys - 1][0][j] = WITH2->genarray[j];
  } else {
    FORLIM = WITH3->nallele;
    for (i = 0; i < FORLIM; i++) {
      FORLIM1 = WITH3->nallele;
      for (j = i; j < FORLIM1; j++) {
	WITH1->possible[whichsys - 1][i]
	  [j] = WITH2->genarray[genenumber[i][j] - 1];
	WITH1->possible[whichsys - 1][j][i] = WITH1->possible[whichsys - 1][i]
	  [j];
      }
    }
  }
  cleanup(proband);
}  /*likelihood*/


static void iterpeds()
{
  long i, j;
  boolean compattest, compat;
  long FORLIM, FORLIM1;

  if (loop1 != NULL || loop2 != NULL)
    return;
  FORLIM = totperson;
  /* This means that this part of unknown is not active for pedigrees with loops! */
  for (i = 1; i <= FORLIM; i++)
    getgene(whichsys, person[i], person[i]->phen);
  first_incompat = false;
  compattest = false;
  compat = false;
  FORLIM = totperson;
  for (i = 1; i <= FORLIM; i++) {
    if (!compattest || person[i]->thisunknown[whichsys - 1]) {
      incompat_traversal=false;
      FORLIM1 = totperson;
      for (j = 1; j <= FORLIM1; j++)
	person[j]->done = false;
      FORLIM1 = totperson;
      for (j = 1; j <= FORLIM1; j++)
	person[j]->up = false;
      likelihood(&person[i]);
      if (!compattest) { /*Only do the overall compatibility test once per
                           locus-pedigree pair*/
	compattest = true;
	FORLIM1 = fgeno;
	for (j = 0; j < FORLIM1; j++)
	  compat = (compat || person[i]->gen->genarray[j]);
	if (!compat) {
	  printf("ERROR: Incompatibility detected in this family for locus %12ld\n",
		 whichsys);
	  respond();
	  first_incompat = true;
	}
      }
    }
  }
  FORLIM = totperson;
  for (i = 1; i <= FORLIM; i++) {
    if (person[i]->unknown)
      cleanup(&person[i]);
  }
}  /*iterpeds*/


static void reinit()
{
  long i, j, FORLIM, FORLIM1;

  FORLIM = totperson;
  for (i = 1; i <= FORLIM; i++) {
    FORLIM1 = nsystem;
    for (j = 0; j < FORLIM1; j++)
      Free(person[i]->phen[j]);
  }
  FORLIM = totperson;
  for (i = 1; i <= FORLIM; i++) {
    if (person[i]->store != NULL)
      Free(person[i]->store);
  }
  FORLIM = totperson;
  for (i = 1; i <= FORLIM; i++) {
    Free(person[i]->gen);
    Free(person[i]);
  }
  FORLIM = totperson;
  for (i = 1; i <= FORLIM; i++)
    person[i] = NULL;
}  /*reinit*/


 boolean testhets()
{
  /*Change: Function added by Joe Terwilliger 7/8/93*/
  double a_, b_, c;
  long prog, d, numl, lc, nall, sexl, nqv, i, j, sexd, int_;
  boolean tmp;
  char fff;

  tmp = true;
  fscanf(datafile, "%ld%lg%ld%ld%*[^\n]", &numl, &b_, &sexl, &prog);
  getc(datafile);
  fscanf(datafile, "%lg%lg%lg%ld%*[^\n]", &a_, &b_, &c, &d);
  getc(datafile);
  if (d == 1) {
    tmp = false;
    goto _L10;
  }
  if (prog != 1 && prog != 3)
    goto _L10;
  fscanf(datafile, "%lg%*[^\n]", &a_);
  getc(datafile);
  for (j = 1; j <= numl; j++) {
    fscanf(datafile, "%ld", &lc);
    switch (lc) {

    case 0:
      fscanf(datafile, "%ld%*[^\n]", &nall);
      getc(datafile);
      fscanf(datafile, "%*[^\n]");
      getc(datafile);
      fscanf(datafile, "%ld%*[^\n]", &nqv);
      getc(datafile);
      for (i = 1; i <= nqv; i++) {
	fscanf(datafile, "%lg%*[^\n]", &a_);
	getc(datafile);
      }
      fscanf(datafile, "%lg%*[^\n]", &a_);
      getc(datafile);
      fscanf(datafile, "%lg%*[^\n]", &a_);
      getc(datafile);
      break;

    case 1:
      fscanf(datafile, "%ld%*[^\n]", &nall);
      getc(datafile);
      fscanf(datafile, "%lg%*[^\n]", &a_);
      getc(datafile);
      fscanf(datafile, "%ld%*[^\n]", &nall);
      getc(datafile);
      if (sexl == 0) {
	for (i = 1; i <= nall; i++) {
	  fscanf(datafile, "%lg%*[^\n]", &a_);
	  getc(datafile);
	}
      } else {
	for (i = 1; i <= nall + nall; i++) {
	  fscanf(datafile, "%lg%*[^\n]", &a_);
	  getc(datafile);
	}
      }
      break;

    case 2:
      fscanf(datafile, "%ld%*[^\n]", &nall);
      getc(datafile);
      for (i = 1; i <= nall + 2; i++) {
	fscanf(datafile, "%lg%*[^\n]", &a_);
	getc(datafile);
      }
      break;

    case 3:
      fscanf(datafile, "%ld%*[^\n]", &nall);
      getc(datafile);
      fscanf(datafile, "%lg%*[^\n]", &a_);
      getc(datafile);
      break;
    }
  }
  fscanf(datafile, "%ld%ld%*[^\n]", &sexd, &int_);
  getc(datafile);
  if (sexd != 0) {
    fscanf(datafile, "%lg%*[^\n]", &a_);
    getc(datafile);
  }
  fscanf(datafile, "%lg%*[^\n]", &a_);
  getc(datafile);
  numl--;
  if (numl == 2 && int_ == 1)
    numl = 3;
  if (sexd == 1)
    numl++;
  if (sexd == 2)
    numl += numl;
  fscanf(datafile, "%lg%*[^\n]", &a_);
  getc(datafile);
  for (i = 1; i <= numl; i++)
    fscanf(datafile, "%ld", &d);
  fff = ' ';
  while (!P_eoln(datafile) && fff != '1') {
    fff = getc(datafile);
    if (fff == '\n')
      fff = ' ';
  }
  if (fff == '1')
    tmp = false;
_L10:
  return tmp;
}  /*testhets*/


static void initunknown()
{
  printf("Program UNKNOWN version %s\n", aVersion);
  printf("The following maximum values are in effect:\n");
  printf("%8ld loci\n", (long)maxlocus);
  printf("%8ld single locus genotypes\n", (long)maxgeno);
  printf("%8ld alleles at a single locus\n", (long)maxall);
  printf("%8ld individuals in one pedigree\n", (long)maxind);
  printf("%8ld marriage(s) for one male\n", (long)maxmarriage);
  printf("%8ld quantitative factor(s) at a single locus\n", (long)maxtrait);
  printf("%8ld liability classes\n", (long)maxliab);
  printf("%8ld binary codes at a single locus\n", (long)maxfact);
  one = 1.00001;   /*change*/

  printf("Opening DATAFILE.DAT\n");
  if (datafile != NULL)
    datafile = freopen("datafile.dat", "r", datafile);
  else
    datafile = fopen("datafile.dat", "r");
  if (datafile == NULL)
    exit(FileNotFound);
  if (P_eof(datafile)) {
    printf(
      "ERROR: File empty or nonexistent. Press <Enter> to continue or <Ctrl-C> to abort\n");
    scanf("%*[^\n]");
    getchar();
  }
  makehomozygous = testhets();   /*Change - 2 lines added*/
  rewind(datafile);

  printf("Opening PEDFILE.DAT\n");
  if (pedfile != NULL)
    pedfile = freopen("pedfile.dat", "r", pedfile);
  else
    pedfile = fopen("pedfile.dat", "r");
  if (pedfile == NULL)
    exit(FileNotFound);
  if (P_eof(pedfile)) {
    printf(
      "ERROR: File empty or nonexistent. Press <Enter> to continue or <Ctrl-C> to abort\n");
    scanf("%*[^\n]");
    getchar();
  }

  if (ipedfile != NULL)
    ipedfile = freopen("ipedfile.dat", "w", ipedfile);
  else
    ipedfile = fopen("ipedfile.dat", "w");
  if (ipedfile == NULL)
    exit(FileNotFound);
  if (speedfile != NULL)
    speedfile = freopen("speedfile.dat", "w", speedfile);
  else
    speedfile = fopen("speedfile.dat", "w");
  if (speedfile == NULL)
    exit(FileNotFound);

  /*changed from speedfile.dat*/
}  /*initunknown*/


int main(argc, argv)
int argc;
char *argv[];
{
  long FORLIM;
  locusvalues *WITH;
  int ind;

  ipedfile = NULL;
  pedfile = NULL;
  datafile = NULL;
  speedfile = NULL;
  initunknown();
  readloci();
  /*CLOSE(datafile);*/
  nsequence = 0;
  if (!P_eof(pedfile))
    fscanf(pedfile, "%ld", &newped);
  while (!P_eof(pedfile)) {
    numind = 0;
    depth = 0;
    readped();
    getunknown();
    FORLIM = nsystem;
    for (whichsys = 1; whichsys <= FORLIM; whichsys++) {
      if (mutsys != whichsys) {
        for(ind = 0; ind <=maxind; ind++)
	  detected_incompat[ind] = false;
	WITH = thislocus[whichsys - 1];
	fgeno = WITH->nallele * (WITH->nallele + 1) / 2;
	if (sexlink)
	  mgeno = WITH->nallele;
	else
	  mgeno = fgeno;
	getlocation(thislocus[whichsys - 1]);
	iterpeds();
      }
    }
    infer();
    writeped();
    writespeed();
    reinit();
  }
  /*CLOSE(pedfile);
  CLOSE(speedfile);
  CLOSE(ipedfile);*/
  if (speedfile != NULL)
    fclose(speedfile);
  if (datafile != NULL)
    fclose(datafile);
  if (pedfile != NULL)
    fclose(pedfile);
  if (ipedfile != NULL)
    fclose(ipedfile);
  exit(EXIT_SUCCESS);
}


