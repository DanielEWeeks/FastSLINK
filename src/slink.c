/* This file is part of a C version of SLINK a simulation program based
on the LINKAGE package. SLINK was originally described in: D. E.
Weeks, J. Ott, G. M. Lathrop (1990), SLINK: a general simulation
program for linkage analysis, American Journal of Human Genetics
47(1990), A204 (abstract). This version of SLINK is adapted to use the
faster pedigree traversal algorithms described in: R. W. Cottingham
Jr., R. M. Idury, A. A. Schaffer (1993), Faster sequential genetic
linkage computations, American Journal of Human Genetics 53(1993), pp.
252-263.*/
/* This file contains the main routines. See editing history below*/


#include <stdio.h>
/* include <p2c/p2c.h>*/
#include "commondefs.h"
#include "sldefs.h"
#include "moddefs.h"

/*29 Dec 1991*/

/*Modified to do general simulation by Daniel E. Weeks with the help of
  Mark Lathrop during May 29 to June 10, 1989.  Updated December 1989.
  Update History:
   5/18/90 - fixed the FUNCTION msegsex.
   5/90    - fixed clipping bug in pedigrees with loops.
   8/90    - Switched codes to be consistent with SIMLINK.
             Added counters for more feedback about input data.
             Corrected msegsex and checked for consistency with
             Version 5.1. Added checking of size of maxfact.
             Combined in features of HLINK, so now SLINK can
             simulate a proportion of unlinked families.
   10/91   - fixed clipping problem (see procedure collapseup)
    7/92   - corrected bug that required sequential numbering of
             families starting with 1.  Used bug fix suggested by
             L.A. Sandkuijl (Implemented Dec 6, 1993)
   12/93   - This file sliced out of entire program by A. A. Schaffer
    3/94   - A.A. Schaffer made syntactic changes in declarations
             to increase portability
   3/2007  - A.A. Schaffer adapted more changes from FASTLINK
             to remove dependence on maxhap and cleanup code
             created by p2c translation
   5/2008  - A.A. Schaffer made changes to allow maxall to be > 32
   8/2008  - A.A. Schaffer added the option -t to print genotypes at the 
             trait locus.
  12/2010  - A.A. Schaffer added more Free() calls to cleanup

  Input Files:
  simdata.dat <= regular LINKAGE data file.
  simped.dat  <= regular LINKAGE pedigree data file with an additional column
  containing the availability code.

  Output Files:
  simout.dat <= contains summary data regarding the simulation
  pedfile.dat <= regular LINKAGE pedigree data file containing the simulated
  data.
  traitfile.dat <= optional file to include trait genotypes; first 9 columns
             of this file should be identical to pedfile.dat

  NOTES:

  1) The last column in the pedigree file must contain the availability codes:
  Meaning
  Code     Markers       Trait
  0      Unavailable   Use original phenotypes
  1      Available     Use simulated phenotypes
  2      Available     Use original phenotypes
  3      Unavailable   Use simulated phenotypes
  So if you want to indicate that someone is unavailable at the trait
  locus, simply make their original phenotype unknown.  Remember that
  phenotypes are simulated CONDITIONAL on the phenotypes appearing in
  the simped.dat file.

  2) Parameters to be input interactively are:
  a) seed for random number generator.
  b) number of replicates desired.
  c) number identifying the trait locus. Use zero
  if there is no trait locus.

  3) This program reads from the simped (No speedfile used) and the
  simulated pedigrees are written to a file named 'pedfile.dat'.

  4) The stream file outputs have been removed, as have the speedfile inputs.

  5) For an affection status locus, the unaffected individuals are
  indicated by the constant integer assigned to unaff.

  6) For a quantitative locus, the simulation only supports ONE trait.
  Unaffected males at a sex-linked locus are indicated by the
  constant integer assigned to unafqu.
  */


/*FUNCTION clock:integer; */
/*SUN*/
/*Vax specific integer function;  returns milliseconds of CPU time
  spent on current process.  In Prospero Pascal, returns milliseconds
  on clock since midnight*/
/*VAR hour,minute,second,hundred:integer;*/
/*Must be 4-byte integers*/
/*BEGIN
  time(hour,minute,second,hundred);
  clock:=hour*3600000+minute*60000+second*1000+hundred*10;
 END;
*/

static void getdatetime(userseconds, seconds)
double *userseconds, *seconds;
{
  /*SUN*/
  struct timeval thistime;
  gettimeofday(&thistime,NULL);   /*SUN*/
  *userseconds = (thistime.tv_sec * 1000000.0 + thistime.tv_usec)/1000000.0;
  gettimeofday(&thistime, NULL);   /*SUN Pascal for elapsed number of seconds*/
  *seconds = (thistime.tv_sec * 1000000.0 + thistime.tv_usec)/1000000.0;
}  /* getdatetime */


double mapfunction(theta1, theta2)
double theta1, theta2;
{
  /*User defined function giving recombination between flanking markers as
    a function of recombination between adjacent markers*/
  return ((theta1 + theta2) / (1 + 4 * theta1 * theta2));
}  /* mapfunction */


double getdist(theta)
double *theta;
{
  if (*theta < 0.5)
    return (log(1.0 - 2.0 * *theta) / -2.0);
  else
    return 10.0;
}  /* getdist */


double invdist(dist)
double *dist;
{
  if (*dist < 10.0)
    return ((1.0 - exp(-2.0 * *dist)) / 2.0);
  else
    return 0.5;
}  /* invdist */


static double randomnum()
{
  /*Global variables used:  seed1,seed2, seed3*/
  double r;

  /*From FORTRAN code:
  THIS FUNCTION GENERATES INDEPENDENT UNIFORM DEVIATES ON
  THE INTERVAL (0.0,1.0).  SEE THE REFERENCE:  WICHMAN B.A.
  AND HILL I.D.(1982). ALGORITHM 183: AN EFFICIENT AND PORTABLE
  PSEUDO-RANDOM NUMBER GENERATOR. APPLIED STATISTICS 31;188-190.
  SEED1, SEED2, AND SEED3 SHOULD BE SET TO INTEGER VALUES
  BETWEEN 1 AND 30000 BEFORE THE FIRST ENTRY.  INTEGER
  ARITHMETIC UP TO 30323 IS NEEDED ON YOUR COMPUTER.

  IMPLICIT DOUBLE PRECISION(A-H,O-Z)
  INTEGER SEED1,SEED2,SEED3
  SAVE SEED2,SEED3
  DATA SEED2,SEED3/2321,18777/
  Note that all three seeds - SEED1, SEED2, and SEED3 - must
  be initialized at the beginning of the main program.
  */
  seed1 = seed1 % 177 * 171 - seed1 / 177 * 2;
  seed2 = seed2 % 176 * 172 - seed2 / 176 * 35;
  seed3 = seed3 % 178 * 170 - seed3 / 178 * 63;
  if (seed1 < 0)
    seed1 += 30269;
  if (seed2 < 0)
    seed2 += 30307;
  if (seed3 < 0)
    seed3 += 30323;
  r = seed1 / 30269.00 + seed2 / 30307.00 + seed3 / 30323.00;
  return (r - (long)r);
  /*  RANDOM:=MOD(R,1.000)  Note: a mod b = a - ((a div b)*b  */
}

/*generate a new seed based on the time of day*/
int generate_new_seed()
{
  int i,n;
  long secs; 
  time_t t1;
  int new_seed; /*value to return*/

  (void) time(&t1);
  srand48((long) t1);

  new_seed = 0;
  do {
    /* Choose new_seed based on the tick count */
    (void) time(&t1);
    secs = (long) t1;
    new_seed = secs % 30000;
    new_seed = ((long) (29999 * randomnum(new_seed,19000,20000))) % 30000;
  } while (new_seed <= 0 || new_seed >= 30000);
  return(new_seed);
}

/* randomnum */

static void selectgeneloc(p, thisped)
thisperson *p;
int thisped;
{
  double sum, unifrm, totprob;
  int loc, i, FORLIM;

  sum = 0.0;
  totprob = 0.0;
  FORLIM = numgen;
  for (i = 0; i < FORLIM; i++)
    totprob += geneprob[i];
  unifrm = randomnum();
  loc = 1;
  do {
    if (loc > numgen) {
      printf("ERROR: Unable to choose genotype for person %12d in pedigree %12d\n",
	     p->id, p->origped);
      fprintf(simout,
	"ERROR: Unable to choose genotype for person %12d in pedigree %12d\n",
	p->id, p->origped);
      printf("Press return to halt...\n");
      scanf("%*[^\n]");
      getchar();
      exit(0);   /*SUN*/
    }
    sum += geneprob[loc - 1] / totprob;
    loc++;
  } while (unifrm >= sum);
  p->geneloc = loc - 1;
  if (p->inloop > 0) {
    /* Need to assign to both looppersons if there is a loop */
    looppers[thisped - 1][p->inloop - 1][0]->geneloc = loc - 1;
    looppers[thisped - 1][p->inloop - 1][1]->geneloc = loc - 1;
  }
  /* QUESTION: Should something be done to holdpoint?  Probably not,
  since holdpoint simply points to the genarray of another */
  /* Now we must zero out the geneprob array for use by the next person */
  for (i = 0; i < maxfemgen; i++)
    geneprob[i] = 0.0;
}


/* selectgeneloc */


static void writeunk(isys, ii)
int isys, ii;
{
  /*Writes unknown genotypes*/
  int k;
  thisperson *WITH;
  locusvalues *WITH1;
  phenotype *WITH2;
  int FORLIM;

  WITH = person[ii];
  WITH1 = thislocus[isys - 1];
  WITH2 = WITH->phen[isys - 1];
  switch (thislocus[isys - 1]->which) {

  case binary_:
    if (WITH1->format == 3)
      if (WITH1->nallele < 100) {
	fprintf(pedfile, "  0  0");
	if (show_trait_genotypes)
	  fprintf(traitfile, "  0  0");
      }
      else {
	fprintf(pedfile, "    0    0");
	if (show_trait_genotypes)
	  fprintf(traitfile, "    0    0");
      }
    else {
      FORLIM = WITH1->UU.U2.numfact;
      for (k = 1; k <= FORLIM; k++) {
	fprintf(pedfile, " 0");
	if (show_trait_genotypes)
	  fprintf(traitfile, " 0");
      }
    }
    break;

  case affection:
    fprintf(pedfile, "%3d", missaff);
    if (show_trait_genotypes)
      fprintf(traitfile, "%3d", missaff);
    if (WITH1->UU.U0.nclass != 1) {
      fprintf(pedfile, "  1");
      if (show_trait_genotypes)
	fprintf(traitfile, "  1");
    }
    if (show_trait_genotypes) {
      fprintf(traitfile, " 0      0");
    }
    break;

  /*Zero liability class is not allowed at an affection locus*/
  case quantitative:
    FORLIM = WITH1->UU.U1.ntrait;
    for (k = 1; k <= FORLIM; k++) {
      fprintf(pedfile, " %9.4f", missval);
      if (show_trait_genotypes)
	fprintf(traitfile, " %9.4f", missval);
    }
    if (show_trait_genotypes) {
      fprintf(traitfile, " 0      0");
    }
    break;
  }
}


/* writeunk */
/*writeunk*/


static void writeorg(isys, ii)
int isys, ii;
{
  /*Echos out original phenotype from phenf found in the pedigree file*/
  int a, b, k;
  thisperson *WITH;
  locusvalues *WITH1;
  phenotype *WITH2;
  int FORLIM;
  int all1, all2;

  WITH = person[ii];
  WITH1 = thislocus[isys - 1];
  WITH2 = WITH->phen[isys - 1];
  /*For the trait locus we wish to write out the original 'phenotypic' data,
  even if the person is unavailable: */
  switch (thislocus[isys - 1]->which) {

  case binary_:  /*this is the case for markers*/
    if (WITH1->format == 3) {  /*This locus is an allele numbers locus */
      a = 0;
      b = 0;
      if (WITH2->alleles[0] <= WITH2->alleles[1]) {
	a = WITH2->alleles[0];
	b = WITH2->alleles[1];
      }
      if (b == 0)
	b = a;
      if (WITH1->nallele < 100) {
	fprintf(pedfile, "%3d%3d", a, b);
	if (show_trait_genotypes)
	  fprintf(traitfile, "%3d%3d", a, b);
      }
      else {
	fprintf(pedfile, "%5d%5d", a, b);
	if (show_trait_genotypes)
	  fprintf(traitfile, "%5d%5d", a, b);
      }
    } else {
      FORLIM = WITH1->UU.U2.numfact;
      for (k = 1; k <= FORLIM; k++) {
	if ((unsigned int)k < 32 && ((1 << k) & WITH2->phenf) != 0) {
	  fprintf(pedfile, " 1");
	  if (show_trait_genotypes) 
	    fprintf(traitfile, " 1");
	}
	else {
	  fprintf(pedfile, " 0");
	  if (show_trait_genotypes) 
	    fprintf(traitfile, " 0");
	}
      }
    }
    break;

  case affection:
    fprintf(pedfile, "%3d", WITH2->aff);
    if (show_trait_genotypes) {
      fprintf(traitfile, "%3d", WITH2->aff);
    }
    if (WITH1->UU.U0.nclass != 1) {
      fprintf(pedfile, "%3d", WITH2->liability);
      if (show_trait_genotypes) {
	fprintf(traitfile, "%3d", WITH2->liability);
      }
    }
    if (show_trait_genotypes) {
      all1 = WITH2->alleles[0];
      all2 = WITH2->alleles[1];
      fprintf(traitfile, " %7d%7d",all1, all2);
    }
    break;

  case quantitative:
    if (!sexlink || !WITH->male) {
      FORLIM = WITH1->UU.U1.ntrait;
      for (k = 0; k < FORLIM; k++) {
	fprintf(pedfile, " %9.4f", WITH2->x[k]);
	if (show_trait_genotypes) {
	  fprintf(traitfile, " %9.4f", WITH2->x[k]);
	}
      }
    } else {
      FORLIM = WITH1->UU.U1.ntrait;
      for (k = 1; k <= FORLIM; k++) {
	fprintf(pedfile, " %9d", WITH2->aff);
	if (show_trait_genotypes) {
	  fprintf(traitfile, " %9d", WITH2->aff);
	}
      }
    }
    if (show_trait_genotypes) {
      all1 = WITH2->alleles[0];
      all2 = WITH2->alleles[1];
      fprintf(traitfile, " %7d%7d",all1, all2);
    }
    break;
  }
  /*Note: for sexlinked, this is not exactly the original,
  but it is equivalent.  i.e. this is what unknown does.*/
}


/* writeorg */
/*writeorig*/


static void writesim(isys, ii, all1, all2)
int isys, ii, all1, all2;
{
  /*all1 = allele number 1 of the simulated genotype;
all2 = allele number 2 of the simulated genotype;
Writes out the simulated genotypes/phenotypes.*/
  int ranaff, k;
  double norm1, ran1, ran2, ranobs;
  thisperson *WITH;
  locusvalues *WITH1;
  phenotype *WITH2;
  int FORLIM;

  WITH = person[ii];
  WITH1 = thislocus[isys - 1];
  WITH2 = WITH->phen[isys - 1];
  switch (thislocus[isys - 1]->which) {

  case binary_:
    if (WITH1->format == 3)   /*This locus is an allele numbers locus */
      if (WITH1->nallele < 100) {
	fprintf(pedfile, "%3d%3d", all1, all2);
	if (show_trait_genotypes) 
	  fprintf(traitfile, "%3d%3d", all1, all2);
      }
      else {
	fprintf(pedfile, "%5d%5d", all1, all2);
	if (show_trait_genotypes)  
	  fprintf(traitfile, "%5d%5d", all1, all2);
      }
    else {
      FORLIM = WITH1->UU.U2.numfact;
      for (k = 1; k <= FORLIM; k++) {
	if (((unsigned int)k < 32 &&
	     ((1 << k) & WITH1->UU.U2.bin_allele[all1 - 1]) != 0) ||
	    ((unsigned int)k < 32 &&
	     ((1 << k) & WITH1->UU.U2.bin_allele[all2 - 1]) != 0)) {
	  fprintf(pedfile, " 1");
	  if (show_trait_genotypes) {
	    fprintf(traitfile, " 1");
	  }
	}
	else {
	  fprintf(pedfile, " 0");
	  if (show_trait_genotypes) {
	    fprintf(traitfile, " 0");
	  }
	}
      }
    }
    break;

  case affection:
    /*Must choose the phenotype at random based on the penetrance*/
    if (randomnum() < WITH1->UU.U0.pen[all1][all2 - 1][2]
	[WITH2->liability - 1])
	  /*Unaffected*/
	    ranaff = affval;   /*Affected => use affected value*/
    else
      ranaff = unaff;
    fprintf(pedfile, "%3d", ranaff);
    if (show_trait_genotypes) {
      fprintf(traitfile, "%3d", ranaff);
    }
    if (WITH1->UU.U0.nclass != 1) {
      fprintf(pedfile, "%3d", WITH2->liability);
      if (show_trait_genotypes) {
	fprintf(traitfile, "%3d", WITH2->liability);
      }
    }
    if (show_trait_genotypes) {
      fprintf(traitfile, " %7d%7d", all1, all2);
    }
    break;

  case quantitative:
    if (!sexlink || !WITH->male) {
      /*Must choose the phenotype at random based on the penetrance*/
      /*Use Box-Muller transformation to obtain an N(0,1) r.v.*/
      ran1 = randomnum();
      ran2 = randomnum();
      norm1 = sqrt(-2.0 * log(ran1)) * cos(6.283185308 * ran2);
      if (all1 != all2)
	ranobs = WITH1->UU.U1.pm[all1][all2 - 1]
		 [0] + sqrt(1.0 / WITH1->UU.U1.conmat) *
		       sqrt(1.0 / WITH1->UU.U1.vmat[0][0]) * norm1;
      else
	ranobs = WITH1->UU.U1.pm[all1][all2 - 1]
		 [0] + sqrt(1.0 / WITH1->UU.U1.vmat[0][0]) * norm1;
      fprintf(pedfile, " %9.4f", ranobs);
      if (show_trait_genotypes) {
	fprintf(traitfile, " %9.4f", ranobs);
      }
    } 
    else {
      if (all1 == affall || all2 == affall) {
	fprintf(pedfile, " %9d", affall);
	if (show_trait_genotypes) {
	  fprintf(traitfile, " %9d", affall);
	}
      }
      else {
	fprintf(pedfile, " %9d", unafqu);
	if (show_trait_genotypes) {
	  fprintf(traitfile, " %9d", unafqu);
	}
      }
    }
    if (show_trait_genotypes)
      fprintf(traitfile, " %7d%7d", all1, all2);
    break;
  }/*quantitative*/

  /*Sexlinked and Male => Either affected or not affected*/
}


/* writesim */


static void writelinkage(ii)
int ii;
{  /* writelinkage - modified from Lodewijk Sandkuyl's code */
  int panr, manr, foffnr, nextpanr, nextmanr, isys;
  thisperson *WITH;
  int FORLIM;

  WITH = person[ii];
  /*  FOR ii:=1 TO totperson DO */
  if (WITH->pa == NULL)
    panr = 0;
  else
    panr = WITH->pa->origid;
  if (WITH->ma == NULL)
    manr = 0;
  else
    manr = WITH->ma->origid;
  if (WITH->foff == NULL)
    foffnr = 0;
  else
    foffnr = WITH->foff->origid;
  if (WITH->nextpa == NULL)
    nextpanr = 0;
  else
    nextpanr = WITH->nextpa->origid;
  if (WITH->nextma == NULL)
    nextmanr = 0;
  else
    nextmanr = WITH->nextma->origid;
  if (totall < 1000) {
    fprintf(pedfile, "%4d%4d%4d%4d%4d%4d%4d",
	  totall, WITH->origid, panr, manr, foffnr, nextpanr, nextmanr);
    if (show_trait_genotypes)
      fprintf(traitfile, "%4d%4d%4d%4d%4d%4d%4d",
	      totall, WITH->origid, panr, manr, foffnr, nextpanr, nextmanr);
  }
  else {
    fprintf(pedfile, "%6d%4d%4d%4d%4d%4d%4d",
	  totall, WITH->origid, panr, manr, foffnr, nextpanr, nextmanr);
    if (show_trait_genotypes)
      fprintf(traitfile, "%6d%4d%4d%4d%4d%4d%4d",
	      totall, WITH->origid, panr, manr, foffnr, nextpanr, nextmanr);
  }
  if (ii == totperson)
    totall++;
  else if (WITH->ped != person[ii + 1]->ped)
    totall++;
  /*Output ped ids must be sequential*/
  if (WITH->male)
    fprintf(pedfile, " 1%3d", WITH->prbnum);
  else
    fprintf(pedfile, " 2%3d", WITH->prbnum);
  if (show_trait_genotypes) {
    if (WITH->male)
      fprintf(traitfile, " 1%3d", WITH->prbnum);
    else
      fprintf(traitfile, " 2%3d", WITH->prbnum);
  }
  /* Since hapstore has been set up in terms of the order of the loci
  within each pedigree, then we must use the correct order for
  printing out the marker data */
  /* BUT hapstore used to be printed out WITHOUT using ORDER[]???? */
  /* In the internal order of the loci, they are numbered 1,2,3...  However,
  if we want to print out in the external order, then we print
  the locus #(order[i]) in position i.  For example, if the external
  order is 2 1 3, and the internal order is 1 2 3, then we want
  to print locus #2 in position 1, locus #1 in position 2, and locus
  #3 in position 3. i.e. internal locus #2 is external locus #1. */
  switch (WITH->availnum) {

  case 0:
    FORLIM = nlocus;
    for (isys = 1; isys <= FORLIM; isys++) {
      if (isys != trait)
	writeunk(order[isys - 1], ii);
      else
	writeorg(order[isys - 1], ii);
    }
    break;

  case 2:
    FORLIM = nlocus;
    for (isys = 0; isys < FORLIM; isys++) {
      if (isys + 1 != trait)
	writesim(order[isys], ii,
	  hapstore[order[isys] - 1][invgenenum1[WITH->geneloc - 1] - 1],
	  hapstore[order[isys] - 1][invgenenum2[WITH->geneloc - 1] - 1]);
      else {
	if (show_trait_genotypes) {
          person[ii]->phen[order[isys] -1]->alleles[0] =
	    hapstore[order[isys] - 1][invgenenum1[WITH->geneloc - 1] - 1];
          person[ii]->phen[order[isys] -1]->alleles[1] =
	    hapstore[order[isys] - 1][invgenenum2[WITH->geneloc - 1] - 1];
	}
	writeorg(order[isys], ii);
      }
    }
    break;

  case 1:
    FORLIM = nlocus;
    for (isys = 0; isys < FORLIM; isys++)
      writesim(order[isys], ii,
	hapstore[order[isys] - 1][invgenenum1[WITH->geneloc - 1] - 1],
	hapstore[order[isys] - 1][invgenenum2[WITH->geneloc - 1] - 1]);
    break;

  case 3:
    FORLIM = nlocus;
    for (isys = 1; isys <= FORLIM; isys++) {
      if (isys != trait)
	writeunk(order[isys - 1], ii);
      else
	writesim(order[isys - 1], ii, hapstore[order[isys - 1] - 1]
				      [invgenenum1[WITH->geneloc - 1] - 1],
		 hapstore[order[isys - 1] - 1]
		 [invgenenum2[WITH->geneloc - 1] - 1]);
    }
    break;
  }
  fprintf(pedfile, "%2d", WITH->input_availnum);
  if (show_trait_genotypes) {
    fprintf(traitfile, "%2d", WITH->input_availnum);
  }
  if (linked)
    fprintf(pedfile, " Linked   ");
  else
    fprintf(pedfile, " Unlinked ");
  if (show_trait_genotypes) {
    if (linked)
      fprintf(traitfile, " Linked   ");
    else
      fprintf(traitfile, " Unlinked ");
  }
  /*      write(pedfile, ' :', geneloc: 7);
  FOR isys:=1 TO nlocus DO
  write(pedfile, ' | ', hapstore[order[isys], invgenenum1[geneloc]]: 2, '/',
  hapstore[order[isys], invgenenum2[geneloc]]: 2); */
  switch (WITH->availnum) {   /*WITH statement*/

  case 0:
    fprintf(pedfile, " Mk Unkno");
    if (show_trait_genotypes) {
      fprintf(traitfile, " Mk Unkno");
    }
    if (trait != 0) {
      fprintf(pedfile, "; Tr orig\n");
      if (show_trait_genotypes) {
	fprintf(traitfile, "; Tr orig\n");
      }
    }
    else {
      putc('\n', pedfile);
      if (show_trait_genotypes) {
	putc('\n', traitfile);
      }
    }
    break;

  case 2:
    fprintf(pedfile, " Mk Avail");
    if (show_trait_genotypes) {
      fprintf(traitfile, " Mk Avail");
    }
    if (trait != 0) {
      fprintf(pedfile, "; Tr orig\n");
      if (show_trait_genotypes) {
	fprintf(traitfile, "; Tr orig\n");
      }
    }
    else {
      putc('\n', pedfile);
      if (show_trait_genotypes) {
	putc('\n', traitfile);
      }
    }
    break;

  case 1:
    fprintf(pedfile, " Mk Avail");
    if (show_trait_genotypes) {
      fprintf(traitfile, " Mk Avail");
    }
    if (trait != 0) {
      fprintf(pedfile, "; Tr simu\n");
      if (show_trait_genotypes) {
	fprintf(traitfile, "; Tr simu\n");
      }
    }
    else {
      putc('\n', pedfile);
      if (show_trait_genotypes) {
	putc('\n', traitfile);
      }
    }
    break;

  case 3:
    fprintf(pedfile, " Mk Unkno");
    if (show_trait_genotypes) {
      fprintf(traitfile, " Mk Unkno");
    }
    if (trait != 0) {
      fprintf(pedfile, "; Tr simu\n");
      if (show_trait_genotypes) {
	fprintf(traitfile, "; Tr simu\n");
      }
    }
    else {
      putc('\n', pedfile);
      if (show_trait_genotypes) {
	putc('\n', traitfile);
      }
    }
    break;
  }
}

/*writelinkage*/


void invertmat(m, n, det)
double (*m)[maxtrait];
int n;
double *det;
{
  covmatrix v;
  double val;
  int i, j, k;

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
}  /* invertmat */


/*static void collapsedown (thisperson *p, struct LOC_likelihood *LINK);*/

static void collapsedown();

void cleanup(p, LINK)
thisperson **p;
struct LOC_likelihood *LINK;
{
  thisperson *WITH;

  WITH = *p;
  if (feedback)
    printf("deleting: %12d\n", WITH->id);
  if (NULL != WITH->gen->genarray) {
    Free(WITH->gen->genarray);
  }
  if (NULL != WITH->gen->sparseflag) {
    Free(WITH->gen->sparseflag);
  }
  Free(WITH->gen);
  /*      DisposPtr(Ptr(gen));*/
  WITH->gen = NULL;
}  /* cleanup */

/*This routine is a simple exit routine, when the sytem is out of memory*/
void malloc_err(message)
char * message;
{

  fprintf(stderr, "\nProblem with malloc, probably not enough space\n");
  fprintf(stderr, "Problem occurred when allocating %s\n", message);
  exit(EXIT_FAILURE);
}

void allocate_thisarray(location, number)
thisarray *location;
int number;
{
   location->genarray = (double *) malloc(number * sizeof(double));
   if (NULL == location->genarray) 
     malloc_err("genarray field");
   location->sparseflag = (unsigned char *) malloc(number * 
                                                     sizeof(unsigned char));
   if (NULL == location->sparseflag) 
     malloc_err("sparseflag field");
}


static void prob(p, LINK)
thisperson **p;
struct LOC_likelihood *LINK;
{
  int i;
  thisperson *WITH;
  unsigned char *tempflag1;
  double *tempwith1; /*R. M. Idury*/

  WITH = *p;
  if (WITH->gen != NULL)
    return;
  WITH->gen = (thisarray *)Malloc(sizeof(thisarray));
  if (WITH->gen == NULL)
    malloc_err("gen field in prob");
  allocate_thisarray(WITH->gen, fgeno);
  /*       gen:=genpoint(NewPtr(SizeOf(thisarray)));*/
  tempwith1 = WITH->gen->genarray;
  tempflag1 = WITH->gen->sparseflag;
  if (feedback)
    printf("making: %12d\n", WITH->id);
  for (i = 0; i < fgeno; i++) {
    tempflag1[i] = 0;
    tempwith1[i] = 0.0;
  }
  if (WITH->inloop > 0) {
    if (looppers[LINK->thisped - 1][WITH->inloop - 1][0] == *p) {
      tempwith1[LINK->loopgen[WITH->inloop - 1] - 1] =
       LINK->holdpoint[WITH->inloop - 1]->
	 genarray[LINK->loopgen[WITH->inloop - 1] - 1];
      if (tempwith1[LINK->loopgen[WITH->inloop - 1] - 1] != 0.0)
	tempflag1[LINK->loopgen[WITH->inloop - 1] - 1] = 1;
    }
    else {
	tempwith1[LINK->loopgen[WITH->inloop - 1] - 1] = 1.0;
	tempflag1[LINK->loopgen[WITH->inloop - 1] - 1] = 1;
      }
    if ((*p)->pa == NULL)
      LINK->nuscale++; /*Alert different scaling from 5.1*/
    return;
  }
  getvect(*p, LINK);
  if ((*p)->pa == NULL)
    LINK->nuscale++;

  /*Changed on basis of ver. 5.03*/
}  /* prob */

Local Void getgenetables()
{
  int locfstart, locfend, locfseg;
  int f1, f2, i, first;

  currentfence = 1;
  for (first = 0; first < fgeno; first++) {
    base[first] = currentfence;
    locfstart = probstart[first];
    locfend = probend[first];
    locfseg = segstart[first];
    for (i = locfstart-1; i < locfend; i++) {
      f1 = invgenenum1[locfseg - 1];
      haps1[currentfence] = f1;
      f2 = invgenenum2[locfseg - 1];
      haps2[currentfence] = f2;
      hind[currentfence++] = i;
      locfseg++;
    }
    fence[first] = currentfence;
  } 
}  /* getgenetables */


/* getgeneindices is a preprocessing routine that compute the encodings
   of all the haplotypes and genotypes and precomputes for each
   genotype a list of the haplotypes it can pass on in one array
   insted of using indirect addressing by genotype */
Void getgeneindices()
{
  int i;

#if (!defined(LESSMEMORY))
  segprob2 = (double *) malloc(nuneed * nuneed * sizeof(double));
  if (segprob2 == NULL) {
    fprintf(stderr, "\nYou probably need to use the slower version of this program");
    malloc_err("segprob2");
  }
#endif

  /*find size of isozygote class and it's square for joint classes*/
  maxclasssize = 2;
  for (i = 3; i<= mlocus; i++)
    maxclasssize *= 2;
  maxisozygclass = maxclasssize * maxclasssize;

  nuprobclass = 2;
  for(i = 2; i <= mlocus; i++)
    nuprobclass = 3 * nuprobclass - 1;

  segval = (double *) malloc(maxisozygclass * sizeof(double));
  if (segval == NULL)
    malloc_err("segval");
  tempseg = (double *) malloc(maxchild *maxisozygclass * sizeof(double));
  if (tempseg == NULL)
    malloc_err("tempseg");
  segindex = (unsigned *) malloc(4 * maxisozygclass * sizeof(unsigned));
  if (segindex == NULL)
    malloc_err("segindex");
#if (defined(LESSMEMORY))
  tempseg2 = (double *) malloc(maxchild *maxisozygclass * sizeof(double));
  if (tempseg2 == NULL)
    malloc_err("tempseg2");
#endif
}  /* getgeneindices */


void initseg(LINK)
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

void exitseg(LINK)
struct LOC_seg *LINK;
{
  LINK->LINK->nuscale++;
  LINK->child = LINK->father->foff;
  do {
    if (LINK->child->ma == LINK->mother && !LINK->child->up)
      cleanup(&LINK->child, LINK->LINK);
    LINK->child = LINK->child->nextpa;
  } while (LINK->child != NULL);   /* exitseg */
}

/*computenumhaps computes the weighted sum over all genotypes of
how many haplotypes can be passed on to a child.
The following explanation assumes some basic understanding of
combinatorics at the undergraduate level. The total number of
haplotypes that can be passed on is
a weighted sum of the number of genotypes, where
different genotypes get different weights,  We first
classify the genotypes by heterozygosity pattern, where 0
means homozygous and 1 means heterozygous. There are mlocus loci
and therefore 2**mlocus different heterozygosity patterns.
The number of genotypes of a given heterozygosity pattern can be computed
by considering one locus at a time and multiplying the contributions
of all the loci, since they are independent.
A homozygous locus with A alleles contributes a factor of A.
The first heterozygous locus contributes a factor of A *(A-1)/2 if
it has A alleles.
Any other heterozygous locus with B alleles, contributes a factor of B*(B-1).

The weight of a genotype is 1 if it is homozygous; otherwise it
is 2**(h-1), where h is the number of heterozygous loci in the genotype.*/

static int computenumhaps()
{
  int numpatterns;
  hapvector *patternmatrix;
  int p; /*index over heterozygosity patterns*/
  int l; /*index over loci*/
  int value; /*accumulator for return value*/
  int term; /*term for this heterozygosity pattern*/
  boolean hetyet; /*have we seen a heterozygous locus yet*/ 

  numpatterns = 2;
  for (l=2; l<=mlocus;l++)
  numpatterns*=2;
  patternmatrix = (hapvector *) malloc(numpatterns*sizeof(hapvector));
  if (patternmatrix == NULL)
    malloc_err("patternmatrix");
  for(l=0; l < mlocus; l++)
    patternmatrix[0][l]=0;
  for(p=1; p< numpatterns; p++) {
    l = 0;
    do{
      patternmatrix[p][l] = !patternmatrix[p-1][l];
      l++;
    } while(!patternmatrix[p][l-1]);
    do{
      patternmatrix[p][l] = patternmatrix[p-1][l];
      l++;
    } while(l < mlocus);
  }

  value = 2; /*offset by 1 because 0 index has special meaning*/
  for(p=0; p < numpatterns; p++){
    term = 1;
    hetyet = 0;
    for(l=0; l < mlocus; l++){
      if (!patternmatrix[p][l]) {
        term *= thislocus[l]->nallele;
      }
      else 
        if (!hetyet) {
          hetyet = 1;
          term *= (thislocus[l]->nallele * (thislocus[l]->nallele -1))/2;
        }
        else 
          term *= thislocus[l]->nallele * (thislocus[l]->nallele - 1) * 2;
    }      
    value+=term;
  }
  free(patternmatrix);
  return(value);
}

Void allocategenetables()
{
  int maxfgeno, maxfinalfence;
  int i;

  maxfemgen = 0;
  maxhaplo = 0;

  maxfgeno = fgeno;
  maxhaplo = nuhap;
  maxfinalfence = computenumhaps();

  maxfemgen = maxfgeno;

#ifndef LESSMEMORY
  indpool = (unsigned *) malloc(2 * maxclasssize * maxfgeno * sizeof (unsigned));
  if (indpool == NULL) {
    fprintf(stderr, "\nYou probably need to use the slower version of this program");
    malloc_err("indpool");
  }
  invpool = (unsigned *) malloc(2 * maxclasssize * maxfgeno * sizeof (unsigned));
  if (invpool == NULL) {
    fprintf(stderr, "\nYou probably need to use the slower version of this program");
    malloc_err("invpool");
  }
  nextpool = (unsigned *) malloc(2 * maxclasssize * maxfgeno * sizeof (unsigned));
  if (nextpool == NULL) {
    fprintf(stderr, "\nYou probably need to use the slower version of this program");
    malloc_err("nextpool");
  }
#endif

  genenumber = (int**) malloc(maxhaplo * sizeof(int*));
  if (genenumber == NULL)
    malloc_err("genenumber");
  for(i = 0; i < maxhaplo; i++) {
    genenumber[i] = (int *) malloc(maxhaplo * sizeof(int));
    if (genenumber[i] == NULL)
      malloc_err("Entry in genenumber");
  }

  if (approximate) {
    approxarray = (unschar **) malloc(nuped * sizeof(unschar *));
    if (approxarray == NULL)
      malloc_err("approxarray");
    for(i = 0; i <= nuped; i++) {
      approxarray[i] = (unschar *) malloc(maxfgeno * sizeof(unschar));
      if (approxarray[i] == NULL)
        malloc_err("Entry in approxarray");
    }
  }

  base = (unsigned *) malloc(maxfgeno * sizeof(unsigned));
  if (base == NULL)
    malloc_err("base");
  fence = (unsigned *) malloc(maxfgeno * sizeof(unsigned));
  if (fence == NULL)
    malloc_err("fence");
  invgenenum1 = (unsigned *) malloc(maxfgeno * sizeof(unsigned));
  if (invgenenum1 == NULL)
    malloc_err("invgenenum1");
  invgenenum2 = (unsigned *) malloc(maxfgeno * sizeof(unsigned));
  if (invgenenum2 == NULL)
    malloc_err("invgenenum2");
  segstart = (int *) malloc(maxfgeno * sizeof(int));
  if (segstart == NULL)
    malloc_err("segstart");
  probstart = (unsigned *) malloc(maxfgeno * sizeof(unsigned));
  if (probstart == NULL)
    malloc_err("probstart");
  probend = (unsigned *) malloc(maxfgeno * sizeof(unsigned));
  if (probend == NULL)
    malloc_err("probend");
  rare = (boolean *) malloc(maxfgeno * sizeof(unsigned));
  if (rare == NULL)
    malloc_err("rare");
  if (risk) {
    risk1 = (boolean *) malloc(maxfgeno * sizeof(boolean));
    if (risk1 == NULL)
      malloc_err("risk1");
    risk2 = (boolean *) malloc(maxfgeno * sizeof(boolean));
    if (risk2 == NULL)
      malloc_err("risk2");
    riskmale = (boolean *) malloc(maxhaplo * sizeof(boolean));
    if (riskmale == NULL)
      malloc_err("riskmale");
  }
  if (mutsys != 0) {
    muthap = (unschar *) malloc(maxhaplo * sizeof(unschar));
    if (muthap == NULL)
      malloc_err("muthap");
  }  

  haps1 = (unsigned short *) malloc(maxfinalfence * sizeof(unsigned short));
  if (haps1 == NULL)
    malloc_err("haps1");
  haps2 = (unsigned short *) malloc(maxfinalfence * sizeof(unsigned short));
  if (haps2 == NULL)
    malloc_err("haps2");
  hind = (unsigned int *) malloc(maxfinalfence * sizeof(unsigned int));
  if (hind == NULL)
    malloc_err("hind");
  nonzgens = (unsigned int *) malloc(maxfgeno * sizeof(unsigned int));
  if (nonzgens == NULL)
    malloc_err("nonzgens");
  gene = (double *) malloc(maxfgeno * sizeof(double));
  if (gene == NULL)
   malloc_err("gene");
  flag = (boolean *) malloc(maxfgeno * sizeof(boolean));
  if (flag == NULL)
    malloc_err("flag");
  psumcache = (double *) malloc(maxhaplo * sizeof(double));
  if (psumcache == NULL)
    malloc_err("psumcache");
  qsumcache = (double *) malloc(maxhaplo * sizeof(double));
  if (qsumcache == NULL)
    malloc_err("qsumcache");
#if !defined(LESSMEMORY)
  phapcache1 = (cache *) malloc(maxhaplo * sizeof(cache));
  if (phapcache1 == NULL)
    malloc_err("phapcache1");
#endif

}  /* allocategenetables */



Void freegenetables()
{
  int i;

#ifndef LESSMEMORY
   free(indpool);
   indpool = NULL;

   free(invpool);
   invpool = NULL;

   free(nextpool);
   nextpool = NULL;
#endif

  for(i = 0; i < maxhaplo; i++) {
    free(genenumber[i]);
    genenumber[i] = NULL;
  }
  free(genenumber);
  genenumber = NULL;

  if (approximate) {
    for(i = 0; i <= nuped; i++) {
      free(approxarray[i]);
      approxarray[i] = NULL;
    }
    free(approxarray);
    approxarray = NULL;
  }

  free(base);
  base = NULL;
  free(fence);
  fence = NULL;
  free(invgenenum1);
  invgenenum1 = NULL;
  free(invgenenum2);
  invgenenum2 = NULL;
  free(segstart);
  segstart = NULL;
  free(probstart);
  probstart = NULL;
  free(probend);
  probend = NULL;
  free(rare);
  rare = NULL;
  if (risk) {
    free(risk1);
    risk1 = NULL;
    free(risk2);
    risk2 = NULL;
    free(riskmale);
    riskmale = NULL;
  }
  if (mutsys != 0) {
    free(muthap);
    muthap = NULL;
  }  

  free(haps1);
  haps1 = NULL;
  free(haps2);
  haps2 = NULL;
  free(hind);
  hind = NULL;

  free(nonzgens);
  nonzgens = NULL;

  free(gene);
  gene = NULL;
  free(flag);
  flag = NULL;

  free(psumcache);
  psumcache = NULL;
  free(qsumcache);
  qsumcache = NULL;
#if !defined(LESSMEMORY)
  free(phapcache1);
  phapcache1 = NULL;
#endif
}  /* freegenetables */


static void seg(p_, q_, r_, peel, LINK)
thisperson **p_, **q_, **r_;
direction peel;
struct LOC_likelihood *LINK;
{
  struct LOC_seg V;

  V.LINK = LINK;
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
  if (peel == peelup) {
    switch (thispath) {

    case auto_:
      segup(&V);
      break;

    case mauto:
      oldsegup(&V);
      break;

    case sex:
      segsexup(&V);
      break;

    case msex:
      oldsegsexup(&V);
      break;
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
}  /* seg */

static void collapseup(p, LINK)
thisperson *p;
struct LOC_likelihood *LINK;
{
  thisperson *q, *child, *nextchild;
  boolean down;

  p->done = true;
  /*10/1/91 Correction to clipping problem */
  /* Step 1: Find the spouse */
  if (p->foff != NULL) {
    if (p->male)
      q = p->foff->ma;
    else
      q = p->foff->pa;
  }
  /* Step 2: Collapseup if there are children AND either spouse */
  /*         has not yet been assigned a genotype */
  if (p->foff == NULL)
    return;
  if (p->geneloc != 0 && q->geneloc != 0 && noloop)
    return;
  /* Old code did not test whether spouse was done or not */
  /*IF ((p^.foff<>NIL) AND (p^.geneloc=0)) OR ((p^.foff<>NIL) AND (NOT noloop))*/
  /*Not already assigned or there are loops*/
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

static void collapsedown(p, LINK)
thisperson *p;
struct LOC_likelihood *LINK;
{
  if ((p->pa == NULL || p->geneloc != 0) && (p->pa == NULL || noloop))
    return;
  /*Not already assigned or there are loops*/
  p->up = true;
  collapseup(p->pa, LINK);
  seg(&p->pa, &p->ma, &p, peeldown, LINK);
}  /* collapsedown */

/*collapsedown*/

static void riskcumul(LINK)
struct LOC_likelihood *LINK;
{
  int i;
  thisarray *WITH;
  int FORLIM;

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

static void riskcalc(LINK)
struct LOC_likelihood *LINK;
{
  double normal;

  LINK->homo /= like;
  LINK->hetero /= like;
  normal = 1 - LINK->homo - LINK->hetero;
  fprintf(simout, "RISK FOR PERSON %6d IN PEDIGREE %7d\n",
	  LINK->proband->id, LINK->proband->origped);
  if (!LINK->proband->male || !sexlink)
    fprintf(simout, "HOMOZYGOTE CARRIER   : %8.5f\n", LINK->homo);
  if (!LINK->proband->male || !sexlink)
    fprintf(simout, "HETEROZYGOTE CARRIER : %8.5f\n", LINK->hetero);
  else
    fprintf(simout, "MALE CARRIER         : %8.5f\n", LINK->hetero);
  fprintf(simout, "NORMAL               : %8.5f\n", normal);
  printf("RISK FOR PERSON %6d IN PEDIGREE %7d\n",
	 LINK->proband->id, LINK->proband->origped);
  if (!LINK->proband->male || !sexlink)
    printf("HOMOZYGOTE CARRIER   : %8.5f\n", LINK->homo);
  if (!LINK->proband->male || !sexlink)
    printf("HETEROZYGOTE CARRIER : %8.5f\n", LINK->hetero);
  else
    printf("MALE CARRIER         : %8.5f\n", LINK->hetero);
  printf("NORMAL               : %8.5f\n", normal);
}  /* riskcalc */


/*seg*/

static void likelihood(thisped_, proband_)
int thisped_;
thisperson *proband_;
{
  int ti, tj; /* Ramana */
  struct LOC_likelihood V;
  int loopmax[maxloop];
  double tmplike;
  int i, j;
  boolean gocalc, alldone;
  thisperson *WITH;
  thisarray *WITH1;
  int FORLIM1;

  /*riskcalc*/

  V.thisped = thisped_;
  V.proband = proband_;
  V.homo = 0.0;
  V.hetero = 0.0;
  tmplike = 0.0;
  alldone = false;
  V.nuscale = 0;
  noloop = true;
  for (i = 0; i < maxloop; i++) {   /*with*/
    V.loopgen[i] = 1;
    loopmax[i] = 1;
    V.holdpoint[i] = NULL;
    if (looppers[V.thisped - 1][i][0] != NULL) {
      WITH = looppers[V.thisped - 1][i][0];
      WITH->gen = (thisarray *)Malloc(sizeof(thisarray));
      if (WITH->gen == NULL)
	malloc_err("gen field in likelihood");
      allocate_thisarray(WITH->gen, fgeno);
      /*        gen:=genpoint(NewPtr(SizeOf(thisarray)));*/
      noloop = false;
      if (feedback)
	printf("making loopperson: %12d\n", WITH->id);
      WITH1 = WITH->gen;
      FORLIM1 = fgeno;
      for (j = 0; j < FORLIM1; j++)
	WITH1->genarray[j] = 0.0;
      /* Version 5.03 goes to maxfem */
      getvect(looppers[V.thisped - 1][i][0], &V);
      /* IF looppers[thisped, i, 1]^.pa=NIL THEN nuscale:=nuscale+1;
      Taken out on basis of version 5.03 */
      /* Alert 5.1 has this in*/
      V.holdpoint[i] = WITH->gen;
      WITH->gen = NULL;
      if (WITH->male)
	loopmax[i] = mgeno;
      else
	loopmax[i] = fgeno;
    }
  }
  /*loop on i*/
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

/* accumulate segprob entries */
/* next test written by R. M. Idury*/
#ifndef LESSMEMORY
  if(V.thisped == 1) {
    for(ti = 0; ti < nuneed; ti++)
      for(tj = 0; tj < nuneed; tj++)
	segprob2[ti*nuneed+tj] = maletheta->segprob[ti] * femaletheta->segprob[tj];
    getprobtable();
  }
#endif
      collapseup(V.proband, &V);
      collapsedown(V.proband, &V);
      if (NULL == V.proband->gen) {
	prob(&V.proband, &V);
      }
      like = 0.0;
      WITH1 = V.proband->gen;
      if (V.proband->male) {
	numgen = mgeno;
	FORLIM1 = mgeno;
	for (i = 0; i < FORLIM1; i++) {
	  like += WITH1->genarray[i];
	  if (WITH1->genarray[i] != 0.0)
	    geneprob[i] = WITH1->genarray[i];
	}
	/*NOTE: with loops, may visit this many times -> don't
	want to wipe out previously assigned values. */
      } else {
	numgen = fgeno;
	FORLIM1 = fgeno;
	for (i = 0; i < FORLIM1; i++) {
	  like += WITH1->genarray[i];
	  if (WITH1->genarray[i] != 0.0)
	    geneprob[i] = WITH1->genarray[i];
	}
      }
      tmplike += like;
      if (risk && like != 0.0)
	riskcumul(&V);
      FORLIM1 = totperson;
      for (i = 1; i <= FORLIM1; i++) {
	if (person[i]->gen != NULL)
	  cleanup(&person[i], &V);
      }
    }
    alldone = true;
    for (i = 0; i < maxloop; i++)
      alldone = (alldone && V.loopgen[i] == loopmax[i]);
  } while (!alldone);
  like = tmplike;
  if (risk && like != 0.0)
    riskcalc(&V);
  if (like != 0.0)
    selectgeneloc(V.proband, V.thisped);
  if (like == 0.0)
    like = zerolike;
  else
    like = log(like) - V.nuscale * log(segscale);
  for (i = 0; i < maxloop; i++) {
    if (V.holdpoint[i] != NULL) {
      Free(V.holdpoint[i]);
      /*      DisposPtr(Ptr(holdpoint[i]));*/
      V.holdpoint[i] = NULL;
    }
  }
}


/*likelihood*/

static void checkzero()
{
  int i;
  thetavalues *WITH;
  int FORLIM;

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


static void iterpeds()
{
  int i, thisone, oldped, nunlink, nfamily, j, numreps;
  double randnum;
  boolean tellneg;
  int FORLIM, FORLIM1, FORLIM2;
  thisperson *WITH;
  int FORLIM3;

  /* numreps = counter of replications */

  nunlink = 0;
  nfamily = 0;
  FORLIM = totperson;
  for (i = 1; i <= FORLIM; i++)
    person[i]->done = false;
  FORLIM = totperson;
  for (i = 1; i <= FORLIM; i++)
    person[i]->firstpass = true;
  thisc = minint;
  /* recombination;
  checkzero;  <= This is not in HLINK. */
  /*Allocation and initialization added by A. A. Schaffer*/
  geneprob = (double *) malloc(fgeno * sizeof(double));
  if (geneprob == NULL)
    malloc_err("geneprob");
  for(i = 0; i < fgeno; i++)
    geneprob[i] = 0.0;

  for (i = 1; i <= 35; i++)
    putc('-', simout);
  putc('\n', simout);
  for (i = 1; i <= 35; i++)
    putc('-', simout);
  putc('\n', simout);
  if (sexdif)
    fprintf(simout, "TRUE MALE THETAS FOR LINKED ORDER ");
  else
    fprintf(simout, "TRUE THETAS FOR LINKED ORDER ");
  FORLIM = nlocus - 2;
  for (i = 0; i <= FORLIM; i++)
    fprintf(simout, "%11.6f", maletheta->theta[i]);
  if (interfer)
    fprintf(simout, "%11.6f", maletheta->theta[nlocus - 1]);
  putc('\n', simout);
  if (sexdif) {
    fprintf(simout, "TRUE FEMALE THETAS FOR LINKED ORDER ");
    FORLIM = nlocus - 2;
    for (i = 0; i <= FORLIM; i++)
      fprintf(simout, "%11.6f", femaletheta->theta[i]);
    if (interfer)
      fprintf(simout, "%11.6f", femaletheta->theta[nlocus - 1]);
    putc('\n', simout);
  }
  for (i = 1; i <= 35; i++)
    putc('-', simout);
  putc('\n', simout);
  for (i = 1; i <= 35; i++)
    putchar('-');
  putchar('\n');
  if (sexdif)
    printf("TRUE MALE THETAS FOR LINKED ORDER   ");
  else
    printf("TRUE THETAS FOR LINKED ORDER ");
  FORLIM = nlocus - 2;
  for (i = 0; i <= FORLIM; i++)
    printf("%11.6f", maletheta->theta[i]);
  if (interfer)
    printf("%11.6f", maletheta->theta[nlocus - 1]);
  putchar('\n');
  if (sexdif) {
    printf("TRUE FEMALE THETAS FOR LINKED ORDER ");
    FORLIM = nlocus - 2;
    for (i = 0; i <= FORLIM; i++)
      printf("%11.6f", femaletheta->theta[i]);
    if (interfer)
      printf("%11.6f", femaletheta->theta[nlocus - 1]);
    putchar('\n');
  }
  for (i = 1; i <= 35; i++)
    putchar('-');
  putchar('\n');
  for (i = 1; i <= 35; i++)
    putchar('-');
  printf("\n UNLINKED ORDER OF LOCI: ");
  FORLIM = nlocus;
  for (i = 1; i <= FORLIM; i++) {
    thissystem = 1;
    while (order2[thissystem - 1] != i)
      thissystem++;
    printf("%3d", thissystem);
  }
  putchar('\n');
  for (i = 1; i <= 35; i++)
    putchar('-');
  putchar('\n');
  for (i = 1; i <= 35; i++)
    putc('-', simout);
  fprintf(simout, "\n UNLINKED ORDER OF LOCI: ");
  FORLIM = nlocus;
  for (i = 1; i <= FORLIM; i++) {
    thissystem = 1;
    while (order2[thissystem - 1] != i)
      thissystem++;
    fprintf(simout, "%3d", thissystem);
  }
  putc('\n', simout);
  for (i = 1; i <= 35; i++)
    putc('-', simout);
  putc('\n', simout);
  for (i = 1; i <= 35; i++)
    putc('-', simout);
  putc('\n', simout);
  if (sexdif)
    fprintf(simout, "TRUE MALE THETAS FOR UNLINKED ORDER ");
  else
    fprintf(simout, "TRUE THETAS FOR UNLINKED ORDER ");
  FORLIM = nlocus - 2;
  for (i = 0; i <= FORLIM; i++)
    fprintf(simout, "%11.6f", maletheta->theta2[i]);
  if (interfer)
    fprintf(simout, "%11.6f", maletheta->theta2[nlocus - 1]);
  putc('\n', simout);
  if (sexdif) {
    fprintf(simout, "TRUE FEMALE THETAS FOR UNLINKED ORDER ");
    FORLIM = nlocus - 2;
    for (i = 0; i <= FORLIM; i++)
      fprintf(simout, "%11.6f", femaletheta->theta2[i]);
    if (interfer)
      fprintf(simout, "%11.6f", femaletheta->theta2[nlocus - 1]);
    putc('\n', simout);
  }
  for (i = 1; i <= 35; i++)
    putc('-', simout);
  putc('\n', simout);
  for (i = 1; i <= 35; i++)
    putchar('-');
  putchar('\n');
  if (sexdif)
    printf("TRUE MALE THETAS FOR UNLINKED ORDER   ");
  else
    printf("TRUE THETAS FOR UNLINKED ORDER ");
  FORLIM = nlocus - 2;
  for (i = 0; i <= FORLIM; i++)
    printf("%11.6f", maletheta->theta2[i]);
  if (interfer)
    printf("%11.6f", maletheta->theta2[nlocus - 1]);
  putchar('\n');
  if (sexdif) {
    printf("TRUE FEMALE THETAS FOR UNLINKED ORDER ");
    FORLIM = nlocus - 2;
    for (i = 0; i <= FORLIM; i++)
      printf("%11.6f", femaletheta->theta2[i]);
    if (interfer)
      printf("%11.6f", femaletheta->theta2[nlocus - 1]);
    putchar('\n');
  }
  for (i = 1; i <= 35; i++)
    putchar('-');
  printf("\n Beginning simulations... \n");


  getdatetime(&isec, &wsec);   /*SUN*/
  if (print)
    printf("PEDIGREE |  PERSON   |   LN LIKE  \n");
  FORLIM = replicates;
  for (numreps = 1; numreps <= FORLIM; numreps++) {
    if (numreps % 10 == 0)
      putchar('*');
    else
      putchar('.');
    fflush(stdout);
/*    P_ioresult = 0;   * old code for SUN*/
    /*PLFlush(output);*/
    if (numreps % 50 == 0)
      printf(" %d replicates generated\n", numreps);
    if (print)
      printf("   Replicate #%7d\n", numreps);
    oldped = 0;
    FORLIM1 = totperson;
    for (i = 1; i <= FORLIM1; i++)
      person[i]->geneloc = 0;
    FORLIM1 = totperson;
    for (thisone = 1; thisone <= FORLIM1; thisone++)
    {   /* Write out pedigree file a person at a time */
      FORLIM2 = totperson;
      for (i = 1; i <= FORLIM2; i++) {
	person[i]->done = false;
	person[i]->firstpass = true;
      }
      if (person[thisone]->ped != oldped)
      {   /* if person[thisone]^.ped<>oldped */
	nfamily++;
	oldped = person[thisone]->ped;
	randnum = randomnum();
	if (randnum < hetvalue) {
	  nunlink++;
	  if (unlinked == false) {
	    linked = false;
	    unlinked = true;
	    memcpy(maletheta->theta, maletheta->theta2, sizeof(thetarray));
	    memcpy(femaletheta->theta, femaletheta->theta2, sizeof(thetarray));
	    memcpy(order, order2, maxlocus * sizeof(int));
	    FORLIM2 = totperson;
	    for (i = 1; i <= FORLIM2; i++)
	    {   /*totperson includes all the pedigrees*/
	      WITH = person[i];
	      FORLIM3 = nlocus;
	      for (j = 0; j < FORLIM3; j++)
		WITH->phen[j] = WITH->phen2[j];
	    }
	    memcpy(thislocus, thislocus2, maxlocus * sizeof(locusvalues *));
	    increment[nlocus - 1] = 1;
	    for (i = nlocus - 1; i >= 1; i--)
	      increment[i - 1] = increment[i] * thislocus[i]->nallele;
	    getlocations();
	    getgenetables(); /*A. A. Schaffer inserted this line*/
	    recombination();
	  }
	} else if (linked == false) {
	  linked = true;
	  unlinked = false;
	  memcpy(maletheta->theta, maletheta->theta1, sizeof(thetarray));
	  memcpy(femaletheta->theta, femaletheta->theta1, sizeof(thetarray));
	  memcpy(order, order1, maxlocus * sizeof(int));
	  FORLIM2 = totperson;
	  for (i = 1; i <= FORLIM2; i++) {
	    WITH = person[i];
	    FORLIM3 = nlocus;
	    for (j = 0; j < FORLIM3; j++)
	      WITH->phen[j] = WITH->phen1[j];
	  }
	  memcpy(thislocus, thislocus1, maxlocus * sizeof(locusvalues *));
	  increment[nlocus - 1] = 1;
	  for (i = nlocus - 1; i >= 1; i--)
	    increment[i - 1] = increment[i] * thislocus[i]->nallele;
	  getlocations();
          getgenetables(); /*A. A. Schaffer inserted this line*/
	  recombination();
	}
      }  /*if linked = false*/
      /*          if linked then
      status[person[thisone]^.ped] := TRUE
      else
      status[person[thisone]^.ped] := FALSE; */
      /* If there is a loop, then the second loopperson will have
      a genotype assigned before doing a risk calculation for them */
      if (person[thisone]->geneloc == 0) {
	likelihood(person[thisone]->ped, person[thisone]);
	if (print)
	  printf("%9d %9d %12.6f \n",
		 person[thisone]->origped, person[thisone]->id, like);
      }  /*if person[thisone]^.geneloc=0*/
      writelinkage(thisone);
    }
    /* for thisone := 1 to totperson */
    if (firsttime) {
      if (thisc < maxcensor)
	printf("Maxcensor can be reduced to %12d\n", thisc);
      else if (thisc > maxcensor)
	printf("you may gain efficiency by increasing maxcensor\n");
      getdatetime(&fsec, &wfsec);
	  /*SUN - first is user + sys time, second is wall time*/
      fsec -= isec;
      /*     if fsec<0.0 then fsec:=fsec+86400.0;*/
      /*SUN */
      fprintf(simout, " Elapsed Time for one replicate = %8.2f seconds\n",
	      fsec);
      printf(" Elapsed Time for one replicate = %8.2f seconds\n", fsec);
      printf(" Estimated Time to finish = %8.2f minutes = %6.2f hours\n",
	     replicates * fsec / 60.0, replicates * fsec / 3600.0);
    }
    firsttime = false;
  }  /*loop on numreps*/

  putc('\n', simout);
  putchar('\n');
  getdatetime(&fsec, &wfsec);   /*SUN*/
  fsec = (wfsec - wsec) / 60.0;   /*Elapsed time in minutes*/
  /*SUN - use wall clock*/
  tellneg = (fsec < 0.0);
  if (tellneg)   /*allow for time across one midnight*/
    fsec += 1440.0;
  fprintf(simout, " Elapsed Time =%9.2f min. or %8.2f hours",
	  fsec, fsec / 60.0);
  if (tellneg)
    fprintf(simout, " or more - ran over midnight");
  putc('\n', simout);
  printf(" Elapsed Time =%9.2f min. or %8.2f hours", fsec, fsec / 60.0);
  if (tellneg)
    printf(" or more - ran over midnight");
  putchar('\n');

  firsttime = false;   /* Should firsttime be here???*/
  for (i = 1; i <= 35; i++)
    putc('-', simout);
  putc('\n', simout);
  for (i = 1; i <= 35; i++)
    putchar('-');
  printf("\nActual proportion of unlinked families: %6.3f\n",
	 (double)nunlink / nfamily);
  fprintf(simout, "Actual proportion of unlinked families: %6.3f\n",
	  (double)nunlink / nfamily);

  /*if unlinked=FALSE*/
  /*if randnum<hetvalue */
}


/* iterpeds */


static void initialize()
{
  int i;

  firsttime = true;
  for (i = 0; i <= 3; i++)   /* Initialize code counters */
    code[i] = 0;
  men = 0;
  women = 0;   /* Initialize counters */
  printf("\nProgram SLINK version %s\n\n", version);
  printf("The program constants are set to the following maxima:\n");
  printf("%6d loci in mapping problem (maxlocus)\n", maxlocus);
  printf("%6d alleles at a single locus (maxall)\n", maxall);
  printf("%6d maximum of censoring array (maxcensor)\n", maxcensor);
  printf("%6d individuals in all pedigrees combined (maxind)\n",
	 maxind);
  printf("%6d pedigrees (maxped)\n", maxped);
  printf("%6d binary codes at a single locus (maxfact)\n", maxfact);
  printf("%6d quantitative factor(s) at a single locus\n", maxtrait);
  printf("%6d liability classes\n", maxliab);
  printf("%6d binary codes at a single locus\n", maxfact);
  printf("%8.2f base scaling factor for likelihood (scale)\n", scale);
  printf("%8.2f scale multiplier for each locus (scalemult)\n", scalemult);
  printf("%8.5f frequency for elimination of heterozygotes (minfreq)\n",
	 minfreq);
  if (minfreq != 0.0) {
    printf("IMPORTANT : RECOMPILE THIS PROGRAM WITH MINFREQ=0.0\n");
    printf("FOR THE ANALYSIS OF RECESSIVE TRAITS\n");
  }

  for (i = 1; i <= 55; i++)
    putchar('*');
  printf("\n This is SLINK, which is a general simulation program.\n");
  printf(" The output file pedfile.dat may be analyzed with the\n");
  printf(" companion programs MSIM, LSIM, and ISIM.\n");
  for (i = 1; i <= 55; i++)
    putchar('*');
  putchar('\n');

  /* assign(simdata,'simdata.dat');
     assign(simped,'simped.dat');
     assign(simout,'simout.dat');
     assign(pedfile,'pedfile.dat'); */
  if (simout != NULL)
    simout = freopen("simout.dat", "w", simout);
  else
    simout = fopen("simout.dat", "w");
  if (simout == NULL)
    exit(FileNotFound);
  if (pedfile != NULL) {
    pedfile = freopen("pedfile.dat", "w", pedfile);
  } else
    pedfile = fopen("pedfile.dat", "w");
  if (pedfile == NULL)
    exit(FileNotFound);
  if (show_trait_genotypes) {
    if (traitfile != NULL) {
      traitfile = freopen("traitfile.dat", "w", traitfile);
    } else
      traitfile = fopen("traitfile.dat", "w");
    if (traitfile == NULL)
      exit(FileNotFound);
  }
  totall = 1;
}


/* initialize */


static void statistics(fout)
FILE **fout;
{
  int i;

  fprintf(*fout, " The random number seed is:%6d\n", seed1);
  fprintf(*fout, " The number of replications is:%6d\n", replicates);
  fprintf(*fout, " The requested proportion of unlinked families is: %6.3f\n",
	  hetvalue);
  if (trait > 0)
    fprintf(*fout, " The trait locus is locus number:%4d\n", trait);
  else
    fprintf(*fout,
      " There is no trait locus -> all loci will be treated as marker loci\n");
  fprintf(*fout, "    Summary Statistics about simped.dat \n");
  fprintf(*fout, " Number of pedigrees %7d\n", nuped);
  fprintf(*fout, " Number of people    %7d\n", totperson);
  fprintf(*fout, " Number of females   %7d\n", women);
  fprintf(*fout, " Number of males     %7d\n", men);
  for (i = 0; i <= 3; i++) {   /* case */
    fprintf(*fout, " There were%4d in category: ", code[i]);
    switch (i) {

    case 0:
      fprintf(*fout, " Marker Unknown");
      if (trait != 0)
	fprintf(*fout, "; Trait original\n");
      else
	putc('\n', *fout);
      break;

    case 2:
      fprintf(*fout, " Marker Available");
      if (trait != 0)
	fprintf(*fout, "; Trait original\n");
      else
	putc('\n', *fout);
      break;

    case 1:
      fprintf(*fout, " Marker Available");
      if (trait != 0)
	fprintf(*fout, "; Trait simulated\n");
      else
	putc('\n', *fout);
      break;

    case 3:
      fprintf(*fout, " Marker Unknown");
      if (trait != 0)
	fprintf(*fout, "; Trait simulated\n");
      else
	putc('\n', *fout);
      break;
    }
  }
  /* loop on i */
}


/* statistics */

/*read parameters of the run from slinkin.dat;*/
static void getparam()
{
  FILE *slinkin;

  /* assign(slinkin,'SLINKIN.DAT');*/
  /*SUN*/
  slinkin = NULL;
  printf("Reading SLINKIN.DAT\n");
  if (slinkin != NULL)
    slinkin = freopen("slinkin.dat", "r", slinkin);
  else
    slinkin = fopen("slinkin.dat", "r");
  if (slinkin == NULL)
    exit(FileNotFound);
  /*SUN*/
  /*Initialize seed2 and seed3*/
  seed2 = 28007;
  seed3 = 24001;
  fscanf(slinkin, "%d%*[^\n]", &seed1);
  getc(slinkin);   /*Random seed*/
  /*SUN*/
  if (seed1 < 1 || seed1 > 30323) {
    printf("ERROR in SLINKIN.DAT: Seed <1 or >30323\n");
    scanf("%*[^\n]");
    getchar();
    exit(0);   /*SUN*/
  }
  fscanf(slinkin, "%d%*[^\n]", &replicates);
  getc(slinkin);   /*Number of replicates*/
  /*SUN*/
  if (replicates < 1) {
    printf("ERROR in SLINKIN.DAT: Number of replicates <1\n");
    scanf("%*[^\n]");
    getchar();
    exit(0);   /*SUN*/
  }
  fscanf(slinkin, "%d%*[^\n]", &trait);
  getc(slinkin);   /*Number of the "trait" locus*/
  /*SUN*/
  if ((unsigned int)trait > nlocus) {
    printf("ERROR in SLINKIN.DAT: Trait number not between 0 and %d\n",
	   nlocus);
    scanf("%*[^\n]");
    getchar();
    exit(0);   /*SUN*/
  }
  fscanf(slinkin, "%lg%*[^\n]", &hetvalue);
  getc(slinkin);   /*Proportion of unlinked families*/
  /*SUN*/
  if ((unsigned)hetvalue > 1.0) {
    printf("ERROR in SLINKIN.DAT: Proportion of het. families not between 0 and 1\n");
    scanf("%*[^\n]");
    getchar();
    exit(0);   /*SUN*/
  }
  if (slinkin != NULL)
    fclose(slinkin);
}

/* update_slinkin  changes slinkin.dat
  to have a new random seed and move the
  old file to slinkin.dat.bak*/
void update_slinkin()
{
  FILE *slinkin;
  int new_seed;

  system("mv slinkin.dat slinkin.dat.bak");
  new_seed = generate_new_seed();
  slinkin = fopen("slinkin.dat", "w");
  fprintf(slinkin,"%d  <= random seed\n", new_seed);
  fprintf(slinkin,"%d <= replicates\n", replicates);
  fprintf(slinkin,"%d <= trait locus\n", trait);
  fprintf(slinkin,"%lf <= proportion of unlinked families\n", hetvalue);
  fclose(slinkin);
  printf("Finished changing seed in slinkin.dat, old version moved to slinkin.dat.bak\n");
}

/* getparam */


int main(argc, argv)
int argc;
char *argv[];
{
  FILE *TEMP;
  int FORLIM;
  thetavalues *WITH;
  thisperson *WITH1;
  int FORLIM1;
  int c;
  int i, j, k;

/*  PASCAL_MAIN(argc, argv);*/
  pedfile = NULL;
  simdata = NULL;
  simped = NULL;
  simout = NULL;
  rewrite_slinkin = false;
  all_available = false;
  show_trait_genotypes = false;
#if !defined(MS_VISUAL_C)
  while ((c = getopt(argc, argv, "ast")) != -1) {
    switch (c) {
    case 'a':
      all_available = true;
      break;
    case 's':
      rewrite_slinkin = true;
      break;
    case 't':
      show_trait_genotypes = true;
      break;
    default:
    break;
    }
  }
#endif
  initialize();
  inputdata();
  getparam();

  statistics(&simout);
  TEMP = stdout;
  statistics(&TEMP);

  firsttime = true;
  lasttime = false;
  dolod = false;
  censorstruct = (censorrec *)Malloc(sizeof(censorrec));
/*  for (k = minint; k <= maxcensor; k++)
    censorstruct->censor[k - minint] = false; */
  getgeneindices(); /*Ramana*/
  allocategenetables(); /*A. A. Schaffer*/
  getlocations();
  getgenetables();
  if (simped != NULL)
    fclose(simped);
  simped = NULL;   /*  close(simped,true);*/
  /*SUN*/
  fprintf(simout, "LINKAGE/SLINK (V%s) WITH%3d-POINT", version, nlocus);
  if (sexlink)
    fprintf(simout, " SEXLINKED DATA\n");
  else
    fprintf(simout, " AUTOSOMAL DATA\n");
  for (i = 1; i <= 35; i++)
    putc('-', simout);
  fprintf(simout, "\n LINKED ORDER OF LOCI: ");
  FORLIM = nlocus;
  for (i = 1; i <= FORLIM; i++) {
    thissystem = 1;
    while (order[thissystem - 1] != i)
      thissystem++;
    fprintf(simout, "%3d", thissystem);
  }
  putc('\n', simout);
  for (i = 1; i <= 35; i++)
    putc('-', simout);
  putc('\n', simout);
  for (i = 1; i <= 35; i++)
    putchar('-');
  printf("\n LINKED ORDER OF LOCI: ");
  FORLIM = nlocus;
  for (i = 1; i <= FORLIM; i++) {
    thissystem = 1;
    while (order[thissystem - 1] != i)
      thissystem++;
    printf("%3d", thissystem);
  }
  putchar('\n');
  for (i = 1; i <= 35; i++)
    putchar('-');
  putchar('\n');
  /*Setup unlinked trait locus*/
  if (trait != 0)   /*Corrected 4/1/91*/
    j = order[trait - 1];
  else
    j = 1;
  FORLIM = nlocus;
  for (i = 1; i <= FORLIM; i++) {
    if (i != trait) {
      if (order[i - 1] < j)
	order2[i - 1] = order[i - 1] + 1;
      else
	order2[i - 1] = order[i - 1];
    }
  }
  if (trait != 0)
    order2[trait - 1] = 1;
  memcpy(order1, order, maxlocus * sizeof(int));
  WITH = maletheta;
  FORLIM = nlocus;
  for (i = j; i < FORLIM; i++)
    WITH->theta2[i - 1] = WITH->theta[i - 1];
  WITH->theta2[0] = 0.5;
  FORLIM = j;
  for (i = 2; i < FORLIM; i++)
    WITH->theta2[i - 1] = WITH->theta[i - 2];
  if (j != 1 && j != nlocus)
    WITH->theta2[j - 1] = WITH->theta[j - 2] * (1.0 - WITH->theta[j - 1]) +
			  WITH->theta[j - 1] * (1.0 - WITH->theta[j - 2]);
  memcpy(WITH->theta1, WITH->theta, sizeof(thetarray));
  WITH = femaletheta;
  FORLIM = nlocus;
  for (i = j; i < FORLIM; i++)
    WITH->theta2[i - 1] = WITH->theta[i - 1];
  WITH->theta2[0] = 0.5;
  FORLIM = j;
  for (i = 2; i < FORLIM; i++)
    WITH->theta2[i - 1] = WITH->theta[i - 2];
  if (j != 1 && j != nlocus)
    WITH->theta2[j - 1] = WITH->theta[j - 2] * (1.0 - WITH->theta[j - 1]) +
			  WITH->theta[j - 1] * (1.0 - WITH->theta[j - 2]);
  memcpy(WITH->theta1, WITH->theta, sizeof(thetarray));
  FORLIM = nlocus;
  for (i = 1; i <= FORLIM; i++) {
    if (order[i - 1] < j)
      thislocus2[order[i - 1]] = thislocus[order[i - 1] - 1];
    else if (order[i - 1] > j)
      thislocus2[order[i - 1] - 1] = thislocus[order[i - 1] - 1];
  }
  thislocus2[0] = thislocus[j - 1];
  memcpy(thislocus1, thislocus, maxlocus * sizeof(locusvalues *));
  FORLIM = totperson;
  for (k = 1; k <= FORLIM; k++) {
    WITH1 = person[k];
    FORLIM1 = nlocus;
    for (i = 1; i <= FORLIM1; i++)
      WITH1->phen1[i - 1] = WITH1->phen[i - 1];
    FORLIM1 = nlocus;
    for (i = 1; i <= FORLIM1; i++) {
      if (order[i - 1] < j)
	WITH1->phen2[order[i - 1]] = WITH1->phen[order[i - 1] - 1];
      else if (order[i - 1] > j)
	WITH1->phen2[order[i - 1] - 1] = WITH1->phen[order[i - 1] - 1];
    }
    WITH1->phen2[0] = WITH1->phen[j - 1];
  }
  linked = false;
  unlinked = false;
  /*Need only call once at the theta value given in simdata*/
  iterpeds();
  printf(" SLINK has completed %12d replicates.\n", replicates);
  if (rewrite_slinkin)
    update_slinkin();
  freegenetables();
  if (simout != NULL)
    fclose(simout);
  if (simped != NULL)
    fclose(simped);
  if (simdata != NULL)
    fclose(simdata);
  if (pedfile != NULL)
    fclose(pedfile);
  if (show_trait_genotypes) {
    if (traitfile != NULL)
      fclose(traitfile);
  }
  exit(EXIT_SUCCESS);
}  /* slink */



/* End. */



