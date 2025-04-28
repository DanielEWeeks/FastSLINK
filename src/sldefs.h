/*This file is part of a C version of SLINK a simulation program based
on the LINKAGE package. SLINK was originally described in: D. E.
Weeks, J. Ott, G. M. Lathrop (1990), SLINK: a general simulation
program for linkage analysis, American Journal of Human Genetics
47(1990), A204 (abstract). This version of SLINK is adapted to use the
faster pedigree traversal algorithms described in: R. W. Cottingham
Jr., R. M. Idury, A. A. Schaffer (1993), Faster sequential genetic
linkage computations, American Journal of Human Genetics 53(1993), pp.
252-263.*/
/*This file contains definitions and declaration, which differ from those
used in the LINKAGE programs. In some cases the declarations are totally
new and in other cases they are changed*/
/* 8 December 1993, A. A. Schaffer pulled these out of slink.c and adapted
for faster algorithms */

#define print           false
    /*FALSE turns off some output for quicker running*/
#define feedback        false   /*DESCRIPTION OF PATH THROUGH FAMILIES*/

#define unafqu          1
/*CODE FOR UNAFFECTED MALES AT A SEX-LINKED
  QUANTITATIVE TRAIT; This must not equal affall*/
/* AFFECTION STATUS */
#define unaff           1   /*CODE FOR UNAFFECTED INDIVIDUAL*/
/*Note: unaff must not equal affval */
    /*MAXIMUM NUMBER OF BINARY CODES AT A SINGLE LOCUS*/
/* OTHERS */

#define scale           1.0   /*SCALE FACTOR*/
#define zerolike        (-1.0e20)

typedef int hapvector[maxlocus];

typedef struct locusvalues {
  int nallele; /*number of alleles at the locus*/
  int format; /*which type of locus, 3 is for numbered alleles*/
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
    struct {
      num_phenarray num_allele; /* new for allowing maxall > 32*/
      bin_phenarray bin_allele; /* new for allowing maxall > 32*/
      int numfact;
    } U2;
  } UU;
} locusvalues;


typedef struct thisperson {
  int origped;   /* L. Sandkuijl's suggestion */
  int id, ped, inloop, geneloc, prbnum, availnum, origid;
  /*geneloc = location of the simulated genotype in the genarray;
  prbnum = proband number; availnum = availability code;
  origid = original id found in pedigree file*/
  int input_availnum; /*use to store the input, when -a
                        is used to make all available*/
  struct thisperson *pa, *ma, *foff, *nextpa, *nextma;
  thisarray *gen;
  indphen phen, phen1, phen2;
  phenotype *privphen;
  boolean unknown, multi, done, up, male, firstpass;
  information *store;
} thisperson;


double wsec, wfsec, isec, fsec, hetvalue, finish, inc;
/* Used for timing replicates */
/*SUN*/
int code[4];   /* Counter of availnum codes */
int men, women;   /* Counters of men and women in simped.dat */
boolean noloop;
genotype geneprob;
    /*Array of genotype probabilities taken from genarray*/
int totall, trait, replicates;
/*trait = number of the trait locus;
replicates = desired number of replications */
int numgen, seed1, seed2, seed3;
/*numgen = number of genotypes; random number seeds*/
int **hapstore;
/*Array of alleles corresponding to the position geneloc in genarray*/
boolean zeromale[maxlocus], zerofemale[maxlocus];
locusvalues *thislocus[maxlocus], *thislocus1[maxlocus], *thislocus2[maxlocus];
int order1[maxlocus], order2[maxlocus];
/*PEOPLE*/
thisperson *person[maxind + 1];
thisperson *proband[maxped];
thisperson *looppers[maxped][maxloop][2];
/*OTHERS*/
int  nlocus, which, thissystem;
boolean  unlinked, linked;
FILE *simout, *simped, *simdata, *pedfile;
/*pedfile = where the simulated pedigrees are written to*/
FILE *traitfile;
/*traitfile is an additional optional file for writing genotypes at the
  trait locus*/
boolean show_trait_genotypes; /*if true, print trait genotypes in traitfile*/
boolean rewrite_slinkin; /*if true, slinkin.dat gets rewritten with a new seed*/
boolean all_available; /*if true, all individuals are treated as if they have availability code 2*/ 


/* Local variables for seg: */
struct LOC_seg {
  struct LOC_likelihood *LINK;
  thisperson **p, **q, **r, *child, *father, *mother;
  int fseg, sseg, sstart, send, fstart, fend, nfirst, nsecond, firstseg,
       secondseg;
  double pf, ps;
  thetavalues *firstsex, *secondsex;
} ;

/*Local variables for likelihood*/
struct LOC_likelihood {
  int thisped;
  thisperson *proband;
  int loopgen[maxloop];
  double homo, hetero;
  int nuscale;
  thisarray *holdpoint[maxloop];
} ;


