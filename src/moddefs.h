/*This file is part of a C version of SLINK a simulation program based
on the LINKAGE package. SLINK was originally described in: D. E.
Weeks, J. Ott, G. M. Lathrop (1990), SLINK: a general simulation
program for linkage analysis, American Journal of Human Genetics
47(1990), A204 (abstract). This version of SLINK is adapted to use the
faster pedigree traversal algorithms described in: R. W. Cottingham
Jr., R. M. Idury, A. A. Schaffer (1993), Faster sequential genetic
linkage computations, American Journal of Human Genetics 53(1993), pp.
252-263. */
/*These are some declarations adapted from the faster LINKAGE programs,
version 1.1*/
/* The versions in this file use a lot of memory */
/* Definitions for somewhat slower versions that use less memory */
/* are in slowmoddefs.h */

#ifndef  _MODDEFS_H

#define  _MODDEFS_H  1



#ifndef AUTOSOMAL_RUN
#define AUTOSOMAL_RUN  1  /*1 is safe; can be set to 0 to make
                            sexlinked runs space efficient*/
#endif

#ifndef SEXDIF_RUN
#define SEXDIF_RUN     1  /*1 is safe; can be set to 0 on sexlinked runs or
                            runs in which maletheta and femaletheta will be
                            the same*/
#endif


/*segprob2 stores products of recombination or
nonrecombination probabilities; one probability
comes from each parent */
double *segprob2;

/*used in the autosomal case when sexdif=1*/
double *probtabledif;
double *probtable;
unsigned *probtableindex;

unsigned *classsize;
unsigned *classbase;

typedef double childprob[maxchild]; 
childprob **partialprob;

unsigned *invpool, *nextpool, *indpool;


/* _MODDEFS_H */

#endif
