/*This file is part of a C version of SLINK a simulation program based
on the LINKAGE package. SLINK was originally described in: D. E.
Weeks, J. Ott, G. M. Lathrop (1990), SLINK: a general simulation
program for linkage analysis, American Journal of Human Genetics
47(1990), A204 (abstract). This version of SLINK is adapted to use the
faster pedigree traversal algorithms described in: R. W. Cottingham
Jr., R. M. Idury, A. A. Schaffer (1993), Faster sequential genetic
linkage computations, American Journal of Human Genetics 53(1993), pp.
252-263.*/
/*This file contains new versions of the routines segup, segdown, */
/*and some auxiliary routines for use with */
/*LODSCORE, ILINK, and LINKMAP */
/* The versions in this file use a lot of memory */
/* Somewhat slower versions that use less memory */
/* are in slowautomodified.c */
/* Most of the code in this file was written by R. M. Idury*/
/* 9 December 1993 A. A. Schaffer put in code to clip pedigree
traversal using the variables parunk and noloop*/
/* 2 August 1994, A. A. Schaffer fixed problem with floating point
   exception, by always zeroing out gene[]*/

#include "commondefs.h"
#include "sldefs.h"
#include "moddefs.h"

Void getprobtable()  /* compute joint recombination probabilities */
{
  static int size = 0, firsttime = 1;
  int index, i, j, ii, jj;
  char *recomb_flag;

  if(firsttime && !sexlink) {
    if (!AUTOSOMAL_RUN)
      fprintf(stderr, "\nWARNING: You are doing an autosomal run but have AUTOSOMAL_RUN set to 0");
    if ((!SEXDIF_RUN) && sexdif)
      fprintf(stderr, "\nWARNING: You are doing a sexdif run but have SEXDIF_RUN set to 0");
    firsttime = 0;
    classbase = (unsigned int *) malloc(nuprobclass * sizeof(unsigned int));
    if (classbase == NULL)
      malloc_err("classbase");
    classsize = (unsigned int *) malloc(nuprobclass * sizeof(unsigned int));
    if (classsize == NULL)
      malloc_err("classsize");
    recomb_flag = (char *) malloc(nuneed * sizeof(char));
    if (recomb_flag == NULL)
      malloc_err("recomb_flag");
    for(i = 0; i < nuneed; i++)
      recomb_flag[i] = 0;
    for(i = 0; i < fgeno; i++)
      recomb_flag[probstart[i]-1] = 1;
    size = 0;
    for(i = 0; i < nuneed; i++)
      if(recomb_flag[i] == 1) classbase[size++] = i;
    classsize[size-1] = nuneed - classbase[size-1];
    for(i = size-2; i >= 0; i--)
      classsize[i] = classbase[i+1] - classbase[i];
    if (!sexlink) {
      probtableindex = (unsigned *) malloc(nuneed * nuneed * sizeof (unsigned));
      if (probtableindex == NULL) {
	fprintf(stderr, "\nYou probably need to run the slower version of this program");
	malloc_err("probtableindex");
      }
      probtable = (double *) malloc(maxisozygclass * nuneed * nuneed 
				 * sizeof(double));
      if (probtable == NULL) {
	fprintf(stderr, "\n You probably need to use the slower version of this program");
	malloc_err("probtable");
      }
    }
    if (sexdif) {
      probtabledif = (double *) malloc(maxisozygclass * nuneed * nuneed 
				 * sizeof(double));
      if (probtabledif == NULL) {
	fprintf(stderr, "\n You probably need to use the slower version of this program");
	malloc_err("probtabledif");
      }
    }
    free(recomb_flag);
  }
			 
  index = 0;
  if (!sexlink) {
    for(i = 0; i < size; i++) {
      for(j = 0; j < size; j++) {
	probtableindex[classbase[i]*nuneed+classbase[j]] = index;
	for(ii = classbase[i]; ii < classbase[i]+classsize[i]; ii++) {
	  for(jj = classbase[j]; jj < classbase[j]+classsize[j]; jj++) {
	    probtable[index++] = segprob2[ii*nuneed+jj];
	  }
	}
      }
    }
  }
  if (sexdif) {
    index = 0;
    for(i = 0; i < size; i++) {
      for(j = 0; j < size; j++) {
        for(ii = classbase[i]; ii < classbase[i]+classsize[i]; ii++) {
          for(jj = classbase[j]; jj < classbase[j]+classsize[j]; jj++) {
	    probtabledif[index++] = segprob2[jj*nuneed+ii];
	  }
	}
      }
    }
  }
}


/*segsum2 is used in segup to compute some common subexpressions in
the probability updates.
first and second are the joint genotypes of the two parents.
fslength is the number of probabilities needed for the combined
isozygote class of both parental genotypes; fs stands for the product
of first and second. LINK is used to pass the genetic information for
parents and children */

/* cgh -- made this void for gcc */

Local void segsum2(first,second,fslength,LINK)
struct LOC_seg *LINK;
{
  int g1, g2, g3, g4; /*indices to store gene numbers*/
  int f1, f2, s1, s2; /*indices to store haplotype numbers*/
  int index; /*counter for combined isozygote equivalence class,
                will range from 0 to fslength-1*/
  int i, j, k, l; /*loop indices*/
  int FORLIM; /*limit on a for loop*/
  int FORLIM1, FORLIM2; /*lower and upper limits on a for loop*/
  int *TEMPGENE1, *TEMPGENE2; /*temporary pointers into genenumber table*/
  double *tempwith3; /*temporarily stores genarray for current child*/

  FORLIM = fence[first]; /*start of the isozygote class for first*/
  FORLIM1 = base[second]; /*start of the isozygote class for second*/
  FORLIM2 = fence[second]; /*end of the isozygote class for second*/
  index = 0; /*start with first isozygote in the class*/
  /*iterate over isozygotes of i*/
  for (i = base[first]; i < FORLIM; i++) {
    f1 = haps1[i]; /*retrieve first haplotype of joint genotype i*/
    f2 = haps2[i]; /*retrieve second haplotype*/
    TEMPGENE1 = genenumber[f1 - 1]; /*lookup partial information
					    in genenumber table*/
    TEMPGENE2 = genenumber[f2 - 1];
    /*iterate over isozygotes of j*/
    for (j = FORLIM1; j < FORLIM2; j++) {
      s1 = haps1[j]; /*retrieve first haplotype of joint genotype j*/
      s2 = haps2[j]; /*retrieve second haplotype*/
     /*lookup the four ways to combine one haplotype from i and one from j*/
      g1 = TEMPGENE1[s1 - 1];
      g2 = TEMPGENE1[s2 - 1];
      g3 = TEMPGENE2[s1 - 1];
      g4 = TEMPGENE2[s2 - 1];

    /*iterate over children; update partial computations of
      new probabilities and store in array tempseg
      note that tempseg has nchild*fslength relevant entries
      The fslength entries for child 1 come first, then
      the fslength entries, for child 2, etc. This is why
      the increment on l is fslength for each change of child*/

      for (l = index, k = 0 ; k < nchild; l += fslength, k++) {
	tempwith3 = thischild[k]->genarray; /*retrieve genarray*/

       /* sum the probabilities for the four joint genotypes that
        the child can get from the current parental joint genotypes*/

	tempseg[l] = (tempwith3[g1 - 1] + tempwith3[g2 - 1] +
			       tempwith3[g3 - 1] + tempwith3[g4 - 1]);
      }
      index++; /*increment isozygote class counter*/
    }
  }
}

/*segsumdown2 is used in segdown to compute some common subexpressions in
the probability updates.
first and second are the joint genotypes of the two parents.
fslength is the number of probabilities needed for the combined
isozygote class for both parental genotypes; fs stands for the
product of first and second. LINK is used to pass the genetic information for
parents and children */

/* cgh -- made this void for gcc */

Local void segsumdown2(first,second,fslength,LINK)
struct LOC_seg *LINK;
{
  int g1, g2, g3, g4; /*indices to store gene numbers*/
  int f1, f2, s1, s2; /*indices to store haplotype numbers*/
  int index; /*counter for combined isozygote equivalence class,
                will range from 0 to fslength-1*/
  int index2; /*used as index for isozygote class*/
  int i, j, k, l; /*loop indices*/
  int FORLIM; /*limit on a for loop*/
  int FORLIM1, FORLIM2; /*lower and upper limits on a for loop*/
  int *TEMPGENE1, *TEMPGENE2; /*temporary pointers into genenumber table*/
  double *tempwith3; /*temporarily stores genarray for current child*/

  FORLIM = fence[first]; /*start of the isozygote class for first*/
  FORLIM1 = base[second]; /*start of the isozygote class for second*/
  FORLIM2 = fence[second]; /*end of the isozygote class for second*/
  index = 0; /*start with first isozygote in the class*/
  index2 = 0;
  /*iterate over isozygotes of i*/
  for (i = base[first]; i < FORLIM; i++) {
    f1 = haps1[i]; /*retrieve first haplotype of joint genotype i*/
    f2 = haps2[i]; /*retrieve second haplotype*/
    TEMPGENE1 = genenumber[f1 - 1]; /*lookup partial information
					    in genenumber table*/
    TEMPGENE2 = genenumber[f2 - 1];
    /*iterate over isozygotes of j*/
    for (j = FORLIM1; j < FORLIM2; j++) {
      s1 = haps1[j]; /*retrieve first haplotype of joint genotype j*/
      s2 = haps2[j]; /*retrieve second haplotype*/
     /*lookup the four ways to combine one haplotype from i and one from j*/
      g1 = TEMPGENE1[s1 - 1];
      g2 = TEMPGENE1[s2 - 1];
      g3 = TEMPGENE2[s1 - 1];
      g4 = TEMPGENE2[s2 - 1];
      /*store these gene numbers for later use; the
        key point is that we will want these numbers
        consecutively later, so we store them consecutively
        in segindex*/
      segindex[index2++] = g1;
      segindex[index2++] = g2;
      segindex[index2++] = g3;
      segindex[index2++] = g4; 

    /*iterate over children; update partial computations of
      new probabilities and store in array tempseg
      note that tempseg has nchild*fslength relevant entries
      The fslength entries for child 1 come first, then
      the fslength entries, for child 2, etc. This is why
      the increment on l is fslength for each change of child*/

      for (l = index, k = 0 ; k < nchild; l += fslength, k++) {
	tempwith3 = thischild[k]->genarray; /*retrieve genarray*/

       /* sum the probabilities for the four joint genotypes that
        the child can get from the current parental joint genotypes*/

	tempseg[l] = (tempwith3[g1 - 1] + tempwith3[g2 - 1] +
			       tempwith3[g3 - 1] + tempwith3[g4 - 1]);
      }
      index++; /*increment isozygote class counter*/
    }
  }
}

/*lsegfun2 does a logical test similar to the computation
done in the original segfun to determine whether segfun
would return 0.0. If segfun would return 0.0, then lsegfun2
return 0 (FALSE), while if segfun would not return 0.0.,
lsegfun2 returns 1 (TRUE). Given, a combined isozygote class,
we want to know whether any elements of that isozygote
class are possible joint genotypes for the parents.
This will be the case if and only if each child has
a nonzero probability for at least one of the joint genotypes
that the parents can produce given their isozygote class.
first and second are the joint genotypes for the parents.
LINK is used to pass the genetic data. */


unsigned lsegfun2(first,second,LINK)
struct LOC_seg *LINK;
{
  int g1, g2, g3, g4; /*four gene numbers*/
  int i, j, k; /* loop indices*/
  int f1, f2, s1, s2; /*haplotype numbers*/
  int FORLIM; /*loop bound*/
  int FORLIM1, FORLIM2; /*loop bounds*/
  int *TEMPGENE1, *TEMPGENE2; /* store pointers into genenumber*/
  unsigned char *tempflag3; /*stores sparsity pattern for child's genarray*/

  FORLIM = fence[first]; /*find end of isozygote class for first*/
  FORLIM1 = base[second];/*find beginning and end of isozygote class
			   for second*/
  FORLIM2 = fence[second];
/*try to find a non-zero value for each child*/
for (k = 0; k < nchild; k++) {
  /*code for the non-boolean version is shown in comments
  tempwith3 = thischild[k]->genarray;*/
  /*retrieve sparsity pattern for child k*/
  tempflag3 = thischild[k]->sparseflag;
  
  /*iterate over all recombined isozygotes of first*/
  for (i = base[first]; i < FORLIM; i++) {
    f1 = haps1[i]; /*retrieve the haplotypes of this genotype*/
    f2 = haps2[i];
    TEMPGENE1 = genenumber[f1 - 1];/*get pointer into genenumber for
					   this haplotype*/
    TEMPGENE2 = genenumber[f2 - 1];

    /*iterate over all the recombined isozygotes of second*/
    for (j = FORLIM1; j < FORLIM2; j++) {
      s1 = haps1[j]; /*get haplotypes of this genotype*/
      s2 = haps2[j];
/*retrieve the four genes that this combination of joint haplotypes
  can produce*/
      g1 = TEMPGENE1[s1 - 1];
      g2 = TEMPGENE1[s2 - 1];
      g3 = TEMPGENE2[s1 - 1];
      g4 = TEMPGENE2[s2 - 1];
      /*if(tempwith3[g1 - 1] != 0.0 || tempwith3[g2 - 1] != 0.0 ||
	 tempwith3[g3 - 1] != 0.0 || tempwith3[g4 - 1] != 0.0) goto notzero;*/

   /*if any of the flags is TRUE, then this child can have this genotype,
     and we move immediately to testing the next child*/
      if(tempflag3[g1 - 1] != 0 || tempflag3[g2 - 1] != 0 ||
	 tempflag3[g3 - 1] != 0 || tempflag3[g4 - 1] != 0) goto notzero;
    }
  }
  return 0; /*this child failed to have any of the genotypes for
	      the isozygote class, so this isozygote class is not possible*/
  notzero: continue;
}
return 1; /*all children passed the test and have a possible joint genotype*/
}



/*segup updates the genotype probabilities for a parent (p) based
on the probabilities of children. The parameter LINK includes
all the genetic information about p, p's spouse q, and their
children*/

Void segup(LINK)
struct LOC_seg *LINK;
{
  int  FORLIM; /*loop bound*/
  int nonzindex, nonzcount;    /*loop index and counter for nonzero values*/
  int  step1, step2; /*size of isozygote classes*/
  unsigned int *tempstart; /*temporary index to probability array*/
  static int current; /*index for haplotype pools*/
  double valtemp; /*intermediate values in probability updates*/
  double val, temp1, tempr; /*temporaries to store intermediate values*/
  unsigned int f1, f2; /*two haplotypes from parents*/
  unsigned int here, i, j, first, second; /*genotype indices*/
  unsigned int fslength; /*size of product isozygote class for p and q*/
  int jointisoindex, fisoindex; /*indices to work within isozygote classes*/
  unsigned l; /*increment used to manage offset into probtableindex
                l is the number of probability values needed for the current
		isozygote class. This is obtained by multplying the
		number of probabilities (maxneed) times the size of the
		class*/
  unsigned k; /*index to loop over children*/
  unsigned char skip; /*used to skip iterations of a loop*/
  thisarray *WITH2; /*stores genetic information about p*/
  thisarray *WITH3; /*stores gnetic information about q*/
  double *newwith2, *newwith3, *newwithr; /*store genarrays for
						       p,q, and children*/
  thetavalues *WITH4; /*stores theta values for male parent*/
  thetavalues *WITH5; /*store theta values for female parent*/

  unsigned int c1, c2; /*haplotypes*/
  unsigned char *newflag2, *newflag3, *newflagr; /*store sparsity patterns for
					p and q genarrays*/


  /* newsegprob is  used to hold a segprob
    array, which contains the probabilities for various patterns of
    recombination events*/
  double *newsegprob; 
  double *tempprob; /*temporary holder for probability array*/
  int ind;      /*used to store offset for probability array*/

  initseg(LINK); /*get data about this p,q,children triple*/

  /*clear out pools, debugging
  for(i=0; i<= current; i++){
    indpool[i]=0;
    invpool[i]=0;
    nextpool[i]=0;
  } */

  /*get sparsity patterns for p and q genarrays*/
  newflag2 = (*LINK->p)->gen->sparseflag; 
  newflag3 = (*LINK->q)->gen->sparseflag;

  WITH2 = (*LINK->p)->gen; /*get genetic data for p*/
  WITH3 = (*LINK->q)->gen; /*get genetic data for q*/
  newwith2=WITH2->genarray; /*get genarray for p*/
  newwith3=WITH3->genarray; /*get genarray for q*/
  newwithr=thischild[0]->genarray; /*get genarray for first child*/
  WITH4 = LINK->firstsex; /*get p probabilities*/
  WITH5 = LINK->secondsex; /*get q probabilities*/

/*The case of 1 child is handled specially because the subcomputations
are much simpler. In particular we do not need to multiply the
probabilities among the different children. In a typical pedigree
many pairs will have only one child about which information is known
(in essence this is the child that keeps the pedigree connected), so
this is an important and common special case. */

if(nchild == 1) {
/*initialize cache data structures*/
  for(i=0;i<maxhaplo;i++) {
    flag[i] = 0;
    phapcache1[i].first = 0;
    qsumcache[i] = 0.0;
  }

/*initialize gene array and set up flag array*/
  for(i = 0; i < fgeno; i++) {
    gene[i] = 0.0;
    if(newwithr[i] == 0.0) continue;
    flag[invgenenum1[i]-1] = 1;
    flag[invgenenum2[i]-1] = 1;
  }

  current = 1; /* current should start with 1, since 0 means a null ptr */

  for (first = 0; first < fgeno; first++) {
    if(newflag2[first] == 0) continue; /*use only possible genotypes*/
    FORLIM = fence[first]; /*find end of isozygote class of first*/

  /*iterate over all members of first's isozygote calss*/
    for (i = base[first]; i < FORLIM; i++) {
      f1 = haps1[i]; /*get haplotypes*/
      f2 = haps2[i];
      ind = hind[i];  /*get probability offset for i*/

   /*if f1 is a possible haplotype for child add the
     probability offset and value of first to the cache*/
      if(flag[f1-1] != 0) {
	indpool[current] = ind; /*store values in cache*/
	invpool[current] = first;
	nextpool[current] = 0;
        
       /* increment cache indices*/
	if(phapcache1[f1-1].first == 0) {  
	  phapcache1[f1-1].first = current;
	  phapcache1[f1-1].last = current;
	}
	else {
	  nextpool[phapcache1[f1-1].last] = current;
	  phapcache1[f1-1].last = current;
	}
	current++;
      }
 
      /*similar to f1, if f2 is a possible haplotype for
        a child, cache the ind and first values for it*/
      if(flag[f2-1] != 0) {
	indpool[current] = ind;
	invpool[current] = first;
	nextpool[current] = 0;
	if(phapcache1[f2-1].first == 0) {
	  phapcache1[f2-1].first = current;
	  phapcache1[f2-1].last = current;
	}
	else {
	  nextpool[phapcache1[f2-1].last] = current;
	  phapcache1[f2-1].last = current;
	}
	current++;
      }
    }
  }

  /*This section of the code precomputes for each haplotype the the probability
  that the child will inherit this haplotype from q. Each genotype has
  two haplotypes, but can produce others by recombination. Therefore,
  for each genotype we must sum over the different haplotypes that can
  be produced by its isozygote class. The contributions for each haplotype
  are stored in qsumcache.*/

  newsegprob = WITH5->segprob; /*use q probabilities to work with q*/

 /*iterate over all joint genotypes*/
  for (first = 0; first < fgeno; first++ ) {
    if(newflag3[first] == 0) continue; /*check if joint genotype is possible*/
    valtemp = newwith3[first]; /*get cond. prob. that q has genotype first*/
    FORLIM = fence[first]; /*find bounds for the isozygote class of first*/

 /*iterate over all the recombined genotypes of the isozygote class*/
    for (i = base[first]; i < FORLIM; i++) {
      f1 = haps1[i];
      f2 = haps2[i];
      if((flag[f1-1] != 0) || (flag[f2-1] != 0)) {
/*conditional probability of first as a genotype multiplied by probability
  of passing on f1 (alternatively f2) as a haplotype*/
        val = valtemp*newsegprob[hind[i]];

       /*store in qsumcache*/
        if(flag[f1-1] != 0) {
	  qsumcache[f1 - 1] += val;
        }
        if(flag[f2-1] != 0) {
	  qsumcache[f2 - 1] += val;
        }
      }
    }
  }


 /*In this section of the code we update the probabilities for
  the parent p based on the probabilities for the child*/  
  newsegprob = WITH4->segprob; /*use probabilities to work with p*/
 /*Iterate over all joint genotypes of the child*/
  newflagr=thischild[0]->sparseflag;
  for(here = 0; here < fgeno; here++) {
    if(newflagr[here] == 0) continue;
    tempr = newwithr[here];
    c1 = invgenenum1[here];
    c2 = invgenenum2[here];
    temp1 = tempr * qsumcache[c2 - 1]; /*multiply probabilities for
					 haplotype c2 from q*/

  /*now find all the values of first that as genotypes of the
 father could have passed on the other haplotype c1; multiply
  by their probabilties, which are stored in phapcache, and
  add to the entry gene[first], which will be used to form
 the new genarray for p*/

    for(i=phapcache1[c1-1].first; i != 0; i = nextpool[i]) {
      first = invpool[i];
      gene[first] += temp1*newsegprob[indpool[i]];
    }

  /*if c1 is distinct from c2 we need to do the same computation reversing
    the roles of the two haplotypes*/
    if(c1 != c2) {
      temp1 = tempr * qsumcache[c1 - 1];
      for(i=phapcache1[c2-1].first; i != 0; i = nextpool[i]) {
        first = invpool[i];
	gene[first] += temp1*newsegprob[indpool[i]];
      }
    }
  }

/*set up new genarray for p; it is gene scaled by segscale*/      
  for (first = 0; first < fgeno; first++) {
    if(newflag2[first] == 0.0) continue;
    newwith2[first] *= segscale*gene[first];
  }

}


else  {  /*nchild is bigger than 1*/
  newwith2=(*LINK->p)->gen->genarray; /*get p's genarray*/


/*find nonzero entries in q's genarray and make a list of them
  stored in nonzgens; just get one per isozygote class*/
  nonzcount = 0;
  newwith3=(*LINK->q)->gen->genarray;

  /*iterate over genotypes for q*/
  for(i = 0; i < fgeno; i += step2) {
    /*number of distinct probabilties needed for i's isoz. class*/
    step2 = probend[i] - probstart[i] + 1; 
    for(j = i; j < i+step2; j++)
      if(newflag3[j] != 0) {
        nonzgens[nonzcount++] = i; /*store index of nonzero value*/
	break;                        /*go to next isozygote class*/
      }
  }

  /*iterate over genotypes for p*/
  for (first = 0; first < fgeno; first += step1) {
      /*number of distinct probabilties needed for first's isoz. class*/
    step1 = probend[first] - probstart[first]+1; 
    skip = 1;

    /*work only on those isozygotes that are possible*/
    for(i = first; i < first+step1; i++)
      if(newflag2[i] != 0) {
	skip = 0;
	break; /*go to next isozygote in class*/
      }
    if(skip) continue;

    /*initialize update multiple for each isozygote in class*/
    for(i = 0; i < step1; i++) segval[i] = 0.0;

    /*iterate over the genotypes representing isozygote classes that
     q may have*/
    for (nonzindex = 0; nonzindex < nonzcount; nonzindex++) {
      second = nonzgens[nonzindex];
      /*check if this first, second pair yield a nonzero value
        among children's probabilities*/
      if(lsegfun2(first,second,LINK) == 0) continue;

      /*number of distinct probabilties needed for second's isoz. class*/
      step2 = probend[second] - probstart[second] + 1; 
      fslength = step1 * step2; /*number of probs. for combined isozygote class*/
      fisoindex = (probstart[first] - 1) * nuneed - 1; /*probability array offset*/
      tempstart = probstart + second; /*set tempstart to the start of
					the section of probstart that
					we want*/
      l = step1 * nuneed; 
  
   /*call segsum2 to compute the part of the conditional probability
     update that is common to all the members of the combined isozygote
     class defined by first and second */
      segsum2(first,second,fslength,LINK);

      /*now specialize update for each member of first's class*/
      for(i = 0; i < step1; i++) {
	if(newflag2[first+i] == 0) {
	  fisoindex += l; /*increment by number of probabilities
			   needed for this iso. class*/
	  continue; /*skip if this isozygote is not possible*/
	}
        /*further specialize update for each member of second's isozygote
          class*/
	for(j = 0; j < step2; j++) {
        /*skip if this isozygote not possible*/
	  if(newflag3[second+j] == 0) continue; 

          /*get offset into probtable; the offset depends on
            the isozygote class size and index of each parent
            note that fisozygoteindex gets incremented by the size
            of the joint class each time, so it represents a
            sum of all the numbers of distinct probabilities
            needed for all joint iso. classes  considered before the
            current p isozygote class*/
#if  AUTOSOMAL_RUN
          if ((!sexdif) || ((*LINK->p)->male))
            tempprob = probtable + probtableindex[fisoindex+tempstart[j]];
#if SEXDIF_RUN
          else
            tempprob = probtabledif + probtableindex[fisoindex+tempstart[j]];
#endif
#endif
        /*combine for all children*/
        /*due to the arrangement in segsum all the probability contributions
          for a given child are contiguous in the newwithr array.
          the number of contributions is fslength which is the number of
          probabilities needed for the joint isozygote class of the parents.
          We get the contribution of the first child (index 0) and
          then loop over the rest*/
	  val = tempprob[0] * tempseg[0];
         /*iterate over joint isozygote class*/
	  for(jointisoindex = 1; jointisoindex < fslength; jointisoindex++)
	    val += tempprob[jointisoindex] * tempseg[jointisoindex];
	  newwithr = tempseg + fslength; /*move the base of newwithr*/

         /*iterate over all the remaining children; notice that the
           loop index k is not used inside the loop. this is because the
           last loop statement which moves the offset of newwithr
           has the effect of incrementing the child*/
	  for(k = 1; k < nchild; k++) {
	    temp1 = tempprob[0] * newwithr[0];
	    for(jointisoindex = 1; jointisoindex < fslength; jointisoindex++)
	      temp1 += tempprob[jointisoindex] * newwithr[jointisoindex];
	    val *= temp1;
	    newwithr += fslength;
	  }
         /*update segval entry for this isozygote of first*/
	  segval[i] += newwith3[second+j] * val;
	}
	fisoindex += l; /*go to next isozygote for p*/
      }
    }
    /*update p's genarray for each isozygote in this class*/
    for(i = 0; i < step1; i++) newwith2[first+i] *= segval[i] * segscale;
  }
}

/* If any of the nonzero entries in p's genarray became 0,
   we want to set them to zero to avoid computations on subsequent
   calls*/

  for(i = 0; i < fgeno; i++)
    if((newflag2[i] != 0) && (newwith2[i] == 0.0)) newflag2[i] = 0;
  cleanup(LINK->q, LINK->LINK);
  exitseg(LINK);
}  /*segup*/

/* segdown updates the genotype probabilities for a child  based
on the probabilities of parents and other children. The parameter LINK includes
all the genetic information about p, p's spouse q, and their
children */

Void segdown(LINK)
struct LOC_seg *LINK;
{
  int  FORLIM; /*loop bound*/
  int nonzindex, nonzcount;    /*loop index and counter for nonzero values*/
  int  step1, step2; /*size of isozygote classes*/
  unsigned int *tempstart; /*temporary index to probability array*/
  double valtemp; /*intermediate values in probability updates*/
  double val, temp1; /*temporaries to store intermediate values*/
  unsigned int f1, f2;  /*two haplotypes from parent*/
  unsigned int here, i, j, first, second; /*genotype indices*/
  unsigned int fslength; /*number of probs. for product isoz. class of p,q*/
  int jointisoindex, fisoindex; /*indices to work within isozygote classes*/
  unsigned currentindex; /*index to update genarray within isozygote class*/
  unsigned l; /*increment used to manage offset into probtableindex
                l is the number of probability values needed for the current
		isozygote class. This is obtained by multplying the
		number of probabilities (maxneed) times the size of the
		class*/
  unsigned k; /*index to loop over children*/
  unsigned char skip; /*used to skip iterations of a loop*/
  thisarray *WITH2; /*stores genetic information about p*/
  thisarray *WITH3; /*stores gnetic information about q*/
  double *newwith2, *newwith3, *newwithr, *newwithc; /*store genarrays for
						       p,q, and children*/
  thetavalues *WITH4; /*stores theta values for p*/
  thetavalues *WITH5; /*store theta values for q*/

  unsigned int c1, c2; /*haplotypes*/
  unsigned char *newflag2, *newflag3, *newflagr; /*store sparsity patterns for
					p and q genarrays*/


  /* newsegprob is  used to hold a segprob
    array, which contains the probabilities for various patterns of
    recombination events*/
  double *newsegprob;
  double *tempprob; /*temporary holder for probability array*/
  int ind;      /*used to store offset for probability array*/
  boolean parunk;

  initseg(LINK); /*get data about this p,q,children triple*/

  /*get sparsity patterns for p, q, and child genarrays*/
  newflag2 = (*LINK->p)->gen->sparseflag; 
  newflag3 = (*LINK->q)->gen->sparseflag;
  newflagr = (*LINK->r)->gen->sparseflag;

  WITH2 = (*LINK->p)->gen; /*get genetic data for p*/
  WITH3 = (*LINK->q)->gen; /*get genetic data for q*/
  newwith2 = WITH2->genarray; /*get genarray for p*/
  newwith3 = WITH3->genarray; /*get genarray for q*/
  newwithr = (*LINK->r)->gen->genarray; /*get genarray for first child*/
  WITH4 = LINK->firstsex; /*get p probabilities*/
  WITH5 = LINK->secondsex; /*get q probabilities*/
/* ERROR:  This code corrected on  15 Nov 1996 
  parunk = (((*LINK->p)->geneloc == 0) || ((*LINK->p)->geneloc == 0) || !noloop);  */
  parunk = (((*LINK->p)->geneloc == 0) || ((*LINK->q)->geneloc == 0) || !noloop);

/*The case of 1 child (nchild==0) is handled specially because the
subcomputations are much simpler. In particular we do not need to multiply the
probabilities among the different children. In a typical pedigree
many pairs will have only one child about which information is known
(in essence this is the child that keeps the pedigree connected), so
this is an important and common special case. */

if((nchild == 0) || (!parunk)) {
/*initialize cache data structures*/
  for(i = 0;i < maxhaplo;i++) {
    flag[i] = 0;
    psumcache[i] = 0.0;
    qsumcache[i] = 0.0;
  }

/*initialize gene array and set up flag array*/
  for(i = 0; i < fgeno; i++) {
    gene[i] = 0.0;
    if(newwithr[i] == 0.0) continue;
    flag[invgenenum1[i]-1] = 1;
    flag[invgenenum2[i]-1] = 1;
  }

  /*This section of the code precomputes for each haplotype the the probability
  that the child will inherit this haplotype from p. Each genotype has
  two haplotypes, but can produce others by recombination. Therefore,
  for each genotype we must sum over the different haplotypes that can
  be produced by its isozygote class. The contributions for each haplotype
  are stored in psumcache.
  Afterwards a similar computation is done for the inheritance from q
  with the results stored in qsumcache.*/


  newsegprob = WITH4->segprob; /*get p probabilities for recomb. patterns*/
  for (first = 0; first < fgeno; first++) {
    if(newflag2[first] == 0) continue; /*use only possible genotypes*/
    FORLIM = fence[first]; /*find end of isozygote class of first*/
    valtemp = newwith2[first]; /*probability of getting this genotype*/

  /*iterate over all members of first's isozygote calss*/
    for (i = base[first]; i < FORLIM; i++) {
      f1 = haps1[i]; /*get haplotypes*/
      f2 = haps2[i];
      if ((flag[f1-1] !=0) || (flag[f2-1] != 0)) {
         ind = hind[i];  /*get probability offset for i*/
       /*multiply probability of getting genotype times
	 probability of this recombination pattern and haplo. choice*/
         val = valtemp * newsegprob[ind];

       /*add to psumcache*/
         if(flag[f1-1] != 0) {
	   psumcache[f1-1] += val;
          }
 
         if(flag[f2-1] != 0) {
	   psumcache[f2-1] += val;
         }
       }
    }
  }

  newsegprob = WITH5->segprob; /*get q probabilities*/

 /*iterate over all joint genotypes*/
  for (first = 0; first < fgeno; first++ ) {
    if(newflag3[first] == 0) continue; /*check if joint genotype is possible*/
    valtemp = newwith3[first]; /*get cond. prob. that q has genotype first*/
    FORLIM = fence[first]; /*find bounds for the isozygote class of first*/

 /*iterate over all the recombined genotypes of the isozygote class*/
    for (i = base[first]; i < FORLIM; i++) {
      f1 = haps1[i];
      f2 = haps2[i];

      if((flag[f1-1] !=0) || (flag[f2-1] !=0)) {

/*condition probability of first as a genotype multiplied by probability
  of passing on f1 (alternatively f2) as a haplotype*/
        ind = hind[i];  /*get probability offset for i*/
        val = valtemp*newsegprob[hind[i]];
     /*store in qsumcache*/
        if(flag[f1-1] != 0) {
	  qsumcache[f1 - 1] += val;
        }
        if(flag[f2-1] != 0) {
	  qsumcache[f2 - 1] += val;
        }
      }
    }
  }


 /*In this section of the code we update the probabilities for
  the child based on the probabilities for the parents*/  

 /*Iterate over all joint genotypes of the child*/

  for(here = 0; here < fgeno; here++) {
    if(newflagr[here] == 0) continue;
    c1 = invgenenum1[here];
    c2 = invgenenum2[here];

   /*probability of child getting genotype here as a
     result of p passing on c1 and q passing on c2 is
     summed to gene[here] */

    gene[here] += psumcache[c1-1] * qsumcache[c2-1];

  /*if c1 is distinct from c2 we need to do the same computation reversing
    the roles of the two haplotypes*/
    if(c1 != c2) {
      gene[here] += psumcache[c2-1] * qsumcache[c1-1];
    }
  }

/*set up new genarray for r; it is gene scaled by segscale*/      
  for (first = 0; first < fgeno; first++) {
    if(newflagr[first] == 0) continue;
    if(gene[first] == 0.0) newflagr[first] = 0;
    newwithr[first] *= segscale*gene[first];
  }

}


else  {  /*nchild is bigger than 0*/
  newwith2=(*LINK->p)->gen->genarray; /*get p's genarray*/

 /*initialize genarray entries for child to 0*/
  for(i = 0; i < fgeno; i++) {
    gene[i] = 0.0;
  }

/*find nonzero entries in q's genarray and make a list of them
  stored in nonzgens; just get one per isozygote class*/
  nonzcount = 0;
  newwith3=(*LINK->q)->gen->genarray;

  /*iterate over genotypes for q*/
  for(i = 0; i < fgeno; i += step2) {
    /*number of distinct probabilties needed for i's isoz. class*/
    step2 = probend[i] - probstart[i] + 1; 
    for(j = i; j < i+step2; j++)
      if(newflag3[j] != 0) {
        nonzgens[nonzcount++] = i; /*store index of nonzero value*/
	break;                        /*go to next isozygote class*/
      }
  }

   /*iterate over genotypes for p*/
  for (first = 0; first < fgeno; first += step1) {
      /*number of distinct probabilties needed for first's isoz. class*/
    step1 = probend[first] - probstart[first]+1; 
    skip = 1;

    /*work only on those isozygotes that are possible*/
    for(i = first; i < first+step1; i++)
      if(newflag2[i] != 0) {
	skip = 0;
	break; /*go to next isozygote in class*/
      }
    if(skip) continue;

    /*iterate over the genotypes representing isozygote classes that
     q may have*/
    for (nonzindex = 0; nonzindex < nonzcount; nonzindex++) {
      second = nonzgens[nonzindex];
      /*check if this first, second pair yield a nonzero value
        among children's probabilities*/
      if(lsegfun2(first,second,LINK) == 0) continue;

      /*number of distinct probabilties needed for second's isoz. class*/
      step2 = probend[second] - probstart[second] + 1; 
      fslength = step1 * step2; /*number of probs. of combined isozygote class*/
      fisoindex = (probstart[first] - 1) * nuneed - 1; /*probability array offset*/
      tempstart = probstart + second; /*set tempstart to the start of
					the section of probstart that
					we want*/
      l = step1 * nuneed; /*number of probabilities for this
			     isozygote class*/
  
   /*call segsumdown2 to compute the part of the conditional probability
     update that is common to all the members of the combined isozygote
     class defined by first and second */
      segsumdown2(first,second,fslength,LINK);

      for(jointisoindex = 0; jointisoindex < fslength; jointisoindex++)
	segval[jointisoindex] = 0.0;

      /*now specialize update for each member of first's class*/
      for(i = 0; i < step1; i++) {
	if(newflag2[first+i] == 0) {
	  fisoindex += l; /*increment by number of probabilities
			   needed for this iso. class*/
	  continue; /*skip if this isozygote is not possible*/
	}
        /*further specialize update for each member of second's isozygote
          class*/
	for(j = 0; j < step2; j++) {
        /*skip if this isozygote not possible*/
	  if(newflag3[second+j] == 0) {
	    continue; 
	  }

          /*get offset into probtable; the offset depends on
            the isozygote class size and index of each parent
            note that fisozygoteindex gets incremented by the size
            of the joint class each time, so it represents a
            sum of all the numbers of distinct probabilities
            needed for all joint iso. classes  considered before the
            current p isozygote class*/
#if AUTOSOMAL_RUN
          if ((!sexdif) || ((*LINK->p)->male))
            tempprob = probtable + probtableindex[fisoindex+tempstart[j]];
#if SEXDIF_RUN
          else
            tempprob = probtabledif + probtableindex[fisoindex+tempstart[j]];
#endif
#endif
        /*combine for all children*/
        /*due to the arrangement in segsum all the probability contributions
          for a given child are contiguous in the newwithr array.
          the number of contributions is fslength which is the number of
          probabilties needed for the joint isozygote class for the parents.
          We get the contribution of the first child (index 0) and
          then loop over the rest*/
	  val = tempprob[0] * tempseg[0];
         /*iterate over joint isozygote class*/
	  for(jointisoindex = 1; jointisoindex < fslength; jointisoindex++)
	    val += tempprob[jointisoindex] * tempseg[jointisoindex];
	  newwithc = tempseg + fslength; /*move the base of newwithc*/

         /*iterate over all the remaining children; notice that the
           loop index k is not used inside the loop. this is because the
           last loop statement which moves the offset of newwithc
           has the effect of incrementing the child*/
	  for(k = 1; k < nchild; k++) {
	    temp1 = tempprob[0] * newwithc[0];
	    for(jointisoindex = 1; jointisoindex < fslength; jointisoindex++)
	      temp1 += tempprob[jointisoindex] * newwithc[jointisoindex];
	    val *= temp1;
	    newwithc += fslength;
	  }
          /*probability of this combination of parent genotypes*/
	  val *= newwith2[first+i] * newwith3[second+j];
	  for(jointisoindex = 0; jointisoindex < fslength; jointisoindex++) {
          /*probability of this recombination pattern (based on other children)*/
	    segval[jointisoindex] += tempprob[jointisoindex] * val;
	  }
	}
	fisoindex += l; /*go to next isozygote for p*/
      }
      /*update the probabilities of four joint genotypes the
       child might get; each different choice of recombination
       pattern will lead to a different set of four genotypes*/
      currentindex = 0;
      for(jointisoindex = 0; jointisoindex < fslength; jointisoindex++) {
	temp1 = segval[jointisoindex];
	gene[segindex[currentindex++]-1] += temp1;
	gene[segindex[currentindex++]-1] += temp1;
	gene[segindex[currentindex++]-1] += temp1;
	gene[segindex[currentindex++]-1] += temp1;
      }
    }
  }
  /*finally update child's real genarray by coppy gene multiplied
    by scale factor segscale*/
  for(i = 0; i < fgeno; i++) {
    if(newflagr[i] == 0) continue;
    if(gene[i] == 0.0) newflagr[i] = 0; /*if probability changes from nonzero
					  to 0.0 change flag to 0*/
    newwithr[i] *= segscale * gene[i];
  }
}

  cleanup(LINK->p, LINK->LINK);
  cleanup(LINK->q, LINK->LINK);
  exitseg(LINK);
}  /*segdown*/


