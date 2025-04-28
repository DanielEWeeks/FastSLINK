/*This file is part of a C version of SLINK a simulation program based
on the LINKAGE package. SLINK was originally described in: D. E.
Weeks, J. Ott, G. M. Lathrop (1990), SLINK: a general simulation
program for linkage analysis, American Journal of Human Genetics
47(1990), A204 (abstract). This version of SLINK is adapted to use the
faster pedigree traversal algorithms described in: R. W. Cottingham
Jr., R. M. Idury, A. A. Schaffer (1993), Faster sequential genetic
linkage computations, American Journal of Human Genetics 53(1993), pp.
252-263.*/
/*This file contains the original C versions of segup*/
/* and segsexup also used in LODSCORE, ILINK, and LINKMAP*/
/* They are separated out so as to have the same file structure*/
/* as is used for the new versions*/

void oldsegsexup(LINK)
struct LOC_seg *LINK;
{
  double segval;
  int first, second;
  censorrec *WITH;
  thisarray *WITH1;
  int FORLIM;
  thisarray *WITH2;
  int FORLIM1;

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
  cleanup(LINK->q, LINK->LINK);
  exitseg(LINK);
}  /*segsexup*/

void oldsegup(LINK)
struct LOC_seg *LINK;
{
  double segval;
  int first, second;
  censorrec *WITH;
  thisarray *WITH1;
  int FORLIM;
  thisarray *WITH2;
  int FORLIM1;

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
  cleanup(LINK->q, LINK->LINK);
  exitseg(LINK);
}  /*segup*/

