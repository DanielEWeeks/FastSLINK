/*This file is part of a C version of SLINK a simulation program based
on the LINKAGE package. SLINK was originally described in: D. E.
Weeks, J. Ott, G. M. Lathrop (1990), SLINK: a general simulation
program for linkage analysis, American Journal of Human Genetics
47(1990), A204 (abstract). This version of SLINK is adapted to use the
faster pedigree traversal algorithms described in: R. W. Cottingham
Jr., R. M. Idury, A. A. Schaffer (1993), Faster sequential genetic
linkage computations, American Journal of Human Genetics 53(1993), pp.
252-263.*/
/*This file  contains some low level probability routines*/
/*8 December 1993 sliced out of slink.c by A. A. Schaffer */
/* May 2008 Edited to allow > 32 alleles at a numbered allele locus*/

/* Local variables for setval: */
struct LOC_setval {
  struct LOC_getvect *LINK;
  double val;
  int nhap1, nhap2;
} ;

static void prior(LINK)
struct LOC_setval *LINK;
{
  int i;   /*prior*/
  int FORLIM;
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

static void setval(val_, LINK)
double val_;
struct LOC_getvect *LINK;
{
  struct LOC_setval V;
  int here, count, i, FORLIM;
  thisarray *WITH1;

  /*prior*/
  /*setval*/

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
    FORLIM = nlocus;
    for (i = 0; i < FORLIM; i++) {
      V.nhap1 += increment[i] * (LINK->hap1[i] - 1);
      V.nhap2 += increment[i] * (LINK->hap2[i] - 1);
      if (LINK->hap1[i] != LINK->hap2[i])
	count *= 2;
      here = genenumber[V.nhap1 - 1][V.nhap2 - 1];
    }
  }
  if (LINK->p->pa == NULL)
    prior(&V);
  WITH1 = LINK->p->gen;
  if (disequi) {
    WITH1->genarray[here - 1] = V.val;
    WITH1->sparseflag[here - 1] = 1; /*R. M. Idury, A. A. Schaffer*/
    return;
  }
  if (count != 1)
    count /= 2;
  for (i = 1; i <= count; i++) {
    WITH1->genarray[here - 1] = V.val;
    WITH1->sparseflag[here - 1] = 1; /*R. M. Idury, A. A. Schaffer*/
    here++;
  }
}  /* setval */

/* Local variables for getgene: */
struct LOC_getgene {
  struct LOC_getvect *LINK;
  int syste;
  double val;
  double newval;
} ;

static void facmale_getgene_bin(LINK)
     struct LOC_getgene *LINK;
{
  int j;   
  thisperson *WITH;
  locusvalues *WITH1;
  int FORLIM;

  WITH = LINK->LINK->p;
  WITH1 = thislocus[LINK->syste - 1];
  FORLIM = WITH1->nallele;
  for (j = 1; j <= FORLIM; j++) {
    if (WITH->phen[LINK->syste - 1]->phenf == WITH1->UU.U2.bin_allele[j - 1] ||
        WITH->phen[LINK->syste - 1]->phenf == 0) {
      LINK->LINK->hap1[LINK->syste - 1] = j;
      if (LINK->syste != nlocus)
        getgene(LINK->syste + 1, LINK->val, LINK->LINK);
      else
        setval(LINK->val, LINK->LINK);
    }
  }
}  /* facmale_getgene_bin */


static void facmale_getgene_num(LINK)
struct LOC_getgene *LINK;
{
  int j;   
  thisperson *WITH;
  locusvalues *WITH1;
  int FORLIM;

  WITH = LINK->LINK->p;
  WITH1 = thislocus[LINK->syste - 1];
  FORLIM = WITH1->nallele;
  for (j = 1; j <= FORLIM; j++) {
    if ((WITH->phen[LINK->syste - 1]->alleles[0] == WITH1->UU.U2.num_allele[j - 1][0] &&
         WITH->phen[LINK->syste - 1]->alleles[1] == WITH1->UU.U2.num_allele[j - 1][0]) ||
	(WITH->phen[LINK->syste - 1]->alleles[0] == 0 && 
         WITH->phen[LINK->syste - 1]->alleles[1] == 0)) { 
      LINK->LINK->hap1[LINK->syste - 1] = j;
      if (LINK->syste != nlocus)
	getgene(LINK->syste + 1, LINK->val, LINK->LINK);
      else
	setval(LINK->val, LINK->LINK);
    }
  }
}  /* facmale_getgene_num */



static void affmale(LINK)
struct LOC_getgene *LINK;
{
  int j;   /*affmale*/
  locusvalues *WITH;
  int FORLIM;

  WITH = thislocus[LINK->syste - 1];
  FORLIM = WITH->nallele;
  for (j = 1; j <= FORLIM; j++) {
    LINK->newval = LINK->val;
    getval(LINK->syste, 0, j, &LINK->newval, LINK->LINK);
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

static void quanmale(LINK)
struct LOC_getgene *LINK;
{
  int j;   /*quanmale*/
  thisperson *WITH;
  locusvalues *WITH1;
  int FORLIM;

  WITH = LINK->LINK->p;
  WITH1 = thislocus[LINK->syste - 1];
  if (WITH->phen[LINK->syste - 1]->aff == affall ||
      WITH->phen[LINK->syste - 1]->aff == missaff) {
    LINK->newval = LINK->val;
    LINK->LINK->hap1[LINK->syste - 1] = affall;
    if (LINK->newval != 0.0) {
      if (LINK->syste != nlocus)
	getgene(LINK->syste + 1, LINK->newval, LINK->LINK);
      else
	setval(LINK->newval, LINK->LINK);
    }
  }
  if (WITH->phen[LINK->syste - 1]->aff == affall &&
      WITH->phen[LINK->syste - 1]->aff != missaff)
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

static void fac_getgene_bin(LINK)
struct LOC_getgene *LINK;
{
  int i, j;   /*fac*/
  thisperson *WITH;
  locusvalues *WITH1;
  int FORLIM, FORLIM1;

  WITH = LINK->LINK->p;
  WITH1 = thislocus[LINK->syste - 1];
  FORLIM = WITH1->nallele;
  for (i = 1; i <= FORLIM; i++) {
    LINK->LINK->hap1[LINK->syste - 1] = i;
    FORLIM1 = WITH1->nallele;
    for (j = i; j <= FORLIM1; j++) {
      if (WITH->phen[LINK->syste - 1]->phenf ==
	  (WITH1->UU.U2.bin_allele[i - 1] | WITH1->UU.U2.bin_allele[j - 1]) ||
	  WITH->phen[LINK->syste - 1]->phenf == 0) {
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
      if (WITH->phen[LINK->syste - 1]->phenf ==
	  (WITH1->UU.U2.bin_allele[i - 1] | WITH1->UU.U2.bin_allele[j - 1]) ||
	  WITH->phen[LINK->syste - 1]->phenf == 0) {
	LINK->LINK->hap1[LINK->syste - 1] = j;
	if (LINK->syste != nlocus)
	  getgene(LINK->syste + 1, LINK->val, LINK->LINK);
	else
	  setval(LINK->val, LINK->LINK);
      }
    }
  }
}  /* fac_getgene_bin */

static void fac_getgene_num(LINK)
struct LOC_getgene *LINK;
{
  int i, j;   /*fac*/
  thisperson *WITH;
  locusvalues *WITH1;
  int FORLIM, FORLIM1;

  WITH = LINK->LINK->p;
  WITH1 = thislocus[LINK->syste - 1];
  FORLIM = WITH1->nallele;
  for (i = 1; i <= FORLIM; i++) {
    LINK->LINK->hap1[LINK->syste - 1] = i;
    FORLIM1 = WITH1->nallele;
    for (j = i; j <= FORLIM1; j++) {
      if (((WITH->phen[LINK->syste - 1]->alleles[0] ==
	     WITH1->UU.U2.num_allele[i - 1][0]) &&
            (WITH->phen[LINK->syste - 1]->alleles[1] ==
	     WITH1->UU.U2.num_allele[j - 1][0])) ||
	  ((WITH->phen[LINK->syste - 1]->alleles[1] ==
	    WITH1->UU.U2.num_allele[i - 1][0]) &&
	   (WITH->phen[LINK->syste - 1]->alleles[0] ==
	    WITH1->UU.U2.num_allele[j - 1][0])) ||
          ((WITH->phen[LINK->syste - 1]->alleles[0] == 0) &&
           (WITH->phen[LINK->syste - 1]->alleles[1] == 0))) {
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
      if (((WITH->phen[LINK->syste - 1]->alleles[0] ==
	     WITH1->UU.U2.num_allele[i - 1][0]) &&
            (WITH->phen[LINK->syste - 1]->alleles[1] ==
	     WITH1->UU.U2.num_allele[j - 1][0])) ||
	  ((WITH->phen[LINK->syste - 1]->alleles[1] ==
	    WITH1->UU.U2.num_allele[i - 1][0]) &&
	   (WITH->phen[LINK->syste - 1]->alleles[0] ==
	    WITH1->UU.U2.num_allele[j - 1][0])) ||
          ((WITH->phen[LINK->syste - 1]->alleles[0] == 0) &&
           (WITH->phen[LINK->syste - 1]->alleles[1] == 0))) {
	LINK->LINK->hap1[LINK->syste - 1] = j;
	if (LINK->syste != nlocus)
	  getgene(LINK->syste + 1, LINK->val, LINK->LINK);
	else
	  setval(LINK->val, LINK->LINK);
      }
    }
  }
}  /* fac_getgene_num */



static void aff(LINK)
struct LOC_getgene *LINK;
{
  /*Used with an affection status phenotype or when
  thislocus[syste]^which is null*/
  int i, j;
  locusvalues *WITH;
  int FORLIM, FORLIM1;

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


static void quanval(LINK)
struct LOC_getgene *LINK;
{
  /*Uses this only when thislocus[syste]^.which is not null*/
  int i, j;
  locusvalues *WITH;
  int FORLIM, FORLIM1;

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

void getgene(syste_, val_, LINK)
int syste_;
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
      if (WITH->format == 3)   /* This locus is an allele numbered locus */
        facmale_getgene_num(&V);
      else /* not numbered alleles */
        facmale_getgene_bin(&V);
      break;
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
    if (WITH->format == 3)   /* This locus is an allele numbered locus */
      fac_getgene_num(&V);
    else /* not numbered alleles */
      fac_getgene_bin(&V);
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
  int syste;
  double val;
  double newval;
} ;

static void facmale_ugetgene_bin(LINK)
struct LOC_ugetgene *LINK;
{
  int j;   
  thisperson *WITH;
  information *WITH1;
  locusvalues *WITH2;
  int FORLIM;

  WITH = LINK->LINK->p;
  WITH1 = LINK->LINK->p->store;
  WITH2 = thislocus[LINK->syste - 1];
  FORLIM = WITH2->nallele;
  for (j = 1; j <= FORLIM; j++) {
    if (WITH->phen[LINK->syste - 1]->phenf == WITH2->UU.U2.bin_allele[j - 1] ||
	WITH->phen[LINK->syste - 1]->phenf == 0) {
      if (WITH1->possible[LINK->syste - 1][0][j - 1]) {
	LINK->LINK->hap1[LINK->syste - 1] = j;
	if (LINK->syste != nlocus)
	  ugetgene(LINK->syste + 1, LINK->val, LINK->LINK);
	else
	  setval(LINK->val, LINK->LINK);
      }
    }
  }
}  /* facmale_ugetgene_bin */



static void facmale_ugetgene_num(LINK)
struct LOC_ugetgene *LINK;
{
  int j;   
  thisperson *WITH;
  information *WITH1;
  locusvalues *WITH2;
  int FORLIM;

  WITH = LINK->LINK->p;
  WITH1 = LINK->LINK->p->store;
  WITH2 = thislocus[LINK->syste - 1];
  FORLIM = WITH2->nallele;
  for (j = 1; j <= FORLIM; j++) {
    if ((WITH->phen[LINK->syste - 1]->alleles[0] == WITH2->UU.U2.num_allele[j - 1][0] &&
         WITH->phen[LINK->syste - 1]->alleles[1] == WITH2->UU.U2.num_allele[j - 1][0]) ||
	(WITH->phen[LINK->syste - 1]->alleles[0] == 0 && 
         WITH->phen[LINK->syste - 1]->alleles[1] == 0)) { 
      if (WITH1->possible[LINK->syste - 1][0][j - 1]) {
	LINK->LINK->hap1[LINK->syste - 1] = j;
	if (LINK->syste != nlocus)
	  ugetgene(LINK->syste + 1, LINK->val, LINK->LINK);
	else
	  setval(LINK->val, LINK->LINK);
      }
    }
  }
}  /* facmale_ugetgene_num */



static void affmale_(LINK)
struct LOC_ugetgene *LINK;
{
  int j;
  information *WITH;
  locusvalues *WITH1;
  int FORLIM;

  WITH = LINK->LINK->p->store;
  WITH1 = thislocus[LINK->syste - 1];
  FORLIM = WITH1->nallele;
  for (j = 1; j <= FORLIM; j++) {
    if (WITH->possible[LINK->syste - 1][0][j - 1]) {
      LINK->newval = LINK->val;
      getval(LINK->syste, 0, j, &LINK->newval, LINK->LINK);
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


static void quanmale_(LINK)
struct LOC_ugetgene *LINK;
{
  int j;
  thisperson *WITH;
  information *WITH1;
  locusvalues *WITH2;
  information *WITH3;
  int FORLIM;

  WITH = LINK->LINK->p;
  WITH1 = LINK->LINK->p->store;
  WITH2 = thislocus[LINK->syste - 1];
  if (WITH->phen[LINK->syste - 1]->aff == affall ||
      WITH->phen[LINK->syste - 1]->aff == missaff) {
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
  if (WITH->phen[LINK->syste - 1]->aff == affall &&
      WITH->phen[LINK->syste - 1]->aff != missaff)
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


static void fac_ugetgene_bin(LINK)
struct LOC_ugetgene *LINK;
{
  int i, j;
  thisperson *WITH;
  information *WITH1;
  locusvalues *WITH2;
  int FORLIM, FORLIM1;

  WITH = LINK->LINK->p;
  WITH1 = LINK->LINK->p->store;
  WITH2 = thislocus[LINK->syste - 1];
  FORLIM = WITH2->nallele;
  for (i = 1; i <= FORLIM; i++) {
    LINK->LINK->hap1[LINK->syste - 1] = i;
    FORLIM1 = WITH2->nallele;
    for (j = i; j <= FORLIM1; j++) {
      if (WITH->phen[LINK->syste - 1]->phenf ==
	  (WITH2->UU.U2.bin_allele[i - 1] | WITH2->UU.U2.bin_allele[j - 1]) ||
	  WITH->phen[LINK->syste - 1]->phenf == 0) {
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
      if (WITH->phen[LINK->syste - 1]->phenf ==
	  (WITH2->UU.U2.bin_allele[i - 1] | WITH2->UU.U2.bin_allele[j - 1]) ||
	  WITH->phen[LINK->syste - 1]->phenf == 0) {
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
}  /* fac_ugetgene_bin */

static void fac_ugetgene_num(LINK)
struct LOC_ugetgene *LINK;
{
  int i, j;
  thisperson *WITH;
  information *WITH1;
  locusvalues *WITH2;
  int FORLIM, FORLIM1;

  WITH = LINK->LINK->p;
  WITH1 = LINK->LINK->p->store;
  WITH2 = thislocus[LINK->syste - 1];
  FORLIM = WITH2->nallele;
  for (i = 1; i <= FORLIM; i++) {
    LINK->LINK->hap1[LINK->syste - 1] = i;
    FORLIM1 = WITH2->nallele;
    for (j = i; j <= FORLIM1; j++) {
      if (((WITH->phen[LINK->syste - 1]->alleles[0] ==
	     WITH2->UU.U2.num_allele[i - 1][0]) &&
            (WITH->phen[LINK->syste - 1]->alleles[1] ==
	     WITH2->UU.U2.num_allele[j - 1][0])) ||
	  ((WITH->phen[LINK->syste - 1]->alleles[1] ==
	    WITH2->UU.U2.num_allele[i - 1][0]) &&
	   (WITH->phen[LINK->syste - 1]->alleles[0] ==
	    WITH2->UU.U2.num_allele[j - 1][0])) ||
          ((WITH->phen[LINK->syste - 1]->alleles[0] == 0) &&
           (WITH->phen[LINK->syste - 1]->alleles[1] == 0))) {
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
      if (((WITH->phen[LINK->syste - 1]->alleles[0] ==
	     WITH2->UU.U2.num_allele[i - 1][0]) &&
            (WITH->phen[LINK->syste - 1]->alleles[1] ==
	     WITH2->UU.U2.num_allele[j - 1][0])) ||
	  ((WITH->phen[LINK->syste - 1]->alleles[1] ==
	    WITH2->UU.U2.num_allele[i - 1][0]) &&
	   (WITH->phen[LINK->syste - 1]->alleles[0] ==
	    WITH2->UU.U2.num_allele[j - 1][0])) ||
          ((WITH->phen[LINK->syste - 1]->alleles[0] == 0) &&
           (WITH->phen[LINK->syste - 1]->alleles[1] == 0))) {
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
}  /* fac_ugetgene_num */


static void aff_(LINK)
struct LOC_ugetgene *LINK;
{
  /*Used with an affection status phenotype or when
  thislocus[syste]^which is null*/
  int i, j;
  thisperson *WITH;
  information *WITH1;
  locusvalues *WITH2;
  int FORLIM, FORLIM1;

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


static void quanval_(LINK)
struct LOC_ugetgene *LINK;
{
  /*Uses this only when thislocus[syste]^.which is not null*/
  int i, j;
  thisperson *WITH;
  information *WITH1;
  locusvalues *WITH2;
  int FORLIM, FORLIM1;

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


void ugetgene(syste_, val_, LINK)
int syste_;
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
      if (WITH->format == 3)   /* This locus is an allele numbered locus */
        facmale_ugetgene_num(&V);
      else /* not numbered alleles */
        facmale_ugetgene_bin(&V);
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
    if (WITH->format == 3)   /* This locus is an allele numbered locus */
      fac_ugetgene_num(&V);
    else /* not numbered alleles */
      fac_ugetgene_bin(&V);
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


void getvect(p_, LINK)
thisperson *p_;
struct LOC_likelihood *LINK;
{
  struct LOC_getvect V;


  V.LINK = LINK;
  V.p = p_;
  /*When a genotype is assigned to the person p, then geneloc is non-zero*/
  if (V.p->geneloc != 0) {
    /*Setting this element of genarray to 1 is equivalent to fixing the
    genotype of person p to be the one in the position geneloc.*/
    V.p->gen->genarray[V.p->geneloc - 1] = 1.0;
    V.p->gen->sparseflag[V.p->geneloc - 1] = 1;
    return;
  }
  if (V.p->unknown)
    ugetgene(1, 1.0, &V);
  else
    getgene(1, 1.0, &V);
}  /* getvect */

/*getvect*/


