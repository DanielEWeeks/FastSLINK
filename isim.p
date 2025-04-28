program isim(input, output, isim, outfile, datafile, speedfile, ipedfile, simout, final);
{3 July 1992}

{ISIM is a modified version of ILINK, modified by Daniel E. Weeks with
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
 }


const
    version = '2.51'; {VERSION OF ISIM, based on version 4.91 of LINKAGE}
    prn = false;				{ Turns off most outfile output for more efficient running }
    finalprn = false; {Turns off output to final.dat for more efficient running }
    { SOME USER DEFINED CONSTANTS }
    {THE PROGRAM WILL TELL YOU IF THE FOLLOWING TWO CONSTANTS CAN BE REDUCED}
    {IF THE PROGRAM TERMINATES WITH AN ERROR IN RECOMBINATION INCREASE MAXNEED}
    maxneed = 800;				{MAXIMUM NUMBER OF RECOMBINATION PROBABILITIES}
    {THE FOLLOWING SHOULD BE LARGER THAN MININT}
    maxcensor = 10000; {MAXIMUM FOR CENSORING ARRAY}
    maxlocus = 4;				{MAXIMUM NUMBER OF LOCI }
    maxseg = 8;				{2 TO THE POWER maxlocus-1 }
    maxall = 6;				{MAX NUMBER OF ALLELES AT A SINGLE LOCUS}
    maxhap = 16;				{MAXIMUM NUMBER OF HAPLOTYPES}
    { = n1 x n2 x ... ETC. WHERE ni = NUMBER OF ALLELES LOCUS I}
    maxfem = maxhap * (maxhap + 1) div 2; {MAX. NO. OF JOINT GENOTYPES FOR A FEMALE}
    maxmal = maxfem;				{MAXIMUM NUMBER OF JOINT GENOTYPES FOR A MALE}
    { = maxfem (AUTOSOMAL) OR maxhap (SEXLINKED)}
    maxind = 600;				{MAXIMUM NUMBER OF INDIVIDUALS}
    maxped = 65;				{MAXIMUM NUMBER OF PEDIGREES}
    maxchild = 20;				{MAXIMUM NUMBER OF FULLSIBS IN A SIBSHIP}
    maxloop = 3;				{MAXIMUM NUMBER OF LOOPS PER PEDIGREE}
    minfreq = 0.0;
    affall = 2;				{DISEASE ALLELE FOR QUANTITATIVE TRAITS OR AFFECTION STATUS}
    { QUANTITATIVE TRAIT }
    maxtrait = 1;				{MAXIMUM NUMBER OF QUANTITATIVE FACTORS AT A SINGLE LOCUS}
    missval = 0.0;				{MISSING VALUES FOR QUANTITATIVE TRAITS }
    { AFFECTION STATUS }
    missaff = 0;				{MISSING VALUE FOR AFFECTION STATUS }
    affval = 2;				{CODE FOR AFFECTED INDIVIDUAL}
    maxliab = 20;				{MAXIMUM NUMBER OF LIABILITY CLASSES }
    { BINARY (FACTOR UNION) SYSTEM }
    maxfact = 10;				{MAXIMUM NUMBER OF BINARY CODES AT A SINGLE LOCUS}
    { OTHERS }
    scale = 3.0;				{SCALE FACTOR}
    scalemult = 2.0;				{SCALE WEIGHT FOR EACH LOCUS INCLUDED}
    fitmodel = false; {TRUE IF ESTIMATING PARAMETERS OTHER THAN REC}
    byfamily = true;				{GIVE LOD SCORES BY FAMILY}
    zerolike = -1.0e20; {FOR INCONSISTENT DATA OR RECOMBINATION } {SUN}
    log10 = 2.30259;
    minint = -32767;				{MINIMUM ALLOWED INTEGER}
    {GRADIENT APPROXIMATIONS}
    approximate = false; {do not change or else change approxrec below}
    epsilon = 1.0e-4; {SUN}
    {GEMINI}
    maxn = 20;				{ MAXIMUM NUMBER OF ITERATED PARAMETERS }
    nbit = 50;				{NUMBER OF BITS OF MACHINE PRECISION}
    tolconst = 1.0e-3; {TOLERANCE ON LIKELIHOOD FOR STOP} {SUN}
    {DO NOT CHANGE THE FOLLOWING CONSTANTS}
    tolg = 1.0e-35; {SUN}
    xpmcon = 1.0e30; {SUN}
    clbcon = 1.0e10; {SUN}
    tbndcon = 1.0e5; {SUN}
    curv = 0.75;
    small = 1.0e-10; {SUN}

type
    gennuptr = ^ gennurec;
    gennurec = 
	record 
	    genenumber: array [1..maxhap, 1..maxhap] of integer
	end;
    approxpnt = ^ approxrec;
    approxrec = 
	record 
	    approxarray: array [1..1, 1..1] of boolean
	end; {Was packed}
    censorpnt = ^ censorrec;
    censorrec = 
	record 
	    censor: array [minint..maxcensor] of boolean
	end; {Was packed}
    genotype = array [1..maxfem] of double;
    covmatrix = array [1..maxtrait, 1..maxtrait] of double;
    hapvector = array [1..maxlocus] of 0..maxall; {Was packed}
    mutarray = array [0..2, 1..2] of integer;
    thesemeans = array [1..maxtrait] of double;
    means = array [0..maxall, 1..maxall] of thesemeans;
    pathway = (auto, mauto, sex, msex);
    direction = (peelup, peeldown);

    binset = set of 1..maxfact;
    phenarray = array [1..maxall] of binset;

    locuspoint = ^ locusvalues;
    phenpoint = ^ phenotype;
    locustype = (affection, quantitative, binary, null);
    locusvalues = 
	record 
	    nallele, format: integer;
	    freq: array [1..maxall] of double;
	    privlocus: locuspoint;
	    case which: locustype of
		affection: (
		    pen: array [0..maxall, 1..maxall, 0..2, 1..maxliab] of double;
		    nclass: integer
		);
		quantitative: (
		    ntrait: integer;
		    pm: means;
		    vmat: covmatrix;
		    det, contrait, conmat: double
		);
		binary: (
		    allele: phenarray
		);
		null: (
		)
	end;
    phenotype = 
	record 
	    case which: locustype of
		binary: (
		    phenf: binset
		);
		quantitative: (
		    x: array [1..maxtrait] of double;
		    missing: boolean
		);
		affection: (
		    aff, liability: integer
		);
		null: (
		)
	end;

    ind = ^ thisperson;
    genpoint = ^ thisarray;
    indphen = array [1..maxlocus] of phenpoint;
    thisarray = 
	record 
	    genarray: genotype
	end;

    possvect = array [1..maxall, 1..maxall] of boolean;
    possarray = array [1..maxlocus] of possvect;

    infoptr = ^ information;
    information = 
	record 
	    possible: possarray
	end;

    thisperson = 
	record 
	    id, ped, inloop: integer;
	    pa, ma, foff, nextpa, nextma: ind;
	    gen: genpoint;
	    phen: indphen;
	    privphen: phenpoint;
	    unknown, multi, done, up, male, firstpass: boolean;
	    store: infoptr
	end;

    subhap = array [1..maxseg] of integer;
    thetarray = array [1..maxlocus] of double;
    thetapoint = ^ thetavalues;
    happrob = array [1..maxneed] of double;
    thetavalues = 
	record 
	    theta: thetarray;
	    segprob: happrob
	end;
    {GEMINI}
    vector = array [1..maxn] of double;
    matrix = array [1..maxn] of vector;

var
    itmax: boolean; {Flag indicating iterated parameters}
    nrep, j: integer; {nrep is a counter of number of replications}
    rep, ltemp, aver, variance, big, smll: double;
    clods: array [1..3] of integer;
    lodlimit: array [1..3] of double;
    {These are where various statistics are stored}
    nfactor: array [1..maxlocus] of integer;
    opeds: integer; {Number of original pedigrees}
    infile: boolean; {TRUE if currently in the middle of the pedigree file}
    oldped: integer; {oldped stores number of next pedigree to be read.}
    ddist: double;
    thispath: pathway;
    informative: array [1..maxped] of boolean;
    rare, risk1, risk2: array [1..maxfem] of boolean; {Was packed}
    riskmale: array [1..maxhap] of boolean; {Was packed}
    approxstruct: approxpnt;
    censorstruct: censorpnt;
    thisc: integer;
    malechild: array [1..maxchild] of boolean;
    thischild: array [1..maxchild] of genpoint;
    nchild: integer;
    seghap1, seghap2: array [1..maxfem] of 1..maxhap;
    segstart: array [1..maxfem] of 1..maxfem;
    probstart, probend: array [1..maxfem] of 1..maxneed;
    nohom: array [1..maxlocus] of boolean;
    thislocus: array [1..maxlocus] of locuspoint;
    increment, order: array [1..maxlocus] of integer;
    {PEOPLE}
    person: array [0..maxind] of ind;
    proband: array [1..maxped] of ind;
    looppers: array [1..maxped, 1..maxloop, 1..2] of ind;
    {MUTATION }
    muthap: array [1..maxhap] of 1..maxhap;
    gennustruct: gennuptr;
    {RECOMBINATION}
    maletheta, femaletheta: thetapoint;
    {FREQUENCIES }
    hapfreq: genpoint;
    {RISK}
    riskall: 1..maxall;
    {OTHERS}
    risksys, mutsys, nlocus, which, lastpriv, i, lastspeed, lastseg, segperson: integer;
    {speedfile modification above; added variable lastspeed,lastseg,segperson}
    nuhap: 1..maxhap;
    fgeno, mgeno: 1..maxfem;
    nuped, totperson: integer;
    segscale, mutmale, mutfemale, like, alike, distratio: double;
    interfer, disequi, sexlink, risk, sexdif, readfemale, inconsistent, dolod, firstapprox, firsttime, lasttime: boolean;
    isim, outfile, datafile, speedfile, ipedfile, simout: text;
    chtemp: char;
    {isim}
    final: text;
    mapping, penalty: boolean;
    penlike: double;
    {GEMINI}
    active: boolean;
    continue: (go, quit, restart);

    {VARIABLES FROM COMMON MIN1}
    tl: matrix;
    d, g, gsave, y, p: vector;

    {VARIABLE FROM COMMON MIN2}
    nit, nfe, idg, idif, isw, iret, ibnd, ihess, ivar, ihx, maxit: integer;
    tol, tmin, gh, trupb, ptg: double;

    {DIMENSIONED IN GEMINI}
    holdxall, xall, x, gv, se: vector;
    holditp, itp: array [1..maxn] of integer;
    bnd: array [1..maxn, 1..2] of double;
    nall, n, icall: integer;
    tbnd, f, fsmf, fsav2, t, hx, xsave, fxph, fxmh, xp, xpm, ytp: double;

    {ARRAYS IN UPDATE}
    wtil, ztil, w, z, wtjp1, ztjp1, s, dp: vector;

    nu, muj, lambj, lambj2, sb, sc, sd, fbcd, alpha, sa, thet1, thet2, aa, bb, cc, del2, alph1, alph2, rdiv, eta, aj, thj, bj, gamlj, betlj, del: double;

    iswup: integer;

    {FROM STEP}
    xt: vector;

    fsave, sumt, twot, ft, f2t, ft2, scal: double;

    {FROM CHKBND }
    clb, check, eh, teq: double;

    {OTHERS}
    itsys: integer;
    {USUALLY IN INIB}
    bmat: matrix;


function min(a, b: double): double;
begin
    if a < b then 
	min := a
    else 
	min := b
end; { min }

function max(a, b: double): double;
begin
    if a > b then 
	max := a
    else 
	max := b
end; { max }


procedure outcontrol(var z: vector);

var
    i, j: integer;

begin
    j := 0;
    for i := 1 to n do begin
	if prn then 
	    write(outfile, z[i]);
	j := j + 1;
	if j = 4 then begin
	    if prn then 
		writeln(outfile);
	    j := 0
	end
    end;
    if (j <> 0) and prn then 
	writeln(outfile)
end; { outcontrol }

function mapfunction(theta1, theta2: double): double;
{User defined function giving recombination between
flanking markers as a function of recombination
between adjacent markers}

begin
    mapfunction := (theta1 + theta2) / (1 + 4 * theta1 * theta2)
end; { mapfunction }

function getdist(var theta: double): double;

begin
    if theta < 0.5 then 
	getdist := -ln(1.0 - 2.0 * theta) / 2.0
    else 
	getdist := 10.0
end; { getdist }

function invdist(var dist: double): double;

begin
    if dist <> 10.0 then 
	invdist := (1.0 - exp(-2 * dist)) / 2.0
    else 
	invdist := 0.5
end; { invdist }

procedure invert(var m: covmatrix; n: integer; var det: double);

var
    v: covmatrix;
    val: double;
    i, j, k: integer;

begin
    det := 1.0;
    for i := 1 to n do begin
	val := m[i, i];
	for k := 1 to i - 1 do 
	    val := val - v[k, i] * v[k, i];
	det := det * val;
	v[i, i] := sqrt(val);
	for j := i + 1 to n do begin
	    val := m[i, j];
	    for k := 1 to i - 1 do 
		val := val - v[k, i] * v[k, j];
	    v[i, j] := val / v[i, i];
	    v[j, i] := 0
	end
    end;
    for i := 1 to n do begin
	m[i, i] := 1 / v[i, i];
	for j := i + 1 to n do begin
	    val := 0;
	    for k := i to j - 1 do 
		val := val - m[k, i] * v[k, j];
	    m[j, i] := val / v[j, j];
	    m[i, j] := 0
	end
    end;
    for i := 1 to n do begin
	for j := 1 to i do begin
	    val := 0;
	    for k := j to n do 
		val := val + m[k, i] * m[k, j];
	    v[i, j] := val;
	    v[j, i] := val
	end
    end;
    m := v
end; { invert }


procedure recombination;

var
    i, here: integer;
    p1, p2, p3, p4: double;
    oldtheta: thetarray;

    procedure recombine(var theta: thetarray; var segprob: happrob);

    var
	there, nhap, system: integer;
	thishet: array [1..maxlocus] of boolean;
	hap1, hap2: hapvector;
	thishap1, thishap2: array [1..maxseg] of hapvector;

	procedure scramble;

	var
	    whichhap, start, length, i, j, k: integer;
	    recval, val: double;

	begin
	    start := 0;
	    repeat
		start := start + 1
	    until thishet[start];
	    length := there - here;
	    for i := 2 to length do begin
		hap1 := thishap1[i];
		for j := 1 to length do begin
		    val := 0.5;
		    whichhap := 1;
		    recval := theta[start];
		    for k := start + 1 to nlocus do 
			if not thishet[k] then 
			    recval := recval * (1.0 - theta[k]) + theta[k] * (1.0 - recval)
			else begin
			    if whichhap = 1 then begin
				if thishap1[j, k] = hap1[k] then 
				    val := val * (1 - recval)
				else begin
				    val := val * recval;
				    whichhap := 2
				end
			    end else begin
				if thishap2[j, k] = hap1[k] then 
				    val := val * (1 - recval)
				else begin
				    val := val * recval;
				    whichhap := 1
				end
			    end;
			    recval := theta[k]
			end;
		    there := there + 1;
		    segprob[there] := val
		end
	    end
	end; { scramble }


	procedure setrec(val: double);

	begin
	    nhap := nhap + 1;
	    thishap1[nhap] := hap1;
	    thishap2[nhap] := hap2;
	    there := there + 1;
	    segprob[there] := val
	end; { setrec }


	procedure dointer;

	var
	    i: integer;
	    temphet: array [1..3] of boolean;

	begin
	    for i := 1 to nlocus do 
		temphet[i] := thishet[i];

	    if temphet[1] and temphet[2] and not temphet[3] then begin
		setrec(0.5 - 0.5 * theta[1]);
		setrec(0.5 * theta[1]);
		setrec(0.5 * theta[1]);
		setrec(0.5 - 0.5 * theta[1])
	    end else if temphet[3] and temphet[2] and not temphet[1] then begin
		setrec(0.5 - 0.5 * theta[nlocus - 1]);
		setrec(0.5 * theta[nlocus - 1]);
		setrec(0.5 * theta[nlocus - 1]);
		setrec(0.5 - 0.5 * theta[nlocus - 1])
	    end else if temphet[3] and temphet[1] and not temphet[2] then begin
		setrec(0.5 - 0.5 * theta[nlocus]);
		setrec(0.5 * theta[nlocus]);
		setrec(0.5 * theta[nlocus]);
		setrec(0.5 - 0.5 * theta[nlocus])
	    end else if temphet[1] and temphet[2] and temphet[3] then begin
		setrec(0.5 * p4);
		setrec(0.5 * p2);
		setrec(0.5 * p3);
		setrec(0.5 * p1);
		setrec(0.5 * p2);
		setrec(0.5 * p4);
		setrec(0.5 * p1);
		setrec(0.5 * p3);
		setrec(0.5 * p3);
		setrec(0.5 * p1);
		setrec(0.5 * p4);
		setrec(0.5 * p2);
		setrec(0.5 * p1);
		setrec(0.5 * p3);
		setrec(0.5 * p2);
		setrec(0.5 * p4)
	    end else  {not informative}
		setrec(0.5)
	end; { dointer }


	procedure nexthet(i: integer; val: double; inphase: boolean);

	var
	    newval, recval: double;

	begin
	    recval := theta[i];
	    repeat
		i := i + 1;
		hap1[i] := 0;
		hap2[i] := 0;
		if not thishet[i] and not (i = nlocus) then 
		    recval := recval * (1 - theta[i]) + (1 - recval) * theta[i]
	    until (i = nlocus) or thishet[i];
	    if i <> nlocus then begin
		if inphase then 
		    newval := val * (1 - recval)
		else 
		    newval := val * recval;
		hap1[i] := 1;
		hap2[i] := 2;
		nexthet(i, newval, true);
		hap2[i] := 1;
		hap1[i] := 2;
		if not inphase then 
		    newval := val * (1 - recval)
		else 
		    newval := val * recval;
		nexthet(i, newval, false)
	    end else if not thishet[i] then 
		setrec(val)
	    else begin
		if inphase then 
		    newval := val * (1 - recval)
		else 
		    newval := val * recval;
		hap1[i] := 1;
		hap2[i] := 2;
		setrec(newval);
		if not inphase then 
		    newval := val * (1 - recval)
		else 
		    newval := val * recval;
		hap2[i] := 1;
		hap1[i] := 2;
		setrec(newval)
	    end
	end; { nexthet }


	procedure getrecprob;

	var
	    i: integer;

	begin
	    nhap := 0;
	    there := here;
	    i := 0;
	    repeat
		i := i + 1;
		if thishet[i] then begin
		    hap1[i] := 1;
		    hap2[i] := 2
		end else begin
		    hap1[i] := 0;
		    hap2[i] := 0
		end
	    until thishet[i] or (i = nlocus);
	    if i = nlocus then 
		setrec(0.5)
	    else if interfer then 
		dointer
	    else 
		nexthet(i, 0.5, true);
	    if (nhap > 1) and not interfer then 
		scramble;
	    here := there
	end; { getrecprob }


	procedure gethet(var system: integer);

	var
	    newsystem: integer;

	begin
	    newsystem := system + 1;
	    thishet[system] := false;
	    if system <> nlocus then 
		gethet(newsystem)
	    else 
		getrecprob;
	    thishet[system] := true;
	    if system <> nlocus then 
		gethet(newsystem)
	    else 
		getrecprob
	end; { gethet }


    begin {RECOMBINE}
	here := 0;
	system := 1;
	gethet(system)
    end; { recombine }


    procedure getfemaletheta;

    var
	dist: double;
	ntheta, i: integer;

    begin
	if interfer then 
	    ntheta := nlocus
	else 
	    ntheta := nlocus - 1;
	for i := 1 to ntheta do begin
	    dist := getdist(maletheta^.theta[i]) * distratio;
	    femaletheta^.theta[i] := invdist(dist)
	end
    end; { getfemaletheta }


begin
    if interfer then 
	with maletheta^ do begin
	    if mapping then 
		theta[nlocus] := mapfunction(theta[nlocus - 1], theta[1]);
	    oldtheta := theta;
	    if not mapping and not dolod then begin
		for i := 1 to nlocus do 
		    oldtheta[i] := 1 / (1 + exp(oldtheta[i]));
		theta[1] := oldtheta[1] + oldtheta[nlocus];
		theta[nlocus - 1] := oldtheta[nlocus - 1] + oldtheta[nlocus];
		theta[nlocus] := oldtheta[1] + oldtheta[nlocus - 1];
		p1 := oldtheta[1];
		p2 := oldtheta[nlocus - 1];
		p3 := oldtheta[nlocus];
		p4 := 1.0 - p1 - p2 - p3
	    end else begin
		p1 := (theta[1] + theta[nlocus] - theta[nlocus - 1]) / 2.0;
		p2 := (theta[nlocus - 1] + theta[nlocus] - theta[1]) / 2.0;
		p3 := (theta[nlocus - 1] + theta[1] - theta[nlocus]) / 2.0;
		p4 := 1.0 - p1 - p2 - p3
	    end;
	    recombine(theta, segprob)
	end
    else 
	recombine(maletheta^.theta, maletheta^.segprob);
    if sexdif then begin
	if not readfemale then begin
	    if interfer and not dolod then 
		with maletheta^ do begin
		    theta[1] := oldtheta[1] + oldtheta[nlocus];
		    theta[nlocus - 1] := oldtheta[nlocus - 1] + oldtheta[nlocus];
		    theta[nlocus] := oldtheta[1] + oldtheta[nlocus - 1]
		end;
	    getfemaletheta
	end;
	if interfer then 
	    with femaletheta^ do begin
		if mapping then 
		    theta[nlocus] := mapfunction(theta[nlocus - 1], theta[1]);
		if readfemale and not mapping and not dolod then begin
		    oldtheta := theta;
		    for i := 1 to nlocus do 
			oldtheta[i] := 1 / (1 + exp(oldtheta[i]));
		    theta[1] := oldtheta[1] + oldtheta[nlocus];
		    theta[nlocus - 1] := oldtheta[nlocus - 1] + oldtheta[nlocus];
		    theta[nlocus] := oldtheta[1] + oldtheta[nlocus - 1];
		    p1 := oldtheta[1];
		    p2 := oldtheta[nlocus - 1];
		    p3 := oldtheta[nlocus];
		    p4 := 1.0 - p1 - p2 - p3
		end else begin
		    p1 := (theta[1] + theta[nlocus] - theta[nlocus - 1]) / 2.0;
		    p2 := (theta[nlocus - 1] + theta[nlocus] - theta[1]) / 2.0;
		    p3 := (theta[nlocus - 1] + theta[1] - theta[nlocus]) / 2.0;
		    p4 := 1.0 - p1 - p2 - p3
		end;
		recombine(theta, segprob)
	    end
	else 
	    recombine(femaletheta^.theta, femaletheta^.segprob)
    end;
    if prn then 
	if firsttime then 
	    if here < maxneed then 
		writeln(output, 'Maxneed can be reduced to ', here)
end; { recombination }

procedure setparam;

var
    i, j, k, l, m: integer;

begin
    penalty := false;
    k := 0;
    with maletheta^ do 
	for i := 1 to nlocus - 1 do begin
	    k := k + 1;
	    theta[i] := xall[k]
	end;
    if interfer and not mapping then begin
	k := k + 1;
	maletheta^.theta[nlocus] := xall[k]
    end;
    if sexdif then 
	if readfemale then begin
	    with femaletheta^ do 
		for i := 1 to nlocus - 1 do begin
		    k := k + 1;
		    theta[i] := xall[k]
		end;
	    if interfer and not mapping then begin
		k := k + 1;
		femaletheta^.theta[nlocus] := xall[k]
	    end
	end else begin
	    k := k + 1;
	    distratio := xall[k]
	end;
    if itsys <> 0 then 
	with thislocus[itsys]^ do begin
	    if disequi then 
		with hapfreq^ do begin
		    for i := 1 to nuhap do 
			genarray[i] := xall[i + k];
		    genarray[nuhap] := 1.0;
		    for i := 1 to nuhap - 1 do 
			genarray[nuhap] := genarray[nuhap] - genarray[i];
		    penalty := genarray[nuhap] < 0.0;
		    k := nuhap - 1 + k
		end
	    else begin
		for i := 1 to nallele - 1 do 
		    freq[i] := xall[i + k];
		freq[nallele] := 1.0;
		for i := 1 to nallele - 1 do 
		    freq[nallele] := freq[nallele] - freq[i];
		penalty := freq[nallele] < 0.0;
		k := nallele - 1 + k
	    end;
	    if which = affection then 
		for l := 1 to nclass do begin
		    for i := 1 to nallele do 
			for j := i to nallele do begin
			    k := k + 1;
			    pen[i, j, 2, l] := xall[k];
			    pen[i, j, 1, l] := 1.0 - xall[k];
			    for m := 0 to 2 do 
				pen[j, i, m, l] := pen[i, j, m, l]
			end;
		    if sexlink then 
			for i := 1 to nallele do begin
			    k := k + 1;
			    pen[0, i, 2, l] := xall[k];
			    pen[0, i, 1, l] := 1.0 - xall[k]
			end
		end
	    else if which = quantitative then begin
		{USER MUST CHANGE FOR SYSTEMS WITH MORE THAN 2 ALLELES}
		for i := 1 to ntrait do begin
		    pm[1, 1, i] := xall[k + 1];
		    pm[2, 2, i] := xall[k + 2] + xall[k + 1];
		    pm[1, 2, i] := xall[k + 3] * xall[k + 2] + xall[k + 1];
		    pm[0, 1, i] := xall[k + 1];
		    pm[0, 2, i] := pm[1, 2, i];
		    pm[2, 1, i] := pm[1, 2, i];
		    k := k + 3
		end;
		for i := 1 to ntrait do 
		    for j := i to ntrait do begin
			k := k + 1;
			vmat[i, j] := xall[k];
			vmat[j, i] := xall[k]
		    end;
		invert(vmat, ntrait, det);
		det := 1 / sqrt(det)
	    end
	end
end; { setparam }


procedure getlocations;

var
    ngene, nseg, here, there, start, nhet, thisseg: integer;
    rarepresent, riskhom, riskhet: boolean;
    hap1, hap2: hapvector;
    thishet: array [1..maxlocus] of boolean;

    function checkrare: boolean;

    var
	i: integer;
	check: boolean;

    begin
	check := false;
	for i := 1 to nlocus do 
	    if nohom[i] then 
		with thislocus[i]^ do 
		    if (freq[hap1[i]] < minfreq) or (freq[hap2[i]] < minfreq) then 
			check := true;
	checkrare := check
    end; { checkrare }


    procedure checkrisk(var riskhet, riskhom: boolean);

    begin
	riskhet := false;
	riskhom := false;
	if (hap1[risksys] = riskall) and (hap2[risksys] = riskall) then 
	    riskhom := true
	else if (hap1[risksys] <> riskall) and (hap2[risksys] = riskall) or (hap2[risksys] <> riskall) and (hap1[risksys] = riskall) then 
	    riskhet := true
    end; { checkrisk }


    function gethapn(var hap: hapvector): integer;

    var
	i, n: integer;

    begin
	n := 1;
	for i := 1 to nlocus do 
	    n := n + increment[i] * (hap[i] - 1);
	gethapn := n
    end; { gethapn }


    procedure domalerisk;

	procedure setrisk;

	var
	    n: integer;

	begin
	    n := gethapn(hap1);
	    if hap1[risksys] = riskall then 
		riskmale[n] := true
	    else 
		riskmale[n] := false
	end; { setrisk }


	procedure getriskhap(system: integer);

	var
	    i: integer;

	begin
	    with thislocus[system]^ do 
		for i := 1 to nallele do begin
		    hap1[system] := i;
		    if system <> nlocus then 
			getriskhap(system + 1)
		    else 
			setrisk
		end
	end; { getriskhap }


    begin
	getriskhap(1)
    end; { domalerisk }


    procedure domutation;

	procedure setmutation;

	var
	    i, n: integer;

	begin
	    n := gethapn(hap1);
	    if hap1[mutsys] = thislocus[mutsys]^.nallele then 
		muthap[n] := n
	    else begin
		i := hap1[mutsys];
		hap1[mutsys] := thislocus[mutsys]^.nallele;
		muthap[n] := gethapn(hap1);
		hap1[mutsys] := i
	    end
	end; { setmutation }


	procedure getmuthap(system: integer);

	var
	    i: integer;

	begin
	    with thislocus[system]^ do 
		for i := 1 to nallele do begin
		    hap1[system] := i;
		    if system <> nlocus then 
			getmuthap(system + 1)
		    else 
			setmutation
		end
	end; { getmuthap }


    begin
	getmuthap(1)
    end; { domutation }


    procedure setnumbers;

    var
	nhap1, nhap2: integer;

    begin
	ngene := ngene + 1;

	segstart[ngene] := here + 1;
	probstart[ngene] := there + 1;
	probend[ngene] := there + nseg;

	there := there + nseg;

	nhap1 := gethapn(hap1);
	nhap2 := gethapn(hap2);
	gennustruct^.genenumber[nhap1, nhap2] := ngene;
	gennustruct^.genenumber[nhap2, nhap1] := ngene;

	if minfreq <> 0.0 then 
	    if rarepresent then 
		rare[ngene] := true
	    else 
		rare[ngene] := false
	else 
	    rare[ngene] := false;
	if risk then begin
	    risk1[ngene] := riskhet;
	    risk2[ngene] := riskhom
	end;

	thisseg := thisseg + 1;
	seghap1[thisseg] := nhap1;
	seghap2[thisseg] := nhap2
    end; { setnumbers }


    procedure hapscr(system, nscramble: integer);

    var
	i, j: integer;

    begin
	if thishet[system] then 
	    nscramble := nscramble + 1;
	if system <> nlocus then 
	    hapscr(system + 1, nscramble)
	else 
	    setnumbers;
	if nscramble > 1 then 
	    if hap1[system] <> hap2[system] then begin
		i := hap1[system];
		j := hap2[system];
		hap1[system] := j;
		hap2[system] := i;
		if system <> nlocus then 
		    hapscr(system + 1, nscramble)
		else 
		    setnumbers;
		hap1[system] := i;
		hap2[system] := j
	    end
    end; { hapscr }


    procedure sethap(system: integer);

    var
	i, j: integer;

    begin
	with thislocus[system]^ do 
	    if thishet[system] then 
		for i := 1 to nallele - 1 do begin
		    hap1[system] := i;
		    for j := i + 1 to nallele do begin
			hap2[system] := j;
			if system <> nlocus then 
			    sethap(system + 1)
			else begin
			    rarepresent := checkrare;
			    if risk then 
				checkrisk(riskhet, riskhom);
			    there := start;
			    thisseg := here;
			    hapscr(1, 0);
			    here := here + nseg
			end
		    end
		end
	    else 
		for i := 1 to nallele do begin
		    hap1[system] := i;
		    hap2[system] := i;
		    if system <> nlocus then 
			sethap(system + 1)
		    else begin
			rarepresent := checkrare;
			if risk then 
			    checkrisk(riskhet, riskhom);
			thisseg := here;
			there := start;
			hapscr(1, 0);
			here := here + nseg
		    end
		end
    end; { sethap }


    procedure starthap;

    var
	i: integer;

    begin
	nseg := 1;
	for i := 2 to nhet do 
	    nseg := nseg * 2;
	sethap(1);
	start := there
    end; { starthap }


    procedure gethet1(system: integer);

    begin
	thishet[system] := false;
	if system <> nlocus then 
	    gethet1(system + 1)
	else 
	    starthap;
	thishet[system] := true;
	nhet := nhet + 1;
	if system <> nlocus then 
	    gethet1(system + 1)
	else 
	    starthap;
	nhet := nhet - 1
    end; { gethet1 }


begin				{getlocations??}
    nhet := 0;
    here := 0;
    there := 0;
    ngene := 0;
    start := 0;
    gethet1(1);
    if mutsys <> 0 then 
	domutation;
    if sexlink and risk then 
	domalerisk
end; { getlocations }				{getlocations??}


procedure inputdata;

    procedure inputerror(nerror, par1, par2: integer);

    begin
	writeln('Fatal error detected in procedure inputdata');
	case nerror of
	    0:
		writeln('Number of loci ', par1: 2, ' exceeds the constant maxlocus');
	    1:
		writeln('Number of loci read ', par1: 2, '. Less than minimum of 1');
	    2:
		writeln('Error detected reading loci order. Locus number ', par2: 2, ' in position ', par1: 2, ' exceeds number of loci');
	    3:
		writeln('Error detected reading loci order. Illegal locus number ', par2: 2, ' in position ', par1: 2);
	    4:
		writeln('Error detected reading loci order. Locus number repeated in positions ', par1: 2, ' and ', par2: 2);
	    5:
		writeln('Error detected reading locus description. Illegal locus type ', par2: 2, ' for locus ', par1: 2);
	    6:
		writeln('Error detected reading locus description for system ', par1: 2, '. Number of alleles  ', par1: 2, ' exceeds maxall');
	    7:
		writeln('Error detected reading locus description for system ', par1: 2, '. Illegal number of alleles  ', par2: 2);
	    8:
		writeln('Error detected reading locus description for system ', par1: 2, '. Number of factors  ', par2: 2, ' exceeds maxfact');
	    9:
		writeln('Error detected reading locus description for system ', par1: 2, '. Illegal number of factors  ', par2: 2);
	    10:
		writeln('Error detected reading locus description for system ', par1: 2, '. Alleles not codominant');
	    11:
		writeln('Error detected reading pedigree record ', par1: 2, '. Illegal code for sex ', par2: 2);
	    12:
		writeln('Error detected reading pedigree record at pedigree', par1: 2, '. Maximum number of pedigree records exceeded');
	    13:
		writeln('Error detected reading pedigree record ', par1: 2, '. Maximum number of individuals exceeded');
	    14:
		writeln('Error detected reading pedigree record ', par1: 2, '. Illegal binary factor code ', par2: 2);
	    15:
		writeln('Error detected reading pedigree record ', par1: 2, '. No allelic pair for genotype');
	    16:
		writeln('Error detected reading pedigree record ', par1: 2, '. Allele number ', par2: 2, ' exceeds maxall');
	    17:
		writeln('Error detected reading pedigree record ', par1: 2, '. Illegal allele number ', par2: 2);
	    18:
		writeln('Number of systems after factorization (', par1: 3, ') exceeds maxsystem');
	    19:
		writeln('Number of systems after factorization (', par1: 3, ') less than minimum of 1');
	    20:
		writeln('Number of recombination types (', par1: 3, ') exceeds maxrectype');
	    21:
		writeln('Number of recombination types (', par1: 3, ') less than minimum of 1');
	    22:
		writeln('End of file detected in tempdat by procedure readthg before all data found');
	    23:
		writeln('Error detected reading iterated locus in datafile. Value (', par1: 3, ') greater than nlocus');
	    24:
		writeln('Error detected reading iterated locus in datafile. Illegal value (', par1: 3, ')');
	    25:
		writeln('Number of iterated parameters greater then maxn');
	    26:
		writeln('Error detected reading pedigree record ', par1: 2, '. Liability class (', par2: 2, ') exceeds nclass');
	    27:
		writeln('Error detected reading pedigree record ', par1: 2, '. Illegal liability class (', par2: 2, ')');
	    28:
		writeln('Error detected reading locus description for system', par1: 2, '. Liability classes (', par2: 3, ') exceed maxliab');
	    29:
		writeln('Error detected reading locus description for system', par1: 2, '. Illegal number of liability classes (', par2: 3, ')');
	    30:
		writeln('Error detected reading locus description for system', par1: 2, '. Penetrance out of range');
	    31:
		writeln('Error detected reading locus description for system', par1: 2, '. Number of traits (', par2: 3, ') exceeds maxtrait');
	    32:
		writeln('Error detected reading locus description for system', par1: 2, '. Number of traits out of range (', par2: 3, ')');
	    33:
		writeln('Error detected reading locus description for system', par1: 2, '. Variance must be positive');
	    34:
		writeln('Error detected reading locus description for system', par1: 2, '. Variance multiplier must be positive');
	    35:
		writeln('Error detected reading locus description for system', par1: 2, '. Risk allele ', par2: 3, ') exceeds nallele');
	    36:
		writeln('Error detected reading locus description for system', par1: 2, '. Illegal risk allele (', par2: 3, ')');
	    37:
		writeln('Error detected reading datafile. Risk locus ', par2: 3, ') exceeds nlocus');
	    38:
		writeln('Error detected reading datafile. Illegal value for risk locus ', par2: 3, ')');
	    39:
		writeln('Error detected reading datafile. Mutation locus ', par2: 3, ') exceeds nlocus');
	    40:
		writeln('Error detected reading datafile. Illegal value for mutation locus ', par2: 3, ')');
	    41:
		writeln('Error in datafile, line 1, item 4: Program code not for ISIM/ILINK')
	end;
	writeln('Press return to halt...', chr(7));
	readln;
	halt; {SUN}
    end; { inputerror }


    procedure inputwarning(nwarning, par1, par2: integer);

    begin
	writeln('Warning number from procedure inputdata');
	case nwarning of
	    0:
		writeln('Illegal sex difference parameter ', par1: 2, ' Parameter should be 0, 1, or 2');
	    1:
		writeln('Illegal interference parameter ', par1: 2, ' Lack of interference assumed');
	    2:
		writeln('Illegal sex difference parameter ', par1: 2, ' Parameter must be 0 with sex-linked data');
	    3:
		writeln('Non-standard affection status', par2: 4, ' interpreted as normal in pedigree record', par1: 5)
	end;
	writeln('Press return to continue...', chr(7));
	readln
    end; { inputwarning }


    procedure setiterate;

    type
	boundary = (nobnd, probnd, zerobnd, logbnd);

    var
	bndtype: array [1..maxn] of boundary;
	i, j, k, l: integer;

    begin
	readln(datafile, i);
	if i > nlocus then 
	    inputerror(23, nlocus, i);
	if i < 0 then 
	    inputerror(24, i, i);
	if i = 0 then 
	    itsys := 0
	else begin
	    itsys := 1;
	    while order[itsys] <> i do 
		itsys := itsys + 1
	end;
	k := 0;
	with maletheta^ do 
	    for i := 1 to nlocus - 1 do begin
		k := k + 1;
		if k > maxn then 
		    inputerror(25, maxn, k);
		xall[k] := theta[i];
		if interfer and not mapping then 
		    bndtype[k] := logbnd
		else 
		    bndtype[k] := probnd
	    end;
	if interfer and not mapping then begin
	    k := k + 1;
	    if k > maxn then 
		inputerror(25, maxn, k);
	    xall[k] := maletheta^.theta[nlocus];
	    bndtype[k] := logbnd
	end;
	if sexdif then 
	    if readfemale then begin
		with femaletheta^ do 
		    for i := 1 to nlocus - 1 do begin
			k := k + 1;
			if k > maxn then 
			    inputerror(25, maxn, k);
			xall[k] := theta[i];
			if interfer and not mapping then 
			    bndtype[k] := logbnd
			else 
			    bndtype[k] := probnd
		    end;
		if interfer and not mapping then begin
		    k := k + 1;
		    if k > maxn then 
			inputerror(25, maxn, k);
		    xall[k] := femaletheta^.theta[nlocus];
		    bndtype[k] := logbnd
		end
	    end else begin
		k := k + 1;
		if k > maxn then 
		    inputerror(25, maxn, k);
		xall[k] := distratio;
		bndtype[k] := zerobnd
	    end;
	nall := k;
	if itsys <> 0 then 
	    with thislocus[itsys]^ do begin
		if disequi then begin
		    for i := 1 to nuhap - 1 do 
			xall[i + k] := hapfreq^.genarray[i];
		    nall := nuhap + k - 1
		end else begin
		    for i := 1 to nallele - 1 do 
			xall[i + k] := freq[i];
		    nall := nallele + k - 1
		end;
		if nall > maxn then 
		    inputerror(25, maxn, nall);
		for i := k + 1 to nall do 
		    bndtype[i] := probnd;
		k := nall;

		if which = quantitative then 
		    nall := nall + ntrait * (ntrait + 1) div 2 + ntrait * 3

		else if which = affection then begin
		    nall := nall + nclass * nallele * (nallele + 1) div 2;
		    if sexlink then 
			nall := nall + nallele
		end;
		if nall > maxn then 
		    inputerror(25, maxn, nall);

		if which = affection then 

		    for l := 1 to nclass do begin
			for i := 1 to nallele do 
			    for j := i to nallele do begin
				k := k + 1;
				if k > maxn then 
				    inputerror(25, maxn, k);
				xall[k] := pen[i, j, 2, l];
				bndtype[k] := probnd
			    end;
			if sexlink then 
			    for i := 1 to nallele do begin
				k := k + 1;
				if k > maxn then 
				    inputerror(25, maxn, k);
				xall[k] := pen[0, i, 2, l];
				bndtype[k] := probnd
			    end
		    end
		else if which = quantitative then begin

		    invert(vmat, ntrait, det);
		    {USER MUST MODIFY FOR A SYSTEM WILL MORE THAN 2 ALLELES }
		    for l := 1 to ntrait do 
			for i := 1 to nallele do 
			    for j := i to nallele do begin
				k := k + 1;
				if (i = 1) and (j = 1) then 
				    xall[k] := pm[i, j, l]
				else if i = 1 then 
				    xall[k] := pm[j, j, l] - xall[k - 1]
				else 
				    xall[k] := (pm[i - 1, j, l] - xall[k - 2]) / xall[k - 1];
				bndtype[k] := nobnd
			    end;
		    for i := 1 to ntrait do 
			for j := i to ntrait do begin
			    k := k + 1;
			    if k > maxn then 
				inputerror(25, maxn, k);
			    xall[k] := vmat[i, j];
			    if i = j then 
				bndtype[k] := zerobnd
			    else 
				bndtype[k] := nobnd
			end
		end
	    end;
	for i := 1 to nall do 
	    itp[i] := 0;
	i := 0;
	while (i <= nall) and not eoln(datafile) do begin
	    i := i + 1;
	    read(datafile, itp[i])
	end;
	n := 0;
	for i := 1 to nall do 
	    if itp[i] = 1 then begin
		n := n + 1;
		x[n] := xall[i]
	    end;
	for i := 1 to nall do begin
	    bnd[i, 1] := -1000.0;
	    bnd[i, 2] := 1000.0;
	    if bndtype[i] = probnd then begin
		bnd[i, 1] := 0.0;
		bnd[i, 2] := 1.0
	    end else if bndtype[i] = logbnd then begin
		bnd[i, 1] := -20.0;
		bnd[i, 2] := 20.0
	    end else if bndtype[i] = zerobnd then 
		bnd[i, 1] := 0.0
	end
    end; { setiterate }
 {setiterate}

    procedure readloci;

    var
	i, j, coupling, autosomal, independent, difference, whichtype, nupriv: integer;

	procedure getlocus(system: integer);

	var
	    j: integer;

	    procedure getbin(var locus: locuspoint; var system: integer);

	    var
		i, j, k: integer;

	    begin
		readln(datafile, nfactor[system]);
		if nfactor[system] > maxfact then 
		    inputerror(8, system, nfactor[system]);
		if nfactor[system] <= 0 then 
		    inputerror(9, system, nfactor[system]);
		with locus^ do begin
		    for i := 1 to nallele do 
			allele[i] := [];
		    for i := 1 to nallele do 
			for j := 1 to nfactor[system] do begin
			    read(datafile, k);
			    if k = 1 then 
				allele[i] := allele[i] + [j]
			end
		end;
		readln(datafile)
	    end; { getbin }


	    procedure getnumber(var locus: locuspoint; var system: integer);

	    var
		i: integer;

	    begin
		with locus^ do 
		    for i := 1 to nallele do 
			allele[i] := [i]
	    end; { getnumber }


	    procedure getpen(var locus: locuspoint);

	    var
		i, j, k, l: integer;

	    begin
		with locus^ do begin
		    readln(datafile, nclass);
		    if nclass > maxliab then 
			inputerror(28, system, nclass);
		    if nclass <= 0 then 
			inputerror(29, system, nclass);
		    for l := 1 to nclass do begin
			for i := 1 to nallele do 
			    for j := i to nallele do begin
				read(datafile, pen[i, j, 2, l]);
				if (pen[i, j, 2, l] < 0) or (pen[i, j, 2, l] > 1.0) then 
				    inputerror(30, system, system);
				pen[i, j, 1, l] := 1 - pen[i, j, 2, l];
				pen[i, j, 0, l] := 1.0;
				for k := 0 to 2 do 
				    pen[j, i, k, l] := pen[i, j, k, l]
			    end;
			readln(datafile);
			for i := 1 to nallele do 
			    pen[0, i, 0, l] := 1.0;
			if sexlink then begin
			    for i := 1 to nallele do begin
				read(datafile, pen[0, i, 2, l]);
				if (pen[0, i, 2, l] < 0) or (pen[0, i, 2, l] > 1.0) then 
				    inputerror(30, system, system)
			    end;
			    for i := 1 to nallele do 
				pen[0, i, 1, l] := 1.0 - pen[0, i, 2, l];
			    readln(datafile)
			end
		    end
		end
	    end; { getpen }


	    procedure getquan(var locus: locuspoint; privelege: boolean);
	    {Get information of a quantitative locus. Privelege says whether it is
	    a priveleged locus or not}

	    var
		i, j, k: integer;

	    begin
		with locus^ do begin
		    readln(datafile, ntrait);
		    if ntrait > maxtrait then 
			inputerror(31, system, ntrait);
		    if ntrait <= 0 then 
			inputerror(32, system, nclass);
		    for i := 1 to ntrait do 
			for j := 1 to nallele do 
			    for k := j to nallele do begin
				read(datafile, pm[j, k, i]);
				pm[k, j, i] := pm[j, k, i]
			    end;
		    readln(datafile);
		    if not privelege or (nupriv = lastpriv) then begin
			for i := 1 to ntrait do 
			    for j := i to ntrait do begin
				read(datafile, vmat[i, j]);
				if (i = j) and (vmat[i, j] <= 0.0) then 
				    inputerror(33, system, system);
				vmat[j, i] := vmat[i, j]
			    end;
			readln(datafile);
			invert(vmat, ntrait, det);
			det := 1 / sqrt(det);
			readln(datafile, conmat);
			if conmat <= 0.0 then 
			    inputerror(34, system, system);
			conmat := 1 / conmat;
			contrait := 1.0;
			for i := 1 to ntrait do 
			    contrait := contrait * conmat;
			contrait := sqrt(contrait)
		    end
		end
	    end; { getquan }


	begin
	    new(thislocus[system]);
	    {     thislocus[system]:=locuspoint(NewPtr(SizeOf(locusvalues)));}
	    with thislocus[system]^ do begin
		privlocus := nil;
		read(datafile, whichtype, nallele);
		if (whichtype < 0) and (whichtype > 4) then 
		    inputerror(5, system, whichtype);
		if nallele > maxall then 
		    inputerror(6, system, nallele);
		if nallele <= 0 then 
		    inputerror(7, system, nallele);
		case whichtype of
		    0:
			which := quantitative;
		    1:
			which := affection;
		    2:
			which := binary;
		    3:
			which := binary
		end;
		format := whichtype;
		if lastpriv = 0 then 
		    readln(datafile)
		else begin
		    readln(datafile, whichtype);
		    if (whichtype = 0) or (whichtype = 1) then begin
			nupriv := nupriv + 1;
			new(privlocus);
			{           privlocus:=locuspoint(NewPtr(SizeOf(locusvalues)));}
			privlocus^.nallele := nallele;
			with privlocus^ do 
			    case whichtype of
				0:
				    which := quantitative;
				1:
				    which := affection
			    end
		    end
		end;
		if not disequi then begin
		    for j := 1 to nallele do 
			read(datafile, freq[j]);
		    readln(datafile)
		end;
		case which of
		    binary:
			if format = 3 then 
			    getnumber(thislocus[system], system)
			else 
			    getbin(thislocus[system], system);
		    affection:
			getpen(thislocus[system]);
		    quantitative:
			getquan(thislocus[system], false)
		end;
		if privlocus <> nil then 
		    case privlocus^.which of
			affection:
			    getpen(privlocus);
			quantitative:
			    getquan(privlocus, true)
		    end;
		if (nupriv = lastpriv) and (lastpriv <> 0) then 
		    lastpriv := system;
		if risk then begin
		    if i = risksys then 
			readln(datafile, riskall);
		    if riskall > thislocus[i]^.nallele then 
			inputerror(35, i, riskall);
		    if riskall < 0 then 
			inputerror(36, i, riskall)
		end
	    end
	end; { getlocus }


	procedure gettheta(var sex: thetapoint);

	var
	    oldtheta: thetarray;
	    i: integer;

	begin
	    new(sex);
	    {     sex:=thetapoint(NewPtr(SizeOf(thetavalues)));}
	    with sex^ do begin
		if (sex = maletheta) or readfemale then begin
		    for i := 1 to nlocus - 1 do 
			read(datafile, theta[i]);
		    if interfer and not mapping then 
			read(datafile, theta[nlocus]);
		    readln(datafile)
		end else 
		    readln(datafile, distratio);
		{         FOR j:=1 TO maxneed DO segprob[j]:=0.0;}
		if interfer and not mapping then begin
		    oldtheta := theta;
		    theta[1] := (oldtheta[1] + oldtheta[nlocus] - oldtheta[nlocus - 1]) / 2.0;
		    theta[nlocus - 1] := (oldtheta[nlocus - 1] + oldtheta[nlocus] - oldtheta[1]) / 2.0;
		    theta[nlocus] := (oldtheta[1] + oldtheta[nlocus - 1] - oldtheta[nlocus]) / 2.0;
		    for i := 1 to nlocus do 
			if theta[i] > 0.0 then 
			    if theta[i] < 0.999 then 
				theta[i] := ln(1.0 / theta[i] - 1.0)
			    else 
				theta[i] := -6.9 {=ln(1/0.999-1.0)}
			else 
			    theta[i] := 9.21
		end {=ln(1/0.0001-1.0)}
	    end
	end; { gettheta }


    begin {readloci}
	lastpriv := 0;
	readln(datafile, nlocus, risksys, autosomal, i); {i=prog.code}
	{Replace the above line with the next when using epistasis}
	{readln(datafile,nlocus,risksys,autosomal,lastpriv);}
	if i <> 3 then 
	    inputerror(41, nlocus, nlocus);
	if nlocus > maxlocus then 
	    inputerror(0, nlocus, nlocus);
	if nlocus <= 0 then 
	    inputerror(1, nlocus, nlocus);
	if risksys > maxlocus then 
	    inputerror(37, risksys, risksys);
	if risksys < 0 then 
	    inputerror(38, risksys, risksys);
	risk := risksys <> 0;
	sexlink := autosomal = 1;
	write('YOU ARE USING LINKAGE/ISIM (V', version, ') WITH', nlocus: 3, '-POINT');
	if sexlink then 
	    writeln(' SEXLINKED DATA')
	else 
	    writeln(' AUTOSOMAL DATA');
	readln(datafile, mutsys, mutmale, mutfemale, coupling);
	if mutsys > maxlocus then 
	    inputerror(39, mutsys, mutsys);
	if mutsys < 0 then 
	    inputerror(40, mutsys, mutsys);
	if coupling = 1 then 
	    disequi := true
	else 
	    disequi := false;
	if disequi then 
	    new(hapfreq);
	{     hapfreq:=genpoint(NewPtr(SizeOf(thisarray)));}
	for i := 1 to nlocus do begin
	    read(datafile, j);
	    if j > nlocus then 
		inputerror(2, i, j);
	    if j <= 0 then 
		inputerror(3, i, j);
	    order[j] := i
	end;
	for i := 1 to nlocus do 
	    for j := 1 to i - 1 do 
		if order[i] = order[j] then 
		    inputerror(4, i, j);
	readln(datafile);
	if mutsys <> 0 then 
	    mutsys := order[mutsys];
	nupriv := 0;
	for i := 1 to nlocus do 
	    getlocus(order[i]);
	if risksys <> 0 then 
	    risksys := order[risksys];
	increment[nlocus] := 1;
	for i := nlocus - 1 downto 1 do 
	    increment[i] := increment[i + 1] * thislocus[i + 1]^.nallele;
	fgeno := 1;
	for j := 1 to nlocus do 
	    fgeno := thislocus[j]^.nallele * fgeno;
	mgeno := fgeno;
	nuhap := fgeno;
	for i := 1 to nlocus do 
	    nohom[i] := false;
	if disequi then begin
	    for i := 1 to mgeno do 
		read(datafile, hapfreq^.genarray[i]);
	    readln(datafile)
	end else begin
	    for i := 1 to nlocus do 
		with thislocus[i]^ do 
		    if (which = affection) or (which = quantitative) then 
			if freq[affall] < minfreq then 
			    nohom[i] := true
	end;
	fgeno := fgeno * (fgeno + 1) div 2;
	if not sexlink then 
	    mgeno := fgeno;
	read(datafile, difference);
	if (difference < 0) or (difference > 2) then begin
	    inputwarning(0, difference, difference);
	    difference := 0
	end;
	sexdif := difference <> 0;
	readfemale := difference = 2;
	readln(datafile, independent);
	if (independent < 0) or (independent > 2) then begin
	    inputwarning(1, independent, independent);
	    independent := 0
	end;
	interfer := independent <> 0;
	mapping := independent = 2;
	gettheta(maletheta);
	if sexdif then 
	    gettheta(femaletheta)
	else 
	    femaletheta := maletheta;
	if sexlink and sexdif then begin
	    inputwarning(2, difference, difference);
	    sexdif := false;
	    readfemale := false
	end;
	if not sexlink then 
	    if mutsys = 0 then 
		thispath := auto
	    else 
		thispath := mauto
	else if mutsys = 0 then 
	    thispath := sex
	else 
	    thispath := msex;
	segscale := scale;
	for i := 1 to nlocus do 
	    segscale := segscale * scalemult
    end; { readloci }
 {readloci}

begin {inputdata}
    writeln('Reading DATAFILE.DAT');
    reset(datafile, 'datafile.dat'); {SUN}
    writeln('Reading SPEEDFILE.DAT'); {SUN}
    reset(speedfile, 'speedfile.dat'); {SUN}
    writeln('Reading IPEDFILE.DAT');
    reset(ipedfile, 'ipedfile.dat'); {SUN}
    rewrite(outfile, 'outfile.dat'); {SUN}
    rewrite(final, 'final.dat'); {SUN}
    readloci;
    setiterate
end; { inputdata }

procedure readspseg;
{speedfile modification}

var
    i, j, a, b, sys: integer;
    ch: char;

begin				{readspeed}
    i := lastspeed;
    j := lastspeed - lastseg;
    if (j > 0) and (i <= segperson) then begin
	person[j]^.unknown := true;
	new(person[j]^.store);
	{     person[j]^.store:=infoptr(NewPtr(SizeOf(information)));}
	with person[j]^.store^ do 
	    for sys := 1 to nlocus do 
		for a := 1 to maxall do 
		    for b := 1 to maxall do 
			possible[sys, a, b] := false
    end;
    while not eof(speedfile) and (i <= segperson) do begin
	read(speedfile, ch);
	if (ch = 'i') or (ch = 'I') then begin
	    readln(speedfile, ch, i);
	    if i <= segperson then begin
		j := i - lastseg;
		person[j]^.unknown := true;
		new(person[j]^.store);
		{         person[j]^.store:=infoptr(NewPtr(SizeOf(information)));}
		with person[j]^.store^ do 
		    for sys := 1 to nlocus do 
			for a := 1 to maxall do 
			    for b := 1 to maxall do 
				possible[sys, a, b] := false
	    end
	end else 
	    with person[j]^.store^ do begin
		readln(speedfile, sys, a, b);
		if sys <= nlocus then 
		    possible[order[sys], a, b] := true
	    end
    end;
    lastspeed := i;
    lastseg := segperson
end; { readspseg }				{readspeed}


procedure readpedseg;

var
    i, newid, sex, profield, newped, sequence, nuperson, thisone, thisped: integer;
    startped, endped: array [1..maxped] of integer;
    holdloop: ind;

label
    10; {Exit from the procedure readpedseg}

    procedure inputerror(nerror, par1, par2: integer);

    begin
	writeln('Fatal error detected in procedure inputdata');
	case nerror of
	    11:
		writeln('Error detected reading pedigree record ', par1: 2, '. Illegal code for sex ', par2: 2);
	    12:
		writeln('Error detected reading pedigree record at pedigree', par1: 2, '. Maximum number of pedigree records exceeded');
	    13:
		writeln('Error detected reading pedigree record ', par1: 2, '. Maximum number of individuals exceeded');
	    14:
		writeln('Error detected reading pedigree record ', par1: 2, '. Illegal binary factor code ', par2: 2);
	    15:
		writeln('Error detected reading pedigree record ', par1: 2, '. No allelic pair for genotype');
	    16:
		writeln('Error detected reading pedigree record ', par1: 2, '. Allele number ', par2: 2, ' exceeds maxall');
	    17:
		writeln('Error detected reading pedigree record ', par1: 2, '. Illegal allele number ', par2: 2);
	    20:
		writeln('Number of recombination types (', par1: 3, ') exceeds maxrectype');
	    21:
		writeln('Number of recombination types (', par1: 3, ') less than minimum of 1');
	    22:
		writeln('End of file detected in tempdat by procedure readthg before all data found');
	    23:
		writeln('Error detected reading iterated locus in datafile. Value (', i: 3, ') greater than nlocus');
	    24:
		writeln('Error detected reading iterated locus in datafile. Illegal value (', par1: 3, ')');
	    25:
		writeln('Number of iterated parameters greater then maxn');
	    26:
		writeln('Error detected reading pedigree record ', par1: 2, '. Liability class (', par2: 2, ') exceeds nclass');
	    27:
		writeln('Error detected reading pedigree record ', par1: 2, '. Illegal liability class (', par2: 2, ')')
	end
    end; { inputerror }


    procedure inputwarning(nwarning, par1, par2: integer);

    begin
	writeln('Warning number from procedure readpedseg');
	case nwarning of
	    3:
		writeln('Non-standard affection status', par2: 4, ' interpreted as normal in pedigree record', par1: 5)
	end
    end; { inputwarning }



    procedure getphenotype(var p: ind);

    var
	thisread, system: integer;

	procedure readbin(var phen: phenpoint; thislocus: locuspoint);

	var
	    i, j: integer;

	begin
	    with phen^ do begin
		which := binary;
		phenf := [];
		for i := 1 to nfactor[system] do begin
		    read(ipedfile, j);
		    if (j <> 0) and (j <> 1) then 
			inputerror(14, p^.id, j);
		    if j = 1 then 
			phenf := phenf + [i]
		end
	    end
	end; { readbin }


	procedure readnumber(var phen: phenpoint; thislocus: locuspoint);

	var
	    i, j: integer;

	begin
	    with phen^ do begin
		which := binary;
		phenf := [];
		for i := 1 to 2 do begin
		    read(ipedfile, j);
		    if j > maxall then 
			inputerror(16, p^.id, j);
		    if j < 0 then 
			inputerror(17, p^.id, j);
		    if j <> 0 then 
			phenf := phenf + [j]
		end
	    end
	end; { readnumber }


	procedure readaff(var phen: phenpoint; thislocus: locuspoint);

	var
	    thisval: integer;

	begin
	    with phen^ do begin
		which := affection;
		read(ipedfile, thisval);
		if thisval = missaff then 
		    aff := 0
		else if thisval = affval then 
		    aff := 2
		else begin
		    if thisval <> 1 then 
			inputwarning(3, p^.id, thisval);
		    aff := 1
		end;
		if thislocus^.nclass = 1 then 
		    liability := 1
		else 
		    read(ipedfile, liability);
		if liability > thislocus^.nclass then 
		    inputerror(26, p^.id, liability);
		if liability <= 0 then 
		    inputerror(27, p^.id, liability)
	    end
	end; { readaff }


	procedure readquan(var phen: phenpoint; thislocus: locuspoint);

	var
	    i: integer;
	    xval: double;

	begin
	    with phen^ do 
		if not sexlink or not p^.male then begin
		    which := quantitative;
		    with thislocus^ do begin
			for i := 1 to ntrait do 
			    read(ipedfile, x[i]);
			missing := true;
			for i := 1 to ntrait do 
			    if x[i] <> missval then 
				missing := false
		    end
		end else begin
		    which := affection;
		    read(ipedfile, xval);
		    with thislocus^ do begin
			if xval = missval then 
			    aff := missaff
			else if xval = affall then 
			    aff := affall
			else 
			    aff := -11;
			liability := 1;
			for i := 2 to ntrait do 
			    read(ipedfile, xval)
		    end
		end
	end; { readquan }


    begin
	with p^ do 
	    for thisread := 1 to nlocus do begin
		system := order[thisread];
		phen[system] := nil;
		if thislocus[system]^.which <> null then 
		    new(phen[system]);
		{        phen[system]:=phenpoint(NewPtr(SizeOf(phenotype)));}
		case thislocus[system]^.which of
		    quantitative:
			readquan(phen[system], thislocus[system]);
		    affection:
			readaff(phen[system], thislocus[system]);
		    binary:
			if thislocus[system]^.format = 3 then 
			    readnumber(phen[system], thislocus[system])
			else 
			    readbin(phen[system], thislocus[system])
		end;
		if lastpriv = system then begin
		    new(privphen);
		    {         privphen:=phenpoint(NewPtr(SizeOf(phenotype)));}
		    case thislocus[system]^.privlocus^.which of
			quantitative:
			    readquan(privphen, thislocus[system]^.privlocus);
			affection:
			    readaff(privphen, thislocus[system]^.privlocus)
		    end
		end
	    end
    end; { getphenotype }


    procedure getinformative;

    var
	i, j, k, l, m, count, nchild: integer;
	child: ind;

    begin
	if fitmodel or risk then 
	    for i := 1 to nuped do 
		informative[i] := true
	else 
	    for i := 1 to nuped do begin
		informative[i] := false;
		for j := startped[i] to endped[i] do 
		    if person[j]^.foff <> nil then begin
			nchild := 0;
			child := person[j]^.foff;
			repeat
			    nchild := nchild + 1;
			    if person[j]^.male then 
				child := child^.nextpa
			    else 
				child := child^.nextma
			until child = nil;
			count := 0;
			if (nchild > 1) or (nchild = 1) and (person[j]^.pa <> nil) then 
			    for k := 1 to nlocus do 
				with person[j]^.phen[k]^ do 
				    with thislocus[k]^ do 
					if which <> binary then 
					    count := count + 1
					else if phenf = [] then 
					    count := count + 1
					else begin
					    l := 0;
					    for m := 1 to nallele do 
						if m in phenf then 
						    l := l + 1;
					    if l > 1 then 
						count := count + 1
					end;
			if count > 1 then 
			    informative[i] := true
		    end
	    end
    end; { getinformative }


    procedure getind(var id: integer);

    begin {getind}
	read(ipedfile, id);
	if id <> 0 then begin
	    id := id + sequence;
	    if id > maxind then 
		inputerror(13, id, id);
	    if person[id] = nil then 
		new(person[id])
		{       person[id]:=ind(NewPtr(SizeOf(thisperson)));}
	end
    end; { getind }
						{getind}

    procedure multimarriage(var p: ind);

    var
	q, child: ind;

    begin {multimarriage}
	if p^.foff <> nil then 
	    with p^ do begin
		if male then 
		    q := foff^.ma
		else 
		    q := foff^.pa;
		child := foff;
		p^.multi := false;
		repeat
		    if male then begin
			multi := q = child^.ma;
			child := child^.nextpa
		    end else begin
			multi := q = child^.pa;
			child := child^.nextma
		    end
		until (child = nil) or multi
	    end
	else 
	    p^.multi := false
    end; { multimarriage }
				{multimarriage}

begin				{readpedseg}
    for i := 0 to maxind do 
	person[i] := nil;
    sequence := 0;
    nuperson := 0;
    nuped := 1;
    for i := 1 to maxloop do begin
	looppers[nuped, i, 1] := nil;
	looppers[nuped, i, 2] := nil
    end;
    proband[nuped] := nil;
    if infile then 
	newped := oldped
    else 
	read(ipedfile, newped);
    thisped := newped;
    startped[1] := 1;
    while not eof(ipedfile) do begin
	nuperson := nuperson + 1;
	getind(thisone);
	if proband[nuped] = nil then 
	    proband[nuped] := person[thisone];
	with person[thisone]^ do begin
	    ped := thisped;
	    id := thisone;
	    getind(newid);
	    pa := person[newid];
	    getind(newid);
	    ma := person[newid];
	    getind(newid);
	    foff := person[newid];
	    getind(newid);
	    nextpa := person[newid];
	    getind(newid);
	    nextma := person[newid];
	    read(ipedfile, sex);
	    if (sex <> 1) and (sex <> 2) then 
		inputerror(11, id, sex);
	    if sex = 1 then 
		male := true
	    else 
		male := false;
	    unknown := false;
	    inloop := 0
	end;
	read(ipedfile, profield);
	if profield = 1 then 
	    proband[nuped] := person[thisone]
	else if (profield > 1) and (profield - 1 <= maxloop) then 
	    if looppers[nuped, profield - 1, 2] = nil then 
		looppers[nuped, profield - 1, 2] := person[thisone]
	    else 
		looppers[nuped, profield - 1, 1] := person[thisone];
	getphenotype(person[thisone]);
	readln(ipedfile);
	if not eof(ipedfile) then 
	    read(ipedfile, newped);
	{We only know that we have reached the end of the current pedigree
	by the fact that the pedigree number changes i.e. newped<>thisped
	So we have to read one more record than I wish to read.
	Test for: newped<>thisped and nuped=oped => save newped in oldped
	and don't read it on the next call to readpedseg}
	if (newped <> thisped) and (nuped = opeds) then begin
	    oldped := newped;
	    infile := true;
	    goto 10
	end {Exit the procedure readpedseg};
	if newped <> thisped then begin
	    sequence := nuperson + sequence;
	    endped[nuped] := sequence;
	    nuperson := 0;
	    nuped := nuped + 1;
	    if nuped > maxped then 
		inputerror(12, newped, nuped);
	    startped[nuped] := sequence + 1;
	    for i := 1 to maxloop do begin
		looppers[nuped, i, 1] := nil;
		looppers[nuped, i, 2] := nil
	    end;
	    proband[nuped] := nil;
	    thisped := newped
	end
    end; {while not eof() loop}
10:
    totperson := sequence + nuperson;
    endped[nuped] := totperson;
    for newped := 1 to nuped do begin
	if (looppers[newped, 1, 2] <> nil) and (looppers[newped, 1, 1] = nil) then 
	    looppers[newped, 1, 1] := proband[newped];
	for i := 1 to maxloop do 
	    if looppers[newped, i, 1] = nil then 
		looppers[newped, i, 2] := nil
	    else begin
		looppers[newped, i, 1]^.inloop := i;
		looppers[newped, i, 2]^.inloop := i;
		if (looppers[newped, i, 1]^.pa = nil) and (looppers[newped, i, 2]^.pa <> nil) then begin
		    holdloop := looppers[newped, i, 1];
		    looppers[newped, i, 1] := looppers[newped, i, 2];
		    looppers[newped, i, 2] := holdloop
		end
	    end
    end;
    {speedfile modification}
    segperson := segperson + totperson;
    for thisone := 1 to totperson do 
	multimarriage(person[thisone]);
    getinformative
end; { readpedseg }				{readpedseg}


procedure cleanup(var p: ind);

begin				{cleanup}
    with p^ do begin
	dispose(gen);
	{     DisposPtr(Ptr(gen));}
	gen := nil
    end
end; { cleanup }				{cleanup}

procedure getvect(p: ind);

var
    hap1, hap2: hapvector;

    function quanfun(phen: phenpoint; thislocus: locuspoint; i, j: integer; var mean: thesemeans): double;

    var
	val: double;
	it, jt: integer;

    begin {quanfun}
	with phen^ do 
	    with thislocus^ do begin
		val := 1.0;
		if not missing then begin
		    val := 0;
		    for it := 1 to ntrait do 
			for jt := 1 to ntrait do 
			    if i = j then 
				val := val + vmat[it, jt] * (x[jt] - mean[jt]) * (x[it] - mean[it])
			    else 
				val := val + conmat * vmat[it, jt] * (x[jt] - mean[jt]) * (x[it] - mean[it]);
		    val := det * exp(-val * 0.5);
		    if i <> j then 
			val := val * contrait
		end;
		quanfun := val
	    end
    end; { quanfun }
						{quanfun}

    procedure getval(syste, i, j: integer; var val: double);

    begin {getval}
	with thislocus[syste]^ do 
	    case which of
		quantitative:
		    val := val * quanfun(p^.phen[syste], thislocus[syste], i, j, pm[i, j]);
		affection:
		    with p^.phen[syste]^ do 
			val := val * pen[i, j, aff, liability]
	    end
    end; { getval }
						{getval}

    procedure setval(val: double);

    var
	here, count, nhap1, nhap2, i: integer;

	procedure prior;

	var
	    i: integer; {prior}

	begin
	    if not disequi then 
		if sexlink and p^.male then 
		    for i := 1 to nlocus do 
			with thislocus[i]^ do 
			    val := val * freq[hap1[i]]
		else begin
		    for i := 1 to nlocus do 
			with thislocus[i]^ do 
			    val := val * freq[hap1[i]] * freq[hap2[i]];
		    if nhap1 <> nhap2 then 
			val := 2.0 * val
		end
	    else 
		with hapfreq^ do 
		    if sexlink and p^.male then 
			val := val * genarray[nhap1]
		    else begin
			val := val * genarray[nhap1] * genarray[nhap2];
			if nhap1 <> nhap2 then 
			    val := 2.0 * val
		    end;
	    val := val * segscale
	end; { prior }
 {prior}

    begin {setval}
	count := 1;
	nhap1 := 1;
	nhap2 := 1;
	if sexlink and p^.male then begin
	    for i := 1 to nlocus do 
		nhap1 := nhap1 + increment[i] * (hap1[i] - 1);
	    here := nhap1
	end else 
	    with gennustruct^ do 
		for i := 1 to nlocus do begin
		    nhap1 := nhap1 + increment[i] * (hap1[i] - 1);
		    nhap2 := nhap2 + increment[i] * (hap2[i] - 1);
		    if hap1[i] <> hap2[i] then 
			count := 2 * count;
		    here := genenumber[nhap1, nhap2]
		end;
	if p^.pa = nil then 
	    prior;
	with p^.gen^ do 
	    if not disequi then begin
		if count <> 1 then 
		    count := count div 2;
		for i := 1 to count do begin
		    genarray[here] := val;
		    here := here + 1
		end
	    end else 
		genarray[here] := val
    end; { setval }
						{setval}

    procedure getgene(syste: integer; val: double);

    var
	newval: double;

	procedure facmale;

	var
	    j: integer; {facmale}

	begin
	    with p^ do 
		with thislocus[syste]^ do 
		    for j := 1 to nallele do 
			if (phen[syste]^.phenf = allele[j]) or (phen[syste]^.phenf = []) then begin
			    hap1[syste] := j;
			    if syste <> nlocus then 
				getgene(syste + 1, val)
			    else 
				setval(val)
			end
	end; { facmale }
 {facmale}

	procedure affmale;

	var
	    j: integer; {affmale}

	begin
	    with thislocus[syste]^ do 
		for j := 1 to nallele do begin
		    newval := val;
		    getval(syste, 0, j, newval);
		    hap1[syste] := j;
		    if newval <> 0.0 then 
			if syste <> nlocus then 
			    getgene(syste + 1, newval)
			else 
			    setval(newval)
		end
	end; { affmale }
 {affmale}

	procedure quanmale;

	var
	    j: integer; {quanmale}

	begin
	    with p^ do 
		with thislocus[syste]^ do begin
		    if (phen[syste]^.aff = affall) or (phen[syste]^.aff = missaff) then begin
			newval := val;
			hap1[syste] := affall;
			if newval <> 0.0 then 
			    if syste <> nlocus then 
				getgene(syste + 1, newval)
			    else 
				setval(newval)
		    end;
		    if (phen[syste]^.aff <> affall) or (phen[syste]^.aff = missaff) then 
			for j := 1 to nallele do 
			    if j <> affall then begin
				newval := val;
				hap1[syste] := j;
				if newval <> 0.0 then 
				    if syste <> nlocus then 
					getgene(syste + 1, newval)
				    else 
					setval(newval)
			    end
		end
	end; { quanmale }
 {quanmale}

	procedure fac;

	var
	    i, j: integer; {fac}

	begin
	    with p^ do 
		with thislocus[syste]^ do 
		    for i := 1 to nallele do begin
			hap1[syste] := i;
			for j := i to nallele do 
			    if (phen[syste]^.phenf = allele[i] + allele[j]) or (phen[syste]^.phenf = []) then begin
				hap2[syste] := j;
				if syste <> nlocus then 
				    getgene(syste + 1, val)
				else 
				    setval(val)
			    end
		    end;
	    if disequi then 
		with p^ do 
		    with thislocus[syste]^ do 
			for i := 1 to nallele do begin
			    hap2[syste] := i;
			    for j := i to nallele do 
				if (phen[syste]^.phenf = allele[i] + allele[j]) or (phen[syste]^.phenf = []) then begin
				    hap1[syste] := j;
				    if syste <> nlocus then 
					getgene(syste + 1, val)
				    else 
					setval(val)
				end
			end
	end; { fac }
 {fac}

	procedure aff;
	{Used with an affection status phenotype or when
	thislocus[syste]^which is null}

	var
	    i, j: integer;

	begin
	    with thislocus[syste]^ do 
		for i := 1 to nallele do begin
		    hap1[syste] := i;
		    for j := i to nallele do 
			if not nohom[syste] or (i <> affall) or (j <> affall) then begin
			    hap2[syste] := j;
			    newval := val;
			    getval(syste, i, j, newval);
			    if newval <> 0.0 then 
				if syste <> nlocus then 
				    getgene(syste + 1, newval)
				else 
				    setval(newval)
			end
		end;
	    if disequi then 
		with thislocus[syste]^ do 
		    for i := 1 to nallele do begin
			hap2[syste] := i;
			for j := i to nallele do 
			    if not nohom[syste] or (i <> affall) or (j <> affall) then begin
				hap1[syste] := j;
				newval := val;
				getval(syste, i, j, newval);
				if newval <> 0.0 then 
				    if syste <> nlocus then 
					getgene(syste + 1, newval)
				    else 
					setval(newval)
			    end
		    end
	end; { aff }


	procedure quanval;
	{Uses this only when thislocus[syste]^.which is not null}

	var
	    i, j: integer;

	begin
	    with thislocus[syste]^ do 
		for i := 1 to nallele do begin
		    hap1[syste] := i;
		    for j := i to nallele do 
			if not nohom[syste] or (i <> affall) or (j <> affall) then begin
			    hap2[syste] := j;
			    newval := val;
			    getval(syste, i, j, newval);
			    if newval <> 0.0 then 
				if syste <> nlocus then 
				    getgene(syste + 1, newval)
				else 
				    setval(newval)
			end
		end;
	    if disequi then 
		with thislocus[syste]^ do 
		    for i := 1 to nallele do begin
			hap2[syste] := i;
			for j := i to nallele do 
			    if not nohom[syste] or (i <> affall) or (j <> affall) then begin
				hap1[syste] := j;
				newval := val;
				getval(syste, i, j, newval);
				if newval <> 0.0 then 
				    if syste <> nlocus then 
					getgene(syste + 1, newval)
				    else 
					setval(newval)
			    end
		    end
	end; { quanval }


    begin
	with thislocus[syste]^ do 
	    if sexlink and p^.male then 
		case which of
		    binary:
			facmale;
		    affection:
			affmale;
		    quantitative:
			quanmale;
		    null:
			if privlocus^.which = affection then 
			    affmale
			else 
			    quanmale
		end
	    else 
		case which of
		    binary:
			fac;
		    affection:
			aff;
		    quantitative:
			quanval;
		    null:
			aff
		end
    end; { getgene }


    procedure ugetgene(syste: integer; val: double);

    var
	newval: double;

	procedure facmale;

	var
	    j: integer; {facmale}

	begin
	    with p^ do 
		with p^.store^ do 
		    with thislocus[syste]^ do 
			for j := 1 to nallele do 
			    if (phen[syste]^.phenf = allele[j]) or (phen[syste]^.phenf = []) then 
				if possible[syste, 1, j] then begin
				    hap1[syste] := j;
				    if syste <> nlocus then 
					ugetgene(syste + 1, val)
				    else 
					setval(val)
				end
	end; { facmale }
 {facmale}

	procedure affmale;

	var
	    j: integer;

	begin
	    with p^.store^ do 
		with thislocus[syste]^ do 
		    for j := 1 to nallele do 
			if possible[syste, 1, j] then begin
			    newval := val;
			    getval(syste, 0, j, newval);
			    hap1[syste] := j;
			    if newval <> 0.0 then 
				if syste <> nlocus then 
				    ugetgene(syste + 1, newval)
				else 
				    setval(newval)
			end
	end; { affmale }


	procedure quanmale;

	var
	    j: integer;

	begin
	    with p^ do 
		with p^.store^ do 
		    with thislocus[syste]^ do begin
			if (phen[syste]^.aff = affall) or (phen[syste]^.aff = missaff) then 
			    if possible[syste, 1, affall] then begin
				newval := val;
				hap1[syste] := affall;
				if newval <> 0.0 then 
				    if syste <> nlocus then 
					ugetgene(syste + 1, newval)
				    else 
					setval(newval)
			    end;
			with p^.store^ do 
			    if (phen[syste]^.aff <> affall) or (phen[syste]^.aff = missaff) then 
				for j := 1 to nallele do 
				    if j <> affall then 
					if possible[syste, 1, j] then begin
					    newval := val;
					    hap1[syste] := j;
					    if newval <> 0.0 then 
						if syste <> nlocus then 
						    ugetgene(syste + 1, newval)
						else 
						    setval(newval)
					end
		    end
	end; { quanmale }


	procedure fac;

	var
	    i, j: integer;

	begin
	    with p^ do 
		with p^.store^ do 
		    with thislocus[syste]^ do 
			for i := 1 to nallele do begin
			    hap1[syste] := i;
			    for j := i to nallele do 
				if (phen[syste]^.phenf = allele[i] + allele[j]) or (phen[syste]^.phenf = []) then 
				    if possible[syste, i, j] then begin
					hap2[syste] := j;
					if syste <> nlocus then 
					    ugetgene(syste + 1, val)
					else 
					    setval(val)
				    end
			end;
	    if disequi then 
		with p^ do 
		    with p^.store^ do 
			with thislocus[syste]^ do 
			    for i := 1 to nallele do begin
				hap2[syste] := i;
				for j := i to nallele do 
				    if (phen[syste]^.phenf = allele[i] + allele[j]) or (phen[syste]^.phenf = []) then 
					if possible[syste, i, j] then begin
					    hap1[syste] := j;
					    if syste <> nlocus then 
						ugetgene(syste + 1, val)
					    else 
						setval(val)
					end
			    end
	end; { fac }


	procedure aff;
	{Used with an affection status phenotype or when
	thislocus[syste]^which is null}

	var
	    i, j: integer;

	begin
	    with p^ do 
		with p^.store^ do 
		    with thislocus[syste]^ do 
			for i := 1 to nallele do begin
			    hap1[syste] := i;
			    for j := i to nallele do 
				if not nohom[syste] or (i <> affall) or (j <> affall) then 
				    if possible[syste, i, j] then begin
					hap2[syste] := j;
					newval := val;
					getval(syste, i, j, newval);
					if newval <> 0.0 then 
					    if syste <> nlocus then 
						ugetgene(syste + 1, newval)
					    else 
						setval(newval)
				    end
			end;
	    if disequi then 
		with p^ do 
		    with p^.store^ do 
			with thislocus[syste]^ do 
			    for i := 1 to nallele do begin
				hap2[syste] := i;
				for j := i to nallele do 
				    if not nohom[syste] or (i <> affall) or (j <> affall) then 
					if possible[syste, i, j] then begin
					    hap1[syste] := j;
					    newval := val;
					    getval(syste, i, j, newval);
					    if newval <> 0.0 then 
						if syste <> nlocus then 
						    ugetgene(syste + 1, newval)
						else 
						    setval(newval)
					end
			    end
	end; { aff }


	procedure quanval;
	{Uses this only when thislocus[syste]^.which is not null}

	var
	    i, j: integer;

	begin
	    with p^ do 
		with p^.store^ do 
		    with thislocus[syste]^ do 
			for i := 1 to nallele do begin
			    hap1[syste] := i;
			    for j := i to nallele do 
				if not nohom[syste] or (i <> affall) or (j <> affall) then 
				    if possible[syste, i, j] then begin
					hap2[syste] := j;
					newval := val;
					getval(syste, i, j, newval);
					if newval <> 0.0 then 
					    if syste <> nlocus then 
						ugetgene(syste + 1, newval)
					    else 
						setval(newval)
				    end
			end;
	    if disequi then 
		with p^ do 
		    with p^.store^ do 
			with thislocus[syste]^ do 
			    for i := 1 to nallele do begin
				hap2[syste] := i;
				for j := i to nallele do 
				    if not nohom[syste] or (i <> affall) or (j <> affall) then 
					if possible[syste, i, j] then begin
					    hap1[syste] := j;
					    newval := val;
					    getval(syste, i, j, newval);
					    if newval <> 0.0 then 
						if syste <> nlocus then 
						    ugetgene(syste + 1, newval)
						else 
						    setval(newval)
					end
			    end
	end; { quanval }


    begin
	with thislocus[syste]^ do 
	    if sexlink and p^.male then 
		case which of
		    binary:
			facmale;
		    affection:
			affmale;
		    quantitative:
			quanmale;
		    null:
			if privlocus^.which = affection then 
			    affmale
			else 
			    quanmale
		end
	    else 
		case which of
		    binary:
			fac;
		    affection:
			aff;
		    quantitative:
			quanval;
		    null:
			aff
		end
    end; { ugetgene }


begin
    if p^.unknown then 
	ugetgene(1, 1.0)
    else 
	getgene(1, 1.0)
end; { getvect }


procedure likelihood(thisped: integer; proband: ind);

var
    loopmax, loopgen: array [1..maxloop] of integer;
    homo, hetero, tmplike: double;
    i, j, nuscale: integer;
    holdpoint: array [1..maxloop] of genpoint;
    gocalc, alldone: boolean;

    procedure prob(var p: ind);

    var
	i: integer;

    begin {prob}
	with p^ do 
	    if gen = nil then begin
		new(gen);
		{       gen:=genpoint(NewPtr(SizeOf(thisarray)));}
		with gen^ do 
		    for i := 1 to fgeno do 
			genarray[i] := 0.0;
		if inloop > 0 then 
		    if looppers[thisped, inloop, 1] = p then 
			gen^.genarray[loopgen[inloop]] := holdpoint[inloop]^.genarray[loopgen[inloop]]
		    else 
			gen^.genarray[loopgen[inloop]] := 1.0
		else begin
		    getvect(p);
		    if p^.pa = nil then 
			nuscale := nuscale + 1
		end
	    end
    end; { prob }
						{prob}

    procedure seg(var p, q, r: ind; peel: direction);

    var
	child, father, mother: ind;
	phaseunkn: boolean;
	fseg, sseg, sstart, send, fstart, fend, nfirst, nsecond, firstseg, secondseg: integer;
	pf, ps: double;
	firstsex, secondsex: thetapoint;

	procedure getapprox;

	var
	    first: integer;
	    maxval: double; {getapprox}

	begin
	    maxval := p^.gen^.genarray[1];
	    with p^.gen^ do 
		for first := 1 to fgeno do 
		    if genarray[first] > maxval then 
			maxval := genarray[first];
	    with approxstruct^ do 
		with p^.gen^ do 
		    for first := 1 to fgeno do 
			approxarray[thisped, first] := genarray[first] > maxval * epsilon;
	    if not lasttime then 
		with approxstruct^ do 
		    with p^.gen^ do 
			for first := 1 to fgeno do 
			    if not approxarray[thisped, first] then 
				genarray[first] := 0.0
	end; { getapprox }
 {getapprox}

	function msegsex: double;

	var
	    g1, g2, g3, g4, g5, g6, g7, g8, ms, mf, ms1, ms2, mf1, mf2, j, k, f1, f2, s1, s2: integer;
	    val, temp2: double;
	    temp: array [1..maxchild] of double; {msegsex}

	begin
	    for k := 1 to nchild do 
		temp[k] := 0.0;
	    with gennustruct^ do 
		if p^.male then begin
		    mf := muthap[fseg];
		    secondseg := sseg;
		    with secondsex^ do 
			for j := sstart to send do begin
			    temp2 := segprob[j];
			    s1 := seghap1[secondseg];
			    s2 := seghap2[secondseg];
			    ms1 := muthap[s1];
			    ms2 := muthap[s2];
			    if s1 <> s2 then begin
				g1 := genenumber[fseg, s1];
				g2 := genenumber[fseg, s2];
				g3 := genenumber[fseg, ms1];
				g4 := genenumber[fseg, ms2];
				g5 := genenumber[mf, s1];
				g6 := genenumber[mf, s2];
				g7 := genenumber[mf, ms1];
				g8 := genenumber[mf, ms2];
				for k := 1 to nchild do 
				    with thischild[k]^ do 
					if malechild[k] then 
					    temp[k] := temp[k] + temp2 * ((1 - ps) * (genarray[s1] + genarray[s2]) + ps * (genarray[ms1] + genarray[ms2]))
					else 
					    temp[k] := temp[k] + temp2 * ((1 - pf) * (1 - ps) * (genarray[g1] + genarray[g2]) + (1 - pf) * ps * (genarray[g3] + genarray[g4]) + pf * (1 - ps) * (genarray[g5] + genarray[g6]) + pf * ps * (genarray[g7] + genarray[g8]))
			    end else begin
				g1 := genenumber[fseg, s1];
				g3 := genenumber[fseg, ms1];
				g5 := genenumber[mf, s1];
				g7 := genenumber[mf, ms1];
				for k := 1 to nchild do 
				    with thischild[k]^ do 
					if malechild[k] then 
					    temp[k] := temp[k] + 2.0 * temp2 * ((1 - ps) * genarray[s1] + ps * genarray[ms1])
					else 
					    temp[k] := temp[k] + 2.0 * temp2 * ((1 - pf) * (1 - ps) * genarray[g1] + ps * (1 - pf) * genarray[g3] + pf * (1 - ps) * genarray[g5] + pf * ps * genarray[g7])
			    end;
			    secondseg := secondseg + 1
			end
		end else begin
		    firstseg := fseg;
		    ms := muthap[sseg];
		    with firstsex^ do 
			for j := fstart to fend do begin
			    temp2 := segprob[j];
			    f1 := seghap1[firstseg];
			    f2 := seghap2[firstseg];
			    mf1 := muthap[f1];
			    mf2 := muthap[f2];
			    if f1 <> f2 then begin
				g1 := genenumber[sseg, f1];
				g2 := genenumber[sseg, f2];
				g3 := genenumber[sseg, mf1];
				g4 := genenumber[sseg, mf2];
				g5 := genenumber[ms, f1];
				g6 := genenumber[ms, f2];
				g7 := genenumber[ms, mf1];
				g8 := genenumber[ms, mf2];
				for k := 1 to nchild do 
				    with thischild[k]^ do 
					if malechild[k] then 
					    temp[k] := temp[k] + temp2 * ((1 - pf) * (genarray[f1] + genarray[f2]) + pf * (genarray[mf1] + genarray[mf2]))
					else 
					    temp[k] := temp[k] + temp2 * ((1 - pf) * (1 - ps) * (genarray[g1] + genarray[g2]) + (1 - ps) * pf * (genarray[g3] + genarray[g4]) + ps * (1 - pf) * (genarray[g5] + genarray[g6]) + pf * ps * (genarray[g7] + genarray[g8]))
			    end else begin
				g1 := genenumber[sseg, f1];
				g3 := genenumber[sseg, mf1];
				g5 := genenumber[ms, f1];
				g7 := genenumber[ms, mf1];
				for k := 1 to nchild do 
				    with thischild[k]^ do 
					if malechild[k] then 
					    temp[k] := temp[k] + 2.0 * temp2 * ((1 - pf) * genarray[f1] + pf * genarray[mf1])
					else 
					    temp[k] := temp[k] + 2.0 * temp2 * ((1 - pf) * (1 - ps) * genarray[g1] + pf * (1 - ps) * genarray[g3] + ps * (1 - pf) * genarray[g5] + pf * ps * genarray[g7])
			    end;
			    firstseg := firstseg + 1
			end
		end;
	    val := 1.0;
	    for k := 1 to nchild do 
		val := val * temp[k];
	    msegsex := val
	end; { msegsex }
 {msegsex}

	function msegsexf: double;

	var
	    g1, g2, g3, g4, g5, g6, g7, g8, mf, ms1, ms2, j, k, l, s1, s2, slength: integer;
	    val, temp1: double;
	    temp: array [1..maxchild, 1..maxseg] of double;
	    temp2: array [1..maxseg] of double; {msegsexf}

	begin
	    mf := muthap[fseg];
	    slength := send - sstart + 1;
	    for k := 1 to nchild do 
		for l := 1 to slength do 
		    temp[k, l] := 0.0;
	    secondseg := sseg;
	    with gennustruct^ do 
		with secondsex^ do 
		    for j := sstart to send do begin
			for l := 1 to slength do 
			    temp2[l] := segprob[j + (l - 1) * slength];
			s1 := seghap1[secondseg];
			s2 := seghap2[secondseg];
			ms1 := muthap[s1];
			ms2 := muthap[s2];
			if s1 <> s2 then begin
			    g1 := genenumber[fseg, s1];
			    g2 := genenumber[fseg, s2];
			    g3 := genenumber[fseg, ms1];
			    g4 := genenumber[fseg, ms2];
			    g5 := genenumber[mf, s1];
			    g6 := genenumber[mf, s2];
			    g7 := genenumber[mf, ms1];
			    g8 := genenumber[mf, ms2];
			    for k := 1 to nchild do begin
				with thischild[k]^ do 
				    if malechild[k] then 
					val := (1 - ps) * (genarray[s1] + genarray[s2]) + ps * (genarray[ms1] + genarray[ms2])
				    else 
					val := (1 - pf) * (1 - ps) * (genarray[g1] + genarray[g2]) + (1 - pf) * ps * (genarray[g3] + genarray[g4]) + pf * (1 - ps) * (genarray[g5] + genarray[g6]) + pf * ps * (genarray[g7] + genarray[g8]);
				for l := 1 to slength do 
				    temp[k, l] := temp[k, l] + temp2[l] * val
			    end
			end else begin
			    g1 := genenumber[fseg, s1];
			    g3 := genenumber[fseg, ms1];
			    g5 := genenumber[mf, s1];
			    g7 := genenumber[mf, ms1];
			    for k := 1 to nchild do begin
				with thischild[k]^ do 
				    if malechild[k] then 
					val := 2.0 * ((1 - ps) * genarray[s1] + ps * genarray[ms1])
				    else 
					val := 2.0 * ((1 - pf) * (1 - ps) * genarray[g1] + ps * (1 - pf) * genarray[g3] + pf * (1 - ps) * genarray[g5] + pf * ps * genarray[g7]);
				for l := 1 to slength do 
				    temp[k, l] := temp[k, l] + temp2[l] * val
			    end
			end;
			secondseg := secondseg + 1
		    end;
	    temp1 := 0.0;
	    for l := 1 to slength do begin
		val := 1.0;
		for k := 1 to nchild do 
		    val := val * temp[k, l];
		temp1 := temp1 + val
	    end;
	    msegsexf := temp1
	end; { msegsexf }
 {msegsexf}

	function segsex: double;

	var
	    g1, g2, j, k, f1, f2, s1, s2: integer;
	    val, temp2: double;
	    temp: array [1..maxchild] of double; {segsex}

	begin
	    for k := 1 to nchild do 
		temp[k] := 0.0;
	    with gennustruct^ do 
		if p^.male then begin
		    secondseg := sseg;
		    with secondsex^ do 
			for j := sstart to send do begin
			    temp2 := segprob[j];
			    s1 := seghap1[secondseg];
			    s2 := seghap2[secondseg];
			    if s1 <> s2 then begin
				g1 := genenumber[fseg, s1];
				g2 := genenumber[fseg, s2];
				for k := 1 to nchild do 
				    with thischild[k]^ do 
					if malechild[k] then 
					    temp[k] := temp[k] + temp2 * (genarray[s1] + genarray[s2])
					else 
					    temp[k] := temp[k] + temp2 * (genarray[g1] + genarray[g2])
			    end else begin
				g1 := genenumber[fseg, s1];
				for k := 1 to nchild do 
				    with thischild[k]^ do 
					if malechild[k] then 
					    temp[k] := temp[k] + 2.0 * temp2 * genarray[s1]
					else 
					    temp[k] := temp[k] + 2.0 * temp2 * genarray[g1]
			    end;
			    secondseg := secondseg + 1
			end
		end else begin
		    firstseg := fseg;
		    with firstsex^ do 
			for j := fstart to fend do begin
			    temp2 := segprob[j];
			    f1 := seghap1[firstseg];
			    f2 := seghap2[firstseg];
			    if f1 <> f2 then begin
				g1 := genenumber[sseg, f1];
				g2 := genenumber[sseg, f2];
				for k := 1 to nchild do 
				    with thischild[k]^ do 
					if malechild[k] then 
					    temp[k] := temp[k] + temp2 * (genarray[f1] + genarray[f2])
					else 
					    temp[k] := temp[k] + temp2 * (genarray[g1] + genarray[g2])
			    end else begin
				g1 := genenumber[sseg, f1];
				for k := 1 to nchild do 
				    with thischild[k]^ do 
					if malechild[k] then 
					    temp[k] := temp[k] + 2.0 * temp2 * genarray[f1]
					else 
					    temp[k] := temp[k] + 2.0 * temp2 * genarray[g1]
			    end;
			    firstseg := firstseg + 1
			end
		end;
	    val := 1.0;
	    for k := 1 to nchild do 
		val := val * temp[k];
	    segsex := val
	end; { segsex }
 {segsex}

	function segsexf: double;

	var
	    g1, g2, j, k, l, s1, s2, slength: integer;
	    val, temp1: double;
	    temp: array [1..maxchild, 1..maxseg] of double;
	    temp2: array [1..maxseg] of double;

	begin
	    slength := send - sstart + 1;
	    for k := 1 to nchild do 
		for l := 1 to slength do 
		    temp[k, l] := 0.0;
	    secondseg := sseg;
	    with gennustruct^ do 
		with secondsex^ do 
		    for j := sstart to send do begin
			for l := 1 to slength do 
			    temp2[l] := segprob[j + (l - 1) * slength];
			s1 := seghap1[secondseg];
			s2 := seghap2[secondseg];
			if s1 <> s2 then begin
			    g1 := genenumber[fseg, s1];
			    g2 := genenumber[fseg, s2];
			    for k := 1 to nchild do 
				with thischild[k]^ do begin
				    if malechild[k] then 
					val := genarray[s1] + genarray[s2]
				    else 
					val := genarray[g1] + genarray[g2];
				    for l := 1 to slength do 
					temp[k, l] := temp[k, l] + temp2[l] * val
				end
			end else begin
			    g1 := genenumber[fseg, s1];
			    for k := 1 to nchild do 
				with thischild[k]^ do begin
				    if malechild[k] then 
					val := 2.0 * genarray[s1]
				    else 
					val := 2.0 * genarray[g1];
				    for l := 1 to slength do 
					temp[k, l] := temp[k, l] + temp2[l] * val
				end
			end;
			secondseg := secondseg + 1
		    end;
	    temp1 := 0.0;
	    for l := 1 to slength do begin
		val := 1.0;
		for k := 1 to nchild do 
		    val := val * temp[k, l];
		temp1 := temp1 + val
	    end;
	    segsexf := temp1
	end; { segsexf }


	function segfun: double;

	var
	    g1, g2, g3, g4, i, j, k, f1, f2, s1, s2: integer;
	    val, temp1, temp2: double;
	    temp: array [1..maxchild] of double;

	begin
	    firstseg := fseg;
	    for k := 1 to nchild do 
		temp[k] := 0.0;
	    with gennustruct^ do 
		with firstsex^ do 
		    for i := fstart to fend do begin
			temp1 := segprob[i];
			f1 := seghap1[firstseg];
			f2 := seghap2[firstseg];
			secondseg := sseg;
			with secondsex^ do 
			    if f1 <> f2 then 
				for j := sstart to send do begin
				    temp2 := temp1 * segprob[j];
				    s1 := seghap1[secondseg];
				    s2 := seghap2[secondseg];
				    if s1 <> s2 then begin
					g1 := genenumber[f1, s1];
					g2 := genenumber[f1, s2];
					g3 := genenumber[f2, s1];
					g4 := genenumber[f2, s2];
					for k := 1 to nchild do 
					    with thischild[k]^ do 
						temp[k] := temp[k] + temp2 * (genarray[g1] + genarray[g2] + genarray[g3] + genarray[g4])
				    end else begin
					g1 := genenumber[f1, s1];
					g3 := genenumber[f2, s1];
					for k := 1 to nchild do 
					    with thischild[k]^ do 
						temp[k] := temp[k] + 2 * temp2 * (genarray[g1] + genarray[g3])
				    end;
				    secondseg := secondseg + 1
				end
			    else 
				for j := sstart to send do begin
				    temp2 := temp1 * segprob[j];
				    s1 := seghap1[secondseg];
				    s2 := seghap2[secondseg];
				    if s1 <> s2 then begin
					g1 := genenumber[f1, s1];
					g2 := genenumber[f1, s2];
					for k := 1 to nchild do 
					    with thischild[k]^ do 
						temp[k] := temp[k] + 2 * temp2 * (genarray[g1] + genarray[g2])
				    end else begin
					g1 := genenumber[f1, s1];
					for k := 1 to nchild do 
					    with thischild[k]^ do 
						temp[k] := temp[k] + 4 * temp2 * genarray[g1]
				    end;
				    secondseg := secondseg + 1
				end;
			firstseg := firstseg + 1
		    end;
	    val := 1.0;
	    for k := 1 to nchild do 
		val := val * temp[k];
	    segfun := val
	end; { segfun }


	function msegfast: double;

	var
	    g1, g2, g3, g4, g5, g6, g7, g8, g9, g10, g11, g12, g13, g14, g15, g16, i, j, k, l, f1, f2, s1, s2, ms1, ms2, mf1, mf2, slength: integer;
	    val, temp1: double;
	    temp: array [1..maxchild, 1..maxseg] of double;
	    temp2: array [1..maxseg] of double; {msegfast}

	begin
	    firstseg := fseg;
	    slength := send - sstart + 1;
	    for k := 1 to nchild do 
		for l := 1 to slength do 
		    temp[k, l] := 0.0;
	    with gennustruct^ do 
		with firstsex^ do 
		    for i := fstart to fend do begin
			temp1 := segprob[i];
			f1 := seghap1[firstseg];
			f2 := seghap2[firstseg];
			mf1 := muthap[f1];
			mf2 := muthap[f2];
			secondseg := sseg;
			with secondsex^ do 
			    if f1 <> f2 then 
				for j := sstart to send do begin
				    for l := 1 to slength do 
					temp2[l] := temp1 * segprob[j + (l - 1) * slength];
				    s1 := seghap1[secondseg];
				    s2 := seghap2[secondseg];
				    ms1 := muthap[s1];
				    ms2 := muthap[s2];
				    if s1 <> s2 then begin
					g1 := genenumber[f1, s1];
					g2 := genenumber[f1, s2];
					g3 := genenumber[f2, s1];
					g4 := genenumber[f2, s2];
					g5 := genenumber[f1, ms1];
					g6 := genenumber[f1, ms2];
					g7 := genenumber[f2, ms1];
					g8 := genenumber[f2, ms2];
					g9 := genenumber[mf1, s1];
					g10 := genenumber[mf1, s2];
					g11 := genenumber[mf2, s1];
					g12 := genenumber[mf2, s2];
					g13 := genenumber[mf1, ms1];
					g14 := genenumber[mf1, ms2];
					g15 := genenumber[mf2, ms1];
					g16 := genenumber[mf2, ms2];
					for k := 1 to nchild do 
					    with thischild[k]^ do begin
						val := (1 - ps) * (1 - pf) * (genarray[g1] + genarray[g2] + genarray[g3] + genarray[g4]) + ps * (1 - pf) * (genarray[g5] + genarray[g6] + genarray[g7] + genarray[g8]) + pf * (1 - ps) * (genarray[g9] + genarray[g10] + genarray[g11] + genarray[g12]) + pf * ps * (genarray[g13] + genarray[g14] + genarray[g15] + genarray[g16]);
						for l := 1 to slength do 
						    temp[k, l] := temp[k, l] + temp2[l] * val
					    end
				    end else begin
					g1 := genenumber[f1, s1];
					g3 := genenumber[f2, s1];
					g5 := genenumber[f1, ms1];
					g7 := genenumber[f2, ms1];
					g9 := genenumber[mf1, s1];
					g11 := genenumber[mf2, s1];
					g13 := genenumber[mf1, ms1];
					g15 := genenumber[mf2, ms1];
					for k := 1 to nchild do 
					    with thischild[k]^ do begin
						val := 2 * ((1 - ps) * (1 - pf) * (genarray[g1] + genarray[g3]) + ps * (1 - pf) * (genarray[g5] + genarray[g7]) + pf * (1 - ps) * (genarray[g9] + genarray[g11]) + pf * ps * (genarray[g13] + genarray[g15]));
						for l := 1 to slength do 
						    temp[k, l] := temp[k, l] + temp2[l] * val
					    end
				    end;
				    secondseg := secondseg + 1
				end
			    else 
				for j := sstart to send do begin
				    for l := 1 to slength do 
					temp2[l] := temp1 * segprob[j + (l - 1) * slength];
				    s1 := seghap1[secondseg];
				    s2 := seghap2[secondseg];
				    ms1 := muthap[s1];
				    ms2 := muthap[s2];
				    if s1 <> s2 then begin
					g1 := genenumber[f1, s1];
					g2 := genenumber[f1, s2];
					g5 := genenumber[f1, ms1];
					g6 := genenumber[f1, ms2];
					g9 := genenumber[mf1, s1];
					g10 := genenumber[mf1, s2];
					g13 := genenumber[mf1, ms1];
					g14 := genenumber[mf1, ms2];
					for k := 1 to nchild do 
					    with thischild[k]^ do begin
						val := 2 * ((1 - ps) * (1 - pf) * (genarray[g1] + genarray[g2]) + ps * (1 - pf) * (genarray[g5] + genarray[g6]) + pf * (1 - ps) * (genarray[g9] + genarray[g10]) + pf * ps * (genarray[g13] + genarray[g14]));
						for l := 1 to slength do 
						    temp[k, l] := temp[k, l] + temp2[l] * val
					    end
				    end else begin
					g1 := genenumber[f1, s1];
					g5 := genenumber[f1, ms1];
					g9 := genenumber[mf1, s1];
					g13 := genenumber[mf1, ms1];
					for k := 1 to nchild do 
					    with thischild[k]^ do begin
						val := 4 * ((1 - ps) * (1 - pf) * genarray[g1] + ps * (1 - pf) * genarray[g5] + pf * (1 - ps) * genarray[g9] + pf * ps * genarray[g13]);
						for l := 1 to slength do 
						    temp[k, l] := temp[k, l] + temp2[l] * val
					    end
				    end;
				    secondseg := secondseg + 1
				end;
			firstseg := firstseg + 1
		    end;
	    temp1 := 0.0;
	    for l := 1 to slength do begin
		val := 1.0;
		for k := 1 to nchild do 
		    val := val * temp[k, l];
		temp1 := temp1 + val
	    end;
	    msegfast := temp1
	end; { msegfast }
 {msegfast}

	function msegfun: double;

	var
	    g1, g2, g3, g4, g5, g6, g7, g8, g9, g10, g11, g12, g13, g14, g15, g16, i, j, k, f1, f2, s1, s2, ms1, ms2, mf1, mf2: integer;
	    val, temp1, temp2: double;
	    temp: array [1..maxchild] of double; {msegfun}

	begin
	    firstseg := fseg;
	    for k := 1 to nchild do 
		temp[k] := 0.0;
	    with gennustruct^ do 
		with firstsex^ do 
		    for i := fstart to fend do begin
			temp1 := segprob[i];
			f1 := seghap1[firstseg];
			f2 := seghap2[firstseg];
			mf1 := muthap[f1];
			mf2 := muthap[f2];
			secondseg := sseg;
			with secondsex^ do 
			    if f1 <> f2 then 
				for j := sstart to send do begin
				    temp2 := temp1 * segprob[j];
				    s1 := seghap1[secondseg];
				    s2 := seghap2[secondseg];
				    ms1 := muthap[s1];
				    ms2 := muthap[s2];
				    if s1 <> s2 then begin
					g1 := genenumber[f1, s1];
					g2 := genenumber[f1, s2];
					g3 := genenumber[f2, s1];
					g4 := genenumber[f2, s2];
					g5 := genenumber[f1, ms1];
					g6 := genenumber[f1, ms2];
					g7 := genenumber[f2, ms1];
					g8 := genenumber[f2, ms2];
					g9 := genenumber[mf1, s1];
					g10 := genenumber[mf1, s2];
					g11 := genenumber[mf2, s1];
					g12 := genenumber[mf2, s2];
					g13 := genenumber[mf1, ms1];
					g14 := genenumber[mf1, ms2];
					g15 := genenumber[mf2, ms1];
					g16 := genenumber[mf2, ms2];
					for k := 1 to nchild do 
					    with thischild[k]^ do 
						temp[k] := temp[k] + temp2 * ((1 - ps) * (1 - pf) * (genarray[g1] + genarray[g2] + genarray[g3] + genarray[g4]) + ps * (1 - pf) * (genarray[g5] + genarray[g6] + genarray[g7] + genarray[g8]) + pf * (1 - ps) * (genarray[g9] + genarray[g10] + genarray[g11] + genarray[g12]) + pf * ps * (genarray[g13] + genarray[g14] + genarray[g15] + genarray[g16]))
				    end else begin
					g1 := genenumber[f1, s1];
					g3 := genenumber[f2, s1];
					g5 := genenumber[f1, ms1];
					g7 := genenumber[f2, ms1];
					g9 := genenumber[mf1, s1];
					g11 := genenumber[mf2, s1];
					g13 := genenumber[mf1, ms1];
					g15 := genenumber[mf2, ms1];
					for k := 1 to nchild do 
					    with thischild[k]^ do 
						temp[k] := temp[k] + 2 * temp2 * ((1 - ps) * (1 - pf) * (genarray[g1] + genarray[g3]) + ps * (1 - pf) * (genarray[g5] + genarray[g7]) + pf * (1 - ps) * (genarray[g9] + genarray[g11]) + pf * ps * (genarray[g13] + genarray[g15]))
				    end;
				    secondseg := secondseg + 1
				end
			    else 
				for j := sstart to send do begin
				    temp2 := temp1 * segprob[j];
				    s1 := seghap1[secondseg];
				    s2 := seghap2[secondseg];
				    ms1 := muthap[s1];
				    ms2 := muthap[s2];
				    if s1 <> s2 then begin
					g1 := genenumber[f1, s1];
					g2 := genenumber[f1, s2];
					g5 := genenumber[f1, ms1];
					g6 := genenumber[f1, ms2];
					g9 := genenumber[mf1, s1];
					g10 := genenumber[mf1, s2];
					g13 := genenumber[mf1, ms1];
					g14 := genenumber[mf1, ms2];
					for k := 1 to nchild do 
					    with thischild[k]^ do 
						temp[k] := temp[k] + 2 * temp2 * ((1 - ps) * (1 - pf) * (genarray[g1] + genarray[g2]) + ps * (1 - pf) * (genarray[g5] + genarray[g6]) + pf * (1 - ps) * (genarray[g9] + genarray[g10]) + pf * ps * (genarray[g13] + genarray[g14]))
				    end else begin
					g1 := genenumber[f1, s1];
					g5 := genenumber[f1, ms1];
					g9 := genenumber[mf1, s1];
					g13 := genenumber[mf1, ms1];
					for k := 1 to nchild do 
					    with thischild[k]^ do 
						temp[k] := temp[k] + 4 * temp2 * ((1 - ps) * (1 - pf) * genarray[g1] + ps * (1 - pf) * genarray[g5] + pf * (1 - ps) * genarray[g9] + pf * ps * genarray[g13])
				    end;
				    secondseg := secondseg + 1
				end;
			firstseg := firstseg + 1
		    end;
	    val := 1.0;
	    for k := 1 to nchild do 
		val := val * temp[k];
	    msegfun := val
	end; { msegfun }
 {msegfun}

	function segfast: double;

	var
	    g1, g2, g3, g4, i, j, k, l, f1, f2, s1, s2, slength: integer;
	    val, temp1: double;
	    temp: array [1..maxchild, 1..maxseg] of double;
	    temp2: array [1..maxseg] of double; {segfast}

	begin
	    firstseg := fseg;
	    slength := send - sstart + 1;
	    for k := 1 to nchild do 
		for l := 1 to slength do 
		    temp[k, l] := 0.0;
	    with gennustruct^ do 
		with firstsex^ do 
		    for i := fstart to fend do begin
			temp1 := segprob[i];
			f1 := seghap1[firstseg];
			f2 := seghap2[firstseg];
			secondseg := sseg;
			with secondsex^ do 
			    if f1 <> f2 then 
				for j := sstart to send do begin
				    for l := 1 to slength do 
					temp2[l] := temp1 * segprob[j + (l - 1) * slength];
				    s1 := seghap1[secondseg];
				    s2 := seghap2[secondseg];
				    if s1 <> s2 then begin
					g1 := genenumber[f1, s1];
					g2 := genenumber[f1, s2];
					g3 := genenumber[f2, s1];
					g4 := genenumber[f2, s2];
					for k := 1 to nchild do 
					    with thischild[k]^ do begin
						val := genarray[g1] + genarray[g2] + genarray[g3] + genarray[g4];
						for l := 1 to slength do 
						    temp[k, l] := temp[k, l] + temp2[l] * val
					    end
				    end else begin
					g1 := genenumber[f1, s1];
					g3 := genenumber[f2, s1];
					for k := 1 to nchild do 
					    with thischild[k]^ do begin
						val := 2.0 * (genarray[g1] + genarray[g3]);
						for l := 1 to slength do 
						    temp[k, l] := temp[k, l] + temp2[l] * val
					    end
				    end;
				    secondseg := secondseg + 1
				end
			    else 
				for j := sstart to send do begin
				    for l := 1 to slength do 
					temp2[l] := temp1 * segprob[j + (l - 1) * slength];
				    s1 := seghap1[secondseg];
				    s2 := seghap2[secondseg];
				    if s1 <> s2 then begin
					g1 := genenumber[f1, s1];
					g2 := genenumber[f1, s2];
					for k := 1 to nchild do 
					    with thischild[k]^ do begin
						val := 2.0 * (genarray[g1] + genarray[g2]);
						for l := 1 to slength do 
						    temp[k, l] := temp[k, l] + temp2[l] * val
					    end
				    end else begin
					g1 := genenumber[f1, s1];
					for k := 1 to nchild do 
					    with thischild[k]^ do begin
						val := 4.0 * genarray[g1];
						for l := 1 to slength do 
						    temp[k, l] := temp[k, l] + temp2[l] * val
					    end
				    end;
				    secondseg := secondseg + 1
				end;
			firstseg := firstseg + 1
		    end;
	    temp1 := 0.0;
	    for l := 1 to slength do begin
		val := 1.0;
		for k := 1 to nchild do 
		    val := val * temp[k, l];
		temp1 := temp1 + val
	    end;
	    segfast := temp1
	end; { segfast }
 {segfast}

	procedure initseg; {initseg}

	begin
	    if p^.male then begin
		nfirst := mgeno;
		nsecond := fgeno;
		firstsex := maletheta;
		secondsex := femaletheta;
		pf := mutmale;
		ps := mutfemale
	    end else begin
		nfirst := fgeno;
		nsecond := mgeno;
		firstsex := femaletheta;
		secondsex := maletheta;
		pf := mutfemale;
		ps := mutmale
	    end;
	    prob(father);
	    prob(mother);
	    child := father^.foff;
	    nchild := 0;
	    repeat
		prob(child);
		if (child^.ma = mother) and not child^.up then begin
		    nchild := nchild + 1;
		    thischild[nchild] := child^.gen;
		    malechild[nchild] := child^.male
		end;
		child := child^.nextpa
	    until child = nil
	end; { initseg }
 {initseg}

	procedure exitseg; {exitseg}

	begin
	    nuscale := nuscale + 1;
	    child := father^.foff;
	    repeat
		if (child^.ma = mother) and not child^.up then 
		    cleanup(child);
		child := child^.nextpa
	    until child = nil
	end; { exitseg }
 {exitseg}

	procedure segsexctop;

	var
	    segval: double;
	    first, second, sstop: integer;
	    thiscensor, thisrare: boolean; {segsexctop}

	begin
	    initseg;
	    with censorstruct^ do 
		if p^.male then begin
		    with p^.gen^ do 
			for first := 1 to nfirst do 
			    if genarray[first] <> 0.0 then begin
				segval := 0.0;
				fseg := first;
				thisrare := rare[first];
				second := 1;
				with q^.gen^ do 
				    repeat
					sstart := probstart[second];
					send := probend[second];
					sseg := segstart[second];
					sstop := second + (send - sstart) + 1;
					if thisc < maxcensor then begin
					    thisc := thisc + 1;
					    thiscensor := censor[thisc]
					end else 
					    thiscensor := (genarray[second] = 0.0) or thisrare and rare[second];
					if not thiscensor then 
					    if mutsys <> 0 then 
						segval := segval + genarray[second] * msegsexf
					    else 
						segval := segval + genarray[second] * segsexf;
					second := sstop
				    until second > nsecond;
				genarray[first] := genarray[first] * segval * segscale
			    end
		end else 
		    with p^.gen^ do 
			for first := 1 to nfirst do 
			    if genarray[first] <> 0.0 then begin
				segval := 0.0;
				fstart := probstart[first];
				fend := probend[first];
				fseg := segstart[first];
				thisrare := rare[first];
				second := 1;
				with q^.gen^ do 
				    repeat
					sseg := second;
					if thisc < maxcensor then begin
					    thisc := thisc + 1;
					    thiscensor := censor[thisc]
					end else 
					    thiscensor := (genarray[second] = 0.0) or thisrare and rare[second];
					if not thiscensor then 
					    if mutsys <> 0 then 
						segval := segval + genarray[second] * msegsex
					    else 
						segval := segval + genarray[second] * segsex;
					second := second + 1
				    until second > nsecond;
				genarray[first] := genarray[first] * segval * segscale
			    end;
	    cleanup(q);
	    exitseg
	end; { segsexctop }
 {segsexctop}

	procedure segsextop;

	var
	    segval, val: double;
	    first, second, sstop: integer; {segsextop}

	begin
	    initseg;
	    with censorstruct^ do 
		if p^.male then begin
		    with p^.gen^ do 
			for first := 1 to nfirst do 
			    if genarray[first] <> 0.0 then begin
				segval := 0.0;
				fseg := first;
				second := 1;
				with q^.gen^ do 
				    repeat
					sstart := probstart[second];
					send := probend[second];
					sseg := segstart[second];
					sstop := second + (send - sstart) + 1;
					if thisc <= maxcensor then 
					    thisc := thisc + 1;
					if genarray[second] = 0.0 then begin
					    second := sstop;
					    if thisc <= maxcensor then 
						censor[thisc] := true
					end else begin
					    if mutsys <> 0 then 
						val := msegsexf
					    else 
						val := segsexf;
					    if thisc <= maxcensor then 
						censor[thisc] := val = 0.0;
					    segval := segval + genarray[second] * val;
					    second := sstop
					end
				    until second > nsecond;
				genarray[first] := genarray[first] * segval * segscale
			    end
		end else 
		    with p^.gen^ do 
			for first := 1 to nfirst do 
			    if genarray[first] <> 0.0 then begin
				segval := 0.0;
				fstart := probstart[first];
				fend := probend[first];
				fseg := segstart[first];
				second := 1;
				with q^.gen^ do 
				    repeat
					sseg := second;
					if thisc <= maxcensor then 
					    thisc := thisc + 1;
					if genarray[second] = 0.0 then begin
					    second := second + 1;
					    if thisc <= maxcensor then 
						censor[thisc] := true
					end else begin
					    if mutsys <> 0 then 
						val := msegsex
					    else 
						val := segsex;
					    if thisc <= maxcensor then 
						censor[thisc] := val = 0.0;
					    segval := segval + genarray[second] * val;
					    second := second + 1
					end
				    until second > nsecond;
				genarray[first] := genarray[first] * segval * segscale
			    end;
	    cleanup(q);
	    exitseg
	end; { segsextop }
 {segsextop}

	procedure segsexup;

	var
	    segval: double;
	    first, second: integer; {segsexup}

	begin
	    initseg;
	    with censorstruct^ do 
		if p^.male then begin
		    with p^.gen^ do 
			for first := 1 to nfirst do 
			    if genarray[first] <> 0.0 then begin
				segval := 0.0;
				fseg := first;
				second := 1;
				with q^.gen^ do 
				    for second := 1 to nsecond do 
					if genarray[second] <> 0.0 then begin
					    sstart := probstart[second];
					    send := probend[second];
					    sseg := segstart[second];
					    if mutsys <> 0 then 
						segval := segval + genarray[second] * msegsex
					    else 
						segval := segval + genarray[second] * segsex
					end;
				genarray[first] := genarray[first] * segval * segscale
			    end
		end else 
		    with p^.gen^ do 
			for first := 1 to nfirst do 
			    if genarray[first] <> 0.0 then begin
				segval := 0.0;
				fstart := probstart[first];
				fend := probend[first];
				fseg := segstart[first];
				with q^.gen^ do 
				    for second := 1 to nsecond do 
					if genarray[second] <> 0.0 then begin
					    sseg := second;
					    if mutsys <> 0 then 
						segval := segval + genarray[second] * msegsex
					    else 
						segval := segval + genarray[second] * segsex
					end;
				genarray[first] := genarray[first] * segval * segscale
			    end;
	    cleanup(q);
	    exitseg
	end; { segsexup }
 {segsexup}

	procedure segsexdown;

	var
	    here: 1..maxfem;
	    gene: genotype;
	    val, temp2: double;
	    j, first, second: integer; {segsexdown}

	begin
	    initseg;
	    for here := 1 to fgeno do 
		gene[here] := 0.0;
	    with gennustruct^ do 
		with censorstruct^ do 
		    with p^.gen^ do 
			for first := 1 to nfirst do 
			    if genarray[first] <> 0.0 then begin
				fseg := first;
				second := 1;
				with q^.gen^ do 
				    for second := 1 to nsecond do 
					if genarray[second] <> 0.0 then begin
					    val := genarray[second] * p^.gen^.genarray[first] * segscale;
					    sstart := probstart[second];
					    send := probend[second];
					    sseg := segstart[second];
					    if nchild <> 0 then 
						val := val * segsex;
					    if val <> 0.0 then begin
						secondseg := sseg;
						with femaletheta^ do 
						    for j := sstart to send do begin
							temp2 := segprob[j];
							if r^.male then begin
							    here := seghap1[secondseg];
							    gene[here] := gene[here] + temp2 * val;
							    here := seghap2[secondseg];
							    gene[here] := gene[here] + temp2 * val
							end else begin
							    here := genenumber[seghap1[secondseg], first];
							    gene[here] := gene[here] + temp2 * val;
							    here := genenumber[seghap2[secondseg], first];
							    gene[here] := gene[here] + temp2 * val
							end;
							secondseg := secondseg + 1
						    end
					    end {second}
					end
			    end; {first}
	    p^.gen^.genarray := r^.gen^.genarray;
	    r^.gen^.genarray := gene;
	    with r^.gen^ do 
		for first := 1 to fgeno do 
		    genarray[first] := genarray[first] * p^.gen^.genarray[first];
	    cleanup(p);
	    cleanup(q);
	    exitseg
	end; { segsexdown }
 {segsexdown}

	procedure segctop;

	var
	    segval: double;
	    first, second, sstop: integer;
	    thiscensor, thisrare: boolean; {segctop}

	begin
	    initseg;
	    with censorstruct^ do 
		with p^.gen^ do 
		    for first := 1 to fgeno do 
			if genarray[first] <> 0.0 then begin
			    segval := 0.0;
			    fstart := probstart[first];
			    fend := probend[first];
			    fseg := segstart[first];
			    thisrare := rare[first];
			    second := 1;
			    with q^.gen^ do 
				repeat
				    sstart := probstart[second];
				    send := probend[second];
				    sseg := segstart[second];
				    sstop := second + (send - sstart) + 1;
				    if thisc < maxcensor then begin
					thisc := thisc + 1;
					thiscensor := censor[thisc]
				    end else 
					thiscensor := (genarray[second] = 0.0) or thisrare and rare[second];
				    if thiscensor then 
					second := sstop
				    else begin
					if mutsys <> 0 then 
					    segval := segval + genarray[second] * msegfast
					else 
					    segval := segval + genarray[second] * segfast;
					second := sstop
				    end
				until second > fgeno;
			    genarray[first] := genarray[first] * segval * segscale
			end;
	    exitseg;
	    if approximate and firstapprox and (p^.pa = nil) then 
		getapprox
	end; { segctop }
 {segctop}

	procedure segtop;

	var
	    segval, val: double;
	    first, second, sstop: integer;
	    thatrare, thisrare: boolean; {segtop}

	begin
	    initseg;
	    with censorstruct^ do 
		with p^.gen^ do 
		    for first := 1 to fgeno do 
			if genarray[first] <> 0.0 then begin
			    thisrare := rare[first];
			    segval := 0.0;
			    fstart := probstart[first];
			    fend := probend[first];
			    fseg := segstart[first];
			    second := 1;
			    with q^.gen^ do 
				repeat
				    sstart := probstart[second];
				    send := probend[second];
				    sseg := segstart[second];
				    sstop := second + (send - sstart) + 1;
				    if thisc <= maxcensor then 
					thisc := thisc + 1;
				    thatrare := thisrare and rare[second];
				    if (genarray[second] = 0.0) or thatrare then begin
					second := sstop;
					if thisc <= maxcensor then 
					    censor[thisc] := true
				    end else begin
					if mutsys <> 0 then 
					    val := msegfast
					else 
					    val := segfast;
					if thisc <= maxcensor then 
					    censor[thisc] := val = 0.0;
					segval := segval + genarray[second] * val;
					second := sstop
				    end
				until second > fgeno;
			    genarray[first] := genarray[first] * segval * segscale
			end;
	    cleanup(q);
	    exitseg;
	    if approximate and firstapprox and (p^.pa = nil) then 
		getapprox
	end; { segtop }
 {segtop}

	procedure segcapprox;

	var
	    segval: double;
	    first, second, sstop: integer;
	    thiscensor: boolean; {segcapprox}

	begin
	    initseg;
	    with censorstruct^ do 
		with approxstruct^ do 
		    with p^.gen^ do 
			for first := 1 to fgeno do 
			    if genarray[first] <> 0.0 then begin
				segval := 0.0;
				fstart := probstart[first];
				fend := probend[first];
				fseg := segstart[first];
				second := 1;
				with q^.gen^ do 
				    repeat
					sstart := probstart[second];
					send := probend[second];
					sseg := segstart[second];
					sstop := second + (send - sstart) + 1;
					if thisc < maxcensor then begin
					    thisc := thisc + 1;
					    thiscensor := censor[thisc]
					end else 
					    thiscensor := genarray[second] = 0.0;
					if thiscensor or not approxarray[thisped, first] then 
					    second := sstop
					else begin
					    if mutsys <> 0 then 
						segval := segval + genarray[second] * msegfast
					    else 
						segval := segval + genarray[second] * segfast;
					    second := sstop
					end
				    until second > fgeno;
				genarray[first] := genarray[first] * segval * segscale
			    end;
	    cleanup(q);
	    exitseg
	end; { segcapprox }
 {segcapprox}

	procedure segup;

	var
	    segval: double;
	    first, second: integer; {segup}

	begin
	    initseg;
	    with censorstruct^ do 
		with p^.gen^ do 
		    for first := 1 to fgeno do 
			if genarray[first] <> 0.0 then begin
			    segval := 0.0;
			    fstart := probstart[first];
			    fend := probend[first];
			    fseg := segstart[first];
			    with q^.gen^ do 
				for second := 1 to fgeno do 
				    if genarray[second] <> 0.0 then begin
					sstart := probstart[second];
					send := probend[second];
					sseg := segstart[second];
					if mutsys <> 0 then 
					    segval := segval + genarray[second] * msegfun
					else 
					    segval := segval + genarray[second] * segfun
				    end;
			    genarray[first] := genarray[first] * segval * segscale
			end;
	    cleanup(q);
	    exitseg
	end; { segup }
 {segup}

	procedure segdown;

	var
	    here: 1..maxfem;
	    gene: genotype;
	    val, temp1, temp2: double;
	    f1, f2, s1, s2, i, j, first, second: integer;

	begin
	    initseg;
	    for here := 1 to fgeno do 
		gene[here] := 0.0;
	    with gennustruct^ do 
		with censorstruct^ do 
		    with p^.gen^ do 
			for first := 1 to mgeno do 
			    if genarray[first] <> 0.0 then begin
				fstart := probstart[first];
				fend := probend[first];
				fseg := segstart[first];
				with q^.gen^ do 
				    for second := 1 to fgeno do 
					if genarray[second] <> 0.0 then begin
					    sstart := probstart[second];
					    send := probend[second];
					    sseg := segstart[second];
					    val := segscale * genarray[second] * p^.gen^.genarray[first];
					    if nchild <> 0 then 
						val := segfun * val;
					    if val <> 0.0 then begin
						firstseg := fseg;
						with maletheta^ do 
						    for i := fstart to fend do begin
							temp1 := segprob[i];
							secondseg := sseg;
							with femaletheta^ do 
							    for j := sstart to send do begin
								temp2 := segprob[j] * temp1 * val;
								f1 := seghap1[firstseg];
								f2 := seghap2[firstseg];
								s1 := seghap1[secondseg];
								s2 := seghap2[secondseg];
								here := genenumber[s1, f1];
								gene[here] := gene[here] + temp2;
								here := genenumber[s1, f2];
								gene[here] := gene[here] + temp2;
								here := genenumber[s2, f1];
								gene[here] := gene[here] + temp2;
								here := genenumber[s2, f2];
								gene[here] := gene[here] + temp2;
								secondseg := secondseg + 1
							    end;
							firstseg := firstseg + 1
						    end
					    end {second}
					end
			    end; {first}
	    p^.gen^.genarray := r^.gen^.genarray;
	    r^.gen^.genarray := gene;
	    with r^.gen^ do 
		for first := 1 to fgeno do 
		    genarray[first] := genarray[first] * p^.gen^.genarray[first];
	    cleanup(p);
	    cleanup(q);
	    exitseg
	end; { segdown }
 {segdown}

	procedure msegsexdown;

	var
	    here: 1..maxfem;
	    gene: genotype;
	    val, temp2: double;
	    ms2, ms1, mf, j, first, second: integer; {msegsexdown}

	begin
	    initseg;
	    for here := 1 to fgeno do 
		gene[here] := 0.0;
	    with gennustruct^ do 
		with censorstruct^ do 
		    with p^.gen^ do 
			for first := 1 to nfirst do 
			    if genarray[first] <> 0.0 then begin
				fseg := first;
				second := 1;
				with q^.gen^ do 
				    for second := 1 to nsecond do 
					if genarray[second] <> 0.0 then begin
					    val := genarray[second] * p^.gen^.genarray[first] * segscale;
					    sstart := probstart[second];
					    send := probend[second];
					    sseg := segstart[second];
					    if nchild <> 0 then 
						val := val * msegsex;
					    if val <> 0.0 then begin
						mf := muthap[seghap1[fseg]];
						secondseg := sseg;
						with femaletheta^ do 
						    for j := sstart to send do begin
							ms1 := muthap[seghap1[secondseg]];
							ms2 := muthap[seghap2[secondseg]];
							temp2 := segprob[j];
							if r^.male then begin
							    here := seghap1[secondseg];
							    gene[here] := gene[here] + (1 - ps) * temp2 * val;
							    here := seghap2[secondseg];
							    gene[here] := gene[here] + (1 - ps) * temp2 * val;
							    here := ms1;
							    gene[here] := gene[here] + ps * temp2 * val;
							    here := ms2;
							    gene[here] := gene[here] + ps * temp2 * val
							end else begin
							    here := genenumber[seghap1[secondseg], first];
							    gene[here] := gene[here] + (1 - pf) * (1 - ps) * temp2 * val;
							    here := genenumber[seghap2[secondseg], first];
							    gene[here] := gene[here] + (1 - pf) * (1 - ps) * temp2 * val;
							    here := genenumber[seghap1[secondseg], mf];
							    gene[here] := gene[here] + pf * (1 - ps) * temp2 * val;
							    here := genenumber[seghap2[secondseg], mf];
							    gene[here] := gene[here] + pf * (1 - ps) * temp2 * val;
							    here := genenumber[ms1, first];
							    gene[here] := gene[here] + (1 - pf) * ps * temp2 * val;
							    here := genenumber[ms2, first];
							    gene[here] := gene[here] + (1 - pf) * ps * temp2 * val;
							    here := genenumber[ms1, mf];
							    gene[here] := gene[here] + pf * ps * temp2 * val;
							    here := genenumber[ms2, mf];
							    gene[here] := gene[here] + pf * ps * temp2 * val
							end;
							secondseg := secondseg + 1
						    end
					    end {second}
					end
			    end; {first}
	    p^.gen^.genarray := r^.gen^.genarray;
	    r^.gen^.genarray := gene;
	    with r^.gen^ do 
		for first := 1 to fgeno do 
		    genarray[first] := genarray[first] * p^.gen^.genarray[first];
	    cleanup(p);
	    cleanup(q);
	    exitseg
	end; { msegsexdown }
 {msegsexdown}

	procedure msegdown;

	var
	    here: 1..maxfem;
	    gene: genotype;
	    val, temp, temp1, temp2: double;
	    i, j, first, second, f1, f2, s1, s2, mf1, mf2, ms1, ms2: integer; {msegdown}

	begin
	    initseg;
	    for here := 1 to fgeno do 
		gene[here] := 0.0;
	    with gennustruct^ do 
		with censorstruct^ do 
		    with p^.gen^ do 
			for first := 1 to fgeno do 
			    if genarray[first] <> 0.0 then begin
				fstart := probstart[first];
				fend := probend[first];
				fseg := segstart[first];
				with q^.gen^ do 
				    for second := 1 to fgeno do 
					if genarray[second] <> 0.0 then begin
					    sstart := probstart[second];
					    send := probend[second];
					    sseg := segstart[second];
					    val := genarray[second] * p^.gen^.genarray[first] * segscale;
					    if nchild <> 0 then 
						val := msegfun * val;
					    if val <> 0.0 then begin
						firstseg := fseg;
						with maletheta^ do 
						    for i := fstart to fend do begin
							temp1 := segprob[i];
							secondseg := sseg;
							with femaletheta^ do 
							    for j := sstart to send do begin
								temp2 := segprob[j];
								f1 := seghap1[firstseg];
								f2 := seghap2[firstseg];
								s1 := seghap1[secondseg];
								s2 := seghap2[secondseg];
								temp := (1 - pf) * (1 - ps) * temp1 * temp2 * val;
								here := genenumber[s1, f1];
								gene[here] := gene[here] + temp;
								here := genenumber[s1, f2];
								gene[here] := gene[here] + temp;
								here := genenumber[s2, f1];
								gene[here] := gene[here] + temp;
								here := genenumber[s2, f2];
								gene[here] := gene[here] + temp;
								ms1 := muthap[s1];
								ms2 := muthap[s2];
								temp := (1 - pf) * ps * temp1 * temp2 * val;
								here := genenumber[ms1, f1];
								gene[here] := gene[here] + temp;
								here := genenumber[ms1, f2];
								gene[here] := gene[here] + temp;
								here := genenumber[ms2, f1];
								gene[here] := gene[here] + temp;
								here := genenumber[ms2, f2];
								gene[here] := gene[here] + temp;
								mf1 := muthap[f1];
								mf2 := muthap[f2];
								temp := pf * (1 - ps) * temp1 * temp2 * val;
								here := genenumber[mf1, s1];
								gene[here] := gene[here] + temp;
								here := genenumber[mf1, s2];
								gene[here] := gene[here] + temp;
								here := genenumber[mf2, s1];
								gene[here] := gene[here] + temp;
								here := genenumber[mf2, s2];
								gene[here] := gene[here] + temp;
								temp := pf * ps * temp1 * temp2 * val;
								here := genenumber[mf1, ms1];
								gene[here] := gene[here] + temp;
								here := genenumber[mf1, ms2];
								gene[here] := gene[here] + temp;
								here := genenumber[mf2, ms1];
								gene[here] := gene[here] + temp;
								here := genenumber[mf2, ms2];
								gene[here] := gene[here] + temp;
								secondseg := secondseg + 1
							    end;
							firstseg := firstseg + 1
						    end
					    end {second}
					end
			    end; {first}
	    p^.gen^.genarray := r^.gen^.genarray;
	    r^.gen^.genarray := gene;
	    with r^.gen^ do 
		for first := 1 to fgeno do 
		    genarray[first] := genarray[first] * p^.gen^.genarray[first];
	    cleanup(p);
	    cleanup(q);
	    exitseg
	end; { msegdown }
 {msegdown}

    begin {seg}
	phaseunkn := (q^.pa = nil) and q^.firstpass and (q^.inloop = 0) and not disequi;
	if p^.male then begin
	    father := p;
	    mother := q
	end else begin
	    father := q;
	    mother := p
	end;
	if peel = peelup then 
	    if (p^.pa = nil) and phaseunkn and approximate then 
		if firstapprox and firsttime then 
		    case thispath of
			auto:
			    segtop;
			mauto:
			    segtop;
			sex:
			    segsextop;
			msex:
			    segsextop
		    end
		else if firstapprox then  {first approximate not first time}
		    case thispath of
			auto:
			    segctop;
			mauto:
			    segctop;
			sex:
			    segsexctop;
			msex:
			    segsexctop
		    end
		else  {approximate}
		    case thispath of
			auto:
			    segcapprox;
			mauto:
			    segcapprox;
			sex:
			    segsexctop;
			msex:
			    segsexctop
		    end
	    else if phaseunkn then  {do not approximate}
		if firsttime then 
		    case thispath of
			auto:
			    segtop;
			mauto:
			    segtop;
			sex:
			    segsextop;
			msex:
			    segsextop
		    end
		else  {not firsttime}
		    case thispath of
			auto:
			    segctop;
			mauto:
			    segctop;
			sex:
			    segsexctop;
			msex:
			    segsexctop
		    end
	    else  {phaseinfo}
		case thispath of
		    auto:
			segup;
		    mauto:
			segup;
		    sex:
			segsexup;
		    msex:
			segsexup
		end
	else  {not peelup}
	    case thispath of
		auto:
		    segdown;
		mauto:
		    msegdown;
		sex:
		    segsexdown;
		msex:
		    msegsexdown
	    end;
	q^.firstpass := false;
	p^.firstpass := false
    end; { seg }
						{seg}

    procedure collapsedown(p: ind);
    forward;

    procedure collapseup(p: ind);

    var
	q, child, nextchild: ind;
	down: boolean;

    begin {collapseup}
	p^.done := true;
	if p^.foff <> nil then begin
	    down := false;
	    child := p^.foff;
	    while child <> nil do begin
		down := false;
		if p^.male then 
		    q := child^.ma
		else 
		    q := child^.pa;
		if not q^.done then begin
		    collapsedown(q);
		    nextchild := child;
		    while nextchild <> nil do begin
			if (nextchild^.pa = q) or (nextchild^.ma = q) then 
			    if not nextchild^.up then 
				collapseup(nextchild)
			    else 
				down := true;
			if p^.male then 
			    nextchild := nextchild^.nextpa
			else 
			    nextchild := nextchild^.nextma
		    end;
		    if q^.multi then 
			collapseup(q);
		    if not down then 
			seg(p, q, child, peelup)
		    else 
			collapsedown(p)
		end;
		if p^.male then 
		    child := child^.nextpa
		else 
		    child := child^.nextma
	    end
	end
    end; { collapseup }
						{collapseup}

    procedure collapsedown;

    begin {collapsedown}
	if p^.pa <> nil then begin
	    p^.up := true;
	    collapseup(p^.pa);
	    seg(p^.pa, p^.ma, p, peeldown)
	end
    end; { collapsedown }
						{collapsedown}

    procedure riskcumul;

    var
	i: integer;

    begin {riskcumul}
	with proband^.gen^ do 
	    if sexlink and proband^.male then begin
		for i := 1 to mgeno do 
		    if riskmale[i] then 
			hetero := hetero + genarray[i]
	    end else begin
		for i := 1 to fgeno do 
		    if risk2[i] then 
			homo := homo + genarray[i]
		    else if risk1[i] then 
			hetero := hetero + genarray[i]
	    end
    end; { riskcumul }
						{riskcumul}

    procedure riskcalc;

    var
	normal: double;

    begin {riskcalc}
	homo := homo / like;
	hetero := hetero / like;
	normal := 1 - homo - hetero;
	writeln(outfile, 'RISK FOR PERSON ', proband^.id: 6, ' IN PEDIGREE ', proband^.ped: 7);
	if not proband^.male or not sexlink then 
	    writeln(outfile, 'HOMOZYGOTE CARRIER   : ', homo: 8: 5);
	if not proband^.male or not sexlink then 
	    writeln(outfile, 'HETEROZYGOTE CARRIER : ', hetero: 8: 5)
	else 
	    writeln(outfile, 'MALE CARRIER         : ', hetero: 8: 5);
	writeln(outfile, 'NORMAL               : ', normal: 8: 5);
	writeln(output, 'RISK FOR PERSON ', proband^.id: 6, ' IN PEDIGREE ', proband^.ped: 7);
	if not proband^.male or not sexlink then 
	    writeln(output, 'HOMOZYGOTE CARRIER   : ', homo: 8: 5);
	if not proband^.male or not sexlink then 
	    writeln(output, 'HETEROZYGOTE CARRIER : ', hetero: 8: 5)
	else 
	    writeln(output, 'MALE CARRIER         : ', hetero: 8: 5);
	writeln(output, 'NORMAL               : ', normal: 8: 5)
    end; { riskcalc }
				{riskcalc}

begin				{likelihood}
    if informative[thisped] then begin
	homo := 0.0;
	hetero := 0.0;
	tmplike := 0.0;
	alldone := false;
	nuscale := 0;
	for i := 1 to maxloop do begin
	    loopgen[i] := 1;
	    loopmax[i] := 1;
	    holdpoint[i] := nil;
	    if looppers[thisped, i, 1] <> nil then 
		with looppers[thisped, i, 1]^ do begin
		    new(gen);
		    {          gen:=genpoint(NewPtr(SizeOf(thisarray)));}
		    with gen^ do 
			for j := 1 to fgeno do 
			    genarray[j] := 0.0;
		    getvect(looppers[thisped, i, 1]);
		    if looppers[thisped, i, 1]^.pa = nil then 
			nuscale := nuscale + 1;
		    holdpoint[i] := gen;
		    gen := nil;
		    if male then 
			loopmax[i] := mgeno
		    else 
			loopmax[i] := fgeno
		end
	end;
	loopgen[1] := 0;
	repeat
	    i := 1;
	    repeat
		loopgen[i] := loopgen[i] + 1;
		if loopgen[i] > loopmax[i] then 
		    loopgen[i] := 1
		else 
		    i := maxloop;
		i := i + 1
	    until i > maxloop;
	    gocalc := true;
	    for i := 1 to maxloop do 
		{ML change}
		if holdpoint[i] <> nil then 
		    if holdpoint[i]^.genarray[loopgen[i]] = 0.0 then 
			gocalc := false;
	    if gocalc then begin
		for i := 1 to totperson do 
		    with person[i]^ do begin
			gen := nil;
			done := false;
			up := false
		    end;
		collapseup(proband);
		collapsedown(proband);
		like := 0.0;
		with proband^.gen^ do 
		    if proband^.male then 
			for i := 1 to mgeno do 
			    like := like + genarray[i]
		    else 
			for i := 1 to fgeno do 
			    like := like + genarray[i];
		tmplike := tmplike + like;
		if risk and (like <> 0.0) then 
		    riskcumul;
		for i := 1 to totperson do 
		    if person[i]^.gen <> nil then 
			cleanup(person[i])
	    end;
	    alldone := true;
	    for i := 1 to maxloop do 
		alldone := alldone and (loopgen[i] = loopmax[i])
	until alldone;
	like := tmplike;
	if risk and (like <> 0.0) then 
	    riskcalc;
	if like = 0.0 then 
	    like := zerolike
	else 
	    like := ln(like) - nuscale * ln(segscale);
	for i := 1 to maxloop do 
	    if holdpoint[i] <> nil then begin
		dispose(holdpoint[i]);
		{        DisposPtr(Ptr(holdpoint[i]));}
		holdpoint[i] := nil
	    end
    end else 
	like := 0.0
end; { likelihood }				{likelihood}


procedure fun(var f: double; var x: vector);

var
    i, k: integer;

begin
    k := 0;
    for i := 1 to nall do 
	if itp[i] = 1 then begin
	    k := k + 1;
	    xall[i] := x[k]
	end;

    setparam;
    if not penalty then begin
	for i := 1 to totperson do 
	    person[i]^.done := false;
	for i := 1 to totperson do 
	    person[i]^.firstpass := true;
	recombination;
	alike := 0.0;
	thisc := minint;
	for i := 1 to nuped do begin
	    likelihood(i, proband[i]);
	    alike := alike + like
	end;
	if prn then 
	    if firsttime then 
		if thisc < maxcensor then 
		    writeln(output, 'Maxcensor can be reduced to ', thisc)
		else if thisc > maxcensor then 
		    writeln(output, 'You may gain efficiency by increasing maxcensor');
	firsttime := false;
	f := -2 * alike;
	penlike := f;
	firsttime := false
    end else 
	f := 1.5 * penlike
end; { fun }


procedure out;

var
    i, j, k, l: integer;
    deriv: array [0..2] of double;
    oldtheta: thetarray;
    lods, dist: double;
    thisval: array [1..maxped] of double;

    procedure getlods(var f: double; var x: vector);

    var
	i, k: integer;

    begin
	firstapprox := true;
	k := 0;
	for i := 1 to nall do 
	    if itp[i] = 1 then begin
		k := k + 1;
		xall[i] := x[k]
	    end;
	for i := 1 to nuped do 
	    thisval[i] := 0.0;
	setparam;
	if not penalty then begin
	    if byfamily then begin
		for i := 1 to 35 do 
		    if finalprn then 
			write(final, '-');
		if finalprn then 
		    writeln(final);
		if finalprn then 
		    writeln(final, 'PEDIGREE |  LN LIKE  | LOG 10 LIKE');
		for i := 1 to 35 do 
		    if finalprn then 
			write(final, '-');
		if finalprn then 
		    writeln(final)
	    end;
	    for i := 1 to totperson do 
		person[i]^.done := false;
	    for i := 1 to totperson do 
		person[i]^.firstpass := true;
	    recombination;
	    thisc := minint;
	    for i := 1 to nuped do begin
		likelihood(i, proband[i]);
		thisval[i] := like
	    end
	end;
	firsttime := true;
	for i := 1 to nlocus do 
	    maletheta^.theta[i] := 0.5;
	for i := 1 to nlocus do 
	    femaletheta^.theta[i] := 0.5;
	if not penalty then begin
	    dolod := true;
	    for i := 1 to totperson do 
		person[i]^.done := false;
	    for i := 1 to totperson do 
		person[i]^.firstpass := true;
	    recombination;
	    alike := 0.0;
	    thisc := minint;
	    for i := 1 to nuped do begin
		likelihood(i, proband[i]);
		alike := alike + like;
		if byfamily then begin
		    if finalprn then 
			write(final, proband[i]^.ped: 9, ' ', thisval[i] - like: 12: 6);
		    if finalprn then 
			write(final, (thisval[i] - like) / log10: 12: 6);
		    if like = zerolike then 
			if finalprn then 
			    writeln(final, ' inconsistent data')
			else if finalprn then 
			    writeln(final);
		    if like <> zerolike then 
			thisval[i] := thisval[i] - like
		end
	    end;
	    f := -2 * alike;
	    penlike := f;
	    dolod := false
	end else 
	    f := 1.5 * penlike
    end; { getlods }


begin
    firstapprox := true;
    lasttime := true;
    setparam;
    recombination;
    fun(f, x);
    if finalprn then 
	writeln(final, 'CHROMOSOME ORDER OF LOCI : ');
    for i := 1 to nlocus do begin
	j := 0;
	repeat
	    j := j + 1
	until order[j] = i;
	if finalprn then 
	    write(final, j: 3)
    end;
    if finalprn then 
	writeln(final);
    if itsys <> 0 then begin
	if finalprn then 
	    writeln(final, '****************** FINAL VALUES **********************');
	if finalprn then 
	    writeln(final, 'PROVIDED FOR LOCUS ', order[itsys]: 3, ' (CHROMOSOME ORDER)');
	if finalprn then 
	    writeln(final, '******************************************************');
	with thislocus[itsys]^ do begin
	    if disequi then begin
		if finalprn then 
		    writeln(final, 'HAPLOTYPE FREQUENCIES:');
		if finalprn then 
		    for i := 1 to nuhap do 
			write(final, hapfreq^.genarray[i]: 9: 6)
	    end else begin
		if finalprn then 
		    writeln(final, 'GENE FREQUENCIES :');
		if finalprn then 
		    for i := 1 to nallele do 
			write(final, freq[i]: 9: 6)
	    end;
	    if finalprn then 
		writeln(final);
	    if which = quantitative then begin
		if finalprn then 
		    writeln(final, 'VALUES FOR GENOTYPE MEANS:');
		invert(vmat, ntrait, det);
		if finalprn then 
		    for i := 1 to ntrait do begin
			for j := 1 to nallele do 
			    for k := j to nallele do 
				write(final, pm[j, k, i]: 6: 3, ' ');
			if finalprn then 
			    writeln(final)
		    end;
		if finalprn then 
		    writeln(final, 'COVARIANCE MATRIX:');
		if finalprn then 
		    for i := 1 to ntrait do begin
			for j := 1 to i do 
			    write(final, vmat[i, j]: 6: 3, ' ');
			writeln(final)
		    end
	    end else if which = affection then begin
		if finalprn then 
		    writeln(final, 'PENETRANCES:');
		for l := 1 to nclass do begin
		    if finalprn then 
			for i := 1 to nallele do 
			    for j := i to nallele do 
				write(final, pen[i, j, 2, l]: 5: 3, ' ');
		    if finalprn then 
			writeln(final);
		    if finalprn then 
			if sexlink then 
			    for i := 1 to nallele do 
				write(final, pen[0, i, 2, l]: 5: 3, ' ');
		    if finalprn then 
			writeln(final)
		end
	    end
	end
    end;
    if finalprn then 
	writeln(final, '******************************************************');
    if interfer and not mapping then begin
	if finalprn then 
	    writeln(final, 'P VALUES:');
	with maletheta^ do begin
	    oldtheta := theta;
	    theta[1] := (oldtheta[1] + oldtheta[nlocus] - oldtheta[nlocus - 1]) / 2.0;
	    theta[nlocus - 1] := (oldtheta[nlocus - 1] + oldtheta[nlocus] - oldtheta[1]) / 2.0;
	    theta[nlocus] := (oldtheta[1] + oldtheta[nlocus - 1] - oldtheta[nlocus]) / 2.0;
	    if finalprn then 
		for i := 1 to nlocus do 
		    write(final, ' ', theta[i]: 5: 3);
	    if finalprn then 
		writeln(final);
	    theta := oldtheta
	end;
	if sexdif then 
	    with femaletheta^ do begin
		oldtheta := theta;
		theta[1] := (oldtheta[1] + oldtheta[nlocus] - oldtheta[nlocus - 1]) / 2.0;
		theta[nlocus - 1] := (oldtheta[nlocus - 1] + oldtheta[nlocus] - oldtheta[1]) / 2.0;
		theta[nlocus] := (oldtheta[1] + oldtheta[nlocus - 1] - oldtheta[nlocus]) / 2.0;
		if finalprn then 
		    writeln(final, 'FEMALE:  ');
		if finalprn then 
		    for i := 1 to nlocus do 
			write(final, ' ', theta[i]: 5: 3);
		if finalprn then 
		    writeln(final);
		theta := oldtheta
	    end;
	if (ivar = 1) and (icall = 0) then begin
	    for i := 0 to 2 do 
		deriv[i] := exp(xall[nall - i]) / sqr(1 + exp(xall[nall - i]));
	    k := -1;
	    for i := 0 to 2 do 
		if itp[nall - i] = 1 then begin
		    k := k + 1;
		    l := -1;
		    for j := 0 to 2 do 
			if itp[nall - j] = 1 then begin
			    l := l + 1;
			    bmat[n - k, n - l] := bmat[n - k, n - l] * deriv[j] * deriv[i]
			end;
		    if l = -1 then 
			l := 0;
		    for j := n - l downto 1 do begin
			bmat[j, n - k] := bmat[j, n - k] * deriv[i];
			bmat[n - k, j] := bmat[j, n - k]
		    end
		end
	end
    end;
    if finalprn then 
	writeln(final, 'THETAS:');
    if finalprn then 
	with maletheta^ do 
	    for i := 1 to nlocus - 1 do 
		write(final, ' ', theta[i]: 5: 3);
    if interfer then 
	if finalprn then 
	    writeln(final, ' ', maletheta^.theta[nlocus]: 5: 3)
	else if finalprn then 
	    writeln(final);
    if sexdif then begin
	if finalprn then 
	    writeln(final, 'FEMALE:');
	if finalprn then 
	    with femaletheta^ do 
		for i := 1 to nlocus - 1 do 
		    write(final, ' ', theta[i]: 5: 3);
	if interfer then 
	    if finalprn then 
		writeln(final, ' ', femaletheta^.theta[nlocus]: 5: 3)
	    else if finalprn then 
		writeln(final);
	if readfemale then begin
	    if finalprn then 
		writeln(final, 'FEMALE/MALE DIST RATIO :');
	    for i := 1 to nlocus - 1 do begin
		dist := getdist(maletheta^.theta[i]);
		if dist <> 0.0 then 
		    dist := getdist(femaletheta^.theta[i]) / dist
		else 
		    dist := 0.0;
		if finalprn then 
		    write(final, ' ', dist: 5: 3)
	    end;
	    if interfer then begin
		dist := getdist(maletheta^.theta[nlocus]);
		if dist <> 0.0 then 
		    dist := getdist(femaletheta^.theta[nlocus]) / dist
		else 
		    dist := 0.0;
		if finalprn then 
		    write(final, ' ', dist: 5: 3)
	    end;
	    if finalprn then 
		writeln(final)
	end else begin
	    if finalprn then 
		writeln(final, 'CONSTANT FEMALE/MALE DIST RATIO :');
	    if finalprn then 
		writeln(final, distratio: 5: 3)
	end
    end;
    if finalprn then 
	writeln(final, '******************************************************');
    if finalprn then 
	writeln(final, '-2 LN(LIKE) = ', f);
    if not itmax and (nrep = 1) then begin
	writeln(isim, 'THETAS:');
	with maletheta^ do 
	    for i := 1 to nlocus - 1 do 
		write(isim, ' ', theta[i]: 5: 3);
	if interfer then 
	    writeln(isim, ' ', maletheta^.theta[nlocus]: 5: 3)
	else 
	    writeln(isim);
	if sexdif then begin
	    writeln(isim, 'FEMALE:');
	    with femaletheta^ do 
		for i := 1 to nlocus - 1 do 
		    write(isim, ' ', theta[i]: 5: 3);
	    if interfer then 
		writeln(isim, ' ', femaletheta^.theta[nlocus]: 5: 3)
	    else 
		writeln(isim);
	    if readfemale then begin
		writeln(isim, 'FEMALE/MALE DIST RATIO :');
		for i := 1 to nlocus - 1 do begin
		    ddist := getdist(maletheta^.theta[i]);
		    if ddist <> 0.0 then 
			ddist := getdist(femaletheta^.theta[i]) / ddist
		    else 
			ddist := 0.0;
		    write(isim, ' ', ddist: 5: 3)
		end;
		if interfer then begin
		    ddist := getdist(maletheta^.theta[nlocus]);
		    if ddist <> 0.0 then 
			ddist := getdist(femaletheta^.theta[nlocus]) / ddist
		    else 
			ddist := 0.0;
		    write(isim, ' ', ddist: 5: 3)
		end;
		writeln(isim)
	    end else begin
		writeln(isim, 'CONSTANT FEMALE/MALE DIST RATIO :');
		writeln(isim, distratio: 5: 3)
	    end
	end
    end;
    getlods(lods, x);
    lods := (lods - f) / (2 * log10);
    if nlocus <> 2 then 
	if finalprn then 
	    writeln(final, 'OTTS GENERALIZED LOD SCORE  =', lods)
	else if finalprn then 
	    writeln(final, 'LOD SCORE  =', lods);
    if lods >= lodlimit[1] then 
	clods[1] := clods[1] + 1;
    if lods >= lodlimit[2] then 
	clods[2] := clods[2] + 1;
    if lods >= lodlimit[3] then 
	clods[3] := clods[3] + 1;
    if nrep = 1 then begin
	aver := lods;
	variance := 0.0;
	smll := aver;
	big := aver
    end else begin
	ltemp := lods;
	rep := nrep;
	{  rep := double(nrep);}
	aver := (rep - 1.0) * aver / rep + ltemp / rep;
	variance := (rep - 1.0) * variance / rep + sqr(ltemp - aver) / (rep - 1.0);
	smll := min(smll, ltemp);
	big := max(big, ltemp)
    end;
    if finalprn then 
	writeln(final, 'NUMBER OF ITERATIONS = ', nit: 5);
    if finalprn then 
	writeln(final, 'NUMBER OF FUNCTION EVALUATIONS = ', nfe: 5);
    if finalprn then 
	writeln(final, 'PTG = ', ptg);
    idg := idg + 1;
    if prn then 
	writeln(outfile, 'EXIT CONDITION', idg);
    if prn then 
	case idg of
	    1:
		writeln(outfile, 'Maximum possible accuracy reached');
	    2:
		writeln(outfile, 'Search direction no longer downhill');
	    3:
		writeln(outfile, 'Accumulation of rounding error prevents further progress');
	    4:
		writeln(outfile, 'All significant differences lost through cancellation in conditioning');
	    5:
		writeln(outfile, 'Specified tolerance on normalized gradient met');
	    6:
		writeln(outfile, 'Specified tolerance on gradient met');
	    7:
		writeln(outfile, 'Maximum number of iterations reached');
	    8:
		writeln(outfile, 'Excessive cancellation in gradient')
	end;
    if (ivar = 1) and (icall = 0) then begin
	if finalprn then 
	    writeln(final, 'VARIANCE-COVARIANCE OF THE ESTIMATES');
	if interfer and not mapping then 
	    if finalprn then 
		writeln(final, 'VALUES GIVEN FOR P VALUES');
	for i := 1 to n do begin
	    for j := 1 to n do 
		if finalprn then 
		    write(final, 2.0 * bmat[i, j], ' ');
	    if finalprn then 
		writeln(final)
	end
    end;
    if finalprn then 
	writeln(final, '******************************************************');
    if finalprn then 
	writeln(final, '******************************************************')
end; { out }

procedure step;

var
    cont: (go, exit1, exit2, exit3);

    procedure firststep;

    var
	i: integer;

    begin
	sumt := 0.0;
	idg := 0;
	for i := 1 to n do 
	    xt[i] := x[i] + t * p[i];
	fun(ft, xt);
	nfe := nfe + 1
    end; { firststep }


    procedure decreaset;

    var
	i: integer;

    begin
	{Modified by M. Lathrop 29/04/86 to trap tmin=0 problem}
	if tmin < small then 
	    tmin := small;
	cont := go;
	while cont = go do begin
	    if (f - fsave = 0.0) and (idif = 2) then 
		idg := 2;
	    if t < tmin then 
		cont := exit3
	    else begin
		t := 0.5 * t;
		for i := 1 to n do 
		    xt[i] := x[i] + t * p[i];
		fun(ft2, xt);
		nfe := nfe + 1;
		if ft2 < f then begin
		    sumt := sumt + t;
		    f := ft2;
		    idg := 0;
		    x := xt;
		    if sumt < tmin then 
			cont := exit3
		    else 
			cont := exit2
		end else begin
		    if ft + f - ft2 - ft2 <= 0.0 then 
			scal := 0.1
		    else begin
			scal := 1.0 + 0.5 * (f - ft) / (f + ft - ft2 - ft2);
			if scal < 0.1 then 
			    scal := 0.1
		    end;
		    t := scal * t;
		    for i := 1 to n do 
			xt[i] := x[i] + t * p[i];
		    fun(ft, xt);
		    nfe := nfe + 1;
		    if f > ft then begin
			sumt := sumt + t;
			idg := 0;
			f := ft;
			x := xt;
			if t < tmin then 
			    cont := exit3
			else 
			    cont := exit2
		    end
		end
	    end
	end
    end; { decreaset }


    procedure increaset;

    var
	i: integer;

    begin
	sumt := sumt + t;
	cont := go;
	while cont = go do begin
	    twot := t + t;
	    if (ibnd > 0) and (tbnd >= 0.0) and (twot > tbnd) then begin
		f := ft;
		x := xt;
		if prn then 
		    writeln(outfile, '****** ACTIVE BOUNDARY CONSTRAINT *****');
		cont := exit1
	    end else begin
		x := xt;
		for i := 1 to n do 
		    xt[i] := x[i] + t * p[i];
		fun(f2t, xt);
		nfe := nfe + 1;
		if f2t > ft then begin
		    f := ft;
		    cont := exit2
		end else if (f2t + f - ft - ft < 0.0) or (ft - f2t - curv * (f - ft) >= 0.0) then begin
		    sumt := sumt + t;
		    t := t + t;
		    ft := f2t
		end else begin
		    sumt := sumt + t;
		    f := f2t;
		    x := xt;
		    cont := exit2
		end
	    end
	end
    end; { increaset }


begin
    firstapprox := true;
    firststep;
    if f > ft then 
	increaset
    else 
	decreaset;
    if cont = exit2 then 
	t := sumt;
    if cont = exit3 then 
	iret := 1
    else begin
	gsave := g;
	isw := 0;
	iret := 2
    end
end; { step }

{UPDATE}

procedure update;

var
    j, jm1: integer;
    cont, dfp: boolean;

    procedure prep;

    var
	i, k: integer;

    begin
	iswup := 1;
	wtil[1] := -gsave[1];
	ztil[1] := y[1];
	for k := 2 to n do begin
	    wtil[k] := -gsave[k];
	    ztil[k] := y[k];
	    for i := 1 to k - 1 do begin
		tl[i, k] := tl[k, i];
		wtil[k] := wtil[k] - tl[k, i] * wtil[i];
		ztil[k] := ztil[k] - tl[k, i] * ztil[i]
	    end
	end;
	sb := 0.0;
	sc := 0.0;
	sd := 0.0;
	for i := 1 to n do begin
	    sb := sb + t * ztil[i] * wtil[i] / d[i];
	    sc := sc + t * t * wtil[i] * wtil[i] / d[i];
	    sd := sd + ztil[i] * ztil[i] / d[i]
	end
    end; { prep }


    procedure sr1update;

    begin
	alpha := -1.0;
	if (sc - sb = 0.0) or (sb - sd = 0.0) or (sd - 2.0 * sb + sc = 0.0) then 
	    idg := 3
	else begin
	    sa := 1.0 / (sqrt(abs(sc - sb)) * sqrt(abs(sb - sd)));
	    thet1 := -((sd - sb) * sa + 1.0) / (sd - 2.0 * sb + sc);
	    thet2 := sa + (sa * (sb - sc) + 1.0) / (sd - 2.0 * sb + sc)
	end
    end; { sr1update }


    procedure sr2update;

    begin
	aa := sb / sc - 2.0 * (sd / sb) + sd / sc;
	bb := sb / sc - 1.0;
	cc := 1.0 - sb / sd;
	del2 := bb * bb - aa * cc;
	dfp := true;
	if del2 > 0.00000001 then begin
	    dfp := false;
	    del := sqrt(del2);
	    alph1 := (-bb + del) / aa;
	    alph2 := (-bb - del) / aa;
	    if abs(alph1) < abs(alph2) then 
		alpha := alph1
	    else 
		alpha := alph2;
	    if abs(alpha) < 0.00001 then 
		dfp := true
	    else begin
		sa := (alpha + 1.0) * (alpha + 1.0) + sc / sb - alpha * alpha * (sc / sb) * (sd / sb) - 1.0 + alpha * alpha * sd / sb;
		if sa <= 0.0 then 
		    sa := 0.0
		else begin
		    sa := sqrt(sa);
		    sa := 1.0 / (sa * sb)
		end;
		rdiv := 1.0 / (alpha * alpha * sd + 2.0 * alpha * sb + sc);
		thet1 := -(sa * alpha * (alpha * sd + sb) + 1.0) * rdiv;
		thet2 := sa + (alpha * sa * (sc + alpha * sb) - alpha) * rdiv
	    end
	end;
	if dfp then begin
	    alpha := 0.0;
	    sa := 1.0 / (sqrt(sb) * sqrt(sc));
	    thet1 := -1.0 / sc;
	    thet2 := sa
	end
    end; { sr2update }


    procedure getwzs;

    var
	i, k: integer;

    begin
	for i := 1 to n do begin
	    w[i] := t * wtil[i] + alpha * ztil[i];
	    z[i] := t * thet1 * wtil[i] + thet2 * ztil[i]
	end;
	s[n] := 0.0;
	for k := 1 to n - 1 do begin
	    j := n + 1 - k;
	    jm1 := j - 1;
	    s[jm1] := s[j] + w[j] * w[j] / d[j]
	end;
	nu := 1.0;
	eta := 0.0
    end; { getwzs }


    procedure recur;

    var
	i, k: integer;

    begin
	cont := true;
	while cont do begin
	    cont := false;
	    if iswup < 2 then begin
		for i := 1 to n do 
		    wtjp1[i] := -gsave[i];
		ztjp1 := y
	    end else begin
		for i := 1 to n do begin
		    wtil[i] := -t * g[i] + alpha * y[i];
		    ztil[i] := -t * thet1 * g[i] + thet2 * y[i]
		end
	    end;
	    j := 0;
	    lambj2 := 0.0;
	    while j < n - 1 do begin
		j := j + 1;
		if iswup < 2 then 
		    for k := j + 1 to n do begin
			wtjp1[k] := wtjp1[k] - wtil[j] * tl[k, j];
			ztjp1[k] := ztjp1[k] - ztil[j] * tl[k, j]
		    end
		else 
		    for k := j + 1 to n do begin
			wtjp1[k] := wtil[k] - w[j] * tl[k, j];
			ztjp1[k] := ztil[k] - z[j] * tl[k, j]
		    end;
		aj := nu * z[j] - eta * w[j];
		thj := 1.0 + aj * w[j] / d[j];
		lambj2 := thj * thj + aj * aj * s[j] / d[j];
		if iswup < 2 then begin
		    if lambj2 > 10.0 then begin
			cont := true;
			iswup := 2;
			for k := 2 to n do 
			    for i := 1 to k - 1 do 
				tl[k, i] := tl[i, k];
			j := n
		    end
		end;
		if not cont then begin
		    dp[j] := d[j] * lambj2;
		    lambj := sqrt(lambj2);
		    if thj > 0.0 then 
			lambj := -lambj;
		    muj := thj - lambj;
		    bj := thj * w[j] + aj * s[j];
		    gamlj := bj * nu / (lambj2 * d[j]);
		    betlj := (aj - bj * eta) / (lambj2 * d[j]);
		    nu := -nu / lambj;
		    eta := -(eta + aj * aj / (muj * d[j])) / lambj;
		    if iswup < 2 then 
			for k := j + 1 to n do 
			    tl[k, j] := tl[k, j] + t * (betlj + thet1 * gamlj) * wtjp1[k] + (alpha * betlj + thet2 * gamlj) * ztjp1[k]
		    else 
			for k := j + 1 to n do begin
			    tl[k, j] := tl[k, j] / lambj2 + betlj * wtil[k] + gamlj * ztil[k];
			    wtil[k] := wtjp1[k];
			    ztil[k] := ztjp1[k]
			end
		end
	    end
	end;
	aj := nu * z[n] - eta * w[n];
	lambj := 1.0 + aj * w[n] / d[n];
	dp[n] := d[n] * lambj * lambj;
	d := dp
    end; { recur }


begin
    prep;
    if sb < sc then 
	fbcd := sb
    else 
	fbcd := sc;
    if sd < fbcd then 
	fbcd := sd;
    if fbcd > small then begin
	fbcd := 2.0 * sc * (sd / sb) / (sc + sd);
	if fbcd < 1.0 then 
	    sr1update
	else 
	    sr2update
    end else 
	sr1update;
    if idg <> 3 then begin
	getwzs;
	recur
    end
end; { update }

{START GEMINI}

procedure gemini;

    procedure bldlt(var b: matrix);

    var
	i, j, ic, k: integer;
	su, tlic, ff, hh, temp, temp1: double;
	s: matrix;
	tu, xvec: vector;

    begin
	firstapprox := true;
	ff := f;
	if (icall > 0) and (ihess <= 0) then begin
	    xvec := xall;
	    hh := 10 * gh;
	    for i := 1 to n do begin
		temp := x[i];
		x[i] := x[i] + hh;
		fun(f, x);
		tu[i] := f;
		for j := 1 to i do begin
		    temp1 := x[j];
		    x[j] := x[j] + hh;
		    fun(f, x);
		    b[i, j] := f;
		    b[j, i] := b[i, j];
		    x[j] := temp1
		end;
		x[i] := temp
	    end;
	    xall := xvec;
	    for i := 1 to n do 
		for j := 1 to i do begin
		    b[i, j] := (ff + b[i, j] - tu[i] - tu[j]) / (hh * hh);
		    b[j, i] := b[i, j]
		end
	end;
	ic := 1;
	while (b[ic, ic] > 0.0) and (ic <= n) do begin
	    temp := b[ic, ic];
	    d[ic] := b[ic, ic];
	    tl[ic, ic] := 1.0;
	    if ic <> n then begin
		for k := ic + 1 to n do 
		    tl[k, ic] := b[k, ic] / temp;
		for i := ic + 1 to n do begin
		    tlic := tl[i, ic];
		    for k := i to n do 
			b[k, i] := b[k, i] - tlic * tl[k, ic] * temp
		end
	    end;
	    ic := ic + 1
	end;
	if ic > n then begin
	    icall := icall - 1;
	    if prn then 
		writeln(outfile, 'FACTORIZATION SUCCEEDED')
	end else begin
	    icall := icall + 1;
	    if prn then 
		writeln(outfile, 'FACTORIZATION FAILED')
	end;
	if icall = 0 then begin
	    s[1, 1] := 1.0;
	    for i := 2 to n do begin
		for k := 1 to i - 1 do begin
		    su := 0.0;
		    for j := k to i - 1 do 
			su := su + tl[i, j] * s[j, k];
		    s[i, k] := -su
		end;
		s[i, i] := 1.0
	    end;
	    for i := 1 to n do 
		for j := i to n do begin
		    su := 0.0;
		    for k := j to n do 
			su := su + s[k, i] * s[k, j] / d[k];
		    b[i, j] := su;
		    b[j, i] := su
		end;
	    if prn then 
		writeln(outfile, 'B-MATRIX');
	    if prn then 
		for i := 1 to n do begin
		    for j := 1 to i do 
			write(outfile, b[i, j]);
		    writeln(outfile)
		end;
	    if ivar = 1 then begin
		k := 0;
		for i := 1 to nall do begin
		    se[i] := 0.0;
		    if itp[i] = 1 then begin
			k := k + 1;
			se[i] := sqrt(b[k, k])
		    end
		end
	    end
	end
    end; { bldlt }


    procedure inib;

    var
	i, j: integer;
	in1: text;

    begin
	if (icall <> 1) and (ihess > 1) then begin
{    assign(in1,'IN1.DAT'); }
	    {SUN}
	    writeln('Reading IN1.DAT (covariance matrix)');
	    reset(in1, 'in1.dat'); {SUN}
	    for i := 1 to n do 
		for j := 1 to i do 
		    read(in1, bmat[i, j]);
	    for i := 1 to n do 
		for j := 1 to i do 
		    bmat[j, i] := bmat[i, j]
	end;
	bldlt(bmat);
	if icall = 1 then begin
	    for i := 1 to n do begin
		for j := 1 to i do 
		    tl[i, j] := 0.0;
		tl[i, i] := 1.0;
		d[i] := 1.0
	    end
	end
    end; { inib }


    procedure initialize;

    var
	i, j: integer;

    begin
	firstapprox := true;
	if prn then 
	    writeln(outfile, 'DIFFER INTER = ', gh, ' TRUNC UPPER = ', trupb);
	nit := 0;
	idg := 0;
	idif := 1;
	t := 0.1;
	tmin := 0.0;
	fsmf := 0.0;
	fun(f, x);
	nfe := 1;
	for i := 1 to n do begin
	    for j := 1 to i do 
		tl[i, j] := 0.0;
	    tl[i, i] := 1.0;
	    d[i] := 1.0
	end;
	if ihess > 0 then begin
	    icall := 0;
	    inib;
	    icall := 1;
	    t := 1.0
	end
    end; { initialize }
 {initialize}

    procedure gforward;

    var
	i: integer;

    begin
	firstapprox := true;
	for i := 1 to n do begin
	    hx := gh;
	    if ihx = 1 then 
		if abs(x[i]) < 0.1 then 
		    hx := gh * 0.1
		else 
		    hx := gh * abs(x[i]);
	    xsave := x[i];
	    x[i] := x[i] + hx;
	    fun(fxph, x);
	    firstapprox := false;
	    g[i] := (fxph - f) / hx;
	    x[i] := xsave
	end;
	nfe := nfe + n
    end; { gforward }


    procedure gcentral;

    var
	i: integer;

    begin
	firstapprox := true;
	for i := 1 to n do begin
	    hx := gh;
	    if ihx = 1 then 
		if abs(x[i]) < 0.1 then 
		    hx := gh * 0.1
		else 
		    hx := gh * abs(x[i]);
	    xsave := x[i];
	    x[i] := x[i] + hx;
	    fun(fxph, x);
	    firstapprox := false;
	    x[i] := xsave - hx;
	    fun(fxmh, x);
	    g[i] := (fxph - fxmh) / (hx + hx);
	    x[i] := xsave
	end;
	nfe := nfe + n + n
    end; { gcentral }


    procedure getp;

    var
	i, nmj, j: integer;

    begin
	nit := nit + 1;
	fsav2 := fsave;
	fsave := f;
	if prn then 
	    writeln(outfile);
	if prn then 
	    writeln(outfile, 'ITERATION ', nit: 5, ' T = ', t: 10: 3, ' NFE = ', nfe: 5, ' F = ', f);
	if prn then 
	    writeln(output, 'ITERATION ', nit: 5, ' T = ', t: 10: 3, ' NFE = ', nfe: 5, ' F = ', f);
	if nit > maxit then 
	    idg := 6
	else begin
	    if n < 20 then begin
		if prn then 
		    write(outfile, 'X= ');
		outcontrol(x);
		if prn then 
		    write(outfile, 'G= ');
		outcontrol(g)
	    end;
	    if n = 1 then 
		p[1] := -g[1] / d[1]
	    else begin
		gv[1] := -g[1];
		for i := 2 to n do begin
		    gv[i] := -g[i];
		    for j := 1 to i - 1 do 
			gv[i] := gv[i] - tl[i, j] * gv[j]
		end;
		p[n] := gv[n] / d[n];
		for j := 1 to n - 1 do begin
		    nmj := n - j;
		    p[nmj] := gv[nmj] / d[nmj];
		    for i := nmj + 1 to n do 
			p[nmj] := p[nmj] - tl[i, nmj] * p[i]
		end
	    end;
	    if prn then 
		write(outfile, 'P= ');
	    outcontrol(p);
	    if ibnd = 1 then 
		for i := 1 to n do 
		    if p[i] = 0.0 then 
			idg := 7;
	    if idg <> 7 then begin
		ptg := 0.0;
		for i := 1 to n do 
		    ptg := ptg + p[i] * g[i];
		if ptg >= 0.0 then 
		    idg := 1
		else begin
		    if abs(ptg) < tol then begin
			idg := 4;
			fsmf := fsav2 - f;
			if prn then 
			    writeln(outfile, 'FSMF = ', fsmf, ' PTG = ', ptg, ' TMIN = ', tmin)
		    end else begin
			xpm := xpmcon;
			if nit <> 1 then begin
			    for i := 1 to n do begin
				if ihx = 1 then begin
				    xp := abs(x[i]);
				    if xp < 0.1 then 
					xp := 0.1;
				    xp := xp / abs(p[i]);
				    if xp < xpm then 
					xpm := xp
				end else if xpm > abs(1.0 / p[i]) then 
				    xpm := abs(1.0 / p[i])
			    end;
			    tmin := 0.5 * xpm * gh / trupb;
			    if idif = 2 then 
				t := 1.0
			    else begin
				t := -2.0 * fsmf / ptg;
				if t <= 0.0 then 
				    t := 1.0
				else if 1.0 < t then 
				    t := 1.0
			    end;
			    if prn then 
				writeln(outfile, 'FSMF = ', fsmf, ' PTG =', ptg, ' TMIN=', tmin);
			    if prn then 
				writeln(outfile, 'INITIAL T = ', t)
			end
		    end
		end
	    end
	end
    end; { getp }


    procedure getytp;

    var
	i: integer;

    begin
	ytp := 0.0;
	fsmf := fsave - f;
	for i := 1 to n do begin
	    y[i] := g[i] - gsave[i];
	    ytp := ytp + y[i] * p[i]
	end
    end; { getytp }


    procedure chkbnd;
    {This procedure modified by M. Lathrop 29/04/86
    with the introduction of continue and repeat}

    var
	i, j, k, ik, ii, jk: integer;
	continue: boolean;

    begin
	repeat
	    clb := clbcon;
	    for j := 1 to 2 do begin
		k := 0;
		for i := 1 to nall do begin
		    if itp[i] = 1 then begin
			k := k + 1;
			check := abs(x[k] - bnd[i, j]);
			if check <= clb then begin
			    clb := check;
			    ik := k;
			    ii := i;
			    jk := j
			end
		    end
		end
	    end;
	    if clb < 0.1 * gh then begin
		ihess := 0;
		if jk = 2 then 
		    eh := -gh - gh
		else 
		    eh := gh + gh;
		x[ik] := bnd[ii, jk] + eh;
		xall[ii] := x[ik];
		itp[ii] := 0;
		k := 0;
		for i := 1 to nall do 
		    if itp[i] = 1 then begin
			k := k + 1;
			x[k] := xall[i]
		    end;
		n := k;
		tol := tolconst;
		{            tol:=tolconst*sqrt(n);}
		continue := true;
		if prn then 
		    writeln(outfile, '******* A VARIABLE WAS SET TO A BOUND ***********')
	    end else begin
		tbnd := tbndcon;
		for j := 1 to 2 do begin
		    k := 0;
		    for i := 1 to nall do 
			if itp[i] = 1 then begin
			    k := k + 1;
			    teq := (bnd[i, j] - x[k]) / p[k];
			    if (teq < tbnd) and (teq > 0.0) then 
				tbnd := teq
			end
		end;
		continue := false;
		if t * (2.0 + gh) >= tbnd then 
		    t := tbnd * (0.5 - gh);
		if prn then 
		    writeln(outfile, 'TBND = ', tbnd, ' RESET T = ', t)
	    end
	until not continue
    end; { chkbnd }


    procedure iterate;

    begin
	while continue = go do begin
	    active := false;
	    getp;
	    if idg <> 0 then 
		continue := quit;
	    if idg = 0 then begin
		iret := 2;
		if ibnd = 1 then 
		    chkbnd;
		if n = 0 then 
		    continue := quit;
		if not active then 
		    step
	    end;
	    { THE NEXT IS NOT TRUE IF ACTIVE }
	    if iret <> 2 then begin
		if idif = 1 then begin
		    idif := 2;
		    isw := 1;
		    gcentral;
		    if prn then 
			writeln(outfile, '***** SWITCHING TO CENTRAL DIFFERERENCE*****')
		end else 
		    continue := quit
	    end else if not active then begin
		if idif = 1 then 
		    gforward
		else 
		    gcentral;
		if isw = 0 then begin
		    getytp;
		    if ytp > 0.0 then begin
			if n = 1 then 
			    d[1] := -y[1] * d[1] / (t * gsave[1])
			else 
			    update;
			if idg = 3 then 
			    continue := quit
		    end
		end
	    end else 
		continue := restart
	end
    end; { iterate }


begin				{gemini}
    continue := go;
    while continue <> quit do begin
	initialize;
	isw := 0;
	gforward;
	continue := go;
	iterate
    end;
    if ivar = 1 then 
	inib
end; { gemini }				{gemini}


procedure initisim;
var
    i: integer;
begin
    writeln(output, 'Program ISIM version ', version);
    writeln(output);
    writeln(output, 'The program constants are set to the following maxima:');
    writeln(output, maxlocus: 6, ' loci in mapping problem');
    writeln(output, maxall: 6, ' alleles at a single locus');
    writeln(output, maxneed: 6, ' recombination probabilities (maxneed)');
    writeln(output, maxcensor: 6, ' maximum of censoring array (maxcensor)');
    writeln(output, maxhap: 6, ' haplotypes = n1 x n2 x ... where ni = current # alleles at locus i');
    writeln(output, maxfem: 6, ' joint genotypes for a female');
    writeln(output, maxmal: 6, ' joint genotypes for a male');
    writeln(output, maxind: 6, ' individuals in all pedigrees combined');
    writeln(output, maxped: 6, ' pedigrees (maxped');
    writeln(output, maxtrait: 6, ' quantitative factor(s) at a single locus');
    writeln(output, maxliab: 6, ' liability classes');
    writeln(output, maxfact: 6, ' binary codes at a single locus');
    writeln(output, scale: 8: 2, ' base scaling factor for likelihood (scale)');
    writeln(output, scalemult: 8: 2, ' scale multiplier for each locus (scalemult)');
    writeln(output, minfreq: 8: 5, ' frequency for elimination of heterozygotes (minfreq)');
    if minfreq <> 0.0 then begin
	writeln(output, 'IMPORTANT : RECOMPILE THIS PROGRAM WITH MINFREQ=0.0');
	writeln(output, 'FOR THE ANALYSIS OF RECESSIVE TRAITS')
    end;
    writeln(output);
{ assign(datafile,'datafile.dat');
  assign(speedfile,'speedfil.dat');
  assign(ipedfile,'ipedfile.dat');
  assign(outfile,'outfile.dat');
  assign(final,'final.dat');}
    {SUN}
    for i := 1 to 48 do 
	write(output, '*');
    writeln(output);
    writeln(output, ' This is ISIM, which is a modified version of ILINK');
    writeln(output, ' designed to be used on datafiles generated by the ');
    writeln(output, ' companion program SLINK.');
    writeln(output);
    writeln(output, ' Pedigree data are read from the file ipedfile.dat');
    writeln(output, ' The statistical summary is written to isim.dat ');
    for i := 1 to 48 do 
	write(output, '*');
    writeln(output);

    {Read the number of pedigrees per replicate from the most recent simout.dat}
{ assign(simout,'simout.dat');}
    {SUN}
    writeln('Reading SIMOUT.DAT');
    reset(simout, 'simout.dat'); {SUN}
    for i := 1 to 5 do 
	readln(simout); {The no. of pedigrees is on the sixth line}
    repeat
	read(simout, chtemp)
    until chtemp = 's'; {Read to the 's' in 'pedigrees' }
    readln(simout, opeds); {Read the number of pedigrees }
    if opeds > 1 then 
	writeln(output, ' There are ', opeds: 2, ' pedigrees per replicate.')
    else 
	writeln(output, ' There is ', opeds: 2, ' pedigree per replicate.')
end; { initisim }


procedure acknowl(var ff: text);
begin
    writeln(ff, 'Please use these two references when reporting results based on SLINK:');
    writeln(ff);
    writeln(ff, 'Ott J (1989) Computer-simulation methods in human linkage analysis.');
    writeln(ff, 'Proc Natl Acad Sci USA 86:4175-4178');
    writeln(ff);
    writeln(ff, 'Weeks DE, Ott J, Lathrop GM (1990) SLINK: a general simulation');
    writeln(ff, 'program for linkage analysis.  Am J Hum Genet 47(3):A204 (Supplement)')
end; { acknowl }


procedure getparam;
var
    limit: text;
    i: integer;
begin
{ assign(limit,'LIMIT.DAT');}
    {SUN}
    writeln('Reading LIMIT.DAT');
    reset(limit, 'limit.dat'); {SUN}
    for i := 1 to 3 do 
	lodlimit[i] := i; {Jurg 10/19/91}
    if not (eoln(limit) or eof(limit)) then 
	read(limit, lodlimit[1], lodlimit[2], lodlimit[3])
end; { getparam } {getparam}


begin				{isim}
    initisim;
    inputdata;
    getparam;
    {speedfile modification}
    lastseg := 0;
    segperson := 0;
    lastspeed := 0;
    close(datafile); {SUN}
    gh := sqrt(exp(-nbit * ln(2.0)));
    tol := tolconst;
    { tol:=tolconst*sqrt(n);}
    trupb := sqrt(gh);
    new(censorstruct);
    {  censorstruct:=censorpnt(NewPtr(SizeOf(censorrec)));}
    new(gennustruct);
    {  gennustruct:=gennuptr(NewPtr(SizeOf(gennurec)));}
    if approximate then 
	new(approxstruct);
    {   approxstruct:=approxpnt(NewPtr(SizeOf(approxrec)));}
    getlocations;
    holditp := itp;
    itmax := false; {Flag indicating iterated parameters}
    for i := 1 to nall do 
	if itp[i] = 1 then 
	    itmax := true;
    holdxall := xall;
    for i := 1 to 3 do 
	clods[i] := 0;
    nrep := 1;
{assign(isim,'isim.dat');}
    {SUN}
    rewrite(isim, 'isim.dat'); {SUN}
    { Copy most recent simout.dat to isim.dat }
    writeln(isim, ' ********* Data from most recent SIMOUT.DAT *********');
{assign(simout,'simout.dat');}
    {SUN}
    reset(simout, 'simout.dat'); {SUN}
    while not eof(simout) do begin
	read(simout, chtemp);
	write(isim, chtemp);
	if eoln(simout) then begin
	    readln(simout);
	    writeln(isim)
	end
    end;
    writeln(isim, ' ********* End of most recent SIMOUT.DAT *********');
    writeln(isim);
    close(simout); {SUN}
    writeln(isim, ' ORDER OF LOCI : ');
    for i := 1 to nlocus do begin
	j := 0;
	repeat
	    j := j + 1
	until order[j] = i;
	write(isim, j: 3)
    end;
    writeln(isim);
    if nlocus <> 2 then 
	if itmax then 
	    writeln(isim, ' Average Maximum of Ott''s Generalized Lod Score')
	else 
	    writeln(isim, ' Average of Ott''s Generalized Lod Score')
    else if itmax then 
	writeln(isim, ' Average Maximum Lod Score')
    else 
	writeln(isim, ' Average Lod Score');
    writeln(isim);
    while not eof(ipedfile) do begin
	itp := holditp;
	xall := holdxall;
	n := 0;
	for i := 1 to nall do 
	    if itp[i] = 1 then begin
		n := n + 1;
		x[n] := xall[i]
	    end;
	ihess := 0;
	ibnd := 1;
	icall := 1;
	ivar := 0;
	ihx := 1;
	f := 0;
	ptg := 0.0;
	inconsistent := false;
	firsttime := true;
	lasttime := false;
	dolod := false;
	firstapprox := true;
	maxit := 50 * n;
	readpedseg;
	for i := 1 to totperson do 
	    person[i]^.store := nil;
	{ Set all store pointers to NIL in order to be able to use Dispose later }
	readspseg;
	gemini;
	out;
	nrep := nrep + 1;
	if nrep mod 10 = 0 then 
	    write(output, '*')
	else 
	    write(output, '.');
	flush(output); {SUN}
	if nrep mod 50 = 0 then 
	    writeln(output, ' ', nrep: 1, ' replicates analyzed');
	{ Need to reuse memory once the pedigrees have been processed }
	for i := 1 to totperson do 
	    if person[i] <> nil then begin
		with person[i]^ do 
		    for j := 1 to nlocus do 
			if phen[j] <> nil then 
			    dispose(phen[j]);
		{           DisposPtr(Ptr(phen[j]));}
		if person[i]^.store <> nil then 
		    dispose(person[i]^.store);
		{         DisposPtr(Ptr(person[i]^.store));}
		dispose(person[i]);
		{        DisposPtr(Ptr(person[i]));}
		person[i] := nil
	    end
    end {IF person[i]<>NIL};
    close(ipedfile); {SUN}
    close(speedfile); {SUN}
    close(outfile); {SUN}
    close(final); {SUN}
    writeln(isim, 'The number of replicates is ', nrep - 1: 6);
    for i := 1 to 66 do 
	write(isim, '-');
    writeln(isim);
    if nlocus <> 2 then 
	if itmax then 
	    writeln(output, ' Average Maximum of Ott''s Generalized Lod Score')
	else 
	    writeln(output, ' Average of Ott''s Generalized Lod Score')
    else if itmax then 
	writeln(output, ' Average Maximum Lod Score')
    else 
	writeln(output, ' Average Lod Score');
    writeln(output);
    writeln(output, 'The number of replicates is ', nrep - 1: 6);
    for i := 1 to 66 do 
	write(output, '-');
    writeln(output);
    if itmax then 
	writeln(output, ' Average Maximum | StdDev       | Min          | Max')
    else 
	writeln(output, '      Average    | StdDev       | Min          | Max');
    if itmax then 
	writeln(isim, ' Average Maximum | StdDev       | Min          | Max')
    else 
	writeln(isim, '      Average    | StdDev       | Min          | Max');
    writeln(output, aver: 15: 6, sqrt(variance): 15: 6, smll: 15: 6, big: 15: 6);
    writeln(isim, aver: 15: 6, sqrt(variance): 15: 6, smll: 15: 6, big: 15: 6);
    for i := 1 to 66 do 
	write(isim, '-');
    writeln(isim);
    for i := 1 to 66 do 
	write(output, '-');
    writeln(output);
    for i := 1 to 66 do 
	write(isim, '-');
    writeln(isim);
    for i := 1 to 66 do 
	write(output, '-');
    writeln(output);
    if itmax then 
	writeln(output, '   Number of maximum lod scores greater than ')
    else 
	writeln(output, '   Number of lod scores greater than ');
    writeln(output, '   a given constant ');
    for i := 1 to 66 do 
	write(output, '-');
    writeln(output);
    writeln(output, ' Constant  Number  Percent ');
    for i := 1 to 3 do 
	writeln(output, lodlimit[i]: 8: 3, clods[i]: 8, 100.0 * clods[i] / (nrep - 1): 9: 3);
    if itmax then 
	writeln(isim, '   Number of maximum lod scores greater than ')
    else 
	writeln(isim, '   Number of lod scores greater than ');
    writeln(isim, '   a given constant ');
    for i := 1 to 66 do 
	write(isim, '-');
    writeln(isim);
    writeln(isim, ' Constant  Number  Percent ');
    for i := 1 to 3 do 
	writeln(isim, lodlimit[i]: 8: 3, clods[i]: 8, 100.0 * clods[i] / (nrep - 1): 9: 3);
    for i := 1 to 66 do 
	write(isim, '-');
    writeln(isim);
    for i := 1 to 66 do 
	write(output, '-');
    writeln(output);
    acknowl(isim)
end. { isim }
