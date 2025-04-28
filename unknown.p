PROGRAM unknown(pedfile,datafile,speedfile,ipedfile,output);

CONST
  version=5.2;                  {PRESENT VERSION OF LINKAGE}
  maxlocus=40;                  {MAXIMUM NUMBER OF LOCI}
  maxall=9;                     {MAX NUMBER OF ALLELES AT A SINGLE LOCUS}
  maxgeno=45;                   {MAXIMUM NUMBER OF SINGLE LOCUS GENOTYPES}
  { = maxall(maxall+1)/2.  Different definition than in analysis programs!}
  maxind=500;                   {MAXIMUM NUMBER OF INDIVIDUALS IN A PEDIGREE}
  maxmarriage=3;                {MAXIMUM NUMBER OF MARRIAGES FOR ONE MALE}
  maxfact=9;                    {MAXIMUM NUMBER OF BINARY CODES AT A SINGLE LOCUS}
  maxtrait=3;                   {MAXIMUM NUMBER OF QUANTITATIVE FACTORS AT A SINGLE LOCUS}
  missval=0.0;                  {MISSING VALUES FOR QUANTITATIVE TRAITS}
  affall=2;                     {DISEASE ALLELE FOR QUANT. TRAITS OR AFFECTION STATUS}
  missaff=0;                    {MISSING VALUE FOR AFFECTION STATUS}
  affval=2;                     {CODE FOR AFFECTED INDIVIDUAL}
  maxliab=20;                   {MAXIMUM NUMBER OF LIABILITY CLASSES}

TYPE
  genotype=ARRAY[1..maxgeno] OF BOOLEAN;
  direction=(peelup,peeldown);
  binset=SET OF 1..maxfact;
  phenarray=ARRAY[1..maxall] OF binset;
  possvect=ARRAY[1..maxall,1..maxall] OF BOOLEAN;
  possarray=ARRAY[1..maxlocus] OF possvect;
  locuspoint=^locusvalues;
  phenpoint=^phenotype;
  locustype=(affection,quantitative,binary);
  locusvalues=
  RECORD
    nallele:INTEGER;
    CASE which:locustype OF
      quantitative:(ntrait:INTEGER);
      affection   :(pen:ARRAY[0..maxall,1..maxall,0..2,1..maxliab] OF REAL;nclass:INTEGER);
      binary      :(allele:phenarray;nfactor,format:INTEGER)
  END;
  phenotype=
  RECORD
    CASE which:locustype OF
      quantitative:(x:ARRAY[1..maxtrait] OF REAL;missing:BOOLEAN);
      affection   :(aff,liability:INTEGER);
      binary      :(phenf:binset)
  END;
  ind=^thisperson;
  genpoint=^thisarray;
  indphen=ARRAY[1..maxlocus] OF phenpoint;
  thisarray=
  RECORD
    genarray:genotype
  END;
  infoptr=^information;
  information=
  RECORD
    possible:possarray
  END;
  thisperson=
  RECORD
    id,paid,maid,offid,npaid,nmaid,sex,profield,oldped,nseq:INTEGER;
    pa,ma,foff,nextpa,nextma:ind;
    gen:genpoint;
    phen:indphen;
    store:infoptr;
    thisunknown:ARRAY[1..maxlocus] OF BOOLEAN;
    unknown,multi,done,up,male:BOOLEAN
  END;
  haplotype=(a,b);
  subhap=ARRAY[a..b] OF INTEGER;

VAR
  seghap:ARRAY[1..maxgeno] OF subhap;
  thislocus:ARRAY[1..maxlocus] OF locuspoint;
  person:ARRAY[0..maxind] OF ind;
  proband,loop1,loop2:ind;
  genenumber:ARRAY[1..maxgeno,1..maxgeno] OF INTEGER;
  risksys,mutsys,nsystem:INTEGER;
  fgeno,mgeno:0..maxgeno;
  nsequence,newped,whichsys,totperson:INTEGER;
  sexlink,risk,disequi:BOOLEAN;
  speedfile,datafile,pedfile,ipedfile:TEXT;
  gene:genotype;


  procedure inputerror(nerror,par1,par2 : INTEGER);
  begin
    writeln('Fatal error detected in procedure inputdata');
    CASE nerror OF
      0: writeln('Number of loci ',par1:2,' exceeds the constant maxlocus');
      1: writeln('Number of loci read ',par1:2,'. Less than minimum of 1');
      2: writeln('Error detected reading loci order. Locus number ',par2:2,' in position ',par1:2,' exceeds number of loci');
      3: writeln('Error detected reading loci order. Illegal locus number ',par2:2,' in position ',par1:2);
      4: writeln('Error detected reading loci order. Locus number repeated in positions ',par1:2,' and ',par2:2);
      5: writeln('Error detected reading locus description. Illegal locus type ',par2:2,' for locus ',par1:2);
      6: writeln('Error detected reading locus description for system ',par1:2,'. Number of alleles  ',par1:2,
                 ' exceeds maxall');
      7: writeln('Error detected reading locus description for system ',par1:2,
                 '. Illegal number of alleles  ',par2:2);
      8: writeln('Error detected reading locus description for system ',par1:2,'. Number of factors  ',par2:2,
                 ' exceeds maxfact');
      9: writeln('Error detected reading locus description for system ',par1:2,
                 '. Illegal number of factors  ',par2:2);
      10: writeln('Error detected reading locus description for system ',par1:2,'. Alleles not codominant');
      11: writeln('Error detected reading pedigree record ',par1:2,'. Illegal code for sex ',par2:2);
      12: writeln('Error detected reading pedigree record at pedigree',par1:2,'. Maximum number of pedigree records exceeded');
      13: writeln('Error detected reading pedigree record ',par1:2,'. Maximum number of individuals exceeded');
      14: writeln('Error detected reading pedigree record ',par1:2,'. Illegal binary factor code ',par2:2);
      15: writeln('Error detected reading pedigree record ',par1:2,'. No allelic pair for genotype');
      16: writeln('Error detected reading pedigree record ',par1:2,'. Allele number ',par2:2,' exceeds maxall');
      17: writeln('Error detected reading pedigree record ',par1:2,'. Illegal allele number ',par2:2);
      18: writeln('Number of systems after factorization (',par1:3,') exceeds maxsystem');
      19: writeln('Number of systems after factorization (',par1:3,') less than minimum of 1');
      20: writeln('Number of recombination types (',par1:3,') exceeds maxrectype');
      21: writeln('Number of recombination types (',par1:3,') less than minimum of 1');
      22: writeln('End of file detected in tempdat by procedure readthg before all data found');
      23: writeln('Error detected reading iterated locus in datafile. Value (',par1:3,') greater than nlocus');
      24: writeln('Error detected reading iterated locus in datafile. Illegal value (',par1:3,')');
      25: writeln('Number of iterated parameters greater then maxn');
      26: writeln('Error detected reading pedigree record ',par1:2,'. Liability class (',par2:2,') exceeds nclass');
      27: writeln('Error detected reading pedigree record ',par1:2,'. Illegal liability class (',par2:2,')');
      28: writeln('Error detected reading locus description for system',par1:2,
                  '. Liability classes (',par2:3,') exceed maxliab');
      29: writeln('Error detected reading locus description for system',par1:2,
                  '. Illegal number of liability classes (',par2:3,')');
      30: writeln('Error detected reading locus description for system',par1:2,'. Penetrance out of range');
      31: writeln('Error detected reading locus description for system',par1:2,
                  '. Number of traits (',par2:3,') exceeds maxtrait');
      32: writeln('Error detected reading locus description for system',par1:2,'. Number of traits out of range (',par2:3,')');
      33: writeln('Error detected reading locus description for system',par1:2,'. Variance must be positive');
      34: writeln('Error detected reading locus description for system',par1:2,'. Variance multiplier must be positive');
      35: writeln('Error detected reading locus description for system',par1:2,'. Risk allele ',par2:3,') exceeds nallele');
      36: writeln('Error detected reading locus description for system',par1:2,'. Illegal risk allele (',par2:3,')');
      37: writeln('Error detected reading datafile. Risk locus ',par2:3,') exceeds nlocus');
      38: writeln('Error detected reading datafile. Illegal value for risk locus ',par2:3,')');
      39: writeln('Error detected reading datafile. Mutation locus ',par2:3,') exceeds nlocus');
      40: writeln('Error detected reading datafile. Illegal value for mutation locus ',par2:3,')');
        41: writeln('Error detected reading datafile. Linkage disequilbirium is not allowed with this program');
        42: writeln('Locus ',par1:5,' in lod score list exceeds nlocus ',par2:5);
        43: writeln('Illegal locus number ',par1:5,' in lod score list');
      44: writeln('Error detected reading pedigree record ',par1:2,'. One 0 allele');
    END;
  END;


  procedure inputwarning(nwarning,par1,par2 : INTEGER);
  begin
    writeln('Warning number from procedure inputdata');
    CASE nwarning OF
      0: writeln('Illegal sex difference parameter ',par1:2,' Parameter should be 0, 1, or 2');
      1: writeln('Illegal interference parameter ',par1:2,' Lack of interference assumed');
      2: writeln('Illegal sex difference parameter ',par1:2,' Parameter must be 0 with sex-linked data');
      3: writeln('Non-standard affection status',par2:4,' interpreted as normal in pedigree record',par1:5);
    END;
  END;

  procedure writespeed;

  VAR
    i,j,a,b:INTEGER;

  begin                         {writespeed}
    FOR i:=1 TO totperson DO
      WITH person[i]^ DO
        IF unknown AND (foff<>NIL)
        THEN WITH store^ DO
          begin
            writeln(speedfile,'id',nseq:7);
            FOR j:=1 TO nsystem DO
              FOR a:=1 TO thislocus[j]^.nallele DO
                FOR b:=1 TO thislocus[j]^.nallele DO
                  IF possible[j,a,b]
                  THEN writeln(speedfile,j:3,a:3,b:3)
          END
  END;                          {writespeed}


  procedure writeped;

  VAR
    i,j,k,a,b:INTEGER;

  begin                         {writeped}
    FOR i:=1 TO totperson DO
      WITH person[i]^ DO
        begin
          WRITE(ipedfile,oldped:7,id:5,paid:5,maid:5,offid:5,npaid:5);
          WRITE(ipedfile,nmaid:5,sex:2,profield:2,' ');
          FOR j:=1 TO nsystem DO
            begin
              WITH phen[j]^ DO
                WITH thislocus[j]^ DO
                  IF which=binary
                  THEN IF format=2
                    THEN FOR k:=1 TO nfactor DO
                      IF k IN phenf
                      THEN WRITE(ipedfile,' 1')
                      ELSE WRITE(ipedfile,' 0')
                    ELSE begin
                      a:=0;
                      b:=0;
                      FOR k:=1 TO nallele DO
                        IF k IN phenf
                        THEN IF a=0
                          THEN a:=k
                          ELSE b:=k;
                      IF b=0 THEN b:=a;
                      WRITE(ipedfile,a:3,b:3)
                      END
                  ELSE IF which=quantitative
                  THEN IF (NOT sexlink) OR (NOT male)
                    THEN FOR k:=1 TO ntrait DO
                      WRITE(ipedfile,' ',x[k]:9:4)
                    ELSE FOR k:=1 TO ntrait DO
                      WRITE(ipedfile,' ',aff:9)
                  ELSE begin
                    WRITE(ipedfile,aff:2);
                    IF nclass<>1
                    THEN WRITE(ipedfile,liability:3)
                    END;
              IF j<>nsystem
              THEN WRITE(ipedfile,' ')
            END;
          writeln(ipedfile);
        END;
  END;                          {writeped}


  procedure infer;

  VAR
    i,j,k,l,kposs,lposs,count,pacount,macount:INTEGER;
    someknown:BOOLEAN;

  begin                         {infer}
    FOR i:=1 TO totperson DO
      IF person[i]^.unknown
      THEN WITH person[i]^ DO
        WITH store^ DO
          begin
            FOR j:=1 TO nsystem DO
              IF thislocus[j]^.which=binary
              THEN IF phen[j]^.phenf=[]
                THEN WITH thislocus[j]^ DO
                  begin
                    count:=0;
                    FOR k:=1 TO nallele DO
                      FOR l:=k TO nallele DO
                        IF possible[j,k,l]
                        THEN begin
                          kposs:=k;
                          lposs:=l;
                          count:=count+1;
                          END;
                    IF count=1
                    THEN begin
                      IF sexlink AND male
                      THEN phen[j]^.phenf:=allele[lposs]
                      ELSE phen[j]^.phenf:=allele[kposs]+allele[lposs]
                      END
                  END;
            count:=0;
            FOR j:=1 TO nsystem DO
              IF (thislocus[j]^.which<>binary)
              THEN count:=count+1
              ELSE IF phen[j]^.phenf=[]
              THEN count:=count+1;
            unknown:=count<>0
          END;
    {Infer children when parents are homozygotes}
    FOR i:=1 TO totperson DO
      IF person[i]^.foff=NIL
      THEN WITH person[i]^ DO
        FOR j:=1 TO nsystem DO
          WITH thislocus[j]^ DO
            IF phen[j]^.which=binary
            THEN IF phen[j]^.phenf=[]
              THEN IF pa<>NIL
                THEN begin
                  pacount:=0;
                  macount:=0;
                  FOR k:=1 TO thislocus[j]^.nallele DO
                    IF allele[k]<=pa^.phen[j]^.phenf
                    THEN begin
                      kposs:=k;
                      pacount:=pacount+1;
                      END;
                  FOR l:=1 TO thislocus[j]^.nallele DO
                    IF allele[l]<=ma^.phen[j]^.phenf
                    THEN begin
                      lposs:=l;
                      macount:=macount+1;
                      END;
                  IF ((macount=1) AND (pacount=1)) AND NOT (male AND sexlink)
                  THEN begin
                    phen[j]^.phenf:=allele[kposs]+allele[lposs];
                    END
                  ELSE IF (macount=1) AND (male AND sexlink)
                  THEN begin
                    phen[j]^.phenf:=allele[lposs];
                    END
                  END;
    {Replace by homozygotes if all unknown in a pedigree}
    FOR j:=1 TO nsystem DO
      WITH thislocus[j]^ DO
        IF which=binary
        THEN begin
          someknown:=FALSE;
          FOR i:=1 TO totperson DO
            IF person[i]^.phen[j]^.phenf<>[]
            THEN someknown:=TRUE;
          IF NOT someknown
          THEN FOR i:=1 TO totperson DO
            person[i]^.phen[j]^.phenf:=allele[1];
          END;
  END;                          {infer}


  procedure getunknown;

  VAR
    i,j,n,ahap,bhap:INTEGER;

  begin                         {getunknown}
    FOR i:=1 TO totperson DO
      person[i]^.unknown:=FALSE;
    FOR i:=1 TO totperson DO
      FOR j:=1 TO nsystem DO
        person[i]^.thisunknown[j]:=FALSE;
    FOR i:=1 TO totperson DO
      WITH person[i]^ DO
        FOR j:=1 TO nsystem DO
          begin
            IF thislocus[j]^.which=binary
            THEN IF phen[j]^.phenf=[]
              THEN thisunknown[j]:=TRUE
              ELSE IF thislocus[j]^.which=quantitative
              THEN IF phen[j]^.x[1]=missval
                THEN thisunknown[j]:=TRUE
                ELSE IF phen[j]^.aff=missaff
                THEN thisunknown[j]:=TRUE;
            IF thisunknown[j]
            THEN unknown:=TRUE
          END;
    FOR i:=1 TO totperson DO
      WITH person[i]^ DO
        IF unknown
        THEN begin
          NEW(store);
          WITH store^ DO
            FOR n:=1 TO nsystem DO
              WITH thislocus[n]^ DO
                FOR ahap:=1 TO nallele DO
                  FOR bhap:=1 TO nallele DO
                    possible[n,ahap,bhap]:=TRUE
          END
  END;                          {getunknown}


  procedure getlocation(thislocus:locuspoint);

  VAR
    ahap,bhap,here:INTEGER;

  begin                         {getlocation}
    here:=0;
    WITH thislocus^ DO
      FOR ahap:=1 TO nallele DO
        FOR bhap:=ahap TO nallele DO
          begin
            here:=here+1;
            genenumber[ahap,bhap]:=here;
            genenumber[bhap,ahap]:=here;
            seghap[here,a]:=ahap;
            seghap[here,b]:=bhap
          END
  END;                          {getlocation}


  procedure readped;

  VAR
    i,newid,thisone,thisped,tempid:INTEGER;


    procedure getphenotype(VAR p : ind);
    VAR
      thisread,system : INTEGER;

      procedure readbin(VAR phen : phenpoint; thislocus : locuspoint);
      VAR
        i,j : INTEGER;
      begin
        WITH thislocus^ DO
          WITH phen^ DO
            begin
              which:=binary;
              phenf:=[];
              FOR i:=1 TO nfactor DO
                begin
                  read(pedfile,j);
                  IF (j<>0) AND (j<>1) THEN inputerror(14,p^.id,j);
                  IF j=1
                  THEN phenf:=phenf+[i];
                END;
            END;
      END;


        PROCEDURE readnumber(VAR phen : phenpoint);
        VAR
          j,k : integer;
        BEGIN
          with phen^ DO
            BEGIN
              which:=binary;
              phenf:=[];
                  read(pedfile,j,k);
                  IF j>maxall THEN
                    inputerror(16,p^.id,j);
                  IF j<0 THEN
                    inputerror(17,p^.id,j);
                  IF k>maxall THEN
                    inputerror(16,p^.id,k);
                  IF k<0 THEN
                    inputerror(17,p^.id,k);
		  if ((j=0) or (k=0)) and (j<>k) then
		  inputerror(44,p^.id,j)
		  else
		  begin
                    IF j<>0 THEN phenf:=phenf+[j];
                    IF k<>0 THEN phenf:=phenf+[k];
                  end;
                end;
        end;

      procedure readaff(VAR phen : phenpoint; thislocus : locuspoint);
      VAR
        thisval : INTEGER;

      begin
        WITH phen^ DO
          begin
            which:=affection;
            read(pedfile,thisval);
            IF thisval=missaff
            THEN
              aff:=0
            ELSE
              IF thisval=affval
              THEN aff:=2
              ELSE
                begin
                  IF thisval<>1 THEN inputwarning(3,p^.id,thisval);
                  aff:=1;
                END;
            IF thislocus^.nclass=1
            THEN liability:=1
            ELSE read(pedfile,liability);
            IF liability>thislocus^.nclass THEN inputerror(26,p^.id,liability);
            IF liability<=0 THEN inputerror(27,p^.id,liability);
          END;
      END;

      procedure readquan(VAR phen : phenpoint; thislocus : locuspoint);
      VAR
        i : INTEGER;
        xval : REAL;
      begin
        WITH phen^ DO
          IF (NOT sexlink) OR (NOT p^.male)
          THEN
            begin
              which:=quantitative;
              WITH thislocus^ DO
                begin
                  FOR i:=1 TO ntrait DO read(pedfile,x[i]);
                  missing:=TRUE;
                  FOR i:=1 TO ntrait DO
                    IF x[i]<>missval
                    THEN missing:=FALSE;
                END;
            END
          ELSE
            begin
              which:=affection;
              read(pedfile,xval);
              WITH thislocus^ DO
                begin
                  IF xval=missval
                  THEN
                    aff:=missaff
                  ELSE
                    IF xval=affall
                    THEN aff:=affall
                    ELSE aff:=-11;
                  liability:=1;
                  FOR i:=2 TO ntrait DO read(pedfile,xval);
                END;
            END;
      END;

    begin
      WITH p^ DO
        FOR thisread:=1 TO nsystem DO
          begin
            system:=thisread;
            phen[system]:=NIL;
            NEW(phen[system]);
            CASE thislocus[system]^.which OF
              quantitative: readquan(phen[system],thislocus[system]);
              affection: readaff(phen[system],thislocus[system]);
              binary: IF thislocus[system]^.format=3 THEN
                readnumber(phen[system])
                ELSE readbin(phen[system],thislocus[system]);
            END;
          END;
    END;                        {getphenotype}


    procedure getind(VAR id,seq:INTEGER);

    begin                       {getind}
      id:=0;
      read(pedfile,seq);
      IF seq<>0
      THEN begin
        id:=seq;
        IF id>maxind THEN
          inputerror(13,id,id);
        IF person[id]=NIL
        THEN begin
          NEW(person[id]);
          WITH person[id]^ DO
            begin
              NEW(gen);
              nseq:=seq+nsequence
            END;
          END;
        END;
    END;                        {getind}


    procedure multimarriage(VAR p:ind);

    VAR
      q,child:ind;

    begin                       {multimarriage}
      IF p^.foff<>NIL
      THEN WITH p^ DO
        begin
          IF male
          THEN q:=foff^.ma
          ELSE q:=foff^.pa;
          child:=foff;
          p^.multi:=FALSE;
          REPEAT
            IF male
            THEN begin
              multi:=q=child^.ma;
              child:=child^.nextpa
              END
            ELSE begin
              multi:=q=child^.pa;
              child:=child^.nextma
              END
          UNTIL (child=NIL) OR (multi)
        END
      ELSE p^.multi:=FALSE
    END;                        {multimarriage}


    procedure getrest;

    VAR
      whichchr:CHAR;

    begin                       {getrest}
      whichchr:=' ';
      WHILE NOT (EOLN(pedfile) OR (whichchr=CHR(12))) DO
        read(pedfile,whichchr);
      READLN(pedfile)
    END;                        {getrest}


  begin                         {readped}
    FOR i:=0 TO maxind
    DO person[i]:=NIL;
    totperson:=0;
    loop1:=NIL;
    loop2:=NIL;
    proband:=NIL;
    thisped:=newped;
    WHILE NOT EOF(pedfile) AND (thisped=newped) DO
      begin
        totperson:=totperson+1;
        getind(thisone,tempid);
        IF proband=NIL
        THEN proband:=person[thisone];
        WITH person[thisone]^ DO
          begin
            id:=tempid;
            oldped:=newped;
            getind(newid,paid);
            pa:=person[newid];
            getind(newid,maid);
            ma:=person[newid];
            getind(newid,offid);
            foff:=person[newid];
            getind(newid,npaid);
            nextpa:=person[newid];
            getind(newid,nmaid);
            nextma:=person[newid];
            store:=NIL;
            up:=FALSE;
            read(pedfile,sex);
            IF (sex<>1) AND (sex<>2) THEN
              inputerror(11,id,sex);
            male:=(sex=1);
            read(pedfile,profield);
            IF profield=1
            THEN proband:=person[thisone]
            ELSE IF profield=2
            THEN IF loop2=NIL
              THEN loop2:=person[thisone]
              ELSE loop1:=person[thisone]
          END;
        getphenotype(person[thisone]);
        getrest;
        IF NOT EOF(pedfile)
        THEN read(pedfile,newped)
      END;
    nsequence:=totperson+nsequence;
    IF (loop2<>NIL) AND (loop1=NIL)
    THEN loop1:=proband;
    FOR thisone:=1 TO totperson
    DO multimarriage(person[thisone])
  END;                          {readped}


  procedure readloci;

  VAR
    i,coupling,autosomal,whichtype:INTEGER;
    mu:REAL;


    procedure getlocus(system:INTEGER);


      procedure getquan(VAR locus:locuspoint);

      VAR
        i:INTEGER;

      begin                     {getquan}
        WITH locus^ DO
          begin
            READLN(datafile,ntrait);
            IF ntrait>maxtrait THEN inputerror(31,system,ntrait);
            IF ntrait<=0 THEN inputerror(32,system,nclass);
          END;
        FOR i:=1 TO 3 DO
          READLN(datafile)
      END;                      {getquan}


      procedure getpen(VAR locus:locuspoint);

      VAR
        i,j,k,l:INTEGER;

      begin                     {getpen}
        WITH locus^ DO
          begin
            READLN(datafile,nclass);
            IF nclass>maxliab THEN inputerror(28,system,nclass);
            IF nclass<=0 THEN inputerror(29,system,nclass);
            FOR l:=1 TO nclass DO
              begin
                FOR i:=1 TO nallele DO
                  FOR j:=i TO nallele DO
                    begin
                      read(datafile,pen[i,j,2,l]);
                      IF (pen[i,j,2,l]<0) OR (pen[i,j,2,l]>1.0) THEN inputerror(30,system,system);
                      pen[i,j,1,l]:=1-pen[i,j,2,l];
                      pen[i,j,0,l]:=1.0;
                      FOR k:=0 TO 2 DO
                        pen[j,i,k,l]:=pen[i,j,k,l]
                    END;
                READLN(datafile);
                FOR i:=1 TO nallele DO
                  pen[0,i,0,l]:=1.0;
                IF sexlink
                THEN begin
                  FOR i:=1 TO nallele DO
                    read(datafile,pen[0,i,2,l]);
                  IF (pen[0,j,2,l]<0) OR (pen[0,j,2,l]>1.0) THEN inputerror(30,system,system);
                  FOR i:=1 TO nallele DO
                    pen[0,i,1,l]:=1.0-pen[0,i,2,l];
                  READLN(datafile)
                  END
              END
          END
      END;                      {getpen}


      procedure getbin(VAR locus:locuspoint);

      VAR
        i,j,k:INTEGER;

      begin                     {getbin}
        WITH locus^ DO
          begin
            READLN(datafile,nfactor);
            IF nfactor>maxfact THEN inputerror(8,system,nfactor);
            IF nfactor<=0 THEN inputerror(9,system,nfactor);
            FOR i:=1 TO nallele DO
              allele[i]:=[];
            FOR i:=1 TO nallele DO
              FOR j:=1 TO nfactor DO
                begin
                  read(datafile,k);
                  IF k=1
                  THEN allele[i]:=allele[i]+[j]
                END
          END;
        READLN(datafile)
      END;                      {getbin}


      procedure getnumber(VAR locus:locuspoint);

      VAR
        i:INTEGER;

      begin                     {getnumber}
        WITH locus^ DO
          FOR i:=1 TO nallele DO
            allele[i]:=[i]
      END;                      {getnumber}


    begin                       {getlocus}
      NEW(thislocus[system]);
      WITH thislocus[system]^ DO
        begin
          read(datafile,whichtype,nallele);
          CASE whichtype OF
            0 : which:=quantitative;
            1 : which:=affection;
            2 : begin
              which:=binary;
              format:=2
              END;
            3 : begin
              which:=binary;
              format:=3
              END
          END;
          READLN(datafile);
          IF NOT disequi
          THEN READLN(datafile);
          CASE which OF
            quantitative : getquan(thislocus[system]);
            affection    : getpen(thislocus[system]);
            binary       : IF format=2
              THEN getbin(thislocus[system])
              ELSE getnumber(thislocus[system])
          END;
          IF risk AND (system=risksys)
          THEN READLN(datafile)
        END
    END;                        {getlocus}


  begin                         {readloci}
    READLN(datafile,nsystem,risksys,autosomal);
    IF nsystem>maxlocus THEN inputerror(0,nsystem,nsystem);
    IF nsystem<=0 THEN inputerror(1,nsystem,nsystem);
    risk:=(risksys<>0);
    sexlink:=(autosomal=1);
    WRITE(output,'YOU ARE USING LINKAGE (V',version:3:1,') WITH',nsystem:3,'-POINT');
    IF sexlink
    THEN writeln(output,' SEXLINKED DATA')
    ELSE writeln(output,' AUTOSOMAL DATA');
    READLN(datafile,mutsys,mu,mu,coupling);
    disequi:=(coupling=1);
    READLN(datafile);
    FOR i:=1 TO nsystem DO getlocus(i)
  END;                          {readloci}


  procedure cleanup(VAR p:ind);

  VAR
    i,j:INTEGER;

  begin                         {cleanup}
    WITH p^ DO
      IF unknown
      THEN WITH store^ DO
        WITH gen^ DO
          IF sexlink AND male
          THEN begin
            FOR i:=1 TO mgeno DO
              IF NOT genarray[i]
              THEN possible[whichsys,1,i]:=FALSE;
            FOR i:=2 TO mgeno DO
              FOR j:=1 TO mgeno DO
                possible[whichsys,i,j]:=FALSE;
            END
          ELSE FOR i:=1 TO fgeno DO
            begin
              IF NOT genarray[i]
              THEN possible[whichsys,seghap[i,a],seghap[i,b]]:=FALSE;
              possible[whichsys,seghap[i,b],seghap[i,a]]:=possible[whichsys,seghap[i,a],seghap[i,b]];
            END
  END;                          {cleanup}


  procedure getgene(system:INTEGER;p:ind;phen:indphen);

  VAR
    here,i,j:INTEGER;

  begin                         {getgene}
    here:=0;
    WITH thislocus[system]^ DO
      IF sexlink AND p^.male
      THEN FOR i:=1 TO nallele DO
        begin
          here:=here+1;
          CASE which OF
            quantitative : WITH phen[system]^ DO
              IF i=affall
              THEN p^.gen^.genarray[here]:=(aff=affall) OR (aff=missaff)
            ELSE p^.gen^.genarray[here]:=(aff<>affall) OR (aff=missaff);
            affection    : WITH phen[system]^ DO p^.gen^.genarray[here]:=(pen[0,i,aff,liability]>0.0);
            binary       : WITH phen[system]^ DO p^.gen^.genarray[here]:=(phenf=allele[i]) OR (phenf=[])
          END;
        END
      ELSE FOR i:=1 TO nallele DO
        FOR j:=i TO nallele DO
          begin
            here:=here+1;
            CASE which OF
              quantitative : p^.gen^.genarray[here]:=TRUE;
              affection    : WITH phen[system]^ DO p^.gen^.genarray[here]:=(pen[i,j,aff,liability]>0.0);
              binary       : WITH phen[system]^ DO p^.gen^.genarray[here]:=(phenf=(allele[i]+allele[j])) OR (phenf=[])
            END
          END
  END;                          {getgene}


  procedure likelihood(VAR proband:ind);

  VAR
    i,j:INTEGER;


    procedure seg(VAR p,q,r:ind;peel:direction);

    VAR
      child,father,mother:ind;
      firsthap,secondhap:subhap;
      nfirst,nsecond:INTEGER;


      FUNCTION segfun(VAR child:ind;first,second:INTEGER):BOOLEAN;

      VAR
        temp:BOOLEAN;
        thishap,thathap:haplotype;

      begin                     {segfun}
        temp:=FALSE;
        WITH child^.gen^ DO
          IF NOT sexlink
          THEN FOR thishap:=a TO b DO
            FOR thathap:=a TO b DO
              temp:=temp OR genarray[genenumber[secondhap[thishap],firsthap[thathap]]]
          ELSE IF child^.male
          THEN IF p^.male
            THEN FOR thathap:=a TO b DO
              temp:=temp OR genarray[secondhap[thathap]]
            ELSE FOR thathap:=a TO b DO
              temp:=temp OR genarray[firsthap[thathap]]
          ELSE IF p^.male
          THEN FOR thathap:=a TO b DO
            temp:=temp OR genarray[genenumber[secondhap[thathap],first]]
          ELSE FOR thathap:=a TO b DO
            temp:=temp OR genarray[genenumber[firsthap[thathap],second]];
        segfun:=temp
      END;                      {segfun}


      procedure segup;

      VAR
        first,second:INTEGER;
        segval,val:BOOLEAN;

      begin                     {segup}
        IF p^.male
        THEN begin
          nfirst:=mgeno;
          nsecond:=fgeno
          END
        ELSE begin
          nfirst:=fgeno;
          nsecond:=mgeno
          END;
        WITH p^.gen^ DO
          FOR first:=1 TO nfirst DO
            IF genarray[first]
            THEN begin
              segval:=FALSE;
              firsthap:=seghap[first];
              WITH q^.gen^ DO
                FOR second:=1 TO nsecond DO
                  IF genarray[second]
                  THEN begin
                    secondhap:=seghap[second];
                    val:=TRUE;
                    child:=father^.foff;
                    REPEAT
                      IF mother=child^.ma
                      THEN val:=segfun(child,first,second);
                      child:=child^.nextpa
                    UNTIL (NOT val) OR (child=NIL);
                    segval:=val OR segval
                    END;
              genarray[first]:=segval
              END;
        cleanup(q);
        child:=father^.foff;
        REPEAT
          IF child^.ma=mother
          THEN cleanup(child);
          child:=child^.nextpa
        UNTIL child=NIL;
      END;                      {segup}


      procedure segdown;

      VAR
        first,second,here:INTEGER;
        val:BOOLEAN;
        thishap,thathap:haplotype;

      begin                     {segdown}
        FOR first:=1 TO fgeno DO gene[first]:=FALSE;
        WITH p^.gen^ DO
          FOR first:=1 TO mgeno DO
            IF genarray[first]
            THEN begin
              firsthap:=seghap[first];
              WITH q^.gen^ DO
                FOR second:=1 TO fgeno DO
                  IF genarray[second]
                  THEN begin
                    secondhap:=seghap[second];
                    val:=genarray[second] AND p^.gen^.genarray[first];
                    child:=father^.foff;
                    REPEAT
                      IF child^.ma=mother
                      THEN IF NOT child^.up
                        THEN val:=segfun(child,first,second);
                      child:=child^.nextpa
                    UNTIL (NOT val) OR (child=NIL);
                    IF val
                    THEN begin
                      IF NOT sexlink
                      THEN FOR thishap:=a TO b DO
                        FOR thathap:=a TO b DO
                          begin
                            here:=genenumber[secondhap[thishap],firsthap[thathap]];
                            gene[here]:=gene[here] OR val
                          END
                      ELSE IF r^.male
                      THEN FOR thathap:=a TO b DO
                        begin
                          here:=secondhap[thathap];
                          gene[here]:=gene[here] OR val
                        END
                      ELSE FOR thathap:=a TO b DO
                        begin
                          here:=genenumber[secondhap[thathap],first];
                          gene[here]:=gene[here] OR val
                        END
                      END
                    END
              END;
        WITH r^.gen^ DO
          FOR first:=1 TO fgeno DO
            genarray[first]:=genarray[first] AND gene[first];
        cleanup(p);
        cleanup(q);
        child:=father^.foff;
        REPEAT
          IF child^.ma=mother
          THEN IF NOT child^.up
            THEN cleanup(child);
          child:=child^.nextpa
        UNTIL child=NIL
      END;                      {segdown}


    begin                       {seg}
      IF p^.male
      THEN begin
        father:=p;
        mother:=q
        END
      ELSE begin
        father:=q;
        mother:=p
        END;
      IF peel=peelup
      THEN segup
      ELSE segdown
    END;                        {seg}


    procedure collapsedown(p:ind);FORWARD;


    procedure collapseup(p:ind);

    VAR
      q,child,nextchild:ind;
      down:BOOLEAN;

    begin                       {collapseup}
      p^.done:=TRUE;
      IF p^.foff<>NIL
      THEN begin
        down:=FALSE;
        child:=p^.foff;
        WHILE child<>NIL DO
          begin
            down:=FALSE;
            IF p^.male
            THEN q:=child^.ma
            ELSE q:=child^.pa;
            IF NOT q^.done
            THEN begin
              collapsedown(q);
              nextchild:=child;
              WHILE nextchild<>NIL DO
                begin
                  IF (nextchild^.pa=q) OR (nextchild^.ma=q)
                  THEN IF NOT nextchild^.up
                    THEN collapseup(nextchild)
                    ELSE down:=TRUE;
                  IF p^.male
                  THEN nextchild:=nextchild^.nextpa
                  ELSE nextchild:=nextchild^.nextma
                END;
              IF q^.multi
              THEN collapseup(q);
              IF NOT down
              THEN seg(p,q,child,peelup)
              ELSE collapsedown(p)
              END;
            IF p^.male
            THEN child:=child^.nextpa
            ELSE child:=child^.nextma
          END
        END
    END;                        {collapseup}


    procedure collapsedown;

    begin                       {collapsedown}
      IF p^.pa<>NIL
      THEN begin
        p^.up:=TRUE;
        collapseup(p^.pa);
        seg(p^.pa,p^.ma,p,peeldown)
        END
    END;                        {collapsedown}


  begin                         {likelihood}
    collapsedown(proband);
    collapseup(proband);
    IF proband^.thisunknown[whichsys]
    THEN begin
      WITH proband^ DO
        WITH store^ DO
          WITH gen^ DO
            WITH thislocus[whichsys]^ DO
              IF sexlink AND male
              THEN FOR j:=1 TO nallele DO possible[whichsys,1,j]:=genarray[j]
              ELSE FOR i:=1 TO nallele DO
                FOR j:=i TO nallele DO
                  begin
                    possible[whichsys,i,j]:=genarray[genenumber[i,j]];
                    possible[whichsys,j,i]:=possible[whichsys,i,j]
                  END;
      cleanup(proband)
      END
  END;                          {likelihood}


  procedure iterpeds;

  VAR
    i,j:INTEGER;
    compattest,compat:BOOLEAN;

  begin                         {iterpeds}
    IF (loop1=NIL) AND (loop2=NIL)
    { This means that this part of unknown is not active for pedigrees with loops! }
    THEN begin
      FOR i:=1 TO totperson DO
        getgene(whichsys,person[i],person[i]^.phen);
      compattest:=FALSE;
      compat:=FALSE;
      FOR i:=1 TO totperson DO
        IF (NOT compattest) OR (person[i]^.thisunknown[whichsys] AND compat)
        THEN begin
          FOR j:=1 TO totperson DO person[j]^.done:=FALSE;
          FOR j:=1 TO totperson DO person[j]^.up:=FALSE;
          likelihood(person[i]);
          IF NOT compattest
          THEN begin
            compattest:=TRUE;
            FOR j:=1 TO fgeno DO
              compat:=compat OR person[i]^.gen^.genarray[j];
            IF NOT compat THEN begin
              writeln(output,'ERROR: Incompatibility detected in this family for locus ',whichsys);
              END;
            END;
          END;
      FOR i:=1 TO totperson DO
        IF person[i]^.unknown
        THEN cleanup(person[i]);
      END;
  END;                          {iterpeds}


  procedure reinit;

  VAR
    i,j:INTEGER;

  begin                         {reinit}
    FOR i:=1 TO totperson DO
      FOR j:=1 TO nsystem DO
        DISPOSE(person[i]^.phen[j]);
    FOR i:=1 TO totperson DO
      IF person[i]^.store<>NIL
      THEN DISPOSE(person[i]^.store);
    FOR i:=1 TO totperson DO
      begin
        DISPOSE(person[i]^.gen);
        DISPOSE(person[i])
      END;
    FOR i:=1 TO totperson DO
      person[i]:=NIL
  END;                          {reinit}


  procedure initunknown;

  begin                         {initunknown}
    writeln('Program UNKNOWN version ',version:3:1);
    writeln('The following maximum values are in effect:');
    writeln(maxlocus:8,' loci');
    writeln(maxgeno:8,' single locus genotypes');
    writeln(maxall:8,' alleles at a single locus');
    writeln(maxind:8,' individuals in one pedigree');
    writeln(maxmarriage:8,' marriage(s) for one male');
    writeln(maxtrait:8,' quantitative factor(s) at a single locus');
    writeln(maxliab:8,' liability classes');
    writeln(maxfact:8,' binary codes at a single locus');
  END;                          {initunknown}

begin
  initunknown;
  reset(datafile,'datafile.dat');
  reset(pedfile,'pedfile.dat');
  rewrite(ipedfile,'ipedfile.dat');
  rewrite(speedfile,'speedfile.dat');
  readloci;
  CLOSE(datafile);
  nsequence:=0;
  IF NOT EOF(pedfile)
  THEN read(pedfile,newped);
  WHILE NOT EOF(pedfile) DO
    begin
      readped;
      getunknown;
      FOR whichsys:=1 TO nsystem DO
        IF mutsys<>whichsys
        THEN begin
          WITH thislocus[whichsys]^ DO
            begin
              fgeno:=nallele*(nallele+1) DIV 2;
              IF sexlink
              THEN mgeno:=nallele
              ELSE mgeno:=fgeno
            END;
          getlocation(thislocus[whichsys]);
          iterpeds
          END;
      infer;
      writeped;
      writespeed;
      reinit
    END;
  CLOSE(pedfile);
  CLOSE(speedfile);
  CLOSE(ipedfile);
END.
