program elodhet(input,output,lod,el,msim,limit);

{Computes ELOD and std.dev for given alpha (proportion of linked ped)
 and theta.  Assumes two loci only.
 Input:
   Simulated lod scores expected in file 'lodfile.dat'
   Parameter information in file 'msim.dat'
   Critical limits of lod score in file 'limit.dat'
 Output:
   Writes to file 'elodhet.dat'}

label 88;

const
 version='1.40'; {20 June 1993}
 maxped=300;     {max. no. of pedigrees per replicate}
 maxth=51;       {max. no. of theta values}

var
 numrep,numped,ii,kk,ith,nth:integer;
 xb,sd,rr,rr2,alpha,a1,ln10,invln10:real;
 theta:array[1..maxth] of real;
 sum,sum2:array[1..maxth,1..maxped] of real;
 sumped,sumped2:array[1..maxth] of real;
 limit,lod,el,msim:text;
 lim:array[1..3] of real; {critical limits of max. lod score}
 prop:array[1..3,1..maxth] of real; {proportions of replicates exceeding lim}


 procedure getparam;
 {Reads parameters from file 'msim.dat';
  relies on exact line numbers in that file}
 var
  jj:integer;
  cc:char;
 begin
  writeln('Reading MSIM.DAT');
  reset(msim,'msim.dat');
  if eof(msim) then writeln('File empty or nonexistent');
  for jj:=1 to 2 do readln(msim);
  for jj:=1 to 32 do read(msim,cc); readln(msim,numrep);
  for jj:=1 to 51 do read(msim,cc); readln(msim,a1);
  alpha:=1.0-a1;
  for jj:=1 to 2 do readln(msim);
  for jj:=1 to 22 do read(msim,cc); readln(msim,numped);
 end; {getparam}


begin {elodhet}
 writeln;
 writeln('Program ELODHET  version ',version);
 writeln;
 getparam;
 if numped>maxped then begin
  writeln;
  writeln('ERROR: Number of pedigrees (',numped:1,') larger than constant');
  writeln('       MAXPED=',maxped:1,'.  Recompile and rerun program.');
  goto 88;
 end;
 writeln('Opening LODFILE.DAT');
 reset(lod,'lodfile.dat');
 if eof(lod) then writeln('File empty or nonexistent');
 ln10:=10.0;
 ln10:=LN(ln10);
 invln10:=1.0/ln10;

 {Find theta values}
 readln(lod,theta[1]); {first theta value}
 nth:=1;
 repeat
  nth:=nth+1;
  if nth>maxth then begin
   writeln('ERROR: Number of theta values exceeds MAXTH=',maxth:1);
   goto 88;
  end;
  readln(lod,rr);
  if rr<>theta[1] then theta[nth]:=rr;
 until (rr=theta[1]);
 nth:=nth-1;

 {Calculate ELODs and standard deviations}
 writeln('ANALYSES FOR FIXED THETA VALUES:');
 reset(lod);  {set pointer back to beginning of file}
 for ith:=1 to nth do
 for ii:=1 to numped do begin
  sum[ith,ii]:=0;
  sum2[ith,ii]:=0;
 end;
 for ith:=1 to nth do begin
  sumped[ith]:=0;
  sumped2[ith]:=0;
 end;
 writeln('  Opening LIMIT.DAT file');
 reset(limit,'limit.dat');
 if eof(limit) then writeln('File empty or nonexistent');
 for ii:=1 to 3 do begin
  read(limit,lim[ii]);
  for ith:=1 to nth do prop[ii,ith]:=0;
 end;

 writeln('  Reading lod scores from LODFILE.DAT...');
 for ii:=1 to numrep do
 for ith:=1 to nth do begin
  read(lod,rr);   {theta value, not used here}
  rr2:=0;  {lods under het. for all pedigrees}
  for kk:=1 to numped do begin
   read(lod,rr);  { lod score }
   rr:=EXP(ln10*rr);  { likelihood, L(theta) }
   rr:=invln10*LN(alpha*rr+a1);  { log10[L(alpha,theta)] for which mean }
   sum[ith,kk]:=sum[ith,kk]+rr;  { and std. dev. will be calculated }
   sum2[ith,kk]:=sum2[ith,kk]+SQR(rr);
   rr2:=rr2+rr;
  end; {for kk}
  sumped[ith]:=sumped[ith]+rr2;
  sumped2[ith]:=sumped2[ith]+SQR(rr2);
  for kk:=1 to 3 do if rr2>=lim[kk] then prop[kk,ith]:=prop[kk,ith]+1.0;
  readln(lod); {This line not strictly needed; inserted as a check}
 end; {for ii, ith}

 writeln('  Computing and outputting summary statistics...');
 rewrite(el,'elodhet.dat');
 writeln(el);
 writeln(el,'Program ELODHET  version ',version);
 writeln(el);
 writeln(el,'Alpha =':10,alpha:10:4);
 writeln(el,numrep:8,' replicates');
 writeln(el,numped:8,' pedigree(s) per replicate');
 writeln(el);
 writeln(el,'RESULTS FOR FIXED THETA VALUES');

 for ith:=1 to nth do begin
  writeln(el);
  writeln(el,'  Theta =',theta[ith]:8:4);
  writeln(el,'Prob. Zmax exceeds':56);
  writeln(el,'Ped no':10,'E[Z]':10,'std.dev.':10,
    lim[1]:10:4,lim[2]:10:4,lim[3]:10:4);
  for kk:=1 to numped do begin
  {Mean rr and standard deviation rr2 for each pedigree}
   rr:=sum[ith,kk]/numrep;
   if numrep>1
    then rr2:=(sum2[ith,kk]-numrep*SQR(rr))/(numrep-1.0)
    else rr2:=0;
   if (rr2>0) then rr2:=SQRT(rr2) else rr2:=0;
   writeln(el,kk:10,rr:10:4,rr2:10:4);
  end;
  {Mean and std.dev. for all pedigrees jointly}
  rr:=sumped[ith]/numrep;
  if numrep>1
   then rr2:=(sumped2[ith]-numrep*SQR(rr))/(numrep-1.0)
   else rr2:=0;
  if (rr2>0) then rr2:=SQRT(rr2) else rr2:=0;
  write(el,'All ped''s':10,rr:10:4,rr2:10:4);
  for kk:=1 to 3 do write(el,prop[kk,ith]/numrep:10:4);  writeln(el);
 end; {for ith}

 writeln('LIKELIHOODS MAXIMIZED OVER ALPHA AND THETA:');
 for ii:=1 to 3 do prop[ii,1]:=0;
 writeln('  Reading lod scores from LODFILE.DAT...');
 reset(lod);  {set pointer back to beginning of file}
 xb:=0;  sd:=0;
 for ii:=1 to numrep do begin
  for ith:=1 to nth do begin
   read(lod,rr);   {theta value, not used here}
   for kk:=1 to numped do begin
    read(lod,rr);  { lod score }
    sum[ith,kk]:=EXP(ln10*rr);  { likelihood, L(theta) }
   end; {for kk}
   readln(lod); {This line not strictly needed; inserted as a check}
  end; {for ith}

  {Find max. L(alpha,theta) in each replicate, alpha step size = 0.025}
  alpha:=1.0;
  rr:=0;  {max Ln L(alpha,theta) for all families in given replic}
  repeat
   a1:=1.0-alpha;
   for ith:=1 to nth do begin
    {Given alpha and theta:}
    rr2:=0;
    for kk:=1 to numped do rr2:=rr2+LN(alpha*sum[ith,kk]+a1);
    if (rr2>rr) then rr:=rr2;
   end; {for ith}
   alpha:=alpha-0.025
  until (alpha<0.01);
  rr:=rr*invln10;  {turns max. Ln L into max. log10 L}
  xb:=xb+rr;
  sd:=sd+SQR(rr);
  for kk:=1 to 3 do if rr>=lim[kk] then prop[kk,1]:=prop[kk,1]+1.0;
 end; {for ii}

 writeln('  Computing and outputting summary statistics...');
 writeln(el);  writeln(el);
 writeln(el,'LIKELIHOOD MAXIMIZED OVER ALPHA AND THETA');
 writeln(el);
 {Mean and std.dev. for all pedigrees jointly}
 xb:=xb/numrep;
 if numrep>1
  then sd:=(sd-numrep*SQR(xb))/(numrep-1.0)
  else sd:=0;
 if (sd>0) then sd:=SQRT(sd) else sd:=0;
 writeln(el,'E(max.lod)','Prob. Zmax exceeds':46);
 writeln(el,'Mean':20,'Std.dev':10,lim[1]:10:4,lim[2]:10:4,lim[3]:10:4);
 write(el,xb:20:4,sd:10:4);
 for ii:=1 to 3 do write(el,prop[ii,1]/numrep:10:4); writeln(el);

 writeln('Output is in file ELODHET.DAT');
 if not eof(lod) then begin
  writeln;
  writeln('WARNING:  Results unreliable, LODFILE.DAT file only partially read.');
 end;
 88:
 { close(el); }
end.
