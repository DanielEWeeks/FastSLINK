Last mod: October 28, 2010

INTRODUCTION
------------
This directory contains the fast version of the simulation program
SLINK. Read the file slink.txt for a general introduction to
the programs.

SLINK implements a simulation algorithm developed by Jurg Ott and
described in:

  1) Ott J (1989) Computer-simulation methods in human linkage
analysis.  Proc Natl Acad Sci USA 86:4175-4178

The algorithm was implemented in the original SLINK computer
program package by Weeks, Ott, and Lathrop:

  2) Weeks DE, Ott J, Lathrop GM (1990) SLINK: a general simulation program for linkage analysis. Am J Hum Genet 47:A204 (abstr)

SLINK is based on the LINKAGE programs version 4.9 (Lathrop
et al., 1984) and accepts slightly modified LINKAGE data files as
explained in the file 'slink.txt'.  The code has been updated to
be consistent with LINKAGE version 5.1.

The SLINK simulation program has been modified by Schaffer and
Weeks to use the algorithms developed by Cottingham et al:

  3) Cottingham Jr RW, Idury RM, Schaffer AA (1993) Faster sequential
genetic linkage computations. Am J Hum Genet 53:252-263

Please cite references 1-3 if you use FastSLINK.  Thank you.


Note that SLINK by itself is quite limited in terms of the number of
markers that it can handle.  If you wish to simulate large number of
markers, then the SUP program will also be necessary:

Lemire M. SUP: an extension to SLINK to allow a larger number of marker
loci to be simulated in pedigrees conditional on trait values. BMC
Genet. 2006 Jul 3;7:40. PubMed PMID: 16803631; PubMed Central PMCID:
PMC1524809.

SUP was previously available from this web site:

http://mlemire.freeshell.org/software.html




The C source files and compilation scripts are in the sub-directory
./src.	These produce fast versions of

	slink

For user convenience, the /src subdirectory also includes C code files
    isim.c 
    lsim.c
    msim.c
    unknown.c

It should be noted that to date we have not implemented the mutation
part of the code due to lack of demand.  If there is interest, please
let us know and we will work on these.

In 2007-2010, some changes were made to improve the code and add new
features. One part of the changes are based on analogous changes to
FASTLINK, implemented long ago, that reduce the usage of #defined
constants for quantities such as the number of multi-locus genotypes. One
effect of these changes is that when the true values are much less than
what the #define maximum value used to be, the program slink runs
demonstrably faster.

Please let us know if you have problems with the programs, including if
you are unhappy with the speedup and are willing to share your data to
the extent that we may be able to study the problem. Note that this does
not mean you have to tell us anything about what disease you are
studying. And of course we will respect any request for confidentiality.
We only wish to consider studying problems to see if we can find
improvements.

Also please tell us if your experience is positive. (email
weeks@pitt.edu or schaffer@helix.nih.gov)



1. SEPARATE COMPILATION. This new version has been split up into
multiple definition and code files. The code is still virtually the
same but it is organized differently. The main advantages of separate
compilation are:

a.) Under some circumstances the program can be recompiled much
    faster than before.

b.) The new arrangement is much easier for us to maintain and work
with.  Duplicate code is well-known to be very prone to bugs (usually
when only some copies are updated) and generally hard to maintain.

c.) Instead of using compilation scripts, we are now using the make
utility to compile. Most LINKAGE users have probably used this utility
since it was the mechanism for compiling previous distributions of
LINKAGE.



SOURCE CODE ORGANIZATION
------------------------

The code for the executable slink is now broken up into .h files which have just definitions
and declarations, and .c files that have just C procedures.  Each file
opens with a comment explaining what role that file plays. Here we
summarize which files are needed for which programs:

File Name              
commondefs.h             
sldefs.h                
moddefs.h               
slomoddefs.h           
slautomodified.c          
commongetvect.c         
commonnuclear.c         
iostuff.c
slgetvect.c           
slink.c               
slinputcode.c         
sloldnuclear.c        
oldsegup.c           
slsexmodified.c          
slsloautomodified.c     
slslosexmodified.c      

The distinction between faster and slower versions is explained in
the COMPILATION section below.

Most of the modified code is in the following files:

```
    moddefs.h       definitions and declarations associated
                    with the fast, but memory intensive
                    version of the new code.

    slautomodified.c    contains the fast, but memory intensive 
                        version of the new code.

    slsexmodified.c	contains the sex-linked version of
                        slautomodified.c

    slomoddefs.h 	 versions that contain the slower but
    slsloautomodified.c  space-efficient alternatives of the
    slslosexmodified.c	 above modules

    oldsegup.c      contains the old (p2c version of original
                    programs) versions of segup() and segsexup(),
                    which are needed for handling mutation data.
```

We have changed seg() in the original code so that the modified
routines are called only for mutation-less data.  For data with
mutation, the old routines are called for compatibility.  These have
not been modified yet.  Two routines with the names oldsegup() and
oldsegsexup() are called which are the same as the old segup() and
segsexup().

The code for
    isim.c
    lsim.c
    msim.c
    unknown.c
is unchanged, except that previous releases had 
PASCAL code files isim.p, lsim.p, msim.p.
The C versions were obtained by translating with p2c.
We distribute the C versions because C compilers are
more readily available and produce faster executable code than
current PASCAL compilers. 


COMPILATION
-----------

Part of the distribution is a file called Makefile, which enables you
to compile the programs. You can issue any of the following
commands for compilation:

1. make slink
2. make sloslink
3. make unknown
4. make isim
5. make lsim
6. make msim
7. make all


1 and 2 put an executable version of SLINK in the file slink.
Command 7 does all of 1, 3, 4, 5, 6 in a single command, so you
get executable files for
  slink
  unknown
  isim
  lsim
  msim 


The version of slink produced by command 1 is faster than that
produced by command 2, but may use much more memory
depending on the data and run parameters. The second version is
slower, though still much faster than early versions of slink, and uses little
memory.  If you do not have enough memory to run the faster version,
delete the faster version and make the slower one instead.  The new
compilation structure makes this recompilation much faster than it was
before because most of the code does not have to be
recompiled.

Please read the remark below about the constant MAXALL before
trying to compile.

The Makefile we are distributing uses the gcc compiler distributed by
the Free Software Foundation. If gcc is not available to you or you
would like to use the cc compiler instead, you can do this in either
of two ways:

1. Change the third line in Makefile from

	CC	= gcc 
   to

	CC	= cc

2. When you issue the make command, add the tag CC=cc; this will override
   the setting in Makefile; e.g.,

	make slink CC=cc

In general, we recommend using gcc instead of cc, because gcc produces
faster machine code.

If you are unfamiliar with these concepts and want help, see your
system administrator for information about how your particular system
is organized.



MEMORY REQUIREMENTS
-------------------

These programs can require large amounts of memory.  For instance
slink as configured in this distribution and compiled with

  make slink
  
requires less than 1 Mb.  Of course the amount of memory required is
very dependent on the number of loci and the number of alleles at each
locus.  However larger amounts of memory is not a problem to run under
Sun OS for instance, because this is a virtual memory operating system. 
Ideally one would want to run a program of this size on a machine with
32 Mb of memory, but in our experience it is possible to run on machines
with as little as 12 Mb.

Of course it is necessary to have a swap file with sufficient space to
run the OS and have enough free space to for the program.

To see how much space a program requires, use the Unix command:

	size <program name>

for instance:

	unix> /usr/bin/size slink

	text    data    bss     dec     hex
	131072  8192    830152  969416  ecac8

This value under "dec" is the decimal number of bytes for the whole
program.  So we see in this case that less than 1 Mbyte  is required.

Then compare this with the unix pstat command:

	unix> /etc/pstat -s

	2992k allocated + 688k reserved = 3680k used, 61676k available

This indicates that a total of 3680 Kbytes has been
allocated by programs in execution  on this system for swap space,
and with the current job mix, another 61.7 Mb are available.  So
in this case slink will be able to run.

To enlarge the swap space consult your local system administrator.

Alternatively, use the "slow" versions of the programs.  The term slow
is a little misleading in that these versions will still be
significantly faster than the originals.  In the case of slink, the
version compiled with

  make sloslink 
  
will be much smaller in size.  Almost any Unix system will have a swap
file large enough for this.



NEW FEATURES ADDED in Version 3.02
----------------------------------
Two new features were added in version 3.02.

The number of alleles at a marker can now be greater than 32. This
feature is used in the simulation method described in the paper:

Lemire M. SUP: an extension to SLINK to allow a larger number of marker
loci to be simulated in pedigrees conditional on trait values. BMC
Genet. 2006 Jul 3;7:40. PubMed PMID: 16803631; PubMed Central PMCID:
PMC1524809.

The article states that

"there is a limit on the number of alleles a marker can have ...
for most computers the total number of alleles a marker may have
cannot exceed 32..."

and this was true in 2006, when the article was published. 

Now the restriction on number of alleles is removed. It remains
necessary to set a maximum number of alleles in commondefs.h by changing
the code line that looks like:

```
#define maxall          129   /*MAX NUMBER OF ALLELES AT A SINGLE LOCUS*/
```

before compiling.

By calling

  slink -t
  
instead of the default

  slink
  
one can get a second output file, besides pedfile.dat which is called

  traitfile.dat

The first nine columns of pedfile.dat and traitfile.dat should be
identical. The column that shows the trait value will be identical also.
The two columns in traitfile.dat that follow the trait column do not
occur in pedfile.dat, and they show the alleles at the disease locus.

By calling 

  slink -a

one gets filled-in genotypes for all individuals, regardless of the
availability code. -a is a mnemonic for "all available".  The -a
feature is needed to simplify the usage of slink in conjunction with
SUP. In that situtation, call

  slink -a -t

or equivalently

  slink -t -a

SUP will take care of suppressing the genotypes for individuals whose
availability codes indicate that they are not available for
genotyping. 

Another possible usage of

  slink -a

is to evaluate what would be the benefits of genoptying
all individuals rather than just those currently coded as
available.

QUANTITATIVE TRAITS AND LINKAGE DISEQUILIBRIUM
------------------------------------------------

Mathieu Lemire, developer of SUP, pointed out that SLINK can handle
either quantitative traits or linkage disequilibrium between loci, but
these topics are not adequately covered in the original documentation.
Linkage disequilibrium is especially relevant to the combined usage of
SLINK and SUP.

Quantitative trait loci can be simulated.
Here is an example of what a QTL specification looks like in
simdata.dat

```
 2 0 0 5  << NO. OF LOCI, RISK LOCUS, SEXLINKED (IF 1) PROGRAM
 0 0.0 0.0 0  << MUT LOCUS, MUT RATE, HAPLOTYPE FREQUENCIES (IF 1)
  1  2
0   2  << QUANTITATIVE, NO. OF ALLELES
 0.990000 0.010000   << GENE FREQUENCIES
 1 << NO. OF TRAITS
1.57 2.10 2.10 << GENOTYPE MEANS
0.059 << GENOTYPE VARIANCE
0.029 << MULTIPLIER FOR VARIANCE IN HOMOZYGOTES
3   4  << ALLELE NUMBERS, NO. OF ALLELES
 0.250000 0.250000 0.250000 0.250000 << GENE FREQUENCIES
 0 0  << SEX DIFFERENCE, INTERFERENCE (IF 1 OR 2)
 0.01000 << RECOMBINATION VALUES
 1 0.00 0.50000 << REC VARIED, INCREMENT, FINISHING VALUE
```

Here we have one quantitative trait locus and one marker locus. The
recombination fraction (second value on the last line) is chosen as 0.00
because that is what would be used in SUP.

If one is using only SLINK (no SUP), then datain.dat looks identical to
simdata.dat except for this second value on the last line, which has to
be a non-zero increment of the recombination fraction, typically 0.05 or
0.10.

In the pedigree file (simped.dat), the trait values are specified as
single floating point numbers. Use 0.00 for an unknown value.


For linkage disequilibrium, here is a variant of the simdata.dat in the
example1 subdirectory. The first key difference is that on the second line,
the fourth and last number is 1 instead of the usual 0. The other key
difference is that haplotype frequencies are provided after both loci,
instead of allele frequencies after each locus.

```
 2 0 0 5  << NO. OF LOCI, RISK LOCUS, SEXLINKED (IF 1) PROGRAM
 0 0.0 0.0 1  << MUT LOCUS, MUT RATE, HAPLOTYPE FREQUENCIES (IF 1)
  1  2
1   2  << AFFECTION, NO. OF ALLELES
 2 << NO. OF LIABILITY CLASSES
 0.0000 0.0000 1.0000
 0.0000 1.0000 0.0000 << PENETRANCES
3   4  << ALLELE NUMBERS, NO. OF ALLELES
 0.24750000 0.2475000 0.2475000 0.24750000 0.0025 0.0025 0.0025 0.0025 << HAPLOTYPE FREQUENCIES
 0 0  << SEX DIFFERENCE, INTERFERENCE (IF 1 OR 2)
 0.01000 << RECOMBINATION VALUES
 1 0.00000 0.50000 << REC VARIED, INCREMENT, FINISHING VALUE
```

For a more in-depth explanation of how the haplotype frequencies
are coded, please see the SUP documentation.


CITING SLINK
------------
Please cite references 1-3 if you use FastSLINK in a publication:

  1) Ott J (1989) Computer-simulation methods in human linkage
analysis.  Proc Natl Acad Sci USA 86:4175-4178

  2) Weeks DE, Ott J, Lathrop GM (1990) SLINK: a general simula-
tion program for linkage analysis. Am J Hum Genet 47:A204 (abstr)

  3) Cottingham Jr RW, Idury RM, Schaffer AA (1993) Faster sequential
genetic linkage computations. Am J Hum Genet 53:252-263


CONTACTS
--------

Daniel E. Weeks:  weeks@pitt.edu

Alejandro Schäffer: alejandro.schaffer@nih.gov

For SUP, contact Mathieu Lemire: Mathieu.Lemire@oicr.on.ca
