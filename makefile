# This makefile assumes the shell is Bash.
#########################################################################
# To get to compile with X, you might need to supplement this path.
LIBPATH= -L. -L./accis/ -L/usr/lib/mesa
#########################################################################
# Test whether X libraries are found. Null => yes.
 TESTGL:=$(shell ld  $(LIBPATH) -lGLU -lGL -o /dev/null 2>&1 | grep GL)
 TESTX11:=$(shell ld $(LIBPATH) -lXt -lX11 -o /dev/null 2>&1 | grep X)
##########################################################################
########################################################################
# Decide accis driver choice. Alternatives are vec4014 vecx or vecglx. 
# Automatic choice can be overriden by commandline option e.g. 
# make VECX=vec4014
# or in accis if libraries are unfound. 
# But we need to be able to tell which ACCISDRV to use based upon 
# accis configuration. So we decide here. 
ifeq ("$(VECX)","vec4014")
 ACCISDRV=accis
 LIBRARIES = $(LIBPATH) -lsceptic -l$(ACCISDRV)
 else
 ifeq ("$(VECX)","vecx")
   ifeq ("$(TESTX11)","")
     ACCISDRV=accisX
     LIBRARIES = $(LIBPATH) -lsceptic -l$(ACCISDRV) -lXt -lX11 $(GLULIBS)
   else
# Wanted vecx but could not have it:
     ACCISDRV=accis
     LIBRARIES = $(LIBPATH) -lsceptic -l$(ACCISDRV)
     VECX:=vec4014
   endif
 else
# VECX not vec4014 or vecx
   ifeq ("$(TESTX11)","")
    ACCISDRV=accisX
    GLULIBS=
    LIBRARIES = $(LIBPATH) -lsceptic -l$(ACCISDRV) -lXt -lX11 $(GLULIBS)
    VECX=vecx
   else
    ifeq ("$(TESTGL)","")
     ACCISDRV=accisX
     GLULIBS= -lGL -lGLU
     LIBRARIES = $(LIBPATH) -lsceptic -l$(ACCISDRV) -lXt -lX11 $(GLULIBS)
     VECX=vecglx
    else
     ACCISDRV=accis
     LIBRARIES = $(LIBPATH) -lsceptic -l$(ACCISDRV)
     VECX:=vec4014
    endif
   endif
 endif
endif
ACCISLIB=./accis/lib$(ACCISDRV).a
# For submakes:
export LIBRARIES
export VECX
##########################################################################
# Decide which compiler to use.
ifeq ("$(G77)","")
# I don't know why this has to be overridden. 
# But within this section of code G77 is not set without an override.
	override G77=$(shell cat compiler 2>/dev/null)
	ifeq ("$(G77)","")
# Default compiler. Ought to be used if a strange make target is used 
# on the very first call.
# After that, compiler ought to be set on disk and used.
		override G77=mpif77 -f77=g77
	endif
endif
# In g77 -Wno-globals silences spurious type messages on reduce.f
# This is unrecognized by gfortan. For which no-unused is better.
NGW=-Wno-unused
ifeq ("$(G77)","mpif77 -f77=g77")	
  NGW=-Wno-globals
endif
# export this so it is inherited by sub-makes.
export G77
###################################################################
NOWARN=-Wno-unused-label -Wno-unused-dummy-argument
XLIB=$LIBPATH
MPIexecutable=scepticmpi
###################################################################
#__________________________________________________________________________
#
#     This code is copyright (c)
#              Ian H Hutchinson    hutch@psfc.mit.edu.
#              Leonardo Patacchini patacchi@mit.edu
#
#     It may be used freely with the stipulation that any scientific or
#     scholarly publication concerning work that uses the code must give
#     an acknowledgement referring to the relevant papers
#
#     I.H. Hutchinson, Plasma Physics and Controlled Fusion, vol 44, p
#     1953 (2002), vol 45, p 1477 (2003).
#
#     L. Patacchini and I.H. Hutchinson, Plasma Physics and Controlled
#     Fusion, vol 49, p1193 (2007), vol 49, p 1719 (2007).
#
#     I.H. Hutchinson and L. Patacchini, Physics of Plasmas, vol 14,
#     p013505 (2007)
#
#     The code may not be redistributed except in its original package.
#
#     No warranty, explicit or implied, is given. If you choose to build
#     or run the code, you do so at your own risk.
#___________________________________________________________________________
# Universal Makefile for sceptic

COMPILE-SWITCHES =-Wall $(NOWARN) -O2  -I.
# For debugging.
#  -g  -ffortran-bounds-check
# For profiling
#COMPILE-SWITCHES = -Wall -O2 -pg

REINJECT=fvinject.o orbitinject.o extint.o maxreinject.o ogeninject.o reindiag.o

MPICOMPILE-SWITCHES = -DMPI $(COMPILE-SWITCHES)

OBJECTS = initiate.o advancing.o randc.o randf.o diags.o outputs.o	\
outputlive.o chargefield.o $(REINJECT) damp.o stringsnames.o		\
rhoinfcalc.o shielding.o

OBJECTSO = initiate.o advancingo.o randc.o randf.o diags.o outputs.o	\
chargefield.o $(REINJECT) damp.o stringsnames.o rhoinfcalc.o		\
shielding.o

MPIOBJECTS=mpibbdy.o sor2dmpi.o shielding_par.o 

COMMONS=piccom.f colncom.f distcom.f sorcom.f
###########################################################################
# So make does not do multiple tries. The default target is makefile first.
all : makefile compiler libsceptic.a sceptic scepticmpi 
	./sceptic -s5 -nr20 -nt20 -ni100000 -v1. -f -g 
	echo $(COMPCOM)
# By default run a test case with no graphics.

libsceptic.a : compiler makefile $(OBJECTS) $(MPIOBJECTS) $(COMMONS)
	rm -f libsceptic.a
	@ar -rs libsceptic.a $(OBJECTS) $(MPIOBJECTS) $(UTILITIES)

scepticmpi : sceptic.F  $(COMMONS) $(ACCISLIB) $(OBJECTS) $(MPIOBJECTS) makefile
	$(G77) $(MPICOMPILE-SWITCHES) ${NGW} -o scepticmpi  sceptic.F $(LIBRARIES)

$(ACCISLIB) : ./accis/*.f
#	make VECX=vecx -C accis
	make -C accis

orbitint : orbitint.f coulflux.o $(OBJECTS) $(ACCISLIB) makefile
	$(G77) $(COMPILE-SWITCHES) -o orbitint orbitint.f $(OBJECTS) coulflux.o $(LIBRARIES)

tools : makefile compiler libsceptic.a
	make -C tools clean
	make -C tools

tools/postproc : makefile compiler libsceptic.a

coulflux.o : tools/coulflux.f
	$(G77) -c $(COMPILE-SWITCHES) tools/coulflux.f

fvinjecttest : fvinjecttest.F makefile fvinject.o reinject.o initiate.o advancing.o chargefield.o randf.o reindiag.o fvcom.f
	$(G77)  -o fvinjecttest $(COMPILE-SWITCHES) fvinjecttest.F fvinject.o reinject.o initiate.o advancing.o chargefield.o randf.o reindiag.o $(LIBRARIES)

fvinject.o : fvinject.f fvcom.f $(COMMONS)
	$(G77) -c $(COMPILE-SWITCHES) fvinject.f

sceptic.tar.gz : makefile $(ACCISLIB) sceptic $(MPIexecutable)
	make -C accis mproper
	make -C tools clean
	make makefile
	make clean
	./copyattach.sh
	tar chzf sceptic.tar.gz  --exclude *.tar.gz -C .. sceptic
	./copyremove.sh

tars : sceptic.tar.gz


# Configure compiler.
compiler : makefile
	@echo -n Compiler tests.
	@if which mpif77 >/dev/null;\
 then echo -n " MPI system. ";\
  if which g77 >/dev/null ;\
  then  echo -n "Force g77. ";GHERE="mpif77 -f77=g77";\
  else GHERE=mpif77 ; fi\
 else echo -n "Not MPI System. ";\
  if which g77 >/dev/null ;\
  then GHERE="g77";\
  else GHERE="f77";fi\
 fi;\
 echo "Chosen G77="$${GHERE}; G77=$${GHERE}; echo $${G77} > compiler;\
# To obtain this information, one has to make a second time.
	@echo 'To set compiler manually do: make cleanall; echo "compiler command" >compiler'
	@echo 'Then might need to use:  make NOWARN='
	@echo "*********** Remaking with chosen G77 ****************"
	@export MAKEFLAGS=; make libsceptic.a

clean :
	make -C accis mproper
	rm -f *.o
	rm -f *.ps
	rm -f *.orb
	rm -f *.html
	rm -f Orbits.txt
	rm -f *~
	rm -f *.liv
	rm -f sceptic.tar.gz
	rm -f libsceptic.a
	make -C tools clean

cleanall :
	make clean
	rm -f sceptic scepticmpi sceptichdf scepticmpihdf fvinjecttest fvinittest
	rm -f *.dat *.dst
	rm -f *.frc
	rm -f compiler

ftnchek :
	ftnchek -nocheck -nof77 -calltree=text,no-sort -mkhtml -quiet -brief sceptic.F *.f

help :
	@echo 'To set compiler manually do: echo "<compilercommand>" >compiler'

#######################################################################
# Allow optional hdf-enabled version if have libraries

# To figure out what to use for the hdf includes and libraries
# run the h5fc script with -show (usr/local/hdf5/bin/h5fc)
HDFINCLUDE = -I/usr/local/hdf5/include
HDFLIBRARIES = -L/usr/local/hdf5/lib -lhdf5hl_fortran -lhdf5_hl -lhdf5_fortran -lhdf5 -lz -lm -Wl,-rpath -Wl,/usr/local/hdf5/lib

#Defaults compiler (mpif90 compiler)
ifeq ("$(G90)","")
	G90=mpif90
endif

outputhdf.o : outputhdf.f piccom.f colncom.f
	$(G90) -c -Wall -O2 -I. $(HDFINCLUDE) outputhdf.f

sceptichdf : sceptic.F  $(COMMONS) $(ACCISLIB) $(OBJECTS) outputhdf.o makefile
	$(G77) $(COMPILE-SWITCHES) -DHDF outputhdf.o $(HDFINCLUDE) $(HDFLIBRARIES) ${NGW} -o sceptichdf sceptic.F $(LIBRARIES)

scepticmpihdf : sceptic.F  $(COMMONS) $(ACCISLIB) $(OBJECTS) $(MPIOBJECTS) outputhdf.o makefile
	$(G77) $(MPICOMPILE-SWITCHES) -DHDF outputhdf.o $(HDFINCLUDE) $(HDFLIBRARIES) ${NGW} -o scepticmpihdf sceptic.F $(LIBRARIES)

#######################################################################
#pattern rules need to be at end not to override specific rules
%.o : %.f $(COMMONS) fvcom.f makefile;
	$(G77) -c $(COMPILE-SWITCHES) $*.f

%.o : %.F $(COMMONS) makefile;
	$(G77) -c $(COMPILE-SWITCHES) $*.F

% : %.f makefile $(ACCISLIB)
	$(G77)  -o $* $(COMPILE-SWITCHES) $*.f  $(LIBRARIES)

% : %.F makefile $(ACCISLIB)
	$(G77)  -o $* $(COMPILE-SWITCHES) $*.F  $(LIBRARIES)
