###################################################################
#___________________________________________________________________________
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

# Universal Makefile for sceptic postprocessing

#Defaults compiler (mpif77 compiler)
ifeq ("$(G77)","")
	G77=g77
endif
#Default Xlib (32 bit)
ifeq ("$(XLIB)","")
	XLIB=/usr/X11R6/lib
endif
#Default Accis lib
ifeq ("$(ACCISLIB)","")
	ACCISLIB=-L./../accis -L/usr/lib/mesa -laccisX -lXt -lX11 -lGL -lGLU
endif

# To silence warnings when compiling with g77 uncomment the next statement:
ifeq ("$(G77)","g77")
	NOGLOBALS=-Wno-globals
endif

COMPILE-SWITCHES = -Wall -O2

LIBRARIES=-L../ -lsceptic $(ACCISLIB)

#pattern rule, compile using the external definitions of commons, no backslash.
%.o : %.f piccom.f makefile  ../libsceptic.a;
	$(G77) -c $(COMPILE-SWITCHES) $*.f

%.o : %.F piccom.f makefile ../libsceptic.a;
	$(G77) -c $(COMPILE-SWITCHES) $*.F

% : %.f ../libsceptic.a
	$(G77)  -o $* $(COMPILE-SWITCHES) $*.f  $(LIBRARIES)

% : %.F ../libsceptic.a
	$(G77)  -o $* $(COMPILE-SWITCHES) $*.F  $(LIBRARIES)

all : ../accis/libaccisX.a postprocB postproc parrot fluxVp newsumplot fcol vdist pdiag

postproc : postproc.f piccom.f ../accis/libaccisX.a coulflux.f makefile
	$(G77) $(COMPILE-SWITCHES) -g  -o postproc postproc.f coulflux.f $(LIBRARIES)

postprocnew : postprocnew.f piccom.f ../accis/libaccisX.a coulflux.f makefile
	$(G77) $(COMPILE-SWITCHES) -g  -o postprocnew postprocnew.f coulflux.f $(LIBRARIES)

../accis/libaccisX.a : ../accis/*.f
	make -C ../accis

parrot : parrot.f piccom.f ../accis/libaccisX.a makefile
	$(G77) $(COMPILE-SWITCHES) -o parrot parrot.f -laccisX $(LIBRARIES)

fluxVp : fluxVp.f ../accis/libaccisX.a coulflux.f makefile
	$(G77) $(COMPILE-SWITCHES) -g -o fluxVp fluxVp.f coulflux.f $(LIBRARIES)

postprocB : postprocB.f piccom.f ../accis/libaccisX.a coulflux.f 
	$(G77) $(COMPILE-SWITCHES) -g  -o postprocB postprocB.f  $(LIBRARIES)

fcol : fcol.f kennedyallen.o makefile
	$(G77)  -o fcol $(COMPILE-SWITCHES) fcol.f  kennedyallen.o $(LIBRARIES)

clean :

	rm -f *.o
	rm -f postproc postprocB parrot fluxVp newsumplot fcol phioutplot rhopicplot pdiag vdist
	rm -f *.dat *.frc
	rm -f *.ps *.txt
