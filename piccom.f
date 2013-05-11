
c Version 2.5; Aug 2004   
      integer npartmax,npart,nr,nth,ndim,np
c Number of particles: npartmax, radial and theta mesh size: nr, nth.
c Don't change anything else.
      parameter (npartmax=1000000,np=1,ndim=6)
c Use of particle advance subcycling in inner regions for accuracy.
      logical lsubcycle
c Integrator type. True=old, False=new symplectic schemes
      logical verlet
c CIC definitions
      logical LCIC
      integer NRUSED,NTHUSED,NRFULL,NTHFULL
      parameter (LCIC=.true.)
      integer nrsize,nthsize
c These correspond to nrfull and nthfull.
      parameter (nrsize=200,nthsize=201)
c Positions and velocities of particles (6-d phase-space).
      real xp(ndim,npartmax)
      real vzinit(npartmax)

      real dtprec(npartmax)
c Flag of particle slot status (e.g. in use or not)
      integer ipf(npartmax)
c The particle number
      real psum(0:nrsize,0:nthsize)
c The trapped particle number
      real ptsum(0:nrsize,0:nthsize)
c The sum of particle radial velocities
      real vrsum(0:nrsize,0:nthsize)
c Sum of theta velocities
      real vtsum(0:nrsize,0:nthsize)
c Sum of phi velocities
      real vpsum(0:nrsize,0:nthsize)
c The sum of particle velocities squared
      real v2sum(0:nrsize,0:nthsize)
c The sum of radial particle velocities squared
      real vr2sum(0:nrsize,0:nthsize)
c The sum of non-radial particle velocities squared
      real vtp2sum(0:nrsize,0:nthsize)
c The sum of particle z-velocities
      real vzsum(0:nrsize,0:nthsize)
c The potential normalized to Te/e
      real phi(0:nrsize,0:nthsize)
c Charge density
      real rho(0:nrsize,0:nthsize)
c Injection complement. How many particles to reinject each step
      integer ninjcomp
c Highest occupied particle slot.
      integer iocprev

      real pi
      parameter (pi=3.1415927)
      real cerr,bdyfc,Ti,vd,Bz
      logical diags,lplot,ldist,linsulate,lfloat,lat0,lfext,localinj
      logical lfixedn
      integer nfvdist
      integer myid,numprocs
      real rmtoz
      common /piccom/xp,npart,vzinit,psum,dtprec,ptsum,
     $     vrsum,vtsum,vpsum,v2sum,vr2sum,vtp2sum,vzsum,
     $     phi,rho,cerr,bdyfc,Ti,vd,diags,ninjcomp,
     $     lplot,ldist,linsulate,lfloat,lat0,lfext,localinj,lfixedn,
     $     nfvdist,myid,numprocs,rmtoz,ipf,iocprev,Bz,lsubcycle,verlet
c*********************************************************************
c Radius mesh
      real r(0:nrsize),rcc(0:nrsize)
c Theta angle mesh
      real th(0:nthsize),tcc(0:nthsize)
c Theta mesh radians
      real thang(0:nthsize)
c Inverse of volume of radial shells.
      real volinv(0:nrsize)
c Precalculation functions
      integer nrpre,ntpre
      parameter (nrpre=4*nrsize,ntpre=4*nthsize)
      integer irpre(nrpre),itpre(ntpre)
      real rfac,tfac
c Non-uniform handling quantities.
      real hr(0:nrsize+1),zeta(0:nrsize+1),zetahalf(0:nrsize+1),
     $     cminus(nrsize),cmid(nrsize),cplus(nrsize)
c Lower limit of averaging range. 0.6 by default
      real avelim
c Parallel or serial solving
      logical sorparallel
c Parallel bloc solver arguments
      integer idim1,idim2

      common /meshcom/r,rcc,th,tcc,thang,volinv,irpre,itpre,rfac,tfac,
     $     hr,zeta,zetahalf,cminus,cmid,cplus,avelim
     $     ,nr,NRFULL,NRUSED,nth,NTHFULL,NTHUSED,sorparallel,
     $     idim1,idim2
c********************************************************************
c Random interpolate data.
      integer nvel,nQth
      parameter (nvel=50)
      parameter (nQth=200)
      real Qcom(nQth) 
      real Gcom(nvel,nQth)
      real Vcom(nvel)
      real pu1(nvel),pu2(nvel)
      real Pc(nQth,nvel)
c New BC
      integer bcphi,bcr
      logical infdbl
c Reinjection flux as a function of cos(theta) (line) and chi (column,
c from 0 to 9)
      common /rancom/Gcom,Vcom,Qcom,pu1,pu2,Pc,infdbl,bcphi,bcr
c********************************************************************
c diagnostic data
      integer nvmax,nrein,nreintry,ninner,nstepmax,nrfrcmax
      parameter (nvmax=60,nstepmax=10001,nrfrcmax=10)
      real nvdiag(nvmax),nvdiagave(nvmax),vdiag(nvmax)
      real vrdiagin(nvmax),vtdiagin(nvmax)
      real vrange
      real diagrho(nrsize),diagphi(nrsize)
cIHH      real diagchi(nthsize)
      real diagchi(0:nthsize)
      real diagvr(nrsize,nthsize),diagtrap(nrsize,nthsize)
      integer collision,partz,fieldz,epressz,enccharge
      parameter(enccharge=1,fieldz=2,epressz=3,partz=4,collision=5)
c Total particle flux to probe
      real fluxprobe(nstepmax)
c Total z-momentum flux to probe
      real zmomprobe 
c Total z-momentum at injection for collected particles
      real collmom
c Total energy collected
      real enerprobe
c Z-momentum flux across outer boundary.
      real zmout
c An array of particle zmomentum data:
      real zmparts(nrfrcmax)
      equivalence (zmomprobe,zmparts(1))
      equivalence (zmout,zmparts(2))
c To be obtained at radial mesh numbers
      integer izmrad(nrfrcmax)
c Combined zmom data: charge, fields, electron pressure, ion momentum.
c For inner 1, outer 2, Then possibly more radii.
      real zmom(nstepmax,5,nrfrcmax)
      integer nrfrc
c collmomtot is the reduced collmom for each time-step
      real collmomtot(nstepmax)
c enertot is the reduced enerprobe for each time-step
      real enertot(nstepmax)
c Number of particles striking probe in theta cell
      integer ninthstep(nthsize,0:nstepmax)
      integer ninth(nthsize)
c Sum of zvelocities of particles striking probe.
      real finthave(nthsize)
c Number of particles reinjected per theta cell.
c      integer noutrein(nth),ivoutrein(nth)
c Sum and average of potentials at which particles were reinjected.
      real spotrein,averein
c Flux of reinjection.
      real fluxrein
c Density at infinity.
      real rhoinf
c Number of trapped reinjections
      integer ntrapre
c Coefficient of density deficit, for external solution
      real adeficit
c Cell in which to accumulate distribution functions
      integer ircell,itcell
      common /diagcom/nvdiag,nvdiagave,vdiag,vrange,diagrho,diagphi,
     $     diagchi,nrein,nreintry,ninner,fluxprobe,ninthstep,ninth,
     $     rhoinf,diagvr,vrdiagin,vtdiagin,
     $     spotrein,averein,fluxrein,ntrapre,adeficit,
     $     ircell,itcell,nrfrc,zmparts,izmrad,collmom,finthave,zmom,
     $     collmomtot,enerprobe,enertot,diagtrap
c*********************************************************************
c Poisson coefficients for iterative solution, etc.

      real debyelen,vprobe,Ezext
      real apc(0:nrsize),bpc(0:nrsize)
      real cpc(0:nrsize,0:nthsize),dpc(0:nrsize,0:nthsize)
      real fpc(0:nrsize,0:nthsize)
      real gpc(0:nthsize,1:5)
      common /poisson/debyelen,vprobe,Ezext,apc,bpc,cpc,dpc,fpc,gpc
c*********************************************************************
c Smoothing steps
      integer nstepsave,nsamax
      common /stepave/nstepsave,nsamax
c*********************************************************************
c Orbit plotting storage for tracking the first norbits orbits.
      integer nobsmax,norbits
      parameter (nobsmax=100)
      real xorbit(nstepmax,nobsmax),yorbit(nstepmax,nobsmax),
     $     zorbit(nstepmax,nobsmax)
      real vxorbit(nstepmax,nobsmax),vyorbit(nstepmax,nobsmax),
     $     vzorbit(nstepmax,nobsmax),rorbit(nstepmax,nobsmax)
      integer iorbitlen(nobsmax)
      common /orbits/norbits,iorbitlen,xorbit,yorbit,zorbit,rorbit,
     $     vxorbit,vyorbit,vzorbit

c*********************************************************************
c Data necessary for the orbit tracking
      logical orbinit
      integer maxsteps,trackinit
      common /orbtrack/orbinit,maxsteps,trackinit
c*********************************************************************
c Monte Carlo based reinjection (mcr)
c     Number of injected particles
      integer mcrninjd
c     Positions and velocities of injected particles (6-d phase-space)
      real mcrxpinjd(ndim,npartmax)
c     Mcr common block
      common /mcr/mcrninjd,mcrxpinjd
