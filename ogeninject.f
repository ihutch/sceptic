c___c___c___c___c___c___c___c___c___c___c___c___c___c___c___c___c___c___
c
c This code is copyright(c) (2003-5) Ian H Hutchinson hutch@psfc.mit.edu
c
c  It may be used freely with the stipulation that any scientific or
c scholarly publication concerning work that uses the code must give an
c acknowledgement referring to the papers I.H.Hutchinson, Plasma Physics
c and Controlled Fusion, vol 44, p 1953 (2002), vol 45, p 1477 (2003).
c  The code may not be redistributed except in its original package.
c
c No warranty, explicit or implied, is given. If you choose to build
c or run the code, you do so at your own risk.
c
c Version 2.6 Aug 2005.
c___c___c___c___c___c___c___c___c___c___c___c___c___c___c___c___c___c___
c***********************************************************************
c Orbit injection for general gyrotropic distribution function.
c***********************************************************************
c Other versions are in other source files.
      subroutine ogenreinject(i,dt)
      integer i
      real dt
c Common data:
      include 'piccom.f'
      parameter (eup=1.e-10)
      external pu
      logical istrapped
c Testing
      real vdist(nvel)
      real tdist(nthsize)
      real crdist(nthsize),cidist(nthsize)
      common/rtest/crdist,cidist,tdist,vdist

c In this routine we work in velocity units relative to ion thermal till end.
      vscale=sqrt(2.*Ti)
      idum=1
      ilaunch=0
 1    continue
      ilaunch=ilaunch+1
      if(ilaunch.gt.1000)then
         write(*,*)'ilaunch excessive. averein=',averein,' brcsq=',
     $        brcsq,' ierr=',ierr,' rp=',rp
         stop
      endif
c Pick normal velocity from cumulative Pu
      y=ran0(idum)
      call finvtfunc(pu,nvel,y,u)
      iv=u
      dv=u-iv
      u=dv*Vcom(iv+1)+(1.-dv)*Vcom(iv)
      if(dv.gt.1)write(*,*)'Error in u calculation',iv,dv
      vdist(iv)=vdist(iv)+1.
c Pick angle from cumulative Pc.
      y=ran0(idum)
c Here the drift velocity is scaled to the ion temperature.
      Uc=vd/vscale
      uu2=2.*Uc*u
      if(uu2.gt.50.) then
         crt=1.+alog(y)/uu2
      elseif(uu2.lt.-50.) then
         crt=-1.+alog(1-y)/uu2
      elseif(abs(uu2).lt.1.e-5)then
         crt=2.*y -1.
      else
         expuu2=exp(uu2)
c This expression is evaluated very inaccurately if expuu2 is nearly 1.
c That is why such a large limit on abs(uu2) is adopted.
         crt=alog(y*expuu2 + (1-y)/expuu2)/uu2
c The following do not do any better at solving this problem.
c         crt=alog( (y*expuu2 + (1-y)/expuu2)**(1./uu2))
c         crt=-1.+alog(1.+(expuu2**2-1.)*y)/uu2
      endif
      if(.not. abs(crt).le.1)then
c         write(*,*)'Strange crt:',crt,y,expuu2,uu2
c It seems impossible to avoid occasional strange results when uu2 is small.
         crt=2.*y-1.
      endif
c Testing angular distribution.
      if(LCIC)then
         icr=(1.+crt)*0.5*(NTHUSED-1) + 1.5
      else
         icr=(1.+crt)*0.5*(nth-1) + 1
      endif
      crdist(icr)=crdist(icr)+1.
c End of testing distribution monitor.
      srt=sqrt(1.- crt**2)
c Pick angle zt of poloidal impact and angle eta of impact parameter
      zt=ran0(idum)*2.*pi
      czt=cos(zt)
      szt=sin(zt)
      eta=ran0(idum)*2.*pi
      ceta=cos(eta)
      seta=sin(eta)
c Choose impact parameter, preventing overflow.
      chium2=-averein/Ti/(u+eup)**2
      if(chium2.le.-1.) then
         write(*,*)'Impossible chium2=',chium2,' averein=', averein,
     $        ' u=',u,' iv=',iv
c         stop
      endif
c      if(.not.lfixedn)chium2=0.
      brcsq=ran0(idum)*(1.+chium2)
c Reject a particle that will not reach boundary.
      if(brcsq.lt.0.) then
         goto 1
      endif
      brc=sqrt(brcsq)
c Get cosine and sine of impact angle relative to distant position.
c Based on integration.
      p2=brcsq*2.*Ti*u**2
      ierr=0
      if(debyelen.gt..001)then
c Orbit integration angle calculation.
c There is an overflow with this at zero debyelen. Ought to be properly fixed.
         call alphaint(p2,brcsq,cosal,sinal,ierr)
         if(ierr.ne.0)goto 1
c      write(*,'(4f9.4)')cosal-alcos(brc,chium2),sinal-alsin(brc,chium2)
c     Now ilaunch is the number of launches at infinity it took to get
c     one that reached the boundary.
      else
c Alternative based on analytic orbit calculation.
c Used for low debyelen, but really assumes negligible boundary potential.
         call alcossin(brc,chium2,cosal,sinal)
         cosal=alcos(brc,chium2)
         sinal=alsin(brc,chium2)
      endif
c Install reinjection position
      a1=crt*ceta*sinal+srt*cosal
      rs=r(nr)
      xp(1,i)=rs*(czt*a1+szt*seta*sinal)
      xp(2,i)=rs*(-szt*a1+czt*seta*sinal)
      xp(3,i)=rs*(-srt*ceta*sinal + crt*cosal)

c Obtain angle coordinate and map back to th for phihere.
      ct=xp(3,i)/rs
      call invtfunc(th(1),nth,ct,x)
      ic1=x
      ic2=ic1+1
      dc=x-ic1
c This expression should work for CIC And NGP.
      phihere=(phi(NRUSED,ic1)+phi(NRFULL,ic1))*0.5*(1.-dc)
     $        +(phi(NRUSED,ic2)+phi(NRFULL,ic2))*0.5*dc
c Section to correct the injection velocity and direction (but not the
c position) to account for local potential. 26 July 2004.
      if(localinj)then
         brcsq=(brcsq*(1.-phihere/Ti/(u+eup)**2)/(1.+chium2))
         if(brcsq.lt. 0.) then
c     This launch cannot penetrate at this angle. But it would have done
c     if the potential were equal here to averein. Thus it probably
c     should not be counted as part of the launch effort. So
            ilaunch=ilaunch-1
            goto 1
         endif
         chium2=-phihere/Ti/(u+eup)**2
         brc=sqrt(brcsq)
      endif
c Injection velocity components normalized in the rotated frame:
      ua1=-brc*cosal -sqrt(1.+chium2-brcsq)*sinal
      ua3=brc*sinal - sqrt(1.+chium2-brcsq)*cosal
      ua=crt*ceta*ua1+srt*ua3
c Install reinjection velocity in Te-scaled units
      u=u*vscale
      xp(4,i)=u*(czt*ua+szt*seta*ua1)
      xp(5,i)=u*(-szt*ua+czt*seta*ua1)
      xp(6,i)=u*(-srt*ceta*ua1 + crt*ua3)
c Increment the position by a random amount of the velocity.
c This is equivalent to the particle having started at an appropriately
c random position prior to reentering the domain.
      xinc=ran0(idum)*dt
c      xinc=0.
      vdx=0.
      do j=1,3
         vdx=vdx+xp(j,i)*xp(j+3,i)
         xp(j,i)=xp(j,i)+xp(j+3,i)*xinc
      enddo
      if(vdx.gt.0.)then
         write(*,*)'Positive projection. u,phi=',u,phihere
 601     format(a,5G10.5)
      endif
      rcyl=xp(1,i)**2+xp(2,i)**2
      rp=rcyl+xp(3,i)**2
c Reject particles that are already outside the mesh.
      if(.not.rp.lt.r(nr)*r(nr))then
c      if(.not.rp.le.r(nr)*r(nr))then
c         write(*,*)'Relaunch',rp,xp(1,i),xp(2,i),xp(3,i)
         goto 1
      else

c Do the outer flux accumulation.
c In order to accumulate the number of launches at infinity, rather than
c just the number of reinjections, we weight this by ilaunch
         spotrein=spotrein+phihere*ilaunch
         nrein=nrein+ilaunch
         fluxrein=fluxrein+1.
         if(istrapped(i))then
            ntrapre=ntrapre+ilaunch
c            v=sqrt(xp(4,i)**2+xp(5,i)**2+xp(6,i)**2)
c            write(*,*)'Trapped',vdx/rp,u,v,sqrt(u**2-2.*averein)
c crt,czt,ceta,cosal
         endif
      endif
      end
c********************************************************************
c Initialize the distributions describing reinjected particles
      subroutine ogeninjinit()
c Common data:
      include 'piccom.f'

c Here the drift velocity is scaled to the ion temperature.
c And U's are in units of sqrt(2T/m), unlike vd.
      Uc=abs(vd)/sqrt(2.*Ti)
c Range of velocities (times (Ti/m_i)^(1/2)) permitted for injection.
      vspread=5.+abs(Uc)

c Can't use these formulas for Uc exactly equal to zero.
      if(abs(Uc).lt.1.e-4)then
         if(Uc.lt.1.e-20) Uc=1.e-20
         do i=1,nvel
            u0= vspread*(i-1.)/(nvel-1.)
            Vcom(i)=u0
            expu0=exp(-u0**2)
            pu2(i)=2.*Uc*expu0
            pu1(i)=0.5*(4.*u0**2*Uc + 2.*Uc)*expu0
     $           +(Uc**2 +0.5)*pu2(i)
         enddo
      else
         do i=1,nvel
            u0= vspread*(i-1.)/(nvel-1.)
            Vcom(i)=u0
            uplus=u0+Uc
            uminus=u0-Uc
            pu2(i)=0.5*sqrt(pi)*(erfcc(uminus)-erfcc(uplus))
            pu1(i)=0.5*(-uminus*exp(-uplus**2)+uplus*exp(-uminus**2))
     $           +(Uc**2 +0.5)*pu2(i)
         enddo
      endif
      call srand(myid)
      end
c***********************************************************************
c***********************************************************************
c Calculate the cumulative probability for velocity index iu such that
c         u= vspread*(iu-1.)/(nvel-1.)   as per injinit
      real function pugen(iu)
      integer iu
c     averein is the average potential of reinjected particles, which is
c     used as an estimate of the potential at the reinjection boundary.
c     It is expressed in units of Te so needs to be scaled to Ti.
      include 'piccom.f'
      pudenom=pu1(1)-pu2(1)*averein/Ti
      pugen=1.- (pu1(iu)-pu2(iu)*averein/Ti)/pudenom
      end
