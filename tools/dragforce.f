      real function dragforce(vd,Ti,rp,phip,debyelen,rmtoz)
c Return the drag force, using the fit [Hutchinson, PPCF, 48, 185
c (2005)] to sceptic results for the Coulomb log, ln(\Lambda)
 
c On input:
c   vd is the drift velocity normalized to sqrt(Te/mi)
c   Ti is the ion temperature in units of Te
c   rp is the normalized sphere radius.
c   phip is the potential of the sphere in units of Te/e
c   debyelen is the normalized electron debye length
c   rmtoz is the ratio of mass-number to charge-number.

c Value returned 
c   is the force normalized to ne.Te.

c The ratio of debyelen to rp, together with the potential phip, are the
c controlling factors. One can choose any length as the normalizing
c length, but usually it is either the debyelen or rp.  

c If one normalizes by rp, then the input values passed should be rp=1,
c debyelen=\lambda_De/radius, and the force returned will be in units of
c ne.Te.rp^2. If one chooses debyelen, then the values passed should be
c debyelen=1, rp=radius/\lambda_De, and the force returned will be in
c units of ne.Te.\lambda_De^2.  

c For small radius/lambda_De, the normalized charge is 
c \bar{Q} = e Q/(4\pi\epsilon_0 \lambda_{De} T_e) (dimensional) 
c       = phip.radius/lambda_{De} = phip.rp/debyelen 
c This relationship determines the value of phip to use if the force is
c required for specified \bar{Q}: phip=\bar{Q}.lambda_De/radius.
c Normally phip, like Q, is negative.
c
c There seems to be a misprint 0.6 vs 0.5 in the vt2. But that's certainly
c not a big deal, because it's only in the enhancement of the transition
c to zero ion shielding.

      implicit none
      real vd,Ti,rp,phip,debyelen,rmtoz

      real fourpi
      parameter (fourpi=4.*3.1415926)
      real vt2,dterm,xls,b90,xlnlambda2,fcoul,fcol,debye,tempmul
      real chandrafunc,colforce
      external chandrafunc,colforce

      debye=debyelen/rp
c The calculation is done in units where rp=1.
      tempmul=2.
      vt2=tempmul*Ti +vd**2
c Enhancement of transition to zero ion shielding.
     $     +vd**2*(vd/(0.6+0.05*alog(rmtoz)
     $     +.1 +(debye/5.)**1*(sqrt(Ti)-.1)
     $     ))**3

      dterm=debye**2*(vt2/(tempmul+vt2)) +1.
      xls=sqrt(dterm)
      b90=-phip/vt2
c Square of the fitted Coulomb Logarithm:
      xlnlambda2=2.*max(log((xls+b90)/(1.+b90)),0.)
c Orbital force:
      fcoul=phip**2*(1./(2.*Ti))*fourpi*xlnlambda2
     $           * chandrafunc(vd/sqrt(2*Ti))
c Direct collection force:
      fcol=sqrt(3.14159)*Ti*colforce(vd/sqrt(2.*Ti),-phip/Ti)

      dragforce=fcoul+fcol
c Scale by rp^2 so force is returned in the appropriate normalization.
      dragforce=dragforce*rp**2

      end
c ********************************************************************
c OML Collection force for a drifting Maxwellian.
      real function colforce(u,chi)
      if(u.ne.0.)then
      colforce=(u*(2.*u**2+1+2.*chi)*exp(-u**2) +
     $     (4.*u**4+4.*u**2-1-2.*(1-2.*u**2)*chi)*sqrt(3.14159)*0.5
     $     *(1.-erfcc(u)) )/u**2
      else
         colforce=0.
      endif
      end
c-------------------------------------------------------------------
c Chandrasekhar function G, for drifting Maxwellian.
      real function chandrafunc(x)
      if(x.ne.0.)then
         chandrafunc=((1.-erfcc(x))-2.*x*(exp(-x**2)/sqrt(3.14159)))
     $     /(2.*x**2)
      else
         chandrafunc=0.
      endif
      end
c-------------------------------------------------------------------
c Complementary error function.
      real FUNCTION ERFCC(X)
      Z=ABS(X)      
      T=1./(1.+0.5*Z)
      ERFCC=T*EXP(-Z*Z-1.26551223+T*(1.00002368+T*(.37409196+
     *    T*(.09678418+T*(-.18628806+T*(.27886807+T*(-1.13520398+
     *    T*(1.48851587+T*(-.82215223+T*.17087277)))))))))
      IF (X.LT.0.) ERFCC=2.-ERFCC
      RETURN
      END
c*********************************************************************
      program dragforcetest
      parameter (np=100)
      real v(np),q(np),force(np)
      rmtoz=1.
      Ti=.01
      if(.false.)then
         vd=.3
         debyelen=1.
         rp=.01
         q0=0.005
         do i=1,np
            q(i)=q0*10**(2.*(i-1.)/(np-1.))
            phip=-q(i)*debyelen/rp
            force(i)=dragforce(vd,Ti,rp,phip,debyelen,rmtoz)
            write(*,*)i,q(i),phip,force(i)
         enddo
         call lautoplot(q,force,np,.true.,.true.)
      endif
      rp=1.
      debyelen=20.
      vmax=5.
      phip=-3.
      do i=1,np
         v(i)=i*vmax/np
         force(i)=dragforce(v(i),Ti,rp,phip,debyelen,rmtoz)
         write(*,*)i,v(i),force(i)
      enddo
      call autoplot(v,force,np)
      call pltend()
      end
