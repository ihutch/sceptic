c****************************************************************
c FVCX function for 1-d drifting CX distribution.
      function fvcx(u,ud)
      real u,ud,v,vd,fvcx
c Return the normalized distribution function v_n f(v) for constant 
c cx collision frequency at a value of normalized velocity u=v/v_n,
c when the normalized drift velocity is ud= (a/\nu_c) /v_n,
c with v_n = sqrt(2T_n/m). a is acceleration, nu_c collision freq.
      if(ud.lt.0.) then
         v=-u
         vd=-ud
      else
         v=u
         vd=ud
      endif
      if(vd.eq.0.)then
         carg=20.
         earg=100
      else
         carg=0.5/vd-v
         earg=(0.5/vd)**2-v/vd
      endif
      if(carg.gt.10)then
c asymptotic form for large exp argument (small vd):
c  exp(-v^2)/[sqrt(\pi)(1-2 v_d v)]:
         fvcx=exp(-v**2)/1.77245385/(1.-2.*vd*v)
      elseif(carg.gt.-5.)then
         fvcx=exp(-v**2)*experfcc(carg)*0.5/vd
      else
c         fvcx=exp(earg)*erfcc(carg)*0.5/vd
         fvcx=exp(earg)/vd
      endif
c      write(*,*)'fvcx:vd,v,earg,fvcx',vd,v,earg,fvcx
      if(.not.fvcx.ge.0) then
         write(*,*)'fvcx error. u=',u,' ud=',ud,' f=',fvcx,carg
         fvcx=0.
      endif
      end
c****************************************************************
c (ERFCC is in randf.f) this is exp*erfc
      FUNCTION expERFCC(X)
      Z=ABS(X)      
      T=1./(1.+0.5*Z)
      expERFCC=T*EXP(-1.26551223+T*(1.00002368+T*(.37409196+
     *    T*(.09678418+T*(-.18628806+T*(.27886807+T*(-1.13520398+
     *    T*(1.48851587+T*(-.82215223+T*.17087277)))))))))
      IF (X.LT.0.) expERFCC=2.*exp(z**2)-expERFCC
      END
c*******************************************************************

c********************************************************************
c Given a monotonic (increasing?) 
c function Q(x) on a 1-D grid x=1..nq, solve Q(x)=y for x.
c That is, invert Q to give x=Q^-1(y).
      subroutine f1invtfunc(Q,nq,y,x)
c Comment out this next declaration if you want Q to be a function
      real Q(nq)
      integer nq
      real y,x
c
      integer iqr,iql,iqx
      real Qx,Qr,Ql
      Ql=Q(1)
      Qr=Q(nq)
      iql=1
      iqr=nq
      if((y-Ql)*(y-Qr).gt.0.) then
c Value is outside the range.
         x=0
         return
      endif
 200  if(iqr-iql.eq.1)goto 210
      iqx=(iqr+iql)/2
      Qx=Q(iqx)
c      write(*,*)y,Ql,Qx,Qr,iql,iqr
c Formerly .lt. which is an error.
      if((Qx-y)*(Qr-y).le.0.) then
         Ql=Qx
         iql=iqx
      else
         Qr=Qx
         iqr=iqx
      endif
      goto 200
 210  continue
c Now iql and iqr, Ql and Qr bracket Q
      if(Qr-Ql.ne.0.)then
         x=(y-Ql)/(Qr-Ql)+iql
      else
         x=iql
         write(*,*)'****** Error!: finvtfunc coincident points'
      endif
      end
c**********************************************************************

c********************************************************************
c Given two monotonic functions (arrays) Q1(x) Q2(x)
c  on a 1-D grid x=1..nq, solve s*Q1(x)+t*Q2(x)=y for x.
c That is, invert Q=s Q1+ t Q2 to give x=Q^-1(y).
c Interpolating inside this routine reduces computational effort.
      subroutine f2invtfunc(Q1,Q2,nq,y,x,s,t)
      integer nq
c Comment out the next declaration if you want Q1, Q2 to be functions.
      real Q1(nq),Q2(nq)
      real y,x,s,t
c Statement function
      Q(j)=s*Q1(j)+t*Q2(j)
c Here on is just finvtfunc:
      integer iqr,iql,iqx
      real Qx,Qr,Ql
      Ql=Q(1)
      Qr=Q(nq)
      iql=1
      iqr=nq
      if((y-Ql)*(y-Qr).gt.0.) then
c Value is outside the range.
         x=0
         return
      endif
 200  if(iqr-iql.eq.1)goto 210
      iqx=(iqr+iql)/2
      Qx=Q(iqx)
c      write(*,*)y,Ql,Qx,Qr,iql,iqr
c Formerly .lt. which is an error.
      if((Qx-y)*(Qr-y).le.0.) then
         Ql=Qx
         iql=iqx
      else
         Qr=Qx
         iqr=iqx
      endif
      goto 200
 210  continue

c Now iql and iqr, Ql and Qr bracket Q
c Trap errors caused by flat sections.
      Qd=Qr-Ql
      if(Qd.eq.0.)then
         x=(iql+iqr)/2.
      else
         x=(y-Ql)/(Qr-Ql)+iql
      endif
      end
c**********************************************************************
c**********************************************************************
c Return the value of f(z) interpolated by index zi
      function finterp(f,zi,nf)
      integer nf
      real f(nf)
      real zi
      i=zi
      if(i.lt.1 .or. i.ge.nf)then 
         write(*,*)'***** Finterp error!',i,zi,nf
      endif
      fz=zi-i
      finterp=f(i)*(1.-fz)+f(i+1)*fz
      end
c**********************************************************************
c*****************************************************************
      program fofvtest
      parameter (nlen=200,ndraw=1000000)
      real fd(nlen)
      common /cxud/cx_ud
      external fvcxud
      cx_ud=10.
      vm=4.+5*cx_ud
      do i=1,nlen
         fd(i)=0.
      enddo
      do i=1,ndraw
         call drawfromfv(fvcxud,vm,v)
         j=1+(v/vm+1.)*0.499999*(nlen-1)
c         write(*,*)j,v
         fd(j)=fd(j)+1.
      enddo
      call yautoplot(fd,nlen)
      call pltend()
      end

      include 'randf.f'
