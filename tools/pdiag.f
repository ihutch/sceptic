c Diagnose the final particle distributions.
c 8 Oct 05.
c Results with this program using particle data saved from runs with
c -x20 -d.2 -s500 -z1.e-5 show very good agreement between 
c the analytic form of f(v) and f(v) from the particle data.
c That confirms that the particle injection in being done correctly, since
c -z1.e-5 is essentially a zero field case, and so one ought to arrive
c at steady state where the distribution is the injection distribution.
      program pdiag
      include 'piccom.f'
      logical success
      integer nv
      parameter (nv=50,nfmax=100)
      real fv1(-nv:nv),fv2(-nv:nv),fv3(-nv:nv),v(-nv:nv)
     $     ,fmaxwell(-nv:nv),fz(-nv:nv)
      character*100 string
      character*100 filename(nfmax)

c Defaults
      Ti=1.
      rmax=50.
      rmin=0.
      cmin=-1.
      cmax=1.
      vd=0.
      vneutral=0.
      myid=0
      ifile=0

c Deal with arguments.
      if(iargc().eq.0) goto 51
      do 1 i=1,iargc()
         call getarg(i,string)
         if(string(1:1) .eq. ' ') then
            goto 3
         endif
         if(string(1:2) .eq. '-s') diags=.true.
         if(string(1:2) .eq. '-t') read(string(3:),*)Ti
         if(string(1:3) .eq. '-rx') read(string(4:),*)rmax
         if(string(1:3) .eq. '-rn') read(string(4:),*)rmin
         if(string(1:3) .eq. '-cx') read(string(4:),*)cmax
         if(string(1:3) .eq. '-cn') read(string(4:),*)cmin
         if(string(1:3) .eq. '-vn')then
            read(string(4:),*)vneutral
         elseif(string(1:2) .eq. '-v')then
            read(string(3:),*)vd
         endif
         if(string(1:3) .eq. '-pf') then
            lfloat=.true.
         elseif(string(1:3) .eq. '-pi') then
            linsulate=.true.
         elseif(string(1:2) .eq. '-p') then
            read(string(3:),*)vprobe
         endif
         if(string(1:2) .eq. '-m') read(string(3:),*)myid
         if(string(1:2) .eq. '-?') goto 51
c         if(string(1:1) .ne. '-')read(string(1:),*)filename
         if(string(1:1) .ne. '-')then
            if(ifile.lt.nfmax)then
               ifile=ifile+1
               filename(ifile)=string
            endif
         endif
 1    continue
 3    continue
      if(vneutral.eq.0.)then
         vrange=max(4.,4.+8.*vd/sqrt(2.*Ti))*sqrt(2.*Ti)
      else
         vrange=max(4.,4.+4.*vd/sqrt(2.*Ti))*sqrt(2.*Ti)
      endif
      do j=-nv,nv
         v(j)=j*vrange/nv
         fv1(j)=0.
         fv2(j)=0.
         fv3(j)=0.
         fmaxwell(j)=exp(-v(j)**2)
      enddo

      vmean=0.
      vmeanplus=0.
      vmeanminus=0.
      ncount=0
      nplus=0
      nminus=0

      do myid=0,ifile-1
         write(*,'(a,i4,a,a)')'Reading particles for myid=',myid
     $        ,' from ',filename(myid+1)(1:40)
         call partrd(filename(myid+1),success)
         if(.not.success)then
            write(*,*)'Quitting read attempts'
            goto 10
         else
            write(*,*)'npart=',npart
            do i=1,npart
               rp=sqrt(xp(1,i)**2+xp(2,i)**2+xp(3,i)**2)
               costheta=xp(3,i)/rp
c               write(*,*)i,rp,costheta
               if(rmin.le.rp .and. rp.le.rmax)then
                  if(cmin.le.costheta.and.costheta.le.cmax)then
                  iv1=min(nv,max(-nv,nint(xp(4,i)*nv/vrange)))
                  iv2=min(nv,max(-nv,nint(xp(5,i)*nv/vrange)))
                  iv3=min(nv,max(-nv,nint(xp(6,i)*nv/vrange)))
                  fv1(iv1)=fv1(iv1)+1
                  fv2(iv2)=fv2(iv2)+1
                  fv3(iv3)=fv3(iv3)+1
                  vmean=vmean+xp(6,i)
                  ncount=ncount+1
                  if(costheta.gt.0)then
                     vmeanplus=vmeanplus+xp(6,i)
                     nplus=nplus+1
                  else
                     vmeanminus=vmeanminus+xp(6,i)
                     nminus=nminus+1
                  endif
                  endif
               endif
            enddo
         endif
      enddo
 10   vmean=vmean/ncount
      write(*,*)'Counted particles',ncount,'  Mean vz=',vmean,' cf vd='
     $     ,vd
      vmeanplus=vmeanplus/nplus
      vmeanminus=vmeanminus/nminus
      write(*,*)'Mean vzplus=',vmeanplus,'  vzminus=',vmeanminus,
     $     ' nplus,nminus=',nplus,nminus
c Don't plot the bottom end if there's no data.
      do i=-nv,nv
         nvmin=i
         if(fv3(nvmin).ne.0 .or. fv1(nvmin).ne.0) goto 50
      enddo
 50   continue
      nvtot=nv+1-nvmin

      write(string,81)rmin,rmax,cmin,cmax
 81   format(' "',f6.2,'<r<',f6.2,' ; ',f6.3,'<cos!Aq!@<',f6.3,'"')

      call pfset(3)
      call minmax(fv1(nvmin),nvtot,vmin,vmax)
      do j=-nv,nv
         fmaxwell(j)=exp(-v(j)**2/(2.*Ti))*vmax
      enddo
      call autoplot(v(nvmin),fv1(nvmin),nvtot)
      call axlabels('v','f(v)')
      call legendline(.6,.91,0,'f(v1)')
      call color(2)
      call polyline(v(nvmin),fv2(nvmin),nvtot)
      call legendline(.6,.86,0,'f(v2)')
      call color(3)
      call polyline(v(nvmin),fv3(nvmin),nvtot)
      call legendline(.6,.81,0,'f(v3)')
      call dashset(1)
      call color(4)
      call polyline(v(nvmin),fmaxwell(nvmin),nvtot)
      call legendline(.6,.76,0,'fmaxwell')
      do j=-nv,nv
         if(vneutral.eq.0)then
            fz(j)=1.77*vmax*fvcx(v(j)/sqrt(2.*Ti),vd/sqrt(2.*Ti))
         else
            fz(j)=exp(-(v(j)-vd)**2/(2.*Ti))*vmax
         endif
      enddo
      call color(5)
      call polyline(v(nvmin),fz(nvmin),nvtot)
      call legendline(.6,.71,0,'fz [analytic]')
      call color(15)
      call boxtitle(string(3:lentrim(string)-1))
      call pltend()
     
      write(*,'(a)')'legend:f(v!d1!d)'
      write(*,'(a)')'legend:f(v!d2!d)'
      write(*,'(a)')'legend:f(v!d3!d)'
      write(*,'(a)')'legend:ftheory'
      write(*,*)nvtot,4
      write(*,'(f12.5,4f12.1)')(v(nvmin+k),fv1(nvmin+k),fv2(nvmin+k)
     $     ,fv3(nvmin+k),fz(nvmin+k),k=0,nvtot-1)
      write(*,'(''annotation:'',2f10.4,a)').3,.67,string
      write(*,*)

      call exit(0)
 51   write(*,*)'Usage: pdiag [-t -rn -rx -cn -cx -v vn -m] filename'
      write(*,*)'    [tempr, rmin, rmax, cmin, cmax,'
     $     ,' velocity, vneutral, mass]'
      end

c**********************************************************************
c Read in the particle data.
      subroutine partrd(filename,success)
      logical success
      character*(*) filename
c Common data:
      include 'piccom.f'

c      write(filename,'(''part'',i3.3,''.dat'')')myid
c      write(*,*)'filename=',filename
      success=.false.
      open(11,file=filename,status='old',err=101)
      read(11,*,err=100,end=100)ipartmax,npart,ir,ith,idim,ip
      if(ipartmax.eq.npartmax)then
         write(*,*)'Using saved particle data.'
         read(11,*,err=100,end=100)((xp(i,j),i=1,ndim),j=1,npart)
c         read(11,*,err=100,end=100)xp
         read(11,*,err=100,end=100)rhoinf,spotrein,averein
      write(*,*)'rhoinf,spotrein,averein',rhoinf,spotrein,averein
         success=.true.
      else
         write(*,*)'Particle data mismatch',ipartmax,npartmax,ipart,
     $        npart,ir,nr,ith,nth
      endif
      close(11)
c      write(*,*)'partrd',npart

      return
 100  close(11)
      write(*,*) 'Error reading particle file',filename
      return
 101  write(*,*) 'No particle file to read.',filename
      end
c**********************************************************************
c****************************************************************
c FVCX function for 1-d drifting CX distribution.
      function fvcx(u,ud)
      real u,ud,v,vd,fvcx
c Return the normalized distribution function v_n f(v) for constant 
c cx collision frequency at a value of normalized velocity u=v/v_n,
c when the normalized drift velocity is ud= (a/\nu_c) /v_n,
c with v_n = sqrt(2T_n/m). a is acceleration, nu_c collision freq.
c If the neutrals are drifting, then u and ud should be the drift velocity
c relative to vneutral.
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
c******************************************************************
c Obtain the length of a string omitting trailing blanks.
      function lentrim(string)
      character*(*) string
      do i=len(string),1,-1
         if(string(i:i).ne.' ') goto 101
      enddo
      i=0
 101  lentrim=i
      end
