      include 'piccom.f'
      include 'distcom.f'
      character*(35) filename,string

      logical lgraph,lphip,ldens,ltempc,lconline,larrows
      logical laspect
      real xplot1(2*nthsize),yplot1(2*nthsize)
c      real xplot2(0:nthsize),yplot2(0:nthsize)

      data lgraph,lphip,ldens,ltempc,lconline,larrows/6*.false./

      filename=' '
      isw=0
      jstepth=1
      rhomax=1.
      ir1=1
      it1=1
      ir2=2
      it2=2
c Deal with arguments
      do i=1,iargc()
         call getarg(i,string)
C         write(*,*)'Argument:',string(1:40)
         if(string(1:1) .eq. ' ') then
            goto 3
         endif
         if(string(1:1) .eq. '-') then
            if(string(1:2) .eq. '-x') lgraph=.not.lgraph
            if(string(1:2) .eq. '-f') lphip=.true.
            if(string(1:2) .eq. '-n') ldens=.true.
            if(string(1:2) .eq. '-t') ltempc=.true.
            if(string(1:2) .eq. '-l') lconline=.true.
            if(string(1:2) .eq. '-a') larrows=.true.
            if(string(1:2) .eq. '-i') read(string(3:),*)isw
            if(string(1:2) .eq. '-p')then
               read(string(3:),*,err=201,end=201)ir1,ir2,it1,it2
 201           if(ir2.le.ir1)ir2=ir1+1
               if(it2.le.it1)it2=it1+1
            endif
            if(string(1:2) .eq. '-m') read(string(3:),*)nbsm
            if(string(1:2) .eq. '-?') goto 11
         else
            filename=string
         endif
      enddo
 3    continue

c      write(*,*)'isw=',isw
      if(filename(1:1).ne.' ')then
         call readvdist(filename,isw)
      else
         write(*,*)'No filename.'
         goto 11
      endif

      call minmax2(rhodist(1,1),nrdist+1,nrused-1,nthused-1,rhomin
     $     ,rhomax)
      write(*,*)'rhomax,rhomin',rhomax,rhomin
c Now process. Draw a contour plot of rhodist, and indicate where we are.
      if(nthused.le.10.and.nrused.le.10)then
         write(*,*)'rcc,tcc,th:'
         write(*,100)(rcc(ir),ir=1,nrused)
         write(*,100)(tcc(it),it=1,nthused)
c         write(*,100)(th(it),it=1,nthused)
         write(*,*)'rhodist:'
         do it=1,nthused
            write(*,100)(rhodist(ir,it),ir=1,nrused)
         enddo
 100     format(10F8.3)
      endif

c Fix up grid edges.
      rcc(0)=1.
      rcc(nrused+1)=rcc(nrused)
      tcc(0)=1.
      tcc(nthused+1)=-1.

c Default positions
      ir=1
      v1=1.

 10   continue
c Radial limits
      r1=0.5*(rcc(ir1)+rcc(ir1-1))
      r2=0.5*(rcc(ir2)+rcc(ir2-1))
c      write(*,*)ir1,ir2,it1,it2,r1,r2
      ic=0
c Angular limits
      do it=it1,it2
         ic=ic+1
         cth=0.5*(tcc(it-1)+tcc(it))
         xplot1(ic)=r1*cth
         yplot1(ic)=r1*sqrt(1.-cth**2)
c         write(*,*)it,ic,xplot1(ic),yplot1(ic),cth
      enddo
      do it=it2,it1,-1
         ic=ic+1
         cth=0.5*(tcc(it-1)+tcc(it))
         xplot1(ic)=r2*cth
         yplot1(ic)=r2*sqrt(1.-cth**2)
c         write(*,*)it,ic,xplot1(ic),yplot1(ic),cth
      enddo
      ic=ic+1
      xplot1(ic)=xplot1(1)
      yplot1(ic)=yplot1(1)

c Plotting
      call axregion(.6,.95,.1,.3)
      call multiframe(2,2,3)

      call plotdists(ir1,ir2,it1,it2)
         call conrho(ir,jstepth,rhomax,rhomin,nrused,nthused,v1,larrows
     $        ,lconline,rhodist)
c         write(*,*)ir1,it1,rhodist(ir1,it1),rhomin,rhomax
      if(abs(rhodist(ir1,it1)-rhomax)
     $     .gt.2.*abs(rhodist(ir1,it1)-rhomin))then
         call color(14)
      else
         call color(15)
      endif
      call polyline(xplot1,yplot1,ic)
      call color(15)

      call uinterface(iquit,laspect,jsw,iclipping
     $     ,ir1,it1,ir2,it2,nrused+1,nthused+1,icontour,iweb)
      if(iquit.eq.0)goto 10

      call exit(0)
 11   continue
      write(*,*)'Usage ...'
      write(*,*)'-p<ir1,ir2,it1,it2> specify initial cell.'
      write(*,*)'-i1 give read-back commentary.'
      end
c***************************************************************************
c Contouring of the charge density, rho, on distorted mesh.
      subroutine conrho(ir,it,rhomax,rhomin,nrhere,nthhere,v1,
     $     larrows,lconline,rholocal)
      integer ir,it
      real rhomax,v1
      logical larrows,lconline
c Common data:
      include 'piccom.f'
      real rholocal(0:nrsize,0:nthsize)
c      save
      character*20 cstring
      character*30 tstring
      character cworka(nrsize*(nthsize+1+1))
      integer ncont
      parameter (ncont=12)
      real zclv(ncont)
      real xrho(nrsize+1,0:nthsize+1),zrho(nrsize+1,0:nthsize+1)
c      save xrho,zrho
      real basesize
      parameter (basesize=.02)

      if(rhomin.gt.0)call setconlog(.true.)
      if(nthhere.gt.nthsize.or.nthhere.le.0)then
         write(*,*)' Conrho error. Mesh required:',nrhere,nthhere,
     $        ' Inconsistent with allocated:',nrsize,nthsize
         stop
      endif

      do i=1,nrhere
         do j=2,nthhere-1
            zrho(i,j)=rcc(i)*tcc(j)
            xrho(i,j)=rcc(i)*sqrt(1.-tcc(j)**2)
         enddo
c Fix the angle ends for plotting:
         zrho(i,1)=rcc(i)*0.25*(3.+tcc(2))
         xrho(i,1)=rcc(i)*sqrt(1.-(0.25*(3.+tcc(2)))**2)
         zrho(i,nthhere)=rcc(i)*0.25*(-3.+tcc(nthhere-1))
         xrho(i,nthhere)=rcc(i)*sqrt(1.-(0.25*(-3.+tcc(nthhere-1)))**2)
         zrho(i,0)=rcc(i)
         xrho(i,0)=0.
         zrho(i,nthhere+1)=-rcc(i)
         xrho(i,nthhere+1)=0.         
         rholocal(i,0)=rholocal(i,1)
         rholocal(i,nthhere+1)=rholocal(i,nthhere)
      enddo

c      rpmax=2.*rcc(nrhere)-rcc(nrhere-1)
      rpmax=rcc(nrhere)
c      write(*,*)rcc(nrhere),nrhere,nthhere,rpmax


      do j=1,ncont
c         zclv(j)=rhomin+(rhomax-rhomin)*(0.95*(j-1)/float(ncont-1))
         zclv(j)=(rhomax-rhomin)*(1.*(j-1)/float(ncont-1))+rhomin
      enddo
c      write(*,*)rhomax,rhomin,(zclv(j),j=1,ncont)
      icl=-ncont

c      call multiframe(2,2,3)

      call ticnumset(8)
      call pltinaspect(-rpmax,rpmax,0.,rpmax)
      call accisgradinit(-25000,00000,25000,130000,65000,130000)
      ntype=2+16+32
      call contourl(rholocal(1,0),cworka,nrsize+1,nrhere,nthhere+2,
     $        zclv,icl,zrho,xrho,ntype)
      call color(15)
      call gradlegend(zclv(1),zclv(abs(icl)),
     $     .1,1.25,.9,1.25,0.02,.true.)
c Call a second time for contours, without the highest.
      tstring(1:1)=char(0)
         if(lconline)then
            call fitrange(rhomin,rhomax,ncont,ipow,
     $           fac10,delta,first,xlast)
c               write(*,*)'rhomax,fac10,first,delta',
c     $              rhomax,fac10,first,delta
            do j=1,ncont
               zclv(j)=(first+j*delta)
            enddo
            ntype=2
            icl=(ncont-1)
            write(*,*)nrsize,nrhere,nthhere
            call contourl(rholocal(1,0),cworka,nrsize+1,nrhere,
     $           nthhere+2,zclv,icl,zrho,xrho,ntype)
            write(*,'(a,30f5.2)')'Density Contours=',zclv
            call fwrite(delta,iwd,1,cstring)
            tstring=' contour spacing: '//cstring(1:10)
c            call legendline(-.1,-.22,258,tstring)
         endif
c      endif
c Fit closer than boxtitle
      call legendline(0.47,1.07,258,'n/n!A!d;!d!@'//tstring)
c      call boxtitle('n/n!A!d;!d!@')
      call color(15)
      call axis()
      call axlabels('z','r sin!Aq!@')
  
      
      if(larrows) then
         write(*,*)ir,it,nthhere,nrhere
         call color(12)
         if(ir.le.0.or.ir.ge.100) ir=10
         do j=1,nthhere,it
            do i=1,nrhere,max(nrhere/ir,1)
               vri=vrsum(i,j)/(psum(i,j)+1.e-5)
               vti=vtsum(i,j)/(psum(i,j)+1.e-5)
c               write(*,*)i,j,vri,vti
               size=basesize/v1*sqrt(vri**2+vti**2)
               angle=atan2(vti,vri)+acos(tcc(j))
               call charsize(size,0.3*size)
               call charangl(180.*angle/3.141593)
               call jdrwstr(wx2nx(zrho(i,j)),wy2ny(xrho(i,j)),
     $              '!A_!@',0.)
            enddo
         enddo
         call charangl(0.)
         size=basesize
         call charsize(size,0.3*size)
         call legendline(0.8,0.95,258,'!A_!@'//char(0))
         call charsize(0.,0.)
         write(cstring,'(''  v='',f4.1)') v1
         call color(15)
         call legendline(0.8,0.95,258,cstring(1:8)//char(0))
      endif

c      call pltend()
      call setconlog(.false.)
      end

c******************************************************************
      subroutine uinterface(iquit,laspect,jsw,iclipping
     $     ,if1,if2,nf1,nf2,maxf1,maxf2,icontour,iweb)
c Modified from accis ui3d.
c [Plane-position] n1 is controlled by up/down arrows, within the range
c 1-iuds(idfix), where iuds(3) are the used dimensional lengths, and
c idfix is the one currently being fixed (sliced).
c if1 nf1, if2 nf2 control the clipping positions  
c jsw is to do with contouring. laspect preserves aspect-ratio.
c iquit is returned as non-zero to command an end to the display.

      logical laspect,ltellslice
      data ips/0/
      save ips
c 3d display user interface.
c-----------------------------------
      iquit=0
c Limit framing rate to 30fps.
      call usleep(15000)

      if(ips.ne.0)then
c We called for a local print of plot. Terminate and switch it off.
         call pltend()
         call pfset(0)
         ips=0
      endif
c User interface interpret key-press.
 24   call eye3d(isw)
c      write(*,*)'isw',isw
      if(isw.eq.ichar('f')) goto 24
      if(isw.eq.0) iquit=1
c      if(isw.eq.65364 .and. n1.gt.1) n1=n1-1
c      if(isw.eq.65362 .and. n1.lt.iuds(idfix)) n1=n1+1
      if(isw.eq.ichar('q')) iquit=1
      if(isw.eq.ichar('a')) laspect=.not.laspect
      if(isw.eq.ichar('d')) call noeye3d(0)
      if(isw.eq.ichar('s')) jsw=1 + 256*6 + 256*256*7
c      if(isw.eq.ichar('t')) call togi3trunc()
      if(isw.eq.ichar('p'))then
         call pfset(3)
         ips=3
      endif
c Adjust clipping
      if(isw.eq.ichar('m'))then
c m
         iclipping=1
         nf2=min(nf2+1,maxf2)
      elseif(isw.eq.44)then
c ,
         iclipping=1
         nf2=max(nf2-1,2)
         if2=min(nf2-1,if2)
      elseif(isw.eq.46)then
c .
         iclipping=1
         if2=min(if2+1,maxf2-1)
         nf2=max(if2+1,nf2)
      elseif(isw.eq.47)then
c /
         iclipping=1
         if2=max(if2-1,1)
      elseif(isw.eq.ichar('l'))then
c l
         iclipping=1
         nf1=max(nf1-1,2)
         if1=min(nf1-1,if1)
      elseif(isw.eq.59)then
c ;
         iclipping=1
         nf1=min(nf1+1,maxf1)
      elseif(isw.eq.ichar('j'))then
         iclipping=1
         if1=max(if1-1,1)
      elseif(isw.eq.ichar('k'))then
         iclipping=1
         if1=min(if1+1,maxf1-1)
         nf1=max(if1+1,nf1)
      elseif(isw.eq.ichar('u'))then
         ltellslice=.not.ltellslice
      endif
c      if(isw.eq.ichar('c'))icontour=mod(icontour+1,4)
c      if(isw.eq.ichar('w'))iweb=mod(iweb+1,2)
      if(isw.eq.ichar('h'))then
         write(*,*)' ======== Selection interface:'
         write(*,*)
     $        ' s: rescale. p: print.'
         write(*,*)' (jkl;) (m,./): control plotting extent in 2 axes.'
c         write(*,*)' a: aspect'
         write(*,*)
     $        ' d: disable interface; run continuously.',
     $        ' depress f: to interrupt running.'
      endif
c      call rotatezoom(isw)
c End of user interface.
      end

c*******************************************************************
      subroutine plotdists(ir1,ir2,it1,it2)

      include 'piccom.f'
      include 'distcom.f'
      real vdist(nvdist)

      real fvave(nvdist,ndtype)

c      write(*,*)ndtype,nvdist,it1,it2,ir1,ir2

      do is=1,ndtype
         do iv=1,nvdist
            fvave(iv,is)=0.
         enddo
      enddo

      do it=it1,it2-1
         do ir=ir1,ir2-1
            do is=1,ndtype
               do iv=1,nvdist
c                  write(*,*)it,ir,iv,is
                  fvave(iv,is)=fvave(iv,is)+fvrtdist(iv,is,ir,it)
               enddo
            enddo
         enddo
      enddo
      do is=1,ndtype
         do iv=1,nvdist
c            write(*,*)iv,is
            fvave(iv,is)=fvave(iv,is)/float((it2-it1)*(ir2-ir1))
         enddo
      enddo
      do kk=1,nvmax
c         write(*,*)kk
         vdist(kk)=vrange*(float(kk-1)/(nvdist-1.)-0.495)
      enddo

      call autoplot(vdist,fvave(1,1),nvdist)
      call axlabels('v','f!dr!d(v!dr!d)')
      call autoplot(vdist,fvave(1,2),nvdist) 
      call axlabels('v','f!d!Af!@!d(v!d!Af!@!d)')
      call autoplot(vdist,fvave(1,3),nvdist)
      call axlabels('v','f!dz!d(v!dz!d)')

      end
