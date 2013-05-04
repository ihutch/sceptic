c Compare output from two files. They should have the same mesh
c parameters or else bad things will happen. 
c********************************************************************
      parameter (nfile=2)
      character*100 string,filename(nfile)
      include 'piccom.f'

      real rholocal(0:nrsize,0:nthsize)
      real comp(0:nrsize,0:nthsize,nfile+1)
c      real ctemp(0:nrsize)
c      real thanglocal(0:nthsize)
      common /forces/ charge1,ffield1,felec1,fion1,ftot1,
     $     charge2,ffield2,felec2,fion2,ftot2

      real phipic(1000),rhopic(1000),rhotrap(1000)
      real rpic(1000),rpicleft(1000),phicos(1000)
      integer nti0,nbsm
      parameter (nti0=100)
      integer jstepth,nrend
      logical lpcic,ltempc,lphip,lreaddiag,lgraph,larrows,lconline
      logical lvfwrt,lseptp,lunlabel,ledge,ldens,lnlog
      data lpcic/.false./
      data ltempc/.false./
      data lconline/.false./
      data lphip/.false./
      data lreaddiag/.false./
      data lgraph/.true./
      data larrows/.false./
      data lvfwrt/.true./
      data lseptp/.false./
      data lunlabel/.false./
      data ledge/.false./
      data ldens/.false./
      data lnlog/.false./
      data jstepth/1/nbsm/0/nrend/0/ifile/0/

      iboxcar=2
      idata=1
      expm1=exp(-1.)
c Deal with arguments
      do 1 i=1,iargc()
         call getarg(i,string)
C         write(*,*)'Argument:',string(1:40)
         if(string(1:1) .eq. ' ') then
            goto 3
         endif
         if(string(1:1) .eq. '-') then
         if(string(1:2) .eq. '-r') lreaddiag=.true.
         if(string(1:2) .eq. '-x') lgraph=.false.
c Legacy usage for summarize:         
         if(string(1:2) .eq. '-p') lgraph=.false.
         if(string(1:2) .eq. '-c') lpcic=.true.
         if(string(1:2) .eq. '-f') lphip=.true.
         if(string(1:2) .eq. '-n') ldens=.true.
         if(string(1:2) .eq. '-t') ltempc=.true.
         if(string(1:2) .eq. '-l') lconline=.true.
         if(string(1:2) .eq. '-a') larrows=.true.
         if(string(1:2) .eq. '-v') lvfwrt=.false.
         if(string(1:2) .eq. '-s') lseptp=.true.
         if(string(1:2) .eq. '-u') lunlabel=.true.
         if(string(1:2) .eq. '-e') ledge=.true.
         if(string(1:2) .eq. '-b') read(string(3:),*)iboxcar
         if(string(1:2) .eq. '-o') lnlog=.true.
         if(string(1:2) .eq. '-j') read(string(3:),*)jstepth
         if(string(1:2) .eq. '-k') read(string(3:),*)nrend
         if(string(1:2) .eq. '-m') read(string(3:),*)nbsm
         if(string(1:2) .eq. '-?') goto 51
         if(string(1:7) .eq. '--vrsum') idata=2
         if(string(1:7) .eq. '--vtsum') idata=3
         if(string(1:7) .eq. '--vpsum') idata=4
         if(string(1:7) .eq. '--v2sum') idata=5
         if(string(1:8) .eq. '--vr2sum') idata=6
         if(string(1:7) .eq. '--vtp2sum') idata=7
         if(string(1:6) .eq. '--psum') idata=10
         if(string(1:5) .eq. '--phi') idata=11
         else
            if(ifile.lt.nfile)ifile=ifile+1
            filename(ifile)=string
            write(*,*)ifile,filename(ifile)
         endif
 1    continue
 3    continue
      if(i.eq.1)goto 51


c Read the outputfile
      do i=1,ifile
         call readoutput(lreaddiag,lpcic,ledge,
     $     filename(i),rholocal,nrhere,nthhere,nphere,
     $     phipic,rhopic,rhotrap,rpic,rpicleft,phicos,
     $     rhomax,rhomin,
     $     nrti,phiinf,nastep,nsteps,
     $     dt,rmax,fave,
     $     icolntype,colnwt,Eneutral,vneutral,Tneutral,
     $     ierr)

         if(ierr.eq.101) goto 101

         v1=max(1.,vd)
         if(nrend.eq.0 .or. nrend.gt.nrhere)nrend=nrhere
         if(idata.eq.1)then
            call compare(i,vzsum,comp,nrend,nthhere,nfile,larrows
     $           ,lconline,v1,ir,jstepth,iboxcar)
            call congen(ir,jstepth,.3,-.3,
     $           nrend,nthhere,v1,larrows,lconline,comp(0,0,i),'vz')
         elseif(idata.eq.2)then
            call compare(i,vrsum,comp,nrend,nthhere,nfile,larrows
     $           ,lconline,v1,ir,jstepth,iboxcar)
            call congen(ir,jstepth,.3,-.3,
     $           nrend,nthhere,v1,larrows,lconline,comp(0,0,i),'vr')
         elseif(idata.eq.3)then
            call compare(i,vtsum,comp,nrend,nthhere,nfile,larrows
     $           ,lconline,v1,ir,jstepth,iboxcar)
            call congen(ir,jstepth,.3,-.3,
     $           nrend,nthhere,v1,larrows,lconline,comp(0,0,i),'vt')
         elseif(idata.eq.4)then
            call compare(i,vpsum,comp,nrend,nthhere,nfile,larrows
     $           ,lconline,v1,ir,jstepth,iboxcar)
            call congen(ir,jstepth,.3,-.3,
     $           nrend,nthhere,v1,larrows,lconline,comp(0,0,i),'vp')
         elseif(idata.eq.5)then
            call compare(i,v2sum,comp,nrend,nthhere,nfile,larrows
     $           ,lconline,v1,ir,jstepth,iboxcar)
            call congen(ir,jstepth,.3,0.,
     $           nrend,nthhere,v1,larrows,lconline,comp(0,0,i),'v2')
         elseif(idata.eq.11)then
            write(*,*)phi(10,10),phi(10,11)
            call compare(i,phi,comp,nrend,nthhere,nfile,larrows
     $           ,lconline,v1,ir,-jstepth,iboxcar)
            call congen(ir,jstepth,0.,-2.,
     $           nrend,nthhere,v1,larrows,lconline,comp(0,0,i),'phi')
         endif

      enddo

      i=3
      larrows=.true.
      call congen(ir,jstepth,.02,-.02, nrend,nthhere,v1,larrows,lconline
     $     ,comp(0,0,i),'Difference 2-1')

      

      call exit(0)
 101  continue
      write(*,*)'No file: ',filename(i)(1:50),ierr
 51   continue
      write(*,*)'Usage'
      write(*,*)'-knn specify the upper r-index to plot.'
      write(*,*
     $     )'--vrsum --vtsum --vpsum --v2sum --vr2sum --vtp2sum'
      write(*,*)'  specify value to compare.'  
      write(*,*)'-bnnn specify the radial boxcar averaging width (0-n)'
      write(*,*)'-jnnn specify angle boxcar ave width plus 1, jstepth'
      end

c***************************************************************************
c***************************************************************************
c Data reading subroutine
      subroutine readoutput(lreaddiag,lpcic,ledge,
     $     filename,rholocal,nrhere,nthhere,nphere,
     $     phipic,rhopic,rhotrap,rpic,rpicleft,phicos,
     $     rhomax,rhomin,
     $     nrti,phiinf,nastep,nsteps,
     $     dt,rmax,fave,
     $     icolntype,colnwt,Eneutral,vneutral,Tneutral,
     $     ierr)
      logical lreaddiag,lpcic,ledge
      character*100 string,filename
      real phipic(1000),rhopic(1000),rhotrap(1000)
      real rpic(1000),rpicleft(1000),phicos(1000)
      include 'piccom.f'
      real rholocal(0:nrsize,0:nthsize)
      character*256 charin
      common /forces/ charge1,ffield1,felec1,fion1,ftot1,
     $     charge2,ffield2,felec2,fion2,ftot2

      ierr=0
      nforcelines=2
c Read the data file:
c__________________________________________________________________

      open(10,file=filename,status='old',err=101)
c Line for nothing.
      read(10,*)charin
      read(10,'(a)')charin
c      write(*,*)charin
      read(charin,*,err=201,end=201)
     $     dt,vd,Ti,isteps,rhoinf,phiinf,fave,debyelen,vprobe,damplen,Bz
     $     ,icolntype,colnwt
 201  continue
      write(*,'(a,a)')'  dt    vd     Ti     steps  rhoinf ' ,
     $       'phiinf  fave  debyelen Vp damplen  Bz...'
      write(*,'(2f7.4,f7.3,i5,f8.1,f7.3,f8.4,f8.3,f8.3,f6.2,f7.3,$)')
     $     dt,vd,Ti,isteps,rhoinf,phiinf,fave,debyelen,vprobe,damplen,Bz
      if(icolntype.gt.0)then
         write(*,'(i2,f7.3)',err=212)icolntype,colnwt
      else
         write(*,*)
      endif
 212  read(10,*,err=202)nrTi
      nrhere=nrTi
c      write(*,*)'nrTi=',nrTi
      do i=1,nrTi
         read(10,*,err=203)rpic(i),phipic(i),diagrho(i)
         diagrho(i)=diagrho(i)/rhoinf
      enddo
      read(10,*)charin
c      write(*,*)charin(1:78)
      read(10,*,err=204)nsteps
c      write(*,*)nsteps
      if(nsteps.gt.nstepmax) then
         write(*,*)'Number of steps',nsteps,
     $        ' exceeds allocation',nstepmax
         call exit
      endif
      read(10,*)(fluxprobe(j),j=1,nsteps)
c Read theta cells
      read(10,*)charin
c      write(*,*)charin(1:78)
      read(10,*)nthhere,nsteps
c      write(*,*)nthhere,nsteps
      do i=1,nsteps
         read(10,*)(ninthstep(j,i),j=1,nthhere)
      enddo

      read(10,*)charin
c      write(*,*)charin(1:78)
      read(10,*)nastep
c      write(*,*)'nastep',nastep
      nread=nthhere
c This is not necessary here and indeed breaks the combined read because
c lpcic has not yet been set.
c      if(.not.lpcic)nread=nread+1
      read(10,*)(ninth(j),j=1,nread)
      if(lreaddiag)then
         write(*,*)'nastep',nastep,' ninth:'
         write(*,*)( ninth(j),j=1,nthhere)
      endif

      read(10,*)charin
      do j=1,nrhere
         read(10,*)(phi(j,k),k=1,nthhere)
      enddo
      read(10,*)charin
      do j=1,nrhere
         read(10,'(10f8.3)')(rholocal(j,k),k=1,nthhere)
      enddo
      read(10,*)charin
      read(10,*)(volinv(k),k=1,nrhere)
c We don't use open and close for combined files.
      if(filename(1:2).eq.'Ti')then
c But this must be using the split version.
         close(10)
c Read in  summed results.
         filename(1:2)='Sp'
         open(10,file=filename,err=210,status='old')
      endif
      read(10,'(a)')string
      read(10,*,err=200)
     $     dt,vd,Ti,i,rmax,rhoinf,debyelen,vprobe
 200  continue
c      write(*,*)'Reading Diagnostics'
      if(rhoinf.lt.1.)rhoinf=exp(phiinf)
      read(10,*)nrhere,nthhere,nphere
      if(nrhere.gt.nrsize .or. nthhere.gt.nthsize)then
         write(*,*)'Required dimensions: nr',nrhere,' nth',nthhere
         write(*,*)'are too large for the allocated values:'
     $        ,nrsize,nthsize
         stop
      endif
      read(10,*)string
      read(10,*)((psum(k1,k2),k1=1,nrhere),k2=1,nthhere)
      read(10,*)string
      read(10,*)((vrsum(k1,k2),k1=1,nrhere),k2=1,nthhere)
      read(10,*)string
      read(10,*)((vtsum(k1,k2),k1=1,nrhere),k2=1,nthhere)
      read(10,*)string
      read(10,*)((vpsum(k1,k2),k1=1,nrhere),k2=1,nthhere)
      read(10,*)string
      read(10,*)((v2sum(k1,k2),k1=1,nrhere),k2=1,nthhere)
      read(10,*)string
      read(10,*)((vr2sum(k1,k2),k1=1,nrhere),k2=1,nthhere)
      read(10,*)string
      read(10,*)((vtp2sum(k1,k2),k1=1,nrhere),k2=1,nthhere)
      read(10,*)string
      if(string(1:5).eq.'vzsum')then
         read(10,*)((vzsum(k1,k2),k1=1,nrhere),k2=1,nthhere)
         read(10,*)string
      endif
      read(10,*)((diagvr(k1,k2),k1=1,nrhere),k2=1,nthhere)
      read(10,*)string
      read(10,*)(rcc(k1),k1=1,nrhere)
      read(10,*)string
      read(10,*)(volinv(k1),k1=1,nrhere)
      read(10,*)string
      read(10,*)(tcc(k2),k2=1,nthhere)      
      read(10,*,err=410)nforcelines
      read(10,*,err=402,end=402)string
 410  read(10,*)charge1,ffield1,felec1,fion1,fcol1,ftot1
      read(10,*)charge2,ffield2,felec2,fion2,fcol2,ftot2
      do kk=1,nforcelines-2
         read(10,*)charge2,ffield2,felec2,fion2,fcol2,ftot2
      enddo
      read(10,*,err=402,end=402)string
      read(10,*,err=402,end=402)
     $     icolntype,colwt,Eneutral,vneutral,Tneutral
      read(10,*,err=402,end=402)string
      read(10,*,err=402,end=402)string
      read(10,*,err=402,end=402)pinfty
      read(10,*,err=402,end=402)string
      read(10,*,err=402,end=402)efprobe
      read(10,*,err=402,end=402)string
      read(10,*,err=402,end=402)((diagtrap(k1,k2),k1=1,nrhere),k2=1
     $     ,nthhere)
c      write(*,*)((diagtrap(k1,k2),k1=1,nrhere),k2=1
c     $     ,nthhere)

 402  close(10)
      write(*,*)'nrhere,nthhere,icolntype,colwt'
      write(*,*)nrhere,nthhere,icolntype,colwt,pinfty,efprobe
      write(*,*)'Final read string=',string(1:60)
      if(lreaddiag)then
         write(*,*)'Finished reading'
         write(*,*)'vrsum(1)'
         write(*,501)(vrsum(1,k2), k2=1,nthhere)
         write(*,*)'psum(1)'
         write(*,501)(psum(1,k2), k2=1,nthhere)
         write(*,*)'vr(1)'
         write(*,501)(vrsum(1,k2)/psum(1,k2), k2=1,nthhere)
         write(*,*)'diagvr(1)'
         write(*,501)(diagvr(1,k2), k2=1,nthhere)
         write(*,*)'rcc'
         write(*,*)(rcc(k1),k1=1,nrhere)
         write(*,*)'volinv'
         write(*,*)(volinv(k1),k1=1,nrhere)
         write(*,*)'tcc'
         write(*,*)(tcc(k2),k2=1,nthhere)
      endif
 501  format(10f8.3)
c__________________________________________________________________
c End of reading the data section

c__________________________________________________________________
c Fix up data 
      if(tcc(1).eq.1)lpcic=.true.
c Correct the outside angle centers if necessary.
      if(lpcic)then
         tcc(1)=0.25*(3.+tcc(2))
         tcc(0)=1.
         tcc(nthhere)=0.25*(-3.+tcc(nthhere-1))
         tcc(nthhere+1)=-1.

         thcells=nthhere-1.
      else
         thcells=nthhere
      endif


c In postproc, th is used as the weighting of the cells, essentially 
c the delta of cos theta that corresponds to each cell.
      do i=1,nthhere
         th(i)=2./(thcells)
      enddo
c The ends are half-weighted for cic.
      if(lpcic)then
         th(1)=0.5*th(1)
         th(nthhere)=0.5*th(nthhere)
      endif

c      write(*,*)'      tcc      th      thang','  lpcic=',lpcic
c      write(*, '(3f10.4)')(tcc(kk),th(kk),thang(kk),kk=0,nthhere+1)

c      call autoplot(th,diagvr,nthhere)

c Calculate rho(i,j) from the psum, volinv. Normalized by rhoinf.
      rhomax=0.
      do k1=1,nrhere
         do k2=1,nthhere
            rho(k1,k2)=psum(k1,k2)*volinv(k1)*thcells*nphere/rhoinf
c For now, use rholocal to fix --ds problem, maybe wrong ?
            if(rholocal(k1,k2).gt.rhomax)rhomax=rholocal(k1,k2)
         enddo
         if(lpcic)then
c fix up rho as double on boundary.
            rho(k1,1)=2.*rho(k1,1)
            rho(k1,nthhere)=2.*rho(k1,nthhere)
         endif
c fix angle ends of rho and phi
         rho(k1,0)=rho(k1,1)
         rho(k1,nthhere+1)=rho(k1,nthhere)
         phi(k1,0)=phi(k1,1)
         phi(k1,nthhere+1)=phi(k1,nthhere)
         rholocal(k1,0)=rholocal(k1,1)
         rholocal(k1,nthhere+1)=rholocal(k1,nthhere)
      enddo
      ir=10
      if(ledge)then
         rhomax=min(rhomax,1.5)
         rhomin=.5
      else
         rhomin=0.
      endif

      if(lreaddiag)then
         write(*,*)'rho   ','rholocal',' ratio',
     $     ' ;  rho is from psum, rholocal from Ti file'
         do i=1,nrhere
            write(*,*)rho(i,1),rholocal(i,1),rho(i,1)/rholocal(i,1)
         enddo
         write(*,*)'End of rho comparison'
      endif

      phiinf=0.

      jmin=1
      jmax=nthhere
      if(lreaddiag)write(*,*)'jmin,jmax',jmin,jmax
      do i=1,nrhere
         rhopic(i)=0.
         rhotrap(i)=0.
         phicos(i)=0.
c         write(*,*)'th   tcc   phi'
         do j=jmin,jmax
            rhopic(i)=rhopic(i)+rholocal(i,j)
            rhotrap(i)=rhotrap(i)+diagtrap(i,j)
c            if(diagtrap(i,j).eq.0)write(*,*)i,j
c     rhopic(i)=rhopic(i)+rho(i,j)
c \int cos(\theta) \phi(\theta) d\cos(\theta)
            phicos(i)=phicos(i)+th(j)*tcc(j)*phi(i,j)
c            write(*,'(4f10.4)')th(j),tcc(j),phi(i,j),phicos(i)
         enddo
         rhopic(i)=rhopic(i)/float(jmax-jmin+1)
         rhotrap(i)=rhotrap(i)/float(jmax-jmin+1)+0.01
c         write(*,*)diagtrap(i,3),rhotrap(i)
      enddo
c     rescale rho; but usually this is the identity transformation.
      do i=1,nrhere
         rpicleft(i)=-rpic(i)
      enddo

      return
c End of data fix-up section
c__________________________________________________________________
 202  write(*,*)"nr error"
      call exit
 203  write(*,*)"rpicphipic error"
      call exit
 204  write(*,*)"nsteps error"
      call exit
 210  write(*,*)'Error opening file: ',filename(1:50)
      call exit
 101  ierr=101
      end
c************************************************************************

c***************************************************************************
c Contouring of any quantity, rholocal, on distorted mesh.
c But the piccompost gives us some stuff that we need.
      subroutine congen(ir,it,rhomax,rhomin,nrhere,nthhere,v1,
     $     larrows,lconline,rholocal,labelstring)
      integer ir,it
      real rhomax,v1
      logical larrows,lconline
      character*(*) labelstring
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
      save xrho,zrho
      real basesize
      parameter (basesize=.02)

      write(*,*)'congen',ir,it,rhomax,rhomin,nrhere,nthhere,v1

      if(rhomin.gt.0)call setconlog(.true.)
      if(nthhere.gt.nthsize)then
         write(*,*)' Congen error. Mesh required:',nrhere,nthhere,
     $        ' Exceeds allocated:',nrsize,nthsize
         stop
      endif
c Correct the outside angle centers if necessary.
      if(tcc(1).eq.1)then
         tcc(1)=0.25*(3.+tcc(2))
         tcc(0)=1.
         tcc(nthhere)=0.25*(-3.+tcc(nthhere-1))
         tcc(nthhere+1)=-1.
      endif

      do i=1,nrhere
          do j=1,nthhere
            zrho(i,j)=rcc(i)*tcc(j)
            xrho(i,j)=rcc(i)*sqrt(1.-tcc(j)**2)
         enddo
         zrho(i,0)=rcc(i)
         xrho(i,0)=0.
         zrho(i,nthhere+1)=-rcc(i)
         xrho(i,nthhere+1)=0.         
      enddo

c      rpmax=2.*rcc(nrhere)-rcc(nrhere-1)
      rpmax=rcc(nrhere)
c      write(*,*)rcc(nrhere),nrhere,nthhere,rpmax


      do j=1,ncont
c         zclv(j)=rhomin+(rhomax-rhomin)*(0.95*(j-1)/float(ncont-1))
         zclv(j)=(rhomax-rhomin)*(1.*(j-1)/float(ncont-1))+rhomin
      enddo
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
            call contourl(rholocal(1,0),cworka,nrsize+1,nrhere,
     $           nthhere+2,zclv,icl,zrho,xrho,ntype)
            write(*,'(a,30f5.2)')'Contours=',zclv
            call fwrite(delta,iwd,1,cstring)
            tstring=' contour spacing: '//cstring(1:10)
c            call legendline(-.1,-.22,258,tstring)
         endif
c      endif
c Fit closer than boxtitle
         call legendline(0.4,1.07,258,labelstring)
         call drcstr(tstring)
         call color(15)
         call axis()
         call axlabels('z','r sin!Aq!@')
  
      
      if(larrows) then
         call color(12)
         if(ir.le.0.or.ir.ge.100) ir=10
         do j=1,nthhere,it
            do i=1,nrhere,max(nrhere/ir,1)
               vri=vrsum(i,j)/(psum(i,j)+1.e-5)
               vti=vtsum(i,j)/(psum(i,j)+1.e-5)
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
      
      call pltend()
      call setconlog(.false.)
      end

c***********************************************************************
c Get the data from array into comp(.,.,i) and if i.gt.1 generate
c the next higher comp as being the difference, boxcar averaged.
c Normalize it if jstepth positive.
      subroutine compare(i,array,comp,nrhere,nthhere,nfile,larrows
     $     ,lconline,v1,ir,jstepth,iboxcar)
      include 'piccom.f'
      real ctemp(0:nrsize)
      real ctempt(0:nthsize),ctempt1(0:nthsize)
      real array(0:nrsize,0:nthsize)
      real comp(0:nrsize,0:nthsize,nfile+1)
      logical larrows,lconline
      real v1
      integer i,ir,jstepth


      do k=0,nthhere+1
         do j=0,nrhere
            if(jstepth.ge.0)then
               comp(j,k,i)=array(j,k)/max(psum(j,k),1.)
            else
               comp(j,k,i)=array(j,k)
            endif
            if(i.gt.1)then
               ctemp(j)=comp(j,k,i)-comp(j,k,i-1)
            endif
         enddo
         if(i.gt.1)then
            call boxcarave(nrhere+1,iboxcar,ctemp,comp(0,k,i+1))
         endif
      enddo
      if(i.gt.1.and.abs(jstepth-1).gt.0)then
         do j=1,nrhere
            do k=1,nthhere+1
               ctempt(k)=comp(j,k,i+1)
            enddo
            call boxcarave(nthhere+1,abs(jstepth-1),ctempt,ctempt1)
            do k=1,nthhere+1
               comp(j,k,i+1)=ctempt1(k)
            enddo
         enddo
      endif


      end
