
c***********************************************************************
c Version 2.6 outputs two files. T...frc traces the force evolution.
c 2012 version appears to corrupt MPI somehow.

c*********************************************************************
c     Writes the main output file
      subroutine output(dt,damplen,itotsteps,fave,icolntype,colnwt)
c Common data:
      include 'piccom.f'
      include 'colncom.f'
      character*35 filename
c      integer iti,it2
      real totalf(nrfrcmax)

      call basefilename(icolntype,colnwt,filename)
      idf=nbcat(filename,'.dat')
c Write out averaged results.
      open(10,file=filename)
      write(10,'(a,a)')'  dt    vd     Ti     steps  rhoinf ' ,
     $       'phiinf  fave  debyelen Vp damplen  Bz...'
      write(10,'(2f7.4,f7.3,i5,f8.1,f7.3,f8.4,f8.3,f8.3,f6.2,f7.3,$)')
     $     dt,vd,Ti,i,rhoinf,log(rhoinf),fave,debyelen,vprobe,damplen,Bz
      if(icolntype.gt.0)then
         write(10,'(i2,f8.4)')icolntype,colnwt
      else
         write(10,*)
      endif
      write(10,*)NRUSED
      do j=1,NRUSED
         write(10,*)rcc(j),diagphi(j),diagrho(j)
c -log(rhoinf)
      enddo
      write(10,'(a)')'Number of steps, Particles to probe each step'
      write(10,*)itotsteps
      write(10,*)(fluxprobe(j),j=1,itotsteps)
      write(10,'(a)')'Number of theta cells, Number of steps, Particles'
      write(10,*)NTHUSED,itotsteps
      do j=1,NTHUSED
         ninth(j)=0
      enddo
      nastep=0
      do k=1,itotsteps
         write(10,*)(ninthstep(j,k),j=1,NTHUSED)
c     Just save the last quarter for the average
         if(k.gt.3*itotsteps/4)then
            nastep=nastep+1
            do j=1,NTHUSED
               ninth(j)=ninth(j)+ninthstep(j,k)
            enddo
         endif
      enddo
      write(10,'(a,a)')'Particle angular distrib summed over last'
     $     ,' quarter of steps, numbering:'
      write(10,*)nastep
      write(10,*)(ninth(j),j=1,NTHUSED)
      write(10,'(a,i4,i4)')'Mesh potential. Grid',NRUSED,NTHUSED
      call minmax2(phi(1,1),nrsize+1,NRUSED,NTHUSED,phimin,phimax)
      if(max(abs(phimin),abs(phimax)).lt.1.)then
         do j=1,NRUSED
            write(10,'(10f8.4)')(phi(j,k),k=1,NTHUSED)
         enddo
      else
         do j=1,NRUSED
            write(10,'(10f8.3)')(phi(j,k),k=1,NTHUSED)
         enddo
      endif
      write(10,'(a,i4,i4)')'Mesh density/infinity. Grid',NRUSED,NTHUSED
      do j=1,NRUSED
         write(10,'(10f8.3)')(rho(j,k),k=1,NTHUSED)
      enddo
      write(10,'(a,i4,i4)')'Volinv. Grid',NRUSED
      write(10,'(10f8.3)')(volinv(k),k=1,NRUSED)

      call outsums(dt,itotsteps+1)

c Output time-averages of z-force components stored in zmom(nstepmax,*,*).
c Particle units nTr^2, Electric nT lambda_D^2.
      do k=1,nrfrc
         totalf(k)=zmom(nstepmax,fieldz,k)*debyelen**2 +zmom(nstepmax
     $        ,epressz,k)+zmom(nstepmax,partz,k)+zmom(nstepmax,collision
     $        ,k)
      enddo
      write(10,*)nrfrc
      write(10,*)'Radius  Charge      E-field       Electrons',
     $     '      Ions     Coll     Total'
      do k=1,nrfrc
         write(10,'(f6.2,$)')r(izmrad(k))
c         write(10,'(i4,$)')izmrad(k)
         write(10,'(6g12.4)')(zmom(nstepmax,j,k),j=1,5),totalf(k)
      enddo
      write(10,*)'Collisions: Type,Weight,Eneutral,vneutral,Tneutral'
      write(10,701) icolntype,colnwt,Eneutral ,vneutral,Tneutral
      write(10,'(''rmtoz='',f10.4)')rmtoz
      write(10,*) 'Ion momentum collection at infinity'
      write(10,*) collmomtot(nstepmax)
      write(10,*) 'Energy flux to the probe'
      write(10,*) enertot(nstepmax)
c Trapped array.
      write(10,'(a,i4,i4)')'Trapped density. Grid',NRUSED,NTHUSED
      write(10,'(10f8.4)')((diagtrap(j,k),j=1,NRUSED),k=1,NTHUSED)
      
 701  format(10x,i3,4f10.5)
c     701  format('Collisions: type=',i4,' weight=',f8.4,' Eneutral=',
c     $     f10.5,' vneutral=',f8.4,' Tneutral=',f8.4)

c End of output file.
      close(10)

      innm=lentrim(filename)
      filename(innm-3:innm)='.frc'
      call outforce(filename,itotsteps)
      filename(innm-3:innm)='.vdg'
      if(ldist)call outvdiag(filename,itotsteps)
      filename(innm-3:innm)='.dst'
      if(nfvdist.gt.0)call outvdist(filename)
      end
c************************************************************************
c     Writes a txt file with the orbits of the traced particles
      subroutine orbitoutput()
c     Common data
      include 'piccom.f'

      character*30 filename
      integer iti,it2

c Construct a filename that contains many parameters

      write(filename,'(a)')'T'
      iti=nint(alog10(Ti)-0.49)
      it2=nint(Ti/10.**iti)
      write(filename(2:2),'(i1.1)')it2
      if(iti.lt.0) then
         filename(3:3)='m'
         iti=-iti
      else
         filename(3:3)='e'
      endif
      write(filename(4:4),'(i1.1)')iti

      filename(5:5)='v'
      write(filename(6:8),'(i3.3)')nint(100*vd)
      filename(9:9)='r'
      write(filename(10:11),'(i2.2)')ifix(r(nr))
      filename(12:12)='P'
      write(filename(13:14),'(i2.2)')ifix(abs(Vprobe))

      filename(15:15)='L'
      if(debyelen.gt.1.e-10)then
         iti=nint(alog10(debyelen)-0.49)
         it2=nint(debyelen/10.**iti)
      else
         it2=0
         iti=0
      endif
      write(filename(16:16),'(i1.1)')it2
      if(iti.lt.0) then
         filename(17:17)='m'
         iti=-iti
      else
         filename(17:17)='e'
      endif
      write(filename(18:18),'(i1.1)')iti

      filename(19:19)='B'
      if(Bz.gt.1.e-10)then
         iti=nint(alog10(Bz)-0.49)-1
         it2=nint(Bz/10.**iti)
      else
         it2=0
         iti=0
      endif
      write(filename(20:21),'(i2.2)')it2
      if(iti.lt.0) then
         filename(22:22)='m'
         iti=-iti
      else
         filename(22:22)='e'
      endif
      write(filename(23:23),'(i1.1)')iti

      filename(24:27)='.orb'

      open(15,file=filename)
      write(15,*) 'Number of orbits'
      write(15,*) norbits
      do k=1,norbits
         write(15,*) k,'th orbit'
         write(15,*) iorbitlen(k)
         do i=1,iorbitlen(k)
            write(15,590) xorbit(i,k),yorbit(i,k),zorbit(i,k),
     $      vxorbit(i,k),vyorbit(i,k),vzorbit(i,k)
         enddo
      enddo
 590  format(6f9.4)
      close(15)
      end
c*********************************************************************
c Write a second file with the force data as a function of step
      subroutine outforce(filename,itotsteps)
      character*(*) filename
      include 'piccom.f'
      real zmn(nstepmax,5,nrfrcmax)

c Apply normalization factors but don't change the zmom.
c Perhaps this extra storage is unnecessary.
      do i=1,itotsteps
         do k=1,nrfrc
            zmn(i,enccharge,k)=zmom(i,enccharge,k)
            zmn(i,partz,k)=zmom(i,partz,k)/rhoinf
            zmn(i,epressz,k)=zmom(i,epressz,k)
            zmn(i,fieldz,k)=zmom(i,fieldz,k)*debyelen**2
            zmn(i,collision,k)=zmom(i,collision,k)/rhoinf
         enddo
      enddo
      open(9,file=filename)
      write(9,*)'Forces may be plotted by e.g. plottraces -st6 filename'
      write(9,*)'Step    Charge     E-field      Electrons',
     $        '      Ions      Colns Total Force'
      do k=1,nrfrc
         write(9,'(a,i3,a,f8.3)')'legend: Sphere',k,' Radius'
     $        ,r(izmrad(k))
         write(9,*)itotsteps,6
         write(9,'(i5,6f12.5)')(i,(zmn(i,j,k),j=1,5),zmn(i,partz,k)
     $        +zmn(i,fieldz,k)+zmn(i,epressz,k)+zmn(i,collision,k),i=1
     $        ,itotsteps)
      enddo
      close(9)
      end
c*********************************************************************
c Write a third file with the averaged velocity distribution data
c from the single region.
      subroutine outvdiag(filename,itotsteps)
      character*(*) filename
      include 'piccom.f'

      open(9,file=filename)
      write(9,*)'Velocity Distributions. Use plottraces.'
      write(9,'(a,a)')'  dt    vd     Ti     steps  rhoinf ' ,
     $       'phiinf  fave  debyelen Vp damplen  Bz...'
      write(9,'(2f7.4,f7.3,i5,f8.1,f7.3,f8.4,f8.3,f8.3,f6.2,f7.3,$)')
     $     dt,vd,Ti,i,rhoinf,log(rhoinf),fave,debyelen,vprobe,damplen,Bz
      if(icolntype.gt.0)then
         write(9,'(i2,f8.4)')icolntype,colnwt
      else
         write(9,*)
      endif
      write(9,*)'v, vrdiagin, vtdiagin '
      write(9,'(a)')':-xv'
      write(9,'(a)')':-yf(v)'
      write(9,*)nvmax,2
      do k=1,nvmax
         write(9,*)vdiag(k),vrdiagin(k),vtdiagin(k)
      enddo
      write(9,*)
      close(9)
      end
c**********************************************************************
c V-distributions over the entire used cell array.
      subroutine outvdist(filename)
      character*(*) filename
      include 'piccom.f'
      include 'distcom.f'
      open(22,file=filename,status='unknown',err=101)
      close(22,status='delete')
      open(22,file=filename,status='new',form='unformatted',err=101)
      write(22)nvdist,5,nrused,nthused,vrange
      write(22)(rcc(ir),ir=1,nrused),(tcc(it),it=1,nthused)
c      write(22)(th(it),it=1,nthused)
c      write(22)(volinv(ir),ir=1,nrused)
      write(22)((((fvrtdist(iv,is,ir,it)
     $     ,iv=1,nvdist),is=1,5),ir=1,nrused),it=1,nthused)
      write(22)((rhodist(ir,it),ir=1,nrused),it=1,nthused)
      if(nthused.le.10.and.nrused.le.10)then
         write(*,*)'Written fvrtdist for',nvdist,5,nrused,nthused
         do it=1,nthused
            write(*,100)(rhodist(ir,it),ir=1,nrused)
         enddo
 100     format(10F8.3)
      endif
      return
 101  continue
      write(*,*)'Error opening for writing file:',filename
      end
c**********************************************************************
      subroutine readvdist(filename,isw)
      character*(*) filename
      include 'piccom.f'
      include 'distcom.f'
      open(22,file=filename,status='old',form='unformatted',err=101)
      read(22)nvd,nst,nrused,nthused,vrange
      if(nvd.ne.nvdist.or.nst.ne.5.or.
     $     nrused.gt.nrsize.or.nthused.gt.nthsize)then 
         write(*,*)'incorrect size parameters:',nvd,nst,nrused,nthused
      else
         if((isw-(isw/2)*2).eq.1)then
            write(*,*)'reading fvrtdist, parameters:',nvd,nst,nrused
     $           ,nthused
         endif
         read(22)(rcc(ir),ir=1,nrused),(tcc(it),it=1,nthused)
c         read(22)(th(it),it=1,nthused)
c         read(22)(volinv(ir),ir=1,nrused)
         read(22)((((fvrtdist(iv,is,ir,it)
     $        ,iv=1,nvdist),is=1,nst),ir=1,nrused),it=1,nthused)
c         write(*,*)'finished reading fvrtdist'
         read(22)((rhodist(ir,it),ir=1,nrused),it=1,nthused)
      endif
      return
 101  continue
      write(*,*)'Error opening for reading file:',filename
      end
c**********************************************************************
c Write out the particle data.
      subroutine partwrt(icolntype,colnwt)
c Common data:
      include 'piccom.f'
      include 'colncom.f'
      character*36 filename

c      write(filename,'(''part'',i3.3,''.dat'')')myid
      call basefilename(icolntype,colnwt,filename)
      write(filename(lentrim(filename)+1:),'(''.'',i3.3)')myid
c Delete the file first to help with nfs problems.
      open(11,file=filename,status='unknown')
      close(11,status='delete')
c
      open(11,file=filename,status='unknown')
      write(11,*)npartmax,npart,nr,nth,ndim,np
      write(11,*)((xp(i,j),i=1,ndim),j=1,npart)
      write(11,*)rhoinf,spotrein,averein
c      write(*,*)'rhoinf,spotrein,averein',rhoinf,spotrein,averein
      close(11)
      end

c**********************************************************************
c Read in the particle data.
      subroutine partrd(filename,success)
      logical success
      character*(*) filename
c Common data:
      include 'piccom.f'
c      include 'colncom.f'

c      write(filename,'(''part'',i3.3,''.dat'')')myid
c      write(*,*)'filename=',filename
      success=.false.
      open(11,file=filename,status='old',err=101)
      read(11,*,err=100,end=100)ipartmax,ipart,ir,ith,idim,ip
      if(ipartmax.eq.npartmax .and. ipart.eq.npart)then
c     $ .and. ir.eq.nr .and. ith.eq.nth )then
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
      return
 100  close(11)
      write(*,*) 'Error reading pardata.dat',filename
      return
 101  write(*,*) 'No particle file to read.',filename
      end
c**********************************************************************
c Get the average and slope over the rmesh range i1,i2.
      subroutine slopegen(phi,r,nr,i1,i2,slope,average)
      integer nr
      real phi(nr),r(nr)

c Assume r-mesh is linear
      rmom0=0.
      rmom1=0.
      rmid=(r(i2)+r(i1))/2.
      do i=i1,i2
         rmom0=rmom0+phi(i)
         rmom1=rmom1+(r(i)-rmid)*phi(i)
      enddo
      average=rmom0/(i2-i1+1)
c      rave=rmom1/(i2-i1+1)
      rlen=r(i2)-r(i1)
      slope=12.*(rmom1)/(rlen*rlen)/(i2-i1+1)
c      write(*,*)rmom0,rmom1,r(i1),r(i2),rmid,rlen
      end

c**********************************************************************
      subroutine outsums(dt,i)
c Common data:
      include 'piccom.f'
c Write out summed results.
      nrhere=NRUSED
      nthhere=NTHUSED
      nphere=np
c Combined files. Don't have to open.
c      open(10,file=filename)
      write(10,'(a,a)')
     $     '    dt       vd       Ti  steps  rmax',
     $     ' rhoinf debyelen Vp /nr,nth,np; sums'
      write(10,'(2f8.5,f8.4,i6,f8.3,f12.3,2f14.5)')
     $     dt,vd,Ti,i,r(nr),rhoinf,debyelen,vprobe
      write(10,*)nrhere,nthhere,nphere
      write(10,*)'psum'
      write(10,*)((psum(k1,k2),k1=1,nrhere),k2=1,nthhere)
      write(10,*)'vrsum'
      write(10,*)((vrsum(k1,k2),k1=1,nrhere),k2=1,nthhere)
      write(10,*)'vtsum'
      write(10,*)((vtsum(k1,k2),k1=1,nrhere),k2=1,nthhere)
      write(10,*)'vpsum'
      write(10,*)((vpsum(k1,k2),k1=1,nrhere),k2=1,nthhere)
      write(10,*)'v2sum'
      write(10,*)((v2sum(k1,k2),k1=1,nrhere),k2=1,nthhere)
      write(10,*)'vr2sum'
      write(10,*)((vr2sum(k1,k2),k1=1,nrhere),k2=1,nthhere)
      write(10,*)'vtp2sum'
      write(10,*)((vtp2sum(k1,k2),k1=1,nrhere),k2=1,nthhere)
      write(10,*)'vzsum'
      write(10,*)((vzsum(k1,k2),k1=1,nrhere),k2=1,nthhere)
      write(10,*)'diagvr'
      write(10,*)((diagvr(k1,k2),k1=1,nrhere),k2=1,nthhere)
      write(10,*)'r[cc]'
      write(10,*)(rcc(k1),k1=1,nrhere)
      write(10,*)'volinv'
      write(10,*)(volinv(k1),k1=1,nrhere)
      write(10,*)'t[cc]'
      write(10,*)(tcc(k2),k2=1,nthhere)
      end
c***************************************************************
      subroutine basefilename(icolntype,colnwt,filename)
      character*(*) filename
      include 'piccom.f'
      include 'colncom.f'

c Construct a filename that contains many parameters
c Using the routines in strings_names.f
      filename=' '
      call nameappendexp(filename,'T',Ti,1)
      call nameappendint(filename,'v',nint(100*vd),3)
      if(r(nr).ge.100)then
         call nameappendint(filename,'R',ifix(r(nr)/10.),2)
      else
         call nameappendint(filename,'r',ifix(r(nr)),2)
      endif
      call nameappendint(filename,'P',ifix(abs(Vprobe)),2)
      if (infdbl) then
         call nameappendexp(filename,'LI',debyelen,1)
      else
         call nameappendexp(filename,'L',debyelen,1)
      endif
      if(Bz.ne.0.) call nameappendexp(filename,'B',Bz,2)
      if(icolntype.eq.1.or.icolntype.eq.5)
     $     call nameappendexp(filename,'c',colnwt,1)
      if(icolntype.eq.2.or.icolntype.eq.6)
     $     call nameappendexp(filename,'C',colnwt,1)
     
      if(vneutral.ne.0)
     $     call nameappendint(filename,'N',nint(100*vneutral),3)

      end
c******************************************************************
