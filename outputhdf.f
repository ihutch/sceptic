c*********************************************************************
c Writes the hdf output file
c Note that this version was hacked from the sceptic3D version, so it
c   is not necessarily correct nor complete
      subroutine outputhdf(dt,k,fave,icolntype,colnwt)
c Load the hdf5 module
      use hdf5
c Don't allow implicit definitions
      implicit none
c Input variables
      real dt, fave, colnwt
      integer k, icolntype
c Common data
      include 'piccom.f'
      include 'colncom.f'
c Functions
      integer nbcat
c Local variables
c     File, group, and dataset names
      CHARACTER(LEN=55) :: filename
      CHARACTER(LEN=55) :: groupname
      CHARACTER(LEN=55) :: dsetname
c     File and group ids
      INTEGER(HID_T) :: file_id
      INTEGER(HID_T) :: group_id
c     Error flag
      INTEGER :: error
c     Dimensions of dataset (limit rank to <=10)
      INTEGER(HSIZE_T), DIMENSION(10) :: data_dims, storage_dims
c     Rank of dataset
      INTEGER :: rank
c     Working variables
      INTEGER :: i, j
      integer idf
      REAL, DIMENSION(3,3) :: dset_data1
      INTEGER, DIMENSION(3,3) :: dset_data2

c Construct a filename that contains many parameters
c   using the routines in strings_names.f
      filename=' '
      call nameappendexp(filename,'T',Ti,1)
      call nameappendint(filename,'v',nint(100*vd),3)
      if(r(nr).ge.100)then
         call nameappendint(filename,'R',ifix(r(nr)/10.),2)
      else
         call nameappendint(filename,'r',ifix(r(nr)),2)
      endif
      call nameappendint(filename,'P',ifix(abs(Vprobe)),3)
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
     $    call nameappendint(filename,'N',nint(100*vneutral),3)

      if(bcr.eq.3)
     $    call nameappendint(filename,'bcr',bcr,1)

c     Turn on these when testing convergence with changing grid size
cc      call nameappendint(filename,'Nr',nrused,3)
cc      call nameappendint(filename,'Nt',nthused,3)

      idf=nbcat(filename,'2D.h5')



      do i = 1, 3
           do j = 1, 3
                dset_data1(i,j) = j*1.1;
                dset_data2(i,j) = j;
           end do
      end do

c Initialize FORTRAN interface.
      CALL h5open_f(error)

c Create a new file using default properties.
      CALL h5fcreate_f(filename, H5F_ACC_TRUNC_F, file_id, error)

c Now output each desired common block as an hdf group, using the two
c   custom made functions writehdfintmat() and writehfdrealmat() to
c   output integer and real (respectively) variables and arrays.
c   The logical type is treated as integer for output.
c   For now all variables in each desired common block are output, but it may
c   be desireable to create a more limited output file to minimize size.
c   Note that size of each array must be specified throught the rank,
c   dimensions of the used array (data_dims), and dimensions of the
c   full arary in memory (storage_dims).



c Create a group in the file.
      groupname = 'piccom'
      CALL h5gcreate_f(file_id, groupname, group_id, error)

c Write the data in this group

c     Single value variables
      rank = 1
      data_dims(1) = 1
      storage_dims(1) = 1

      dsetname = 'npart'
      call writehdfintmat(group_id,dsetname,
     $ npart,storage_dims,data_dims,rank)

      dsetname = 'nr'
      call writehdfintmat(group_id,dsetname,
     $ nr,storage_dims,data_dims,rank)

      dsetname = 'nth'
      call writehdfintmat(group_id,dsetname,
     $ nth,storage_dims,data_dims,rank)

      dsetname = 'lsubcycle'
      call writehdfintmat(group_id,dsetname,
     $ lsubcycle,storage_dims,data_dims,rank)

      dsetname = 'verlet'
      call writehdfintmat(group_id,dsetname,
     $ verlet,storage_dims,data_dims,rank)

      dsetname = 'LCIC'
      call writehdfintmat(group_id,dsetname,
     $ LCIC,storage_dims,data_dims,rank)

      dsetname = 'NRUSED'
      call writehdfintmat(group_id,dsetname,
     $ NRUSED,storage_dims,data_dims,rank)

      dsetname = 'NTHUSED'
      call writehdfintmat(group_id,dsetname,
     $ NTHUSED,storage_dims,data_dims,rank)

      dsetname = 'NRFULL'
      call writehdfintmat(group_id,dsetname,
     $ NRFULL,storage_dims,data_dims,rank)

      dsetname = 'NTHFULL'
      call writehdfintmat(group_id,dsetname,
     $ NTHFULL,storage_dims,data_dims,rank)

      dsetname = 'ninjcomp'
      call writehdfintmat(group_id,dsetname,
     $ ninjcomp,storage_dims,data_dims,rank)

      dsetname = 'iocprev'
      call writehdfintmat(group_id,dsetname,
     $ iocprev,storage_dims,data_dims,rank)

      dsetname = 'cerr'
      call writehdfrealmat(group_id,dsetname,
     $ cerr,storage_dims,data_dims,rank)

      dsetname = 'bdyfc'
      call writehdfrealmat(group_id,dsetname,
     $ bdyfc,storage_dims,data_dims,rank)

      dsetname = 'Ti'
      call writehdfrealmat(group_id,dsetname,
     $ Ti,storage_dims,data_dims,rank)

      dsetname = 'vd'
      call writehdfrealmat(group_id,dsetname,
     $ vd,storage_dims,data_dims,rank)

      dsetname = 'Bz'
      call writehdfrealmat(group_id,dsetname,
     $ Bz,storage_dims,data_dims,rank)

      dsetname = 'diags'
      call writehdfintmat(group_id,dsetname,
     $ diags,storage_dims,data_dims,rank)

      dsetname = 'lplot'
      call writehdfintmat(group_id,dsetname,
     $ lplot,storage_dims,data_dims,rank)

      dsetname = 'ldist'
      call writehdfintmat(group_id,dsetname,
     $ ldist,storage_dims,data_dims,rank)

      dsetname = 'linsulate'
      call writehdfintmat(group_id,dsetname,
     $ linsulate,storage_dims,data_dims,rank)

      dsetname = 'lfloat'
      call writehdfintmat(group_id,dsetname,
     $ lfloat,storage_dims,data_dims,rank)

      dsetname = 'lat0'
      call writehdfintmat(group_id,dsetname,
     $ lat0,storage_dims,data_dims,rank)

      dsetname = 'localinj'
      call writehdfintmat(group_id,dsetname,
     $ localinj,storage_dims,data_dims,rank)

      dsetname = 'lfixedn'
      call writehdfintmat(group_id,dsetname,
     $ lfixedn,storage_dims,data_dims,rank)

      dsetname = 'myid'
      call writehdfintmat(group_id,dsetname,
     $ myid,storage_dims,data_dims,rank)

      dsetname = 'numprocs'
      call writehdfintmat(group_id,dsetname,
     $ numprocs,storage_dims,data_dims,rank)

      dsetname = 'rmtoz'
      call writehdfrealmat(group_id,dsetname,
     $ rmtoz,storage_dims,data_dims,rank)

c     Variable arrays
      dsetname = 'xp'
      rank = 2
      data_dims(1) = ndim
      data_dims(2) = npart
      storage_dims(1) = ndim
      storage_dims(2) = npartmax
      call writehdfrealmat(group_id,dsetname,
     $  xp,storage_dims,data_dims,rank)

      dsetname = 'vzinit'
      rank = 1
      data_dims(1) = npart
      storage_dims(1) = npartmax
      call writehdfrealmat(group_id,dsetname,
     $  vzinit,storage_dims,data_dims,rank)

      dsetname = 'dtprec'
      rank = 1
      data_dims(1) = npart
      storage_dims(1) = npartmax
      call writehdfrealmat(group_id,dsetname,
     $  dtprec,storage_dims,data_dims,rank)

      dsetname = 'ipf'
      rank = 1
      data_dims(1) = npart
      storage_dims(1) = npartmax
      call writehdfintmat(group_id,dsetname,
     $  ipf,storage_dims,data_dims,rank)

      dsetname = 'phi'
      rank = 2
      data_dims(1) = nr+1
      data_dims(2) = nth+1
      storage_dims(1) = nrsize+1
      storage_dims(2) = nthsize+1
      call writehdfrealmat(group_id,dsetname,
     $  phi,storage_dims,data_dims,rank)

      dsetname = 'rho'
      rank = 2
      data_dims(1) = nr+1
      data_dims(2) = nth+1
      storage_dims(1) = nrsize+1
      storage_dims(2) = nthsize+1
      call writehdfrealmat(group_id,dsetname,
     $  rho,storage_dims,data_dims,rank)

c Close the group.
      CALL h5gclose_f(group_id, error)


c Create a group in the file.
      groupname = 'momcom'
      CALL h5gcreate_f(file_id, groupname, group_id, error)

c Write the data in this group

c     Variable arrays
      rank = 2
      data_dims(1) = nr-1
      data_dims(2) = nth-1
      storage_dims(1) = nrsize-1
      storage_dims(2) = nthsize-1

      dsetname = 'psum'
      call writehdfrealmat(group_id,dsetname,
     $  psum,storage_dims,data_dims,rank)

      dsetname = 'vrsum'
      call writehdfrealmat(group_id,dsetname,
     $  vrsum,storage_dims,data_dims,rank)

      dsetname = 'vtsum'
      call writehdfrealmat(group_id,dsetname,
     $  vtsum,storage_dims,data_dims,rank)

      dsetname = 'vpsum'
      call writehdfrealmat(group_id,dsetname,
     $  vpsum,storage_dims,data_dims,rank)

      dsetname = 'vr2sum'
      call writehdfrealmat(group_id,dsetname,
     $  vr2sum,storage_dims,data_dims,rank)

      dsetname = 'vzsum'
      call writehdfrealmat(group_id,dsetname,
     $  vzsum,storage_dims,data_dims,rank)

c Close the group.
      CALL h5gclose_f(group_id, error)


c Create a group in the file.
      groupname = 'meshcom'
      CALL h5gcreate_f(file_id, groupname, group_id, error)

c Write the data in this group

c     Single value variables
      rank = 1
      data_dims(1) = 1
      storage_dims(1) = 1

      dsetname = 'rfac'
      call writehdfrealmat(group_id,dsetname,
     $ rfac,storage_dims,data_dims,rank)

      dsetname = 'tfac'
      call writehdfrealmat(group_id,dsetname,
     $ tfac,storage_dims,data_dims,rank)

      dsetname = 'avelim'
      call writehdfrealmat(group_id,dsetname,
     $ avelim,storage_dims,data_dims,rank)

      dsetname = 'idim1'
      call writehdfintmat(group_id,dsetname,
     $ idim1,storage_dims,data_dims,rank)

      dsetname = 'idim2'
      call writehdfintmat(group_id,dsetname,
     $ idim2,storage_dims,data_dims,rank)

c     Variable arrays
      dsetname = 'r'
      rank = 1
      data_dims(1) = nr+1
      storage_dims(1) = nrsize+1
      call writehdfrealmat(group_id,dsetname,
     $  r,storage_dims,data_dims,rank)

      dsetname = 'rcc'
      rank = 1
      data_dims(1) = nr+1
      storage_dims(1) = nrsize+1
      call writehdfrealmat(group_id,dsetname,
     $  rcc,storage_dims,data_dims,rank)

      dsetname = 'th'
      rank = 1
      data_dims(1) = nth+1
      storage_dims(1) = nthsize+1
      call writehdfrealmat(group_id,dsetname,
     $  th,storage_dims,data_dims,rank)

      dsetname = 'tcc'
      rank = 1
      data_dims(1) = nth+1
      storage_dims(1) = nthsize+1
      call writehdfrealmat(group_id,dsetname,
     $  tcc,storage_dims,data_dims,rank)

      dsetname = 'thang'
      rank = 1
      data_dims(1) = nth+1
      storage_dims(1) = nthsize+1
      call writehdfrealmat(group_id,dsetname,
     $  thang,storage_dims,data_dims,rank)

      dsetname = 'volinv'
      rank = 1
      data_dims(1) = nr+1
      storage_dims(1) = nrsize+1
      call writehdfrealmat(group_id,dsetname,
     $  volinv,storage_dims,data_dims,rank)

      dsetname = 'irpre'
      rank = 1
      data_dims(1) = 4*nr
      storage_dims(1) = nrpre
      call writehdfintmat(group_id,dsetname,
     $  irpre,storage_dims,data_dims,rank)

      dsetname = 'itpre'
      rank = 1
      data_dims(1) = 4*nth
      storage_dims(1) = ntpre
      call writehdfintmat(group_id,dsetname,
     $  itpre,storage_dims,data_dims,rank)

      dsetname = 'hr'
      rank = 1
      data_dims(1) = nr+2
      storage_dims(1) = nrsize+2
      call writehdfrealmat(group_id,dsetname,
     $  hr,storage_dims,data_dims,rank)

      dsetname = 'zeta'
      rank = 1
      data_dims(1) = nr+2
      storage_dims(1) = nrsize+2
      call writehdfrealmat(group_id,dsetname,
     $  zeta,storage_dims,data_dims,rank)

      dsetname = 'zetahalf'
      rank = 1
      data_dims(1) = nr+2
      storage_dims(1) = nrsize+2
      call writehdfrealmat(group_id,dsetname,
     $  zetahalf,storage_dims,data_dims,rank)

      dsetname = 'cminus'
      rank = 1
      data_dims(1) = nr
      storage_dims(1) = nrsize
      call writehdfrealmat(group_id,dsetname,
     $  cminus,storage_dims,data_dims,rank)

      dsetname = 'cmid'
      rank = 1
      data_dims(1) = nr
      storage_dims(1) = nrsize
      call writehdfrealmat(group_id,dsetname,
     $  cmid,storage_dims,data_dims,rank)

      dsetname = 'cplus'
      rank = 1
      data_dims(1) = nr
      storage_dims(1) = nrsize
      call writehdfrealmat(group_id,dsetname,
     $  cplus,storage_dims,data_dims,rank)

c Close the group.
      CALL h5gclose_f(group_id, error)


c Create a group in the file.
      groupname = 'rancom'
      CALL h5gcreate_f(file_id, groupname, group_id, error)

c Write the data in this group

c     Single value variables
      rank = 1
      data_dims(1) = 1
      storage_dims(1) = 1

      dsetname = 'bcphi'
      call writehdfintmat(group_id,dsetname,
     $ bcphi,storage_dims,data_dims,rank)

      dsetname = 'bcr'
      call writehdfintmat(group_id,dsetname,
     $ bcr,storage_dims,data_dims,rank)

      dsetname = 'infdbl'
      call writehdfintmat(group_id,dsetname,
     $ infdbl,storage_dims,data_dims,rank)

c     Variable arrays
      dsetname = 'Qcom'
      rank = 1
      data_dims(1) = nQth
      storage_dims(1) = nQth
      call writehdfrealmat(group_id,dsetname,
     $  Qcom,storage_dims,data_dims,rank)

      dsetname = 'Gcom'
      rank = 2
      data_dims(1) = nvel
      data_dims(2) = nQth
      storage_dims(1) = nvel
      storage_dims(2) = nQth
      call writehdfrealmat(group_id,dsetname,
     $  Gcom,storage_dims,data_dims,rank)

      dsetname = 'Vcom'
      rank = 1
      data_dims(1) = nvel
      storage_dims(1) = nvel
      call writehdfrealmat(group_id,dsetname,
     $  Vcom,storage_dims,data_dims,rank)

      dsetname = 'pu1'
      rank = 1
      data_dims(1) = nvel
      storage_dims(1) = nvel
      call writehdfrealmat(group_id,dsetname,
     $  pu1,storage_dims,data_dims,rank)

      dsetname = 'pu2'
      rank = 1
      data_dims(1) = nQth
      storage_dims(1) = nQth
      call writehdfrealmat(group_id,dsetname,
     $  pu2,storage_dims,data_dims,rank)

      dsetname = 'Pc'
      rank = 2
      data_dims(1) = nQth
      data_dims(2) = nvel
      storage_dims(1) = nQth
      storage_dims(2) = nvel
      call writehdfrealmat(group_id,dsetname,
     $  Pc,storage_dims,data_dims,rank)

c Close the group.
      CALL h5gclose_f(group_id, error)


c Create a group in the file.
      groupname = 'diagcom'
      CALL h5gcreate_f(file_id, groupname, group_id, error)

c Write the data in this group

c     Single value variables
      rank = 1
      data_dims(1) = 1
      storage_dims(1) = 1

      dsetname = 'nrein'
      call writehdfintmat(group_id,dsetname,
     $ nrein,storage_dims,data_dims,rank)

      dsetname = 'nreintry'
      call writehdfintmat(group_id,dsetname,
     $ nreintry,storage_dims,data_dims,rank)

      dsetname = 'ninner'
      call writehdfintmat(group_id,dsetname,
     $ ninner,storage_dims,data_dims,rank)

      dsetname = 'vrange'
      call writehdfrealmat(group_id,dsetname,
     $ vrange,storage_dims,data_dims,rank)

      dsetname = 'zmomprobe'
      call writehdfrealmat(group_id,dsetname,
     $ zmomprobe,storage_dims,data_dims,rank)

      dsetname = 'enerprobe'
      call writehdfrealmat(group_id,dsetname,
     $ enerprobe,storage_dims,data_dims,rank)

      dsetname = 'zmout'
      call writehdfrealmat(group_id,dsetname,
     $ zmout,storage_dims,data_dims,rank)

      dsetname = 'spotrein'
      call writehdfrealmat(group_id,dsetname,
     $ spotrein,storage_dims,data_dims,rank)

      dsetname = 'averein'
      call writehdfrealmat(group_id,dsetname,
     $ averein,storage_dims,data_dims,rank)

      dsetname = 'fluxrein'
      call writehdfrealmat(group_id,dsetname,
     $ fluxrein,storage_dims,data_dims,rank)

      dsetname = 'rhoinf'
      call writehdfrealmat(group_id,dsetname,
     $ rhoinf,storage_dims,data_dims,rank)

      dsetname = 'ntrapre'
      call writehdfintmat(group_id,dsetname,
     $ ntrapre,storage_dims,data_dims,rank)

      dsetname = 'adeficit'
      call writehdfrealmat(group_id,dsetname,
     $ adeficit,storage_dims,data_dims,rank)

      dsetname = 'ircell'
      call writehdfintmat(group_id,dsetname,
     $ ircell,storage_dims,data_dims,rank)

      dsetname = 'itcell'
      call writehdfintmat(group_id,dsetname,
     $ itcell,storage_dims,data_dims,rank)

c     Variable arrays
      rank = 1
      data_dims(1) = nvmax
      storage_dims(1) = nvmax
      dsetname = 'nvdiag'
      call writehdfrealmat(group_id,dsetname,
     $  nvdiag,storage_dims,data_dims,rank)

      dsetname = 'nvdiagave'
      call writehdfrealmat(group_id,dsetname,
     $  nvdiagave,storage_dims,data_dims,rank)

      dsetname = 'vdiag'
      call writehdfrealmat(group_id,dsetname,
     $  vdiag,storage_dims,data_dims,rank)

      dsetname = 'vrdiagin'
      call writehdfrealmat(group_id,dsetname,
     $  vrdiagin,storage_dims,data_dims,rank)

      dsetname = 'vtdiagin'
      call writehdfrealmat(group_id,dsetname,
     $  vtdiagin,storage_dims,data_dims,rank)

      rank = 1
      data_dims(1) = nr
      storage_dims(1) = nrsize
      dsetname = 'diagrho'
      call writehdfrealmat(group_id,dsetname,
     $  diagrho,storage_dims,data_dims,rank)

      dsetname = 'diagphi'
      call writehdfrealmat(group_id,dsetname,
     $  diagphi,storage_dims,data_dims,rank)

      dsetname = 'diagchi'
      rank = 1
      data_dims(1) = nth+1
      storage_dims(1) = nthsize+1
      call writehdfrealmat(group_id,dsetname,
     $  diagchi,storage_dims,data_dims,rank)

      dsetname = 'fluxprobe'
      rank = 1
      data_dims(1) = maxsteps
      storage_dims(1) = nstepmax
      call writehdfrealmat(group_id,dsetname,
     $  fluxprobe,storage_dims,data_dims,rank)

      dsetname = 'zmom'
      rank = 3
      data_dims(1) = maxsteps
      data_dims(2) = 5
      data_dims(3) = 2
      storage_dims(1) = nstepmax
      storage_dims(2) = 5
      storage_dims(3) = 2
      call writehdfrealmat(group_id,dsetname,
     $  zmom,storage_dims,data_dims,rank)

      dsetname = 'enertot'
      rank = 1
      data_dims(1) = maxsteps
      storage_dims(1) = nstepmax
      call writehdfrealmat(group_id,dsetname,
     $  enertot,storage_dims,data_dims,rank)

c Close the group.
      CALL h5gclose_f(group_id, error)


c Create a group in the file.
      groupname = 'poisson'
      CALL h5gcreate_f(file_id, groupname, group_id, error)

c Write the data in this group

c     Single value variables
      rank = 1
      data_dims(1) = 1
      storage_dims(1) = 1

      dsetname = 'debyelen'
      call writehdfrealmat(group_id,dsetname,
     $ debyelen,storage_dims,data_dims,rank)

c     Variable arrays
      rank = 1
      data_dims(1) = nr+1
      storage_dims(1) = nrsize+1
      dsetname = 'apc'
      call writehdfrealmat(group_id,dsetname,
     $  apc,storage_dims,data_dims,rank)

      dsetname = 'bpc'
      call writehdfrealmat(group_id,dsetname,
     $  bpc,storage_dims,data_dims,rank)

      rank = 2
      data_dims(1) = nr+1
      data_dims(2) = nth+1
      storage_dims(1) = nrsize+1
      storage_dims(2) = nthsize+1
      dsetname = 'cpc'
      call writehdfrealmat(group_id,dsetname,
     $  cpc,storage_dims,data_dims,rank)

      dsetname = 'dpc'
      call writehdfrealmat(group_id,dsetname,
     $  dpc,storage_dims,data_dims,rank)

      dsetname = 'fpc'
      call writehdfrealmat(group_id,dsetname,
     $  fpc,storage_dims,data_dims,rank)

      dsetname = 'gpc'
      rank = 2
      data_dims(1) = nth+1
      data_dims(3) = 5
      storage_dims(1) = nthsize+1
      storage_dims(3) = 5
      call writehdfrealmat(group_id,dsetname,
     $  gpc,storage_dims,data_dims,rank)

c Close the group.
      CALL h5gclose_f(group_id, error)


c Create a group in the file.
      groupname = 'stepave'
      CALL h5gcreate_f(file_id, groupname, group_id, error)

c Write the data in this group

c     Single value variables
      rank = 1
      data_dims(1) = 1
      storage_dims(1) = 1

      dsetname = 'nstepsave'
      call writehdfintmat(group_id,dsetname,
     $ nstepsave,storage_dims,data_dims,rank)

      dsetname = 'nsamax'
      call writehdfintmat(group_id,dsetname,
     $ nsamax,storage_dims,data_dims,rank)

c Close the group.
      CALL h5gclose_f(group_id, error)


c Create a group in the file.
      groupname = 'orbits'
      CALL h5gcreate_f(file_id, groupname, group_id, error)

c Write the data in this group

c     Single value variables
      rank = 1
      data_dims(1) = 1
      storage_dims(1) = 1

      dsetname = 'norbits'
      call writehdfintmat(group_id,dsetname,
     $ norbits,storage_dims,data_dims,rank)

c     Variable arrays
      if (norbits.gt.0) then

      rank = 2
      data_dims(1) = maxsteps
      data_dims(2) = norbits
      storage_dims(1) = nstepmax
      storage_dims(2) = nobsmax
      dsetname = 'xorbit'
      call writehdfrealmat(group_id,dsetname,
     $  xorbit,storage_dims,data_dims,rank)

      dsetname = 'yorbit'
      call writehdfrealmat(group_id,dsetname,
     $  yorbit,storage_dims,data_dims,rank)

      dsetname = 'zorbit'
      call writehdfrealmat(group_id,dsetname,
     $  zorbit,storage_dims,data_dims,rank)

      dsetname = 'vxorbit'
      call writehdfrealmat(group_id,dsetname,
     $  vxorbit,storage_dims,data_dims,rank)

      dsetname = 'vyorbit'
      call writehdfrealmat(group_id,dsetname,
     $  vyorbit,storage_dims,data_dims,rank)

      dsetname = 'vzorbit'
      call writehdfrealmat(group_id,dsetname,
     $  vzorbit,storage_dims,data_dims,rank)

      dsetname = 'rorbit'
      call writehdfrealmat(group_id,dsetname,
     $  rorbit,storage_dims,data_dims,rank)

      rank = 1
      data_dims(1) = norbits
      storage_dims(1) = nobsmax
      dsetname = 'iorbitlen'
      call writehdfintmat(group_id,dsetname,
     $  iorbitlen,storage_dims,data_dims,rank)

      endif

c Close the group.
      CALL h5gclose_f(group_id, error)


c Create a group in the file.
      groupname = 'orbtrack'
      CALL h5gcreate_f(file_id, groupname, group_id, error)

c Write the data in this group

c     Single value variables
      rank = 1
      data_dims(1) = 1
      storage_dims(1) = 1

      dsetname = 'orbinit'
      call writehdfintmat(group_id,dsetname,
     $ orbinit,storage_dims,data_dims,rank)

      dsetname = 'maxsteps'
      call writehdfintmat(group_id,dsetname,
     $ maxsteps,storage_dims,data_dims,rank)

      dsetname = 'trackinit'
      call writehdfintmat(group_id,dsetname,
     $ trackinit,storage_dims,data_dims,rank)

c Close the group.
      CALL h5gclose_f(group_id, error)


c Create a group in the file.
      groupname = 'mcr'
      CALL h5gcreate_f(file_id, groupname, group_id, error)

c Write the data in this group

c     Single value variables
      rank = 1
      data_dims(1) = 1
      storage_dims(1) = 1

      dsetname = 'mcrninjd'
      call writehdfintmat(group_id,dsetname,
     $ mcrninjd,storage_dims,data_dims,rank)

c     Variable arrays

      dsetname = 'mcrxpinjd'
      rank = 2
      data_dims(1) = ndim
      data_dims(2) = mcrninjd
      storage_dims(1) = ndim
      storage_dims(2) = npartmax
      call writehdfrealmat(group_id,dsetname,
     $ mcrxpinjd,storage_dims,data_dims,rank)

c Close the group.
      CALL h5gclose_f(group_id, error)


c Create a group in the file.
      groupname = 'colncom'
      CALL h5gcreate_f(file_id, groupname, group_id, error)

c Write the data in this group

c     Single value variables
      rank = 1
      data_dims(1) = 1
      storage_dims(1) = 1

      dsetname = 'vneutral'
      call writehdfrealmat(group_id,dsetname,
     $ vneutral,storage_dims,data_dims,rank)

      dsetname = 'Tneutral'
      call writehdfrealmat(group_id,dsetname,
     $ Tneutral,storage_dims,data_dims,rank)

      dsetname = 'Eneutral'
      call writehdfrealmat(group_id,dsetname,
     $ Eneutral,storage_dims,data_dims,rank)

      dsetname = 'NCneutral'
      call writehdfintmat(group_id,dsetname,
     $ NCneutral,storage_dims,data_dims,rank)

c Close the group.
      CALL h5gclose_f(group_id, error)


c Create a group in the file.
      groupname = 'noncommonblock'
      CALL h5gcreate_f(file_id, groupname, group_id, error)

c Write the data in this group

c     Single value variables
      rank = 1
      data_dims(1) = 1
      storage_dims(1) = 1

      dsetname = 'dt'
      call writehdfrealmat(group_id,dsetname,
     $ dt,storage_dims,data_dims,rank)

      dsetname = 'fave'
      call writehdfrealmat(group_id,dsetname,
     $ fave,storage_dims,data_dims,rank)

      dsetname = 'colnwt'
      call writehdfrealmat(group_id,dsetname,
     $ colnwt,storage_dims,data_dims,rank)

      dsetname = 'icolntype'
      call writehdfintmat(group_id,dsetname,
     $ icolntype,storage_dims,data_dims,rank)

c Close the group.
 999  CALL h5gclose_f(group_id, error)


c Terminate access to the file.
      CALL h5fclose_f(file_id, error)

c Close FORTRAN interface.
      CALL h5close_f(error)

      end


c*********************************************************************
c Outputs a real matrix to a specified group in an hdf5 file
      subroutine writehdfrealmat2(group_id,dsetname,dset_data,
     $  storage_dims,data_dims,rank)
c Load the hdf5 module
      use hdf5
c Don't allow implicit definitions
      implicit none
c Input variables
      INTEGER(HID_T) :: group_id
      CHARACTER(LEN=55) :: dsetname
      REAL, DIMENSION(*) :: dset_data
      INTEGER(HSIZE_T), DIMENSION(*) :: data_dims, storage_dims
      INTEGER     ::   rank
c Local variables
      INTEGER(HID_T) :: dataset_id
      INTEGER(HID_T) :: dataspace_id
      INTEGER     ::   error

c Create a data space for the dataset.
      CALL h5screate_simple_f(rank, storage_dims, dataspace_id, error)

c Create a dataset in the group with default properties.
      CALL h5dcreate_f(group_id, dsetname, H5T_NATIVE_REAL,
     $  dataspace_id, dataset_id, error)

c Write the dataset.
      CALL h5dwrite_f(dataset_id, H5T_NATIVE_REAL, dset_data,
     $  data_dims, error)

c Close the dataspace.
      CALL h5sclose_f(dataspace_id, error)

c Close the dataset.
      CALL h5dclose_f(dataset_id, error)

      end

c*********************************************************************
c Outputs a real matrix to a specified group in an hdf5 file
      subroutine writehdfrealmat(group_id,dsetname,dset_data,
     $  storage_dims,data_dims,rank)
c Load the hdf5 module
      use hdf5
c Don't allow implicit definitions
      implicit none
c Input variables
      INTEGER(HID_T) :: group_id
      CHARACTER(LEN=55) :: dsetname
      REAL, DIMENSION(*) :: dset_data
      INTEGER(HSIZE_T), DIMENSION(*) :: data_dims, storage_dims
      INTEGER     ::   rank
c Local variables
      INTEGER(HID_T) :: dataset_id
      INTEGER(HID_T) :: dataspace_id, memspace_id
      INTEGER     ::   error
c     Hyperslab offset in memory
      INTEGER(HSIZE_T), DIMENSION(10) :: offset =
     $  (/0,0,0,0,0,0,0,0,0,0/)

c Create a data space for the dataset.
      CALL h5screate_simple_f(rank,data_dims,dataspace_id,error)

c Create a data space for the memory (since only part of array used).
      CALL h5screate_simple_f(rank,storage_dims,memspace_id,error)

c Select hyperslab in memory.
      CALL h5sselect_hyperslab_f(memspace_id, H5S_SELECT_SET_F,
     $  offset, data_dims, error)

c Create a dataset in the group with default properties.
      CALL h5dcreate_f(group_id, dsetname, H5T_NATIVE_REAL,
     $  dataspace_id, dataset_id, error)

c Write the dataset.
      CALL h5dwrite_f(dataset_id, H5T_NATIVE_REAL, dset_data,
     $  data_dims, error, memspace_id)

c Close the dataspace.
      CALL h5sclose_f(dataspace_id, error)

c Close the dataset.
      CALL h5dclose_f(dataset_id, error)

      end

c*********************************************************************
c Outputs an integer matrix to a specified group in an hdf5 file
      subroutine writehdfintmat(group_id,dsetname,dset_data,
     $  storage_dims,data_dims,rank)
c Load the hdf5 module
      use hdf5
c Don't allow implicit definitions
      implicit none
c Input variables
      INTEGER(HID_T) :: group_id
      CHARACTER(LEN=55) :: dsetname
      INTEGER, DIMENSION(*) :: dset_data
      INTEGER(HSIZE_T), DIMENSION(*) :: data_dims, storage_dims
      INTEGER     ::   rank
c Local variables
      INTEGER(HID_T) :: dataset_id
      INTEGER(HID_T) :: dataspace_id, memspace_id
      INTEGER     ::   error
c     Hyperslab offset in memory
      INTEGER(HSIZE_T), DIMENSION(10) :: offset =
     $  (/0,0,0,0,0,0,0,0,0,0/)

c Create a data space for the dataset.
      CALL h5screate_simple_f(rank,data_dims,dataspace_id,error)

c Create a data space for the memory (since only part of array used).
      CALL h5screate_simple_f(rank,storage_dims,memspace_id,error)

c Select hyperslab in memory.
      CALL h5sselect_hyperslab_f(memspace_id, H5S_SELECT_SET_F,
     $  offset, data_dims, error)

c Create a dataset in the group with default properties.
      CALL h5dcreate_f(group_id, dsetname, H5T_NATIVE_INTEGER,
     $  dataspace_id, dataset_id, error)

c Write the dataset.
      CALL h5dwrite_f(dataset_id, H5T_NATIVE_INTEGER, dset_data,
     $  data_dims, error, memspace_id)

c Close the dataspace.
      CALL h5sclose_f(dataspace_id, error)

c Close the dataset.
      CALL h5dclose_f(dataset_id, error)

      end
