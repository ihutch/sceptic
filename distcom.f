c Requires piccom.f to be input first.
      integer nrdist,nthdist,nvdist
      parameter (nrdist=nrsize,nthdist=nthsize,nvdist=60,ndtype=5)
      real fvrtdist(nvdist,ndtype,0:nrdist,0:nthdist)
      real rhodist(0:nrdist,0:nthdist)

      common /distcom/fvrtdist,rhodist

c fvrtdist is the accumulator of distribution function information for
c the following components of velocity, in this order:
cc Cylindrical Radial velocity
c         vxyr=xp(4,i)*cp+xp(5,i)*sp
cc Spherical and cylindrical azimuthal velocity:
c         vp=-xp(4,i)*sp+xp(5,i)*cp
cc Longitudinal velocity
c         vz=xp(6,i)
cc Spherical Radial velocity
c         vr=vz*ct+vxyr*st
cc Spherical Tangential velocity without azimuthal component.
cc I.e. v_theta.
c         vt=-vz*st+vxyr*ct
cc Cylindrical radial, azimuthal, and longitudinal bins.
c         ivxy=min(1+max(0,nint(nvmax*(vxyr/vrange+.499))),nvmax)
c         ivp=min(1+max(0,nint(nvmax*(vp/vrange + .499))),nvmax)
c         ivz=min(1+max(0,nint(nvmax*(vz/vrange+.499))),nvmax)
cc Spherical Radial and angular direction velocity bins: 
c         ivr=min(1+max(0,nint(nvmax*(vr/vrange + .499))),nvmax)
c         ivt=min(1+max(0,nint(nvmax*(vt/vrange + .499))),nvmax)
c
cc Update accumulators.
c         fvrtdist(ivxy,1,irl,ithl)=fvrtdist(ivxy,1,irl,ithl)+1.
c         fvrtdist(ivp,2,irl,ithl)=fvrtdist(ivp,2,irl,ithl)+1.
c         fvrtdist(ivz,3,irl,ithl)=fvrtdist(ivz,3,irl,ithl)+1.
c         fvrtdist(ivr,4,irl,ithl)=fvrtdist(ivr,4,irl,ithl)+1.
c         fvrtdist(ivt,5,irl,ithl)=fvrtdist(ivt,5,irl,ithl)+1.
