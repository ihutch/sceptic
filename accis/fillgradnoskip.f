
c******************************************************************
c Given a triangle consisting of 3x3D points with coordinates x,y,z,
c fill with a gradient. zg0,zg1 are the lower and upper parameter values
c of the which are mapped to colors in the range ng0,ng1. The switch isw
c determines how d(3) is used for parameter values and whether polyline
c or poly3line is called (2-D x,y or 3-D projection).
c First bit of isw determines 2 or 3D.
c Second bit of isw uses d as absolute quantity at each vertex (0)
c or d as direction vector (not nec normalized) and parameter is
c position of the vertex in that direction (1). Thus
c 0: 2-D d-value  1: 3-D d-value  2: 2-D d-direction  3: 3D d-direction.
c Second byte of isw tells the number of color levels to step by for
c each gradient fill. Default (0 or 1) is 1 level, the finest gradient.

c Examples:
c isw=0, color contour value of d on x,y plane (z-irrelevant).
c isw=1, color by value of d on 3-D perspective surface xyz.
c isw=2, d=(d1,d2,0), z irrelevant, 2-D fill in direction d.
c        d=(0,0,1), color contour the height z. 
c isw=3, color contour by direction d, and 3-D perspect (needs setup).


      subroutine gradtri(x,y,z,d,zg0,zg1,ng0,ng1,isw)
      real x(3),y(3),z(3),d(3),zg0,zg1
      integer ng0,ng1,isw

c There were X-problems using implied closure of polygons. So don't.
      real h(3),hmin,hmax,xp(6),yp(6),zp(6)
      integer id(3),imin,imax,l,ll,lmin,lmax
      logical labs

      include 'plotcom.h'

c Stepping resolution is in second byte.
      istep=isw/256
      if(istep.lt.1)istep=1

      if((isw/2-2*(isw/4)).ne.0)then
c direction
         labs=.false.
         dn=sqrt(d(1)**2+d(2)**2+d(3)**2)
         if(dn.eq.0.) stop 'ERROR gradtri: zero dn'
      else
c absolute
         labs=.true.
      endif
      if(zg1-zg0.eq.0.) stop 'ERROR gradtri: zero dzng'
      hr=(ng1-ng0)/(zg1-zg0)

      hmin= 1.e30
      hmax=-1.e30
      do i=1,3
         if(labs)then
c absolute normalized to color index
            h(i)=(d(i)-zg0)*hr+ng0
         else
c Distance in d-direction normalized to effective color index.
            h(i)=((x(i)*d(1)+y(i)*d(2)+z(i)*d(3))/dn-zg0)*hr+ng0
         endif
         if(h(i).gt.hmax)then 
            hmax=h(i)
            imax=i
         endif
         if(h(i).lt.hmin)then
            hmin=h(i)
            imin=i
         endif
      enddo
c Index the vertices in the standard order:
      id(1)=imin
      id(2)=imax
c Arithmetic trick to find the third index:
      id(3)=6-(imin+imax)
c      write(*,*)id(3),imin,imax
      hmid=h(id(3))

      if((min(ng0,ng1).gt.hmax) .or.
     $     (max(ng0,ng1).lt.hmin) ) then
c         write(*,*)'Color range does not intersect triangle.'
         return
      endif
      ll=hmin-1
      lmin=max(int(hmin/istep),min(ng0,ng1)/istep)*istep
      lmax=min(int(hmax+istep),max(ng0,ng1))
c Ensure we don't divide by zero. Because we are normalized to level 
c number, we can just add a small quantity. 
      hdn=hmid-hmin+1.e-20
      hxd=hmax-hmid+1.e-20
      hxn=hmax-hmin+1.e-20
c Initialize mid point.
      xp(5)=x(id(3))
      yp(5)=y(id(3))
      zp(5)=z(id(3))
c Do over the levels:
      do l=lmin,lmax,istep
         if(hmid.gt.l)then
            fp=(l-hmin)/hdn
            if(fp.gt.1.)fp=1.
            if(fp.lt.0.)fp=0.
            xp(1)=x(id(1))*(1-fp)+x(id(3))*fp
            yp(1)=y(id(1))*(1-fp)+y(id(3))*fp
            zp(1)=z(id(1))*(1-fp)+z(id(3))*fp
         else
            fp=(l-hmid)/hxd
            if(fp.gt.1.)fp=1.
            if(fp.lt.0.)fp=0.
            xp(1)=x(id(3))*(1-fp)+x(id(2))*fp
            yp(1)=y(id(3))*(1-fp)+y(id(2))*fp
            zp(1)=z(id(3))*(1-fp)+z(id(2))*fp
         endif
         fp=(l-hmin)/hxn
         if(fp.gt.1.)fp=1.
         if(fp.lt.0.)fp=0.
         xp(2)=x(id(1))*(1-fp)+x(id(2))*fp
         yp(2)=y(id(1))*(1-fp)+y(id(2))*fp
         zp(2)=z(id(1))*(1-fp)+z(id(2))*fp
         if(l.gt.lmin)then
c Except the first time, fill the polygon.
c If we crossed the mid point.
            if((ll-hmid)*(l-hmid).le.0.)then
               np=5
            else
               np=4
            endif
            call gradcolor(l)
            if((isw-2*(isw/2)).eq.0)then
               call polyline(xp,yp,np)
            else
               call poly3line(xp,yp,zp,np)
            endif
            call pathfill()
c Testing only
c            call color(15)
c            call polyline(xp,yp,np)
         endif
         ll=l
         xp(4)=xp(1)
         yp(4)=yp(1)
         zp(4)=zp(1)
         xp(3)=xp(2)
         yp(3)=yp(2)
         zp(3)=zp(2)
 2       continue
      enddo


      end
c***********************************************************************
c Fill a quadrilateral by dividing it up into 4 triangles at centroid.
c Quadrilateral vertices must be in sequential order. 
c See gradtri for the full significance of isw.
c isw2 set => d is direction vector, else absolute value.
c Second byte of isw: level-skip. 
      subroutine gradquad(xq,yq,zq,dq,zg0,zg1,ng0,ng1,isw)
      real xq(4),yq(4),zq(4),dq(4),zg0,zg1
      integer ng0,ng1

      real x(3),y(3),z(3),d(3)

      include 'plotcom.h'
c Make line width minimal for contouring
       if(abs(pfsw).eq.2 .or. abs(pfsw).eq.3)
     $     call abufwrt(' 0 setlinewidth ',16,12)
      

c centroid point
      x(1)=(xq(1)+xq(2)+xq(3)+xq(4))*.25
      y(1)=(yq(1)+yq(2)+yq(3)+yq(4))*.25
      z(1)=(zq(1)+zq(2)+zq(3)+zq(4))*.25
      d(1)=(dq(1)+dq(2)+dq(3)+dq(4))*.25

c fill each subtriangle in turn.
      if((isw/2-2*(isw/4)).eq.0)then
c Absolute. Swap d as well as xyz.
         do i=1,4
            x(2)=xq(i)
            y(2)=yq(i)
            z(2)=zq(i)
            d(2)=dq(i)
            x(3)=xq(1+mod(i,4))
            y(3)=yq(1+mod(i,4))
            z(3)=zq(1+mod(i,4))
            d(3)=dq(1+mod(i,4))
            call gradtri(x,y,z,d,zg0,zg1,ng0,ng1,isw)
         enddo
      else
c Direction. No d changes needed. dq->d
         do i=1,4
            x(2)=xq(i)
            y(2)=yq(i)
            z(2)=zq(i)
            x(3)=xq(1+mod(i,4))
            y(3)=yq(1+mod(i,4))
            z(3)=zq(1+mod(i,4))
            call gradtri(x,y,z,dq,zg0,zg1,ng0,ng1,isw)
         enddo
      endif
      end
c***********************************************************************
c Fill a quadrilateral by dividing it up into 2 triangles
c Quadrilateral vertices must be in sequential order. 
c See gradtri for the full significance of isw.
c isw2 set => d is direction vector, else absolute value.
c The diagonal along which the quadrilateral is split is the one along
c which the value changes most. 
      subroutine gradquad2(xq,yq,zq,dq,zg0,zg1,ng0,ng1,isw)
      real xq(4),yq(4),zq(4),dq(4),zg0,zg1
      integer ng0,ng1

      real x(3),y(3),z(3),d(3)

      include 'plotcom.h'
c Make line width minimal for contouring
       if(abs(pfsw).eq.2 .or. abs(pfsw).eq.3)
     $     call abufwrt(' 0 setlinewidth ',16,12)
      
       id=0
      if((isw/2-2*(isw/4)).eq.0)then
c Absolute
         if(abs(dq(1)-dq(3)).lt.abs(dq(2)-dq(4))) id=1
      else
c Direction.
         a1=dq(1)*xq(1)+dq(2)*yq(1)+dq(3)*zq(1)
         a2=dq(1)*xq(2)+dq(2)*yq(2)+dq(3)*zq(2)
         a3=dq(1)*xq(3)+dq(2)*yq(3)+dq(3)*zq(3)
         a4=dq(1)*xq(4)+dq(2)*yq(4)+dq(3)*zq(4)
         if(abs(a1-a3).lt.abs(a2-a4)) id=1
      endif
      do i=1,3,2
         ip=i+id
         x(i)=xq(ip)
         y(i)=yq(ip)
         z(i)=zq(ip)
         d(i)=dq(ip)
      enddo
      do i=2,4,2
         ip=i-id
         x(2)=xq(ip)
         y(2)=yq(ip)
         z(2)=zq(ip)
         d(2)=dq(ip)
c         write(*,'(3i2,9f8.4)')i,id,ip,x,y,d
         if((isw/2-2*(isw/4)).eq.0)then
            call gradtri(x,y,z,d,zg0,zg1,ng0,ng1,isw)
         else
            call gradtri(x,y,z,dq,zg0,zg1,ng0,ng1,isw)
         endif
      enddo

      end
c***********************************************************************
