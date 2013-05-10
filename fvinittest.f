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


