      subroutine testaom (aom,tolaom,ncent,nmo,lr)
      implicit none
      real (kind=8) aom(ncent,nmo,nmo)
      real (kind=8) sumaom,valueaom,diffaom,tolaom
      integer (kind=4) ncent,nmo,i,j,m,lr

      do m=1,nmo
        do j=1,m
          sumaom=sum(aom(:,j,m))
          valueaom=1d0
          if (j.ne.m) valueaom=0d0
          diffaom=abs(sumaom-valueaom)
          if (diffaom.gt.tolaom) write (lr,544) m,j,sumaom
        enddo
      enddo
      return
 544  format (' # !! AOM is not Ok !! : SUM_n AOM(n,',2I3,') = ',E15.8)
      end
