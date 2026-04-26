c
c-----------------------------------------------------------------------
c
      subroutine writeaom (aom,ncent,nmo,lw,lr,uhf)
      implicit none
      integer (kind=4) ncent,nmo,lw,lr,i,m,j
      real    (kind=8) aom(ncent,nmo,nmo)
      logical  uhf
c
      write (lw,*) '#'
      write (lw,*) '# AOM elements'
      do i=1,ncent
        write (lw,*) '# Atom ',i
        do m=1,nmo
          write (lw,'(5(1x,F22.15))') (aom(i,m,j),j=1,m)
        enddo
      enddo
      write (lw,*) '# SUMAOM'
      do m=1,nmo
        write (lw,'(5(1x,F22.15))') (sum(aom(1:ncent,m,j)),j=1,m)
      enddo
      return
      end
