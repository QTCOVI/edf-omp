c
c.....Read AOM for atoms in groups 1 to NGROUP-1
c
      subroutine readaomstd (aom,nfugrp,ifugrp,ncent,ngroup,
     &  nmo,naom,lu18,ler)
      include    'constants.inc'
      include    'mline.inc'
      integer(kind=4)  nmo,lu18,ler,i,m,j
      real   (kind=8), allocatable,dimension (:,:,:)  :: aomr,aomi

      integer(kind=4)  nfugrp(ngroup)
      integer(kind=4)  ifugrp(ncent,ngroup)
      real   (kind=8)  aom(ncent,nmo,nmo)
      real   (kind=8)  aom1
      real   (kind=8)  aom1cmplx
      integer(kind=4)  nread(ncent)
      character*(mline) line
c
      ntot = 0
      do
        read (lu18,*,end=100) iat
        if (iat.le.0.or.iat.gt.ncent) then
          write (ler,*) '# readaomstd.f: There are ',ncent,' atoms and'
          write (ler,*) '# you are trying to read AOM of atom',iat
          stop
        else
          do i=1,ntot
            if (iat.eq.nread(i)) then
              write (ler,*) '# readaomstd.f: AOM for atom ',iat
              write (ler,*) '# appears at least twice in the AOM file'
              stop
            endif
          enddo
          ntot = ntot + 1
          nread(ntot) = iat
          read (lu18,80) ((aom(iat,m,j),m=1,j),j=1,nmo) ! promolden format
        endif
      enddo
 100  continue
      naom = ntot
c
c-----Test that all AOMs for groups 1 to NGROUT-1 have been read in
c
      do i=1,ngroup-1
        cyclej: do j=1,nfugrp(i)
          ij=ifugrp(j,i)
          do k=1,ntot
            if (ij.eq.nread(k)) cycle cyclej
          enddo
          write (ler,*) '# readaomstd.f: AOM not read for atom ',ij
          stop 
        enddo cyclej
      enddo
c
c.....Fill in the lower part of AOM
c
      do i=1,ntot
        iat=nread(i)
        do j=1,nmo
          do m=1,j
            aom(iat,j,m)=aom(iat,m,j)
          enddo
        enddo
      enddo
      return
 80   format (6(1x,e16.10))
 101  format ('# readaomstd.f: The number of read in AOMs ( ',I3,
     &  ' ) does not coincide with the number atoms in the first'
     &  I3,' groups ( ',I3,' )')
      end
