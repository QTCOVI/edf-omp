c
c-----------------------------------------------------------------------
c
      subroutine rnprobs (elpol,resnc,numelec,nfrag,nprob,stdout)
c
c-----Recurvive version of the nprobs.f routine.
c
c.......................................................................
c
c     Given 'nfrag' 3D domains (for instance, the QTAM atomic basins), 
c     the total number of electrons 'numelec', and a minimum and a maximum 
c     value for the number of electrons that each 3D domain can accommo-
c     date (elpol), this routine finds all the resonance structures of 
c     the molecule compatible with the 'numelec' and the 'elpol' values.
c
c     example: The CH_4 molecule, each QTAM basin being a 3D domain and
c     each C (H) atom having a minimum number of electrons equal to 4 
c     (0) and a maximum number of electrons equal to 8 (2) has the 
c     following input of this routine.
c
c     nel=10         ngroup=5 
c     elpol(1,1)=4   elpol(2,1)=8
c     elpol(1,2)=0   elpol(2,2)=2
c     elpol(1,3)=0   elpol(2,3)=2
c     elpol(1,4)=0   elpol(2,4)=2
c     elpol(1,5)=0   elpol(2,5)=2
c
c.......................................................................
c
      include   'implicit.inc'
      integer,   allocatable, dimension (:,:) :: ielpol
      integer,   allocatable, dimension (:)   :: ires
c
      integer elpol(2,nfrag)
      integer resnc(nprob,nfrag)
      integer stdout

      integer nel,ngroup
      common /ngnel/ ngroup,nel,lurs
c
      ngroup=nfrag
      allocate (ielpol(1:2,1:ngroup),stat=ier)
      if (ier.ne.0) stop 'rnprobs.f: Cannot allocate ielpol()'
      allocate (ires(1:ngroup),stat=ier)
      if (ier.ne.0) stop 'rnprobs.f: Cannot allocate ires()'

      nel=numelec
      ng=ngroup
      ne=nel
      ielpol(1:2,1:ngroup)=elpol(1:2,1:ngroup)
      np=1
      lurs=71
      open (unit=lurs,file='rsrs.dat')
      call folp (ielpol,ires,ne,ng,np)
      nprob=np-1
      rewind(lurs)
      do i=1,nprob
        read (lurs,1) (resnc(i,k),k=1,nfrag)
      enddo
      close (lurs,status='delete')
      return
 1    format (50(I4))
      end

      recursive subroutine folp (ielpol,ires,ne,ng,np)
      include   'implicit.inc'
      common /ngnel/ ngroup,nel,lurs
      integer ielpol(1:2,1:ngroup),ires(1:ngroup)

      do i=ielpol(1,ng),ielpol(2,ng)
        if (ng.eq.1) then
          if (i.eq.ne) then
            ires(1)=ne
            write (lurs,1) (ires(k),k=1,ngroup)
*           write (6,1) (ires(k),k=1,ngroup)
            np=np+1
          endif
        else
          ires(ng)=i
          call folp (ielpol,ires,ne-i,ng-1,np)
        endif
      enddo
      return
 1    format (50(I4))
      end
