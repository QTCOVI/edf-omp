c
c-----------------------------------------------------------------------
c
      subroutine locastd (sg,nmo,ngroup,ncent,nfugrp,ifugrp,pcut,
     &   lw,epsloc,okw)
c
c.....Rather similar to ruedmis routine
c-----------------------------------------------------------------------
c                         
      include 'implicit.inc'
      include 'constants.inc'
      include 'error.inc'
      include 'stderr.inc'
c
      real   (kind=8), allocatable,dimension (:,:,:)  :: sgn
      real   (kind=8), allocatable,dimension (:)      :: xloci
      real   (kind=8), allocatable,dimension (:,:)    :: c
      integer(kind=4), allocatable,dimension (:)      :: iloc
      integer(kind=4), allocatable,dimension (:)      :: mocore,eignz
      integer(kind=4), allocatable,dimension (:,:)    :: mimael
      integer(kind=4), allocatable,dimension (:,:)    :: resnc

      parameter    (closetone=0.001D+00)
      real(kind=8) aom1
      logical      noloc
      real   (kind=8) sg(ngroup,nmo,nmo)
      integer(kind=4) nfugrp(ngroup),ifugrp(ncent,ngroup)
c
      call timer (2,ilocastd,'_locastd  ',-1)
c
c.....Allocate arrays
c
      allocate (c(nmo,nmo))
      allocate (xloci(nmo))
      allocate (iloc(nmo))

      write (lw,111)
      xloci(1:nmo) = 1.0d+00
      call ruedmis (sg,c,xloci,nmo,ngroup,lw,.false.,.true.)
      do i=1,nmo
        xloci(i)=dot_product(sg(1:ngroup,i,i),sg(1:ngroup,i,i))
      enddo
c
c-----Determine the localized and non localized MOs and
c     determine to which center each localized MO belongs
c
      allocate (eignz(nmo))
      allocate (mocore(ngroup))
      nloc  = 0
      nzeig = 0
      mocore(1:ngroup)=0
      do i=1,nmo
        if (abs(xloci(i)-1.0d+00).lt.epsloc) then 
          nloc=nloc+1
          iloc(nloc)=i

          smax=0.0d+00
          do k=1,ngroup
             if (sg(k,i,i).gt.smax) then
               smax=sg(k,i,i)
               kmaximo = k
             endif
          enddo
          mocore(kmaximo)=mocore(kmaximo)+1
          write (lw,112) i,kmaximo,smax
        else
          nzeig=nzeig+1
          eignz(nzeig)=i
        endif
      enddo
      write (lw,1111) nloc,' localized MOs: ',(iloc(i),i=1,nloc)
      write (lw,1111) nzeig,' NON localized MOs:'
      write (lw,113) (mocore(i),i=1,ngroup)
c
c-----Reconstruct the Group Overlap Matrix.
c
      allocate (sgn(ngroup,nzeig,nzeig))
      do k=1,ngroup-1
        do i=1,nzeig
          i1=eignz(i)
          do j=1,nzeig
            j1=eignz(j)
            sgn(k,i,j)=sg(k,i1,j1)
          enddo
        enddo
      enddo
c
c-----Compute Group Overlap integrals for the last group
c
      do i=1,nzeig
       do j=1,nzeig
         if (i.eq.j) delta=1.0d+00
         if (i.ne.j) delta=0.0d+00
         sgn(ngroup,i,j)=delta-sum(sgn(1:ngroup-1,i,j))
       enddo
      enddo

      allocate (mimael(2,ngroup))
      mimael(1,1:ngroup)=0
      mimael(2,1:ngroup)=nzeig
      write (lw,766) ngroup
      do i=1,ngroup
        write (lw,800) i,nfugrp(i),mimael(1,i),mimael(2,i)
      enddo
c
c-----Compute number of probabilities.
c
      npab=1
      do i=1,ngroup-1
        npab=npab*(mimael(2,i)-mimael(1,i)+1)
      enddo
      nprev=npab
c
      allocate (resnc(nprev,ngroup))
c
c.....Computation of resonance structures.
c
      call rnprobs (mimael,resnc,nzeig,ngroup,npab,lw)
      if (nprev.lt.npab) then
        write (stderr,*) 'locastd.f: NPREV < NPAB returned by rnprobs.f'
        stop
      endif
c
c-----Compute EDF.
c
      call compedfstd (pcut,npab,nprev,ngroup,nzeig,lw,mimael,sgn,
     &                resnc,mocore)
c
c-----Deallocate arrays.
c
      deallocate (c)
      deallocate (xloci)
      deallocate (iloc)
      deallocate (mocore)
      deallocate (eignz)
      deallocate (sgn)
      deallocate (mimael)
      deallocate (resnc)
      call timer (4,ilocastd,'_locastd  ',-1)
      return
c
 111  format (' #',/,1x,'# ',75('-'),/,
     &  25x,'ISOPYCNIC LOCALIZATION MODULE')
 112  format (' # MO ',I3,' : Localized on fragment',I3,' = ',F16.9)
 113  format (' # MOs almost fully-localized in each fragment:',/,
     &    1000(' # ',20I3,/))
 766  format (' # NUMBER OF GROUPS = ',I2)
 764  format (' # GROUP ',I2,
     &  ' MinElec and MaxElec (alpha or beta) = ',2I6)
 767  format (1x,'# Number of atoms in the group = ',I3)
 765  format (' # ',10I6)
 1111 format (' # ',I3,a,1000(30I3,/))
 800  format (' # GROUP ',i4,' has ',I4,' atoms',4x,
     &    'MinElec,MaxElec (ALPHA or BETA) = ',2I4)
      end
