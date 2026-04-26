c
c-----------------------------------------------------------------------
c
c     This routine uses a pure brute force method to compute the p(S)
c     probabilities of all the real space resonance structures (RSRS) 
c     associated to a partitionn of the molecule into NGROUP framents
c     such that the electron population of fragment i is between minp(i)
c     and maxp(i), i.e. minp(i) <= elec(i) <= maxp(i)
c
c
      subroutine dirprob (pcut,na,nb,ngroup,lw,sa,sb,minpopul,maxpopul)
      implicit none
      integer(kind=4)  na,nb,ngroup,lw,i
      real   (kind=8)  sa(ngroup,na,na)
      real   (kind=8)  sb(ngroup,nb,nb)
      real   (kind=8)  pcut
      integer(kind=4)  maxp(ngroup),maxpopul(ngroup)
      integer(kind=4)  minp(ngroup),minpopul(ngroup)
      character (len=100) filea,fileb
c
c-----Compute minumum and maximum population for each spin channel
c
      do i=1,ngroup
        maxp(i)=(maxpopul(i)+1)/2
        minp(i)=max((minpopul(i)-1)/2,0)
      enddo
      write (lw,1)
 1    format (//,1x,'# ',85('+'),/,' #',6x,'DIRPROB Routine',
     &        /,1x,'# ',85('+'))
      write (lw,5) na,nb
 5    format (1x,'# (Alpha,Beta) MOs:',2x,2I4)
c
c-----ALPHA spin channel probabilities
c
      call dirproba (pcut,na,ngroup,lw,1,sa,minp,maxp,filea)
c
c-----BETA spin channel probabilities
c
      call dirproba (pcut,nb,ngroup,lw,2,sb,minp,maxp,fileb)

c
c-----Full Spinless probabilities
c
*     filea='probALPHA'
*     fileb='probBETA'
      call pfromab (pcut,ngroup,lw,minpopul,maxpopul,filea,fileb)
      return
      end
c
c-----------------------------------------------------------------------
c
      subroutine dirproba (pcut,nab,ngroup,lw,isp,over,minp,maxp,fprob)
      implicit none
      integer(kind=4)  nab,ngroup,lw,isp
      real   (kind=8)  over(ngroup,nab,nab)
      real   (kind=8)  pcut
      integer(kind=4)  maxp(ngroup)
      integer(kind=4)  minp(ngroup)
c
c-----Arrays used in the nested symmation symbol (NSS) strategy.
c     of Emili Besalú.
c
      integer(kind=4), allocatable,dimension(:) :: ilo
      integer(kind=4), allocatable,dimension(:) :: jlo
      integer(kind=4), allocatable,dimension(:) :: flo
      integer(kind=4), allocatable,dimension(:) :: slo
c
      integer(kind=4), parameter :: nstack = 5000
      integer(kind=4), allocatable,dimension(:,:) :: resnc
      real   (kind=8), allocatable,dimension(:)   :: prob
      real   (kind=8)  :: s(nab,nab)
      integer(kind=4)  :: ind(nab)
      integer(kind=4)  :: nres(ngroup)
      integer(kind=4)  :: i,j,k,l,np,info,ls,ios
      real   (kind=8)  :: det
      real   (kind=8)  :: sumpcut
      logical          :: first,docalc,inlist,exfil
      character(len=5) :: spin
*     character(len=80):: fprob
      character(len=*) :: fprob
c
c-----The nested symmation symbol (NSS) strategy for do loop's 
c     of Emili Besalú is used
c
      allocate (ilo(nab))
      allocate (jlo(nab))
      allocate (flo(nab))
      allocate (slo(nab))
      do k=1,nab
        ilo(k)=1        ! Initial value
        flo(k)=ngroup   ! Final value
        slo(k)=1        ! Step increment
        jlo(k)=ilo(k)   ! Initialization of NSS variable
      enddo
c
      allocate (resnc(nstack,ngroup))
      allocate (prob(nstack))
      k = nab
      first = .true.
      do while (k.gt.0)
         if ((jlo(k)-flo(k))*slo(k).gt.0) then
            jlo(k)=ilo(k)
            k=k-1
         else
            write (6,'(36I3)') jlo(:)
            nres(1:ngroup)=0
            do i=1,nab
              nres(jlo(i))=nres(jlo(i))+1
            enddo
            docalc=.true.
            do i=1,ngroup
              if (nres(i).gt.maxp(i).or.nres(i).lt.minp(i)) then
                docalc=.false.
                exit
              endif
            enddo
            if (docalc) then
              do i=1,nab
                do j=1,nab
                  s(i,j)=over(jlo(j),i,j)
                enddo
              enddo
              call detlapack (s,ind,nab,info,det)
              if (first) then
                np=1
                resnc(np,1:ngroup)=nres(1:ngroup)
                first=.false.
                prob(np)=det
              else
                do j=1,np
                  inlist=.true.
                  do l=1,ngroup
                    inlist=inlist.and.(nres(l).eq.resnc(j,l))
                  enddo
                  if (inlist) then
                    prob(j)=prob(j)+det
                    goto 1
                  endif
                enddo
                np=np+1
                if (np.gt.nstack) write (0,4) nstack
                if (np.gt.nstack) stop
                prob(np)=det
                resnc(np,1:ngroup)=nres(1:ngroup)
 1              continue
              endif
            endif
            k=nab
         endif
         if (k.gt.0) jlo(k)=jlo(k)+slo(k)
      enddo
c
c-----Print results
c
      continue
      if (isp.eq.1) spin='ALPHA'
      if (isp.eq.2) spin='BETA '
c
c-----Open file to write probabilities
c
      fprob='prob'//trim(spin)
      ls  =  20
      exfil = .true.
      do while (exfil)
         inquire (unit=ls,opened=exfil)
         if (.not.exfil) then
           open (unit=ls,file=trim(fprob),iostat=ios)
           if (ios.ne.0) then
             write (0,*) ' # dirprob.f: Error openning '//trim(fprob)
             stop
           endif
         else
           ls=ls+1
         endif
      enddo

      write (lw,2) np,trim(spin),pcut
      write (ls,'(2I10)') ngroup,np
      sumpcut=0d0
      do j=1,np
        if (prob(j).ge.pcut) then
          write (lw,10) prob(j),(resnc(j,k),k=1,ngroup)
          sumpcut=sumpcut+prob(j)
        endif
        write (ls,'(E25.18,20I6)') prob(j),(resnc(j,k),k=1,ngroup)
      enddo
      write (lw,3) sum(prob(1:np)),sumpcut,pcut
      write (lw,5) trim(fprob), trim(spin)
 5    format (' # File "',a,'" with ',a,' probabilities is written')
      close (ls)
c
c-----Deallocate arrays
c
      deallocate (ilo,jlo,flo,slo)

 10   format (1x,'# ',F22.15,3x,20I6)
 4    format (' # Increase the "NSTACK" value, currently',I8)
 2    format (1x,'# ',85('-'),/,1x,
     & '#',6x,'There are ',I10,' probabilities for spin channel ',a,/,
     & 1x,'#',6x,'Printing probabilities greater than ',F15.8,/,1x,
     & '#',6x,'probability',13x,'RSRS --->',/,1x,'#',6x,80('-'))
 3    format (1x,'# ',F22.15,2x,'<--- SUM',/,
     &        1x,'# ',F22.15,2x,'<--- SUM ( p > ',F15.8,' )')
      return
      end
