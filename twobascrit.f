c
c-----------------------------------------------------------------------
c
      subroutine twobascrit (sga,sgb,nmo,nspin,ngroup,pcut,lr,lw,lerr)
c
c.......................................................................
c
      include    'implicit.inc'
      include    'param.inc'
      include    'constants.inc'
      include    'stderr.inc'
      complex*16  sga(ngroup,nmo,nmo)
      complex*16  sgb(ngroup,nmo,nmo)

      complex*16,      allocatable,dimension (:)     :: alph,beta
      complex*16,      allocatable,dimension (:,:)   :: sdiag,bdiag
      complex*16,      allocatable,dimension (:)     :: wdiag
      complex*16,      allocatable,dimension (:,:,:) :: sg
      complex*16,      allocatable,dimension (:,:)   :: anun
      real   (kind=8), allocatable,dimension (:,:)   :: prob
      real   (kind=8), allocatable,dimension (:)     :: probt
      real   (kind=8), allocatable,dimension (:)     :: seigen
      real   (kind=8), allocatable,dimension (:)     :: p1a,p1b
      integer(kind=4), allocatable,dimension (:)     :: nnonz
      real   (kind=8), allocatable,dimension (:)     :: pnew
      integer(kind=4), allocatable,dimension (:)     :: ioprob
      integer(kind=4), allocatable,dimension (:,:)   :: resncord
      character(len=5) spin(2)
c
      call timer (2,itwobascr,'_twobascr ',-1)
c
c.....Compute complex Group Overlap integrals of all but the last one.
c
      allocate (sg(ngroup,nmo,nmo))
      allocate (anun(0:nmo,0:nmo))
      allocate (alph(nmo))
      allocate (beta(nmo))
      allocate (prob(nmo+1,2))
      allocate (sdiag(nmo,nmo))
      allocate (bdiag(nmo,nmo))
      allocate (seigen(nmo))
      allocate (wdiag(nmo+nmo-1))
      allocate (nnonz(2))
c
      spin(1) = 'ALPHA'
      spin(2) = 'BETA '
      ispin = 1
 22   continue
      if (ispin.eq.1) sg = sga
      if (ispin.eq.2) sg = sgb
c
c.....Compute Group Overlap matrix of group 1 and diagonalize itº
c
      sdiag(:,:) = sg(1,:,:)
      call zjacobi (sdiag,bdiag,wdiag,seigen,nmo)
      write (lw,112) trim(spin(ispin))
      write (lw,222) (seigen(i),i=1,nmo)
      write (lw,113) sum(seigen(1:nmo))
c
c-----Compute anun() array 
c     (See A. Savin et al, Theor. Chem. Acc. 111, 373 (2004)).
c
      do k=1,nmo
        beta(k) = seigen(k)
        alph(k) = 1.0d+00-beta(k)
      enddo
c
      anun(0,0)=cmplx(1.0d+00,0.0d+00)
      do k=1,nmo
        anun(0,k)=anun(0,k-1)*alph(k)
        do j=1,k-1
          anun(j,k)=beta(k)*anun(j-1,k-1)+alph(k)*anun(j,k-1)
        enddo
        anun(k,k)=beta(k)*anun(k-1,k-1)
      enddo
      do nu=0,nmo
        prob(nu+1,ispin) = real(anun(nu,nmo))
      enddo
      write (lw,200) trim(spin(ispin)),ngroup,nmo+1
      sumgt = 0.0d+00
      sumtot = sum(prob(1:nmo+1,ispin))
      nnonz(ispin) = 0
      do n=0,nmo
        if (prob(n+1,ispin).ge.pcut) then
          nnonz(ispin) = nnonz(ispin) + 1
          nnz = nnonz(ispin)
          sumgt = sumgt + prob(n+1,ispin)
          write (lw,100) prob(n+1,ispin),n,nmo-n
        endif
      enddo
      write (lw,1060) sumgt, nnonz(ispin), pcut, sumtot

      if (ispin.eq.1.and.nspin.eq.2) then
        ispin = 2
        goto 22
      else
        prob(1:nmo+1,2) = prob(1:nmo+1,1)
        nnonz(2) = nnonz(1)
      endif
c
c.....Total spin-splitted probabilities 
c
      allocate (probt((nmo+1)*(nmo+1)))
      nprobg = 0
      sumgt = 0.0d+00
      sumtot = 0.0d+00
      npais = (nmo+1)
      npbis = (nmo+1)
      write (lw,2001) ngroup,npais*npbis,nnonz(1)*nnonz(2)
      ij = 0
      nprobt = nnonz(1) * nnonz(2)
      allocate (pnew(nprobt))
      pnew = 0.0d+00
      npnew = nnonz(1)+nnonz(2)-1
      allocate (resncord(npnew,2))
      do i=1,nnonz(1)
        do j=1,nnonz(2)
          ij = ij+1
          probt(ij) = prob(i,1)*prob(j,2)
          sumtot = sumtot+probt(ij)
          ijt = i+j-1
          pnew(ijt) = pnew(ijt) + probt(ij)
          resncord(ijt,1) = ijt-1
          resncord(ijt,2) = nmo+nmo-ijt+1
          if (probt(ij).ge.pcut) then
            nprobg=nprobg+1
            sumgt = sumgt+probt(ij)
            write (lw,100) probt(ij),i-1,nmo-i+1,j-1,nmo-j+1
          endif
        enddo
      enddo
      write (lw,106) sumgt, nprobg, pcut, sumtot
c
c.....Average population of alpha and beta electrons in each group.
c
      allocate (p1a(ngroup))
      allocate (p1b(ngroup))
      p1a = 0.0d+00
      do i=0,nmo
        p1a(1) = p1a(1) + prob(i+1,1) * i  
        p1a(2) = p1a(2) + prob(i+1,1) * (nmo-i) 
        p1b(1) = p1b(1) + prob(i+1,2) * i
        p1b(2) = p1b(2) + prob(i+1,2) * (nmo-i)
      enddo

      write (lw,1110)
      do i=1,ngroup
        write (lw,151) i,'ALPHA',p1a(i)
        write (lw,151) i,'BETA ',p1b(i)
      enddo
      tota = p1a(1)+p1a(2)
      totb = p1b(1)+p1b(2)
      totab = tota+totb
      write (lw,152) tota,totb,totab
c
c.....Obtain multiple-fragment electron population covariances.
c
      dela = 0.0d+00
      delb = 0.0d+00
      do i=0,nmo
        proda = (i-p1a(1))*(nmo-i-p1a(2))
        prodb = (i-p1b(1))*(nmo-i-p1b(2))
        dela = dela + proda * prob(i+1,1)
        delb = delb + prodb * prob(i+1,2)
      enddo
      delab = dela+delb
      dia = -2.0*dela
      dib = -2.0*delb
      diab = dia+dib
      write (lw,988) dela,delb,delab,dia,dib,diab

      write (lw,*) '#'
      write (lw,*) '# SPINLESS ELECTRON NUMBER DISTRIBUTION'
      if (npnew.le.nstack) then
        allocate (ioprob (npnew))
        forall (i=1:npnew) ioprob(i) = i
        call qqsort (pnew(1:npnew),ioprob,1,npnew,npnew)
        write (lw,*) '#     Ordered Probabilities'
        write (lw,*) '# -------------------------------------'
        write (lw,*) '#      Probability          n(G1) n(G2)'
        write (lw,*) '# -------------------------------------'
        do n=npnew,1,-1
          m = ioprob(n)
          write (lw,878) pnew(m),resncord(m,1:2)
        enddo
      else
        write (lw,*) '# Non-Ordered Probabilities'
        write (lw,*) '# -------------------------------------'
        write (lw,*) '#      Probability          n(G1) n(G2)'
        write (lw,*) '# -------------------------------------'
        do i=1,npnew
          write (lw,878) pnew(i),resncord(i,1:2)
        enddo
      endif
      write (lw,106) sumgt, npnew, pcut, sumtot
c
      deallocate (sg)
      deallocate (anun)
      deallocate (alph)
      deallocate (beta)
      deallocate (prob)
      deallocate (sdiag)
      deallocate (bdiag)
      deallocate (seigen)
      deallocate (wdiag)
      deallocate (nnonz)
      deallocate (probt)
      deallocate (pnew)
      deallocate (resncord)
      deallocate (p1a)
      deallocate (p1b)
      deallocate (ioprob)
c
      call timer (4,itwobascr,'_twobascr ',-1)
      return
c
c.....Formats.
c
 112  format (/,' # ',70('+'),/,' # Spin Block ',a,/,
     & ' # Ordered Eigenvalues of Group Overlap Matrix of GROUP 1',/,
     & ' #',70('+'))
 113  format (/,' SUM =',13x,F17.10)
 202  format (/,' # M-BASINS ELECTRON NUMBER PROBABILITY DISTRIBUTION',
     & ' NOT INCLUDING SPIN',
     & /,1x,'#',72('-'),/,
     &  ' # NUMBER OF GROUPS               = ',I8,/,
     &  ' # TOTAL NUMBER OF PROBABILITIES  = ',I8,/,
     &   1x,'#',72('-'),/,
     & ' #     Probability            n1    n2    n3 ...')
 100  format (' # ',F22.15,1x,12I6)
 1060 format (1x,'#',72('-'),/,' # ',F22.15,2x,'<-- SUM,',I7,
     & ' PROBABILITIES >',E16.9,/,' # ',F22.15,2x,
     &   '<-- TOTAL SUM',/,1x,'#',72('-'))
 106  format (1x,'#',72('-'),/,' # ',F22.15,2x,'<-- SUM',I7,
     & ' PROBABILITIES > ',E16.9,/,' # ',F22.15,2x,
     &   '<-- TOTAL SUM (NON NECESSARILY EQUAL TO 1.0)',
     &   /,1x,'#',72('-'))
 11   format (//1x,'Average populations and localization indices')
 15   format (1x,'<n(',I3,')>               = ',F16.9)
 150  format (1x,'Total population       = ',F16.9)
 151  format (1x,'# <n(',I3,')>_{',a,'} = ',F16.9)
 152  format (
     & ' # <n>_{ALPHA}',6x,'= ',F16.9,/,
     & ' # <n>_{BETA }',6x,'= ',F16.9,/,' # <n>_{ALPHA+BETA} = ',F16.9)
 988  format (/,' # COV(1,2) (ALPHA,BETA,SUM) = ',3(1x,F17.10),/,
     &          ' #  DI(1,2) (ALPHA,BETA,SUM) = ',3(1x,F17.10))
 115  format (1x,'#',/,1x,'# Delocalization indices',/,1x,'#')
 222  format (4(1x,F17.10))
 200  format (/,' # M-BASINS EDFs (COMPLEX AOM VERSION) ',
     & 'FOR ',a,' ELECTRONS',/,1x,'#',72('-'),/,
     &  ' # NUMBER OF GROUPS',13x,' = ',I8,/,
     &  ' # TOTAL NUMBER OF PROBABILITIES = ',I8,/,1x,'#',72('-'),/,
     & ' #     Probability            n1    n2    n3 ...')
 2001 format (/,' # M-BASINS ELECTRON NUMBER PROBABILITY DISTRIBUTION',
     & ' INCLUDING SPIN', /,1x, '# ',72('-'),/,
     & ' # NUMBER OF GROUPS                              =  ', I8,/,
     & ' # TOTAL AND RELEVANT NUMBER OF PROBABILITIES    = ',2(1x,I8),/,
     & ' # Gi(a) Gi(b) ARE THE NUMBER OF ALPHA AND BETA ELECTRONS ',
     & ' IN GROUP i',/,1x,'#',72('-'),/,1x,'#',5x,
     & 'Probability',11x,'G1(a) G1(b) G2(a) G2(b)',/,' #',72('-'))
 1110 format (/,' # ',35('+'),/,' # Average electron populations',/,
     & ' # ',35('+'))
 878  format (' # ',F22.15,2x,2I6)
      end
