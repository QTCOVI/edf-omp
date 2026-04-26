c
c-----------------------------------------------------------------------
c
      subroutine sopenrho (aom,nfugrp,ifugrp,nprims,ncent,nmo,ngroup,
     &  mal,mbe,epsneg,lw,lerr,ifilc,udat,wfnfile,largwr)
      USE          space_for_wfncoef
      USE          space_for_wfnbasis
      USE          space_for_rdm1
      include     'implicit.inc'
      include     'constants.inc'
      include     'param.inc'
      include     'corr.inc'
      include     'mline.inc'
      real(kind=8) aom(ncent,nmo,nmo)
      integer(kind=4) isigma(nmo)
      real(kind=8) epsneg
      integer(kind=4) nfugrp(ngroup)
      integer(kind=4) ifugrp(ncent,ngroup)
      logical largwr
      real(kind=8),allocatable,dimension (:,:,:) :: aomab
      real(kind=8),allocatable,dimension (:,:)   :: rhoab
      real(kind=8),allocatable,dimension (:,:)   :: sg,shalf,rhop,dno
      real(kind=8),allocatable,dimension (:,:)   :: vecprop,smhalf
      real(kind=8),allocatable,dimension (:,:)   :: vecan
      real(kind=8),allocatable,dimension (:)     :: valprop,sumeigen
      real(kind=8),allocatable,dimension (:)     :: teigen,percent
      real(kind=8),allocatable,dimension (:)     :: work
      integer(kind=4), allocatable,dimension (:) :: elec,eleca,elecb
      real(kind=8),allocatable,dimension (:)     :: totalz,eos
      integer(kind=4), allocatable,dimension (:) :: ilabeig,iord
      real(kind=8), parameter ::   eps   = 1D-14
      integer(kind=4) udat,udatnw,code
      character*(mline) wfnloc
      character*(*) wfnfile
      character*(4) fourchar
      character*(1) spin
c
c-----------------------------------------------------------------------
c
      allocate (eleca(ngroup))
      allocate (elecb(ngroup))
      allocate (elec(ngroup))
      allocate (totalz(ngroup))
      allocate (eos(ngroup))
      write (lw,100)
      do code=1,2
        if (code.eq.1) then
          spin='a'
          write (lw,*) '# ALPHA BLOCK +++++++++++++++++++++++++++++++++'
          nab=nalpha
          nelec=mal
          allocate (rhoab(nab,nab),stat=ier)
          if (ier.ne.0) stop ' # sopenrho.f: Cannot allocate rhoab()'
          isigma(1:nab)=ialpha(1:nab)
          rhoab(1:nab,1:nab)=c1ea(isigma(1:nab),isigma(1:nab))
        else
          spin='b'
          write (lw,*) '# '
          write (lw,*) '# BETA  BLOCK +++++++++++++++++++++++++++++++++'
          nab=nbeta
          nelec=mbe
          allocate (rhoab(nab,nab),stat=ier)
          if (ier.ne.0) stop ' # sopenrho.f: Cannot allocate rhoab()'
          isigma(1:nab)=ibeta(1:nab)
          rhoab(1:nab,1:nab)=c1eb(isigma(1:nab),isigma(1:nab))
        endif
        allocate (aomab(ncent,nab,nab))
        if (ier.ne.0) stop ' # sopenrho.f: Cannot allocate aomab()'
        allocate (sg(nab,nab))
        if (ier.ne.0) stop ' # sopenrho.f: Cannot allocate sg()'
        allocate (shalf(nab,nab))
        if (ier.ne.0) stop ' # sopenrho.f: Cannot allocate shalf()'
        allocate (smhalf(nab,nab))
        if (ier.ne.0) stop ' # sopenrho.f: Cannot allocate smhalf()'
        allocate (rhop(nab,nab))
        if (ier.ne.0) stop ' # sopenrho.f: Cannot allocate rhop()'
        allocate (vecprop(nab,nab))
        if (ier.ne.0) stop ' # sopenrho.f: Cannot allocate vecprop()'
        allocate (vecan(nab,nab))
        if (ier.ne.0) stop ' # sopenrho.f: Cannot allocate vecan()'
        allocate (dno(nab+nab,nprims))
        if (ier.ne.0) stop ' # sopenrho.f: Cannot allocate dno()'
        allocate (valprop(nab))
        if (ier.ne.0) stop ' # sopenrho.f: Cannot allocate valprop()'
        allocate (sumeigen(nab))
        if (ier.ne.0) stop ' # sopenrho.f: Cannot allocate sumeigen()'
        allocate (teigen(nab*ngroup))
        if (ier.ne.0) stop ' # sopenrho.f: Cannot allocate teigen()'
        allocate (percent(nab))
        if (ier.ne.0) stop ' # sopenrho.f: Cannot allocate percent()'
        allocate (ilabeig(nab*ngroup))
        if (ier.ne.0) stop ' # sopenrho.f: Cannot allocate ilabeig()'
        allocate (iord(nab*ngroup))
        if (ier.ne.0) stop ' # sopenrho.f: Cannot allocate iord()'
        allocate (work(3*nab-1))
        if (ier.ne.0) stop ' # sopenrho.f: Cannot allocate work()'
c
c.......Obtain group overlap integrals.
c
        aomab(:,1:nab,1:nab)=aom(:,isigma(1:nab),isigma(1:nab))
        inmo   = nab+nab
        inpr   = nprims
        udatnw = 2
c
c-------Run over the different groups
c
        sumeigen = 0.0d+00
        init=1
        ifin=nab
        do ig=1,ngroup
          sg = 0.0d+00
          totalz(ig)=0.0d+00
          do j=1,nfugrp(ig)
            sg(:,:)=sg(:,:)+aomab(ifugrp(j,ig),:,:)
            totalz(ig)=totalz(ig)+charge(ifugrp(j,ig))
          enddo
c
c---------Diagonalize group overlap matrix
c
          n3=3*nab-1
c
c---------In this diagonalization, sg() is destroyed but it is not used
c         anymore.
c
          call dsyev ('V','U',nab,sg,nab,valprop,work,n3,info)
          shalf = 0.0d+00
          smhalf = 0.0d+00
          do k=1,nab
            if (valprop(k).gt.0.0d+00) then
            elseif (valprop(k).gt.epsneg) then
              write (lerr,222) k,valprop(k),eps+eps
              valprop(k) = eps+eps
            elseif (abs(valprop(k)).lt.eps) then
              valprop(k) = eps+eps
            else
              write (lerr,*) '# sopenrho.f: Negative eigenvalue of sg()'
              valprop(k) = eps
            endif
            do i=1,nab
              do j=1,i
                pvec=sg(i,k)*sg(j,k)
                shalf(i,j)=shalf(i,j)+pvec*sqrt(valprop(k))
                smhalf(i,j)=smhalf(i,j)+pvec/sqrt(valprop(k))
              enddo
            enddo
          enddo
c
c---------Symmetrize shalf() and smhalf().
c
          do i=1,nab
            do j=1,i
              shalf(j,i)=shalf(i,j)
              smhalf(j,i)=smhalf(i,j)
            enddo
          enddo
c
c---------Built S_A^1/2 * 1-RDM * S_A^1/2
c
          rhop = matmul (shalf,matmul(rhoab,transpose(shalf)))
c
c---------Diagonalize it and obtain the eigenvalues lambda.
c
          vecprop=rhop
          call dsyev ('V','U',nab,vecprop,nab,valprop,work,n3,info)
c
c---------Accumulate eigenvalues
c
          teigen(init:ifin)  = valprop(1:nab)
          ilabeig(init:ifin) = ig
          init = init+nab
          ifin = ifin+nab
c
c---------Add eigenvalues of this group to the total
c
          sumeigen = sumeigen+valprop
          write (lw,500) ig,nfugrp(ig),(ifugrp(j,ig),j=1,nfugrp(ig))
          write (lw,*) '# 1RDM Eigenvalues'
          write (lw,445) (valprop(i),i=1,nab)
          write (lw,446) sum(valprop(1:nab))
c
c---------Matrix to express the FNOs in the canonical basis,
c         S^{-1/2} \tilde{U}.
c
          vecan = matmul(smhalf,vecprop)
c
c---------Compute S^{-1}
c
          smhalf=matmul(smhalf,smhalf)
c
          if (largwr) then
            write (lw,'(a)') ' #'
            write (lw,'(a)') ' # Eigenvectors in the canonical MO basis'
            write (lw,'(a)') ' #'
          endif
          do i=1,nab
            if (largwr) then
              write (lw,'(a,I3)') ' # Eigenvector ',i
              write (lw,444) (vecan(j,i),j=1,nab)
            endif
            do j=1,nprims
              tmp = 0.0d+00
              do k=1,nab
                tmp = tmp + vecan(k,i) * coef(isigma(k)+nmo,j)
              enddo
              dno(i,j) = tmp
            enddo
          enddo
c
c---------Overlaps of FNOs in R^3
c
          vecan = matmul(transpose(vecan),vecan)   ! vecan is now <dno|dno>_R^3
          do i=1,nab
            percent(i)=100.0d+00/vecan(i,i)
          enddo
          if (largwr) then
            write (lw,'(a)') ' #'
            write (lw,'(a)') ' # Overlaps of FNOs in R^3'
            write (lw,'(a)') ' #'
            do i=1,nab
              percent(i)=100.0d+00/vecan(i,i)
              write (lw,444) (vecan(i,j),j=1,i)
            enddo
            write (lw,'(a)') ' #'
            write (lw,'(a)') ' #'
            do i=1,nab
              write (lw,'(a,I3,a)') ' # FNO ',i,' in primitive basis'
              write (lw,444) (dno(i,j),j=1,nprims)
            enddo
          endif
          write (lw,26) (percent(i),i=1,nab)
c
c---------Write a WFN file with FNOs. The i^th FNO is multiplied by 
c         vecan(i,i)^{-1/2} such that it is normalized to 1.0 in R^3
c
CCCC      occ(1:nab)=valprop(1:nab)
CCCC      do i=1,nab
CCCC        xnorm=1.0D0/sqrt(vecan(i,i))
CCCC        forall (j=1:nprims) dno(i,j)=dno(i,j)*xnorm
CCCC        forall (j=1:nprims) dno(i+nab,j)=coef(isigma(i)+nmo,j)
CCCC      enddo
CCCC      wfnloc = trim(wfnfile)//"-frag"//fourchar(ig)//trim(spin)
CCCC      write (lw,'(a)') ' # '
CCCC      write (lw,22) trim(wfnloc)
CCCC      write (lw,'(a)') ' # '
CCCC      open (unit = udatnw,file = trim(wfnloc),status='unknown')
CCCC      call wrtwfn (udat,udatnw,lerr,ifilc,dno,inmo,inpr)
CCCC      close (unit = udatnw)
        enddo
        write (lw,23) sum(sumeigen(1:nab))
c
c-------Order the accumulated eigenvalues
c
        forall (i=1:ngroup*nab) iord(i)=i
        call qcksort (teigen, iord, 1,ngroup*nab)
        write (lw,240) 
c       
c-------Compute electrons that are associated to each group.
c       All the eigenvalues are analyzed by decreasing order.
c
        k=0       
        elec=0
        do i=ngroup*nab,1,-1
          ii=iord(i)
          igr=ilabeig(ii)
          write (lw,24) i,teigen(ii),ilabeig(ii)
          k=k+1
          if (k.le.nelec) elec(igr)=elec(igr)+1
          if (k.eq.nelec) exit
        enddo
c
c-------Reliability EOS value
c
        rsigma = 100.0*min(1d0,max(0.0d+00,0.5d0+teigen(k)-teigen(k+1)))
c
c-------Deallocate arrays
c
        deallocate (rhoab,stat=ier)
        if (ier.ne.0) stop ' # sopenrho.f: Cannot deallocate rhoab()'
        deallocate (aomab,stat=ier)
        if (ier.ne.0) stop ' # sopenrho.f: Cannot deallocate aomab()'
        deallocate (sg,stat=ier)
        if (ier.ne.0) stop ' # sopenrho.f: Cannot deallocate sg()'
        deallocate (shalf,stat=ier)
        if (ier.ne.0) stop ' # sopenrho.f: Cannot deallocate shalf()'
        deallocate (smhalf,stat=ier)
        if (ier.ne.0) stop ' # sopenrho.f: Cannot deallocate shalf()'
        deallocate (rhop,stat=ier)
        if (ier.ne.0) stop ' # sopenrho.f: Cannot deallocate rhop()'
        deallocate (vecprop,stat=ier)
        if (ier.ne.0) stop ' # sopenrho.f: Cannot deallocate vecprop()'
        deallocate (vecan,stat=ier)
        if (ier.ne.0) stop ' # sopenrho.f: Cannot deallocate vecan()'
        deallocate (dno,stat=ier)
        if (ier.ne.0) stop ' # sopenrho.f: Cannot deallocate dno()'
        deallocate (valprop,stat=ier)
        if (ier.ne.0) stop ' # sopenrho.f: Cannot deallocate valprop()'
        deallocate (sumeigen,stat=ier)
        if (ier.ne.0) stop ' # sopenrho.f: Cannot deallocate sumeigen()'
        deallocate (teigen,stat=ier)
        if (ier.ne.0) stop ' # sopenrho.f: Cannot deallocate teigen()'
        deallocate (percent,stat=ier)
        if (ier.ne.0) stop ' # sopenrho.f: Cannot deallocate percent()'
        deallocate (ilabeig,stat=ier)
        if (ier.ne.0) stop ' # sopenrho.f: Cannot deallocate ilabeig()'
        deallocate (iord,stat=ier)
        if (ier.ne.0) stop ' # sopenrho.f: Cannot deallocate iord()'
        deallocate (work,stat=ier)
        if (ier.ne.0) stop ' # sopenrho.f: Cannot deallocate work()'
        if (code.eq.1) then
          eleca=elec
          ralpha=rsigma
        else
          elecb=elec
          rbeta=rsigma
        endif
      enddo
c
c-----Effective oxidation state of the group
c
      rsigma=min(ralpha,rbeta)
      do i=1,ngroup
        eos(i)=totalz(i)-eleca(i)-elecb(i)
      enddo
      write (lw,*) '#  OXIDATION STATES OF ALL THE GROUPS'
      write (lw,51) (i,anint(eos(i)),i=1,ngroup)
      write (lw,52) ralpha,rbeta,rsigma
 51   format (' #  Group ',I2,4x,'EOS = ',F15.8)
 52   format (1x,'#  Reliability indices (R_alpha, R_beta, R)% = ',
     &   3(1x,F9.4))
      write (lw,101)
      return
 100  format (//,
     & ' ',/,
     & ' +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++',/,
     & ' +                                                         +',/,
     & ' + B E G I N   O P E N   S Y S T E M S   A N A L Y S I S   +',/,
     & ' +                                                         +',/,
     & ' +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++',/,
     & ' ')
 101  format (//,
     & ' ',/,
     & ' +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++',/,
     & ' +                                                         +',/,
     & ' +    E N D   O P E N   S Y S T E M S   A N A L Y S I S    +',/,
     & ' +                                                         +',/,
     & ' +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++',/,
     & ' ',4/)
 444  format (5(1x,E15.8))
 445  format (5(1x,F15.8))
 20   format (///' # ',77('-'),/,' # FRAGMENT ',I2)
 21   format (' # ',77('-'))
 22   format (' # ',/,
     & ' # Writing file "',a,'". Each FNO is divided by <FNO|FNO>^{1/2}'
     & /' # such that it is normalized to 1.0 in R^3',/,' # ')
 23   format (' # Total number of electrons, N = ',F25.18,/,' #')
 240  format (' # Ordered Eigenvalues of S^1/2 \rho S^1/2',/,
     &        ' # Number     Eigenvalue   Group',/,
     &        ' # -----------------------------')
 24   format (I7,2x,F15.8,2x,I4)
 500  format (' # FRAGMENT ',I2,' FORMED BY ',I3,' ATOMS:',
     & 15I4,/,1000(' # ',33x,15I4,/))
 26   format (/' # Percentage of localization of FNOs in the group',/,
     &   1000(' # ',10(1x,F6.2),/))
 222  format (' # Negative eigenvalue ',I4,' : ',E15.8,
     &  ' is set to ',E17.10)
 223  format (' # Too low eigenvalue ',I4,' of sg(): ',
     &  E15.8,', is set to',E15.8)
 1124 format (' # Eigenvalue ',I4,' of sg() is negative, ',E17.10)
 446  format (' # ELECTRONS =   ',F15.8)
      end
