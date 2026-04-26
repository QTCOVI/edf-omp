c
c-----------------------------------------------------------------------
c
      subroutine casuhfdrvloc (aominp,mxcent,covx,epsneg,maxb,nloctot,
     & mal,mbe,udat,wfnfile,maxcritov,critov,ixcent,largwr,warn,skiph,
     & lw,lerr)
c
c-----------------------------------------------------------------------
c                         
c-----PARAMETERS--------------------------------------------------------
c
c-----AOMINP()                       INPUT
c
c     Atomic Overlap Matrix (AOM) between Canonical MOs in all centers.
c
c.....MXCENT                         INPUT
c
c     Maximum number of N in the search of (Nc,2e) bonds.
c
c.....COVX                           INPUT
c
c     Factor multiplying the sum of covalent radii of two atoms A and B  
c     used in the connect routine to determined whether A and B are 
c     linked or not from a purely geomtrical point of view.
c
c-----EPSNEG                         INPUT
c
c     Diagonalizations in this routine involve definite positive 
c     matrices that should have all of their eigvenvalues positive.
c     Numerical errors, however, make that this not necessarily happens 
c     and some of the eigenvalues can be marginally negative. The value 
c     of 'epsneg' is the largest value (in absolute value) allowed in order
c     that the method ends nicely. For instance, if epsneg=-1D-6
c     and all the eigenvalues are greater that this number the 
c     process continues up to the end. However, if an eigvenvalue
c     is smaller that -1D-6, for instance -1D-5, this routine stops.
c
c-----MAXB                           INPUT
c
c     Maximum number of links between two atoms A and B for which
c     a (2c,2e) localized MOs between these atoms is searched. Atoms
c     directly connected in a geometrical sence have a 'geometrical
c     distance' or link of 1, if A is connected to C which, in turn,
c     is connected to B the 'geometrical distance' or link between
c     A and B is 2, and so on. 
c
c-----NLOCTOT                        OUTPUT
c
c     Total number of Localized MOs (alpha+beta) found.
c
c-----MAL,MBE                        INPUT
c
c     Number of Localized ALPHA and BETA MOs that should be found.
c     If the methods works properly the final OUTPUT value of NLOCTOT
c     should be equal to MAL+MBE
c
c.....UDAT                           INPUT
c
c     Unit number of original WFN file
c
c-----WFNFILE                        INPUT
c
c     Name of the WFN file
c
c-----MAXCRITOV                      INPUT
c
c     Dimension of the arrays critov() and ixcent()
c
c.....CRITOV(1)...CRITOV(MXCENT)     INPUT
c
c     Tolerance used in localization: If an eigvenvalue of
c     S^(1/2) rho S^{1/2} is greater that critov(i) it means that the
c     associated eigenvector corresponds to an MO assumed to be loca-
c     lized in the 'i' centers involved in the diagonalization.
c
c.....IXCENT(1)...IXCENT(MXCENT)     INPUT
c
c     Logical variable: If ixcent(i)=.true. (i-center,2e) bonds are
c     searched, otherwise (i-center,2e) bonds are skipped.
c
c-----LARGWR                         INPUT
c
c-----WARN                           OUTPUT
c
c     Returned  .false. if the routine ends nicely. Otherwise warn is
c     returned .true.
c
c-----SKIPH                          INPUT
c
c     Logical variable: If .true. hydrogen atoms are skipped in one-
c     center localizations. Otherwise hydrogens are not skipped.
c
c-----LW                             INPUT
c
c     Logical output unit
c
c-----LERR                           INPUT
c
c     Logical output unit of error messages
c
c-----------------------------------------------------------------------
c                         
      USE      space_for_wfncoef
      USE      space_for_wfnbasis
      USE      space_for_rdm1
      include 'implicit.inc'
      include 'constants.inc'
      include 'param.inc'
      include 'corr.inc'
      include 'wfn.inc'
      include 'error.inc'
      include 'mline.inc'
      integer(kind=4) code
      real   (kind=8) rho(nmo,nmo)
      integer(kind=4) isigma(nmo)
      real   (kind=8) aom(ncent,nmo,nmo)
      real   (kind=8) aominp(ncent,nmo,nmo)
      real   (kind=8),allocatable, dimension (:,:)   :: rhosig
      real   (kind=8),allocatable, dimension (:,:)   :: aomab
      real   (kind=8),allocatable, dimension (:,:)   :: aomabo
      real   (kind=8),allocatable, dimension (:,:)   :: ca,cb,cao,cbo
      real   (kind=8),allocatable, dimension (:,:,:) :: aomloc,aomloco
      real   (kind=8),allocatable, dimension (:,:)   :: ornor
      real   (kind=8),allocatable, dimension (:,:)   :: c,co
      real   (kind=8),allocatable, dimension (:,:)   :: cret,coret
      real   (kind=8),allocatable, dimension (:,:)   :: ccprev,ccprevo
      real   (kind=8),allocatable, dimension (:,:)   :: dcoef
      real   (kind=8),allocatable, dimension (:)     :: ocup
      logical warn,skiph
c
      real   (kind=8) :: aomnew(ncent,mal+mbe,mal+mbe)
      real   (kind=8) :: aomnewo(ncent,mal+mbe,mal+mbe)
      real   (kind=8) :: cc(mal+mbe,nprims)
      real   (kind=8) :: cco(mal+mbe,nprims)
c
      real   (kind=8)       epsneg
      real   (kind=8)       critov(maxcritov)
      logical               ixcent(maxcritov)
      logical               largwr
      integer(kind=4)       udat,udatnw
      character(len=*)      wfnfile
      character(len=mline)  rdmfile
      character (len=mline) wfnloc
      character (len=1)     spin
c
c-----------------------------------------------------------------------
c
      n = nmo

      aom = aominp

      allocate (ca(mal,nalpha))
      allocate (cb(mbe,nbeta))
      allocate (cao(mal,nalpha))
      allocate (cbo(mbe,nbeta))
c
c-----Run over alpha (code=1) and beta (code=2) blocks
c
      naddmo=0
      nloctot=0
      do code=1,2
        lu1 = 55
        if (code.eq.1) then
          spin='a'
          nab=nalpha
          nloctarget=mal
          ocupation=+1d0
          isigma(1:nab)=ialpha(1:nab)
          rho(1:nmo,1:nmo)=c1ea(1:nmo,1:nmo)
        else
          spin='b'
          nab=nbeta
          nloctarget=mbe
          ocupation=-1d0
          isigma(1:nab)=ibeta(1:nab)
          rho(1:nmo,1:nmo)=c1eb(1:nmo,1:nmo)
        endif
        allocate (rhosig(nab,nab))
        allocate (aomloc(ncent,nab,nab))
        allocate (aomloco(ncent,nab,nab))
        allocate (cret(nloctarget,nab))
        allocate (coret(nloctarget,nab))

C       if (cciqa) then
C         rhosig(1:nab,1:nab)=c1et(1:nab,1:nab)*0.5d0
C       else
C         lu1 = 55
C         rdmfile=trim(wfnfile)//".1rdm"//spin
C         open (lu1,file=trim(rdmfile),form='unformatted',
C    &          iostat=ios,status='old')
C         if (ios.ne.0) write (lerr,100) trim(rdmfile)
C         if (ios.ne.0) stop
C         read (lu1) nmox
C         if (nmox.ne.nmo) then
C           write (lerr,1123) trim(rdmfile),nmox,nmo
C           stop
C         endif
C         rho=0d0
C         do
C           read (lu1,iostat=ios) i,j,value
C           if (ios.ne.0) exit
C           if (i.gt.nmo.or.j.gt.nmo) then
C             stop '# casuhfdrvloc.f: Bad indices of rho() matrix'
C           endif
C           rho(i,j) = value
C           rho(j,i) = value
C         enddo
C         rhosig(1:nab,1:nab)=rho(isigma(1:nab),isigma(1:nab))
C       endif
        rhosig(1:nab,1:nab)=rho(isigma(1:nab),isigma(1:nab))
        aomloc(:,1:nab,1:nab)=aom(:,isigma(1:nab),isigma(1:nab))
        call casuhfoqsloc (aomloc,aomloco,rhosig,cret,coret,code,
     &    nab,nloctarget,mxcent,covx,epsneg,maxb,nloc,udat,wfnfile,
     &    maxcritov,critov,ixcent,largwr,warn,skiph,lw,lerr)

        if (nloc.ne.nloctarget) then
          stop ' casuhfdrvloc.f: Failed localization of Canonical MOs'
        endif
c
c-----------------------------------------------------------------------
c
        inic=naddmo+1
        nfin=naddmo+nloc
        if (code.eq.1) ca=cret
        if (code.eq.2) cb=cret
        if (code.eq.1) cao=coret
        if (code.eq.2) cbo=coret

        allocate (ccprev (nloc,nprims))
        allocate (ccprevo(nloc,nprims))
        allocate (dcoef(nab,nprims))
        do i=1,nab
          dcoef(i,1:nprims) = coef(nmo+isigma(i),1:nprims)
        enddo
        ccprev  = matmul(cret, dcoef)
        ccprevo = matmul(coret,dcoef)
        do i=1,nloc
          cc (i+naddmo,:)  = ccprev (i,:)
          cco(i+naddmo,:)  = ccprevo(i,:)
        enddo
        deallocate (ccprev )
        deallocate (ccprevo)
        deallocate (dcoef  )
c
c-------Compute the overlaps in R^3 between non-orthonormal and
c       orthonormal localized MOs
c
        allocate (ornor(nloc,nloc))
        ornor=matmul(cret,transpose(coret))
        write (lw,200) 'non-orthonormal and orthonormal MOs'
        do i=1,nloc
          write (lw,101)(ornor(i,j),j=1,nloc)
        enddo
        occ(inic:nfin) = ocupation
        aomnew(:,inic:nfin,inic:nfin)=aomloc(:,1:nloc,1:nloc)
        aomnewo(:,inic:nfin,inic:nfin)=aomloco(:,1:nloc,1:nloc)
        naddmo=nloc
        nloctot=nloctot+nloc
        deallocate (rhosig,aomloc,aomloco,ornor,cret,coret)
      enddo
c
      allocate (aomab(mal,mbe))
      allocate (aomabo(nalpha,nbeta))
c
c-----Compute AOM between the localized (non-ortho and orthonormal) MOs
c
      do i=1,ncent
        aomabo = aom(i,ialpha(1:nalpha),ibeta(1:nbeta))
        aomab = matmul (ca,matmul(aomabo,transpose(cb)))
        aomnew(i,1:mal,mal+1:mal+mbe) = aomab
        aomnew(i,mal+1:mal+mbe,1:mal) = transpose(aomab)
      enddo
      do i=1,ncent
        aomabo = aom(i,ialpha(1:nalpha),ibeta(1:nbeta))
        aomab = matmul (cao,matmul(aomabo,transpose(cbo)))
        aomnewo(i,1:mal,mal+1:mal+mbe) = aomab
        aomnewo(i,mal+1:mal+mbe,1:mal) = transpose(aomab)
      enddo
c
c-----Write AOM between the non-orthonormal localized MOs
c
      if (largwr) then
        write (lw,115) 'Non-orthonormal'
        do k=1,ncent
          write (lw,*) '# Center',k
          do i=1,nloctot
            write (lw,101) (aomnew(k,i,j),j=1,i)
          enddo
        enddo
      endif
c
c-----Write AOM between the orthonormal localized MOs
c
      if (largwr) then
        write (lw,115) 'orthonormal'
        do k=1,ncent
          write (lw,*) '# Center',k
          do i=1,nloctot
            write (lw,101) (aomnewo(k,i,j),j=1,i)
          enddo
        enddo
      endif
      deallocate (aomab, aomabo, ca, cb, cao, cbo)
c
c-----write a WFN file with non-orthonormal localized MOs
c
      inpr = nprims
      wfnloc = trim(wfnfile)//"-open"
      udatnw = udat+10
      open (unit=udatnw,file=trim(wfnloc),status='unknown')
      write (lw,*)
      write (lw,351) trim(wfnloc)

      allocate (ocup(nloctot),stat=ier)
      if (ier.ne.0) stop '# casuhfdrvloc.f: Cannot allocate ocup()'
      ocup(1:nloctot) = occ(1:nloctot)
      call cdafhmos (udat,udatnw,0,ocup,cc,nloctot,inpr)
c
c-----write a WFN file with orthonormal localized MOs
c
      inpr = nprims
      wfnloc = trim(wfnfile)//"-open-ortho"
      udatnw = udat+10
      open (unit=udatnw,file=trim(wfnloc),status='unknown')
      write (lw,352) trim(wfnloc)
      call cdafhmos (udat,udatnw,0,ocup,cco,nloctot,inpr)
      deallocate (ocup,stat=ier)
      if (ier.ne.0) stop '# casuhfdrvloc.f: Cannot deallocate ocup()'
c
c-----Write file with AOM between the non-orthonormal localized MOs
c     and between the orthonormal localized MOs
c
      lu18 = udat+11
      wfnloc = trim(wfnfile)//"-open.aom"
      open (unit=lu18,file=trim(wfnloc),status='unknown')
      write (lw,353) trim(wfnloc)
      do k=1,ncent
        write (lu18,'(I4,a)') k,' <=== AOM within this center'
        write (lu18,80) ((aomnew(k,m,j),m=1,j),j=1,nloctot)
      enddo
      close (unit=lu18)
      wfnloc = trim(wfnfile)//"-open-ortho.aom"
      open (unit=lu18,file=trim(wfnloc),status='unknown')
      write (lw,354) trim(wfnloc)
      do k=1,ncent
        write (lu18,'(I4,a)') k,' <=== AOM within this center'
        write (lu18,80) ((aomnewo(k,m,j),m=1,j),j=1,nloctot)
      enddo
      close (unit=lu18)
c
c-----Analyze the MOs and their localizations in all the centers
c
      critics=0.02d0
      call aomtest (aomnew,nloctot,ncent,critics,lw,largwr)
      write (lw,3340) 'Open quantum-systems localized MOs'
      call newdrvmult (nloctot,nprims,ncent,cc,lw)

      return

 100  format (' # casuhfdrvloc.f: File ',a,' NOT FOUND')
 1123 format (' # NMO in file ',a,'(',I4,') # Current value = ',I4)
 351  format (1x,
     &   "# Writing file '",a,"'",1x,"with non-orthonormal MOs")
 352  format (1x,
     &   "# Writing file '",a,"'",1x,"with orthonormal MOs")
 353  format (1x,"# Writing file '",a,"'",1x,
     &   "with AOM for non-orthonormal MOs")
 354  format (1x,"# Writing file '",a,"'",1x,
     &   "with AOM for orthonormal MOs")
 101  format (1000(' # ',5(1x,F15.8),/))
 115  format (' # AOM between ',a,' localized MOs')
 200  format (' # Overlaps in R^3 between ',a,/,
     & ' # Non-orthonormal(CHANGING DOWN)\Orthonormal(CHANGING RIGHT)')
 80   format (6(1x,e16.10))
 3340 format (/,' # Second & fourth moment orbital spread of ',a)
 30   format (5(1x,E15.8))
 33   format (' # MO ',I3)
      end
