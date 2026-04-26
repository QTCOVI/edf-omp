c
c-----------------------------------------------------------------------
c
      SUBROUTINE lmowfn (crot,epswfnloc,ndets,nmo,ncore,nelact,nel,
     & mal,mbe,ndetnw,largwr,wfnfile,cicoef,ifilinp,ifilout,
     & stdout,stderr)
c
c.....Determines the expansion coefficients of the wave function in 
c     terms of Slater determinants built in using localized MOs.
c
c.......................................................................
c
      USE         space_for_cidet
      include    'implicit.inc'
      include    'param.inc'
      include    'constants.inc'
      include    'lengrec.inc'
      include    'mline.inc'
      integer     stdout,stderr,leng
      real(kind=8)     deter(2)
      logical     firsta,firstb,largwr
*     parameter   (mline=200)
      character*(mline) dcoef,wfnfile,cicoef
      real   (kind=8)  crot(nmo,nmo)
      real   (kind=8), allocatable,dimension(:)    :: uda,wwa,wwb
      real   (kind=8), allocatable,dimension(:,:)  :: unaoa,unaob
      real   (kind=8), allocatable,dimension(:)    :: cj
      integer(kind=4), allocatable,dimension(:)    :: ordia,ordib
      integer(kind=4), allocatable,dimension(:)    :: comba,combb
      integer(kind=4), allocatable,dimension(:)    :: cjcja,iord
      integer(kind=4), allocatable,dimension(:)    :: indxa,indxb
      integer(kind=4), allocatable,dimension(:,:)  :: icj
c
      call timer (2,ipid,'_lmowfn   ',-1)
c
c.....Renormalize the original wavefunction taking into account 
c     the actual number of Slater determinants that will be used.
c
      write (stdout,112)
c
      cicoef=cicoef(1:leng(cicoef))
      xnorm=zero
      do i=1,ndets
         read (ifilinp,rec=i) (cidet(m),m=0,nelact)
         xnorm=xnorm+cidet(0)*cidet(0)
      enddo
      ynorm=one/sqrt(xnorm)
      if (largwr) write (stdout,'(1x,a,E16.10)') 
     & '# Normalization constant of the input wave function = ',ynorm
c
c.....Allocate memory to store det (U_rj (alpha)).
c
      if (.not.allocated(uda)) then
        allocate (uda(ndets),stat=ier)
        if (ier.ne.0) stop 'lmowfn.f: Cannot allocate uda()'
      endif
c
c.....ncoma is the combinatorial number (NMO choose MAL)
c.....ncomb is the combinatorial number (NMO choose MBE)
c
      combi=one
      do i=0,mal-1
         combi=combi*dble(nmo-i)/dble(mal-i)
      enddo
      ncoma=int(combi)
c
      combi=one
      do i=0,mbe-1
         combi=combi*dble(nmo-i)/dble(mbe-i)
      enddo
      ncomb=int(combi)
      ncomab=int(ncoma*ncomb)
c
c.....Allocate arrays.
c
      if (.not.allocated(unaoa)) then
        allocate    (unaoa(mal,mal),stat=ier)
        if (ier.ne.0) stop 'lmowfn.f: Cannot allocate unaoa()'
      endif
      if (.not.allocated(unaob)) then
        allocate    (unaob(mbe,mbe),stat=ier)
        if (ier.ne.0) stop 'lmowfn.f: Cannot allocate unaob()'
      endif
      if (.not.allocated(ordia)) then
        allocate    (ordia(nmo),stat=ier)
        if (ier.ne.0) stop 'lmowfn.f: Cannot allocate ordia()'
      endif
      if (.not.allocated(ordib)) then
        allocate    (ordib(nmo),stat=ier)
        if (ier.ne.0) stop 'lmowfn.f: Cannot allocate ordib()'
      endif
c
c.....Core orbitals are always occupied.
c
      do mm=1,ncore
        ordia(mm)=mm
        ordib(mm)=mm
      enddo
c
c.....Open file containing the wave functions expansion coefficients
c     in the transformed basis of orbitals.
c
      dcoef=wfnfile(1:leng(wfnfile))//".TCIcoef"
      open (ifilout,file=dcoef,access='direct',
     &    recl=RecLength*(nel+1),form='unformatted')
c
      ndetnw=0
      firsta=.true.
c
      iper=0
      shundred=ncomab/hundred
      xhundred=shundred
      indet=0
c
      if (.not.allocated(indxa)) then
        allocate (indxa(mal),stat=ier)
        if (ier.ne.0) stop 'lmowfn.f: Cannot allocate indxa()'
      endif
      if (.not.allocated(wwa)) then
        allocate (wwa(mal),stat=ier)
        if (ier.ne.0) stop 'lmowfn.f: Cannot allocate wwa()'
      endif
      if (.not.allocated(indxb)) then
        allocate (indxb(mbe),stat=ier)
        if (ier.ne.0) stop 'lmowfn.f: Cannot allocate indxb()'
      endif
      if (.not.allocated(wwb)) then
        allocate (wwb(mbe),stat=ier)
        if (ier.ne.0) stop 'lmowfn.f: Cannot allocate wwb()'
      endif
      if (.not.allocated(comba)) then
        allocate (comba(nmo),stat=ier)
        if (ier.ne.0) stop 'lmowfn.f: Cannot allocate comba()'
      endif
      if (.not.allocated(combb)) then
        allocate (combb(nmo),stat=ier)
        if (ier.ne.0) stop 'lmowfn.f: Cannot allocate combb()'
      endif
c
      do ia=1,ncoma
        call ngen (mal,nmo,nmo,comba,firsta)
        do i=1,ndets
          ma=ncore
          read (ifilinp,rec=i) (cidet(m),m=0,nelact)
          do m=1,nelact
            mm=int(cidet(m))
            if (mm.gt.0) then
              ma=ma+1
              ordia(ma)=mm+ncore
            endif
          enddo
          if (ma.ne.mal) stop 'Inconsistent determinant'
c
          unaoa=zero
          do k=1,mal
            ioka=ordia(k)
            do l=1,mal
              iola=comba(l)
              unaoa(k,l)=crot(ioka,iola)
            enddo
          enddo
          if (mal.gt.izero) then
            call dgeco (unaoa,mal,mal,indxa,rcond,wwa)
            job=10
            call dgedi (unaoa,mal,mal,indxa,deter,wwa,job)
            deta = deter(1)*tenp**deter(2)
          else
            deta=one
          endif
          uda(i)=deta
        enddo
        firstb=.true.
        do ib=1,ncomb
          indet=indet+1
c
          if (indet.ge.xhundred) then
             xhundred=xhundred+shundred
             iper=iper+1
             if (ncomab.gt.10.and.ncomab.le.1000.and.largwr) then
               write (stderr,'(a,i4,a)') ' # Configs :', iper,' % done'
             endif
             call flush(stdout)
          endif
          call ngen (mbe,nmo,nmo,combb,firstb)
c
c.........Run over all the input Slater determinants.
c
          coefj=zero
          do i=1,ndets
            deta=uda(i)
            mb=ncore
            read (ifilinp,rec=i) (cidet(m),m=0,nelact)
            cd=cidet(0)
            do m=1,nelact
              mm=int(cidet(m))
              if (mm.lt.0) then
                mb=mb+1
                ordib(mb)=-mm+ncore
              endif
            enddo   
            if (mb.ne.mbe) stop 'Inconsistent determinant'
            unaob=zero
            do k=1,mbe
              iokb=ordib(k)
              do l=1,mbe
                iolb=combb(l)
                unaob(k,l)=crot(iokb,iolb)
              enddo
            enddo
c
            if (mbe.gt.izero) then
              call dgeco (unaob,mbe,mbe,indxb,rcond,wwb)
              job=10
              call dgedi (unaob,mbe,mbe,indxb,deter,wwb,job)
              detb = deter(1)*tenp**deter(2)
            else
              detb = one
            endif
c
c...........det (U_rj) = det (U_rj (alpha)) * det (U_rj (beta)) 
c
            det=deta*detb
            coefj=coefj+det*cd*ynorm
          enddo
c
c.........Write coefficients and indices of QTAM NAOs in 'ifilout' file.
c
          if (abs(coefj).ge.epswfnloc) then
            ndetnw=ndetnw+1
c
c...........Test that the number of dets in the expansion of the WFN
c           in terms of localized MOs does not exceed MAXNDET.
c
            if (ndetnw.gt.maxndet) then
              write (stderr,*)
              write (stderr,*) "# ndetnw, maxndet = ",ndetnw,maxndet
              write (stderr,*) "# Too many Slater dets in 'lmowfn.f'"
              write (stderr,*) "# Increase 'epswfnloc' in the input,"//
     &                         " or increase 'MAXNDET' in 'param.inc'"
              write (stderr,*)
              stop
            endif
c
            cidet(0)=coefj
            do i=1,mal
              cidet(i) = comba(i)
            enddo 
            do i=1,mbe
              cidet(i+mal) = -combb(i)
            enddo
            if (largwr) write (stdout,211) 
     &      ndetnw,cidet(0),(int(cidet(j)),j=1,nel)
            write (ifilout,rec=ndetnw) (cidet(j),j=0,nel)
          endif
        enddo
      enddo
c
      if (allocated (unaoa)) then
        deallocate (unaoa,stat=ier)
        if (ier.ne.0) stop 'lmowfn.f: Cannot deallocate unaoa()'
      endif
      if (allocated (unaob)) then
        deallocate (unaob,stat=ier)
        if (ier.ne.0) stop 'lmowfn.f: Cannot deallocate unaob()'
      endif
      if (allocated (uda)) then
        deallocate (uda,stat=ier)
        if (ier.ne.0) stop 'lmowfn.f: Cannot deallocate uda()'
      endif
      if (allocated (indxa)) then
        deallocate (indxa,stat=ier)
        if (ier.ne.0) stop 'lmowfn.f: Cannot deallocate indxa()'
      endif
      if (allocated (indxa)) then
        deallocate (indxa,stat=ier)
        if (ier.ne.0) stop 'lmowfn.f: Cannot deallocate indxa()'
      endif
      if (allocated (indxb)) then
        deallocate (indxb,stat=ier)
        if (ier.ne.0) stop 'lmowfn.f: Cannot deallocate indxb()'
      endif
      if (allocated (wwb)) then
        deallocate (wwb,stat=ier)
        if (ier.ne.0) stop 'lmowfn.f: Cannot deallocate wwb()'
      endif
c
c.....Allocate other (possibly very big) arrays.
c
      if (.not.allocated(cj)) then
        allocate    (cj(ndetnw),stat=ier)
        if (ier.ne.0) stop 'lmowfn.f: Cannot allocate cj()'
      endif
      if (.not.allocated(cjcja)) then
        allocate    (cjcja(ndetnw),stat=ier)
        if (ier.ne.0) stop 'lmowfn.f: Cannot allocate cjcja()'
      endif
      if (.not.allocated(iord)) then
        allocate    (iord(ndetnw),stat=ier)
        if (ier.ne.0) stop 'lmowfn.f: Cannot allocate iord()'
      endif
      if (.not.allocated(icj)) then
        allocate    (icj(ndetnw,nel),stat=ier)
        if (ier.ne.0) stop 'lmowfn.f: Cannot allocate icj()'
      endif
c
c.....Renormalize
c
      xnormloc=zero
      do i=1,ndetnw
        read (ifilout,rec=i) cidet(0)
        xnormloc=xnormloc+cidet(0)*cidet(0)
      enddo
      xnormloc=one/sqrt(xnormloc)
c
c.....Re-write 'ifilout' if there are not too many configurations.
c
      if (largwr) then
        write (stdout,*) '#' 
        write (stdout,*) '# Re-ordering configurations'
        write (stdout,*) '#' 
      endif
      do i=1,ndetnw
        iord(i)=i
        read (ifilout,rec=i) (cidet(j),j=0,nel)
        cidet(0)=cidet(0)*xnormloc
        cj(i)=abs(cidet(0))
        if (cidet(0).ge.zero) then
          cjcja(i)=+1
        else
          cjcja(i)=-1
        endif
        do j=1,nel
         icj(i,j)=int(cidet(j))
        enddo
      enddo
      ifir=1
      call qcksort (cj, iord, ifir, ndetnw)
      do i=ndetnw,1,-1
        ioi=iord(i)
        ico=ndetnw-i+1
        coefi=cjcja(ioi)*cj(ioi)
        if (largwr) then
           write (stdout,211) ico,coefi,(icj(ioi,j),j=1,nel)
        endif
        cidet(0)=coefi
        do j=1,nel
          cidet(j)=icj(ioi,j)
        enddo
        if (largwr) write (ifilout,rec=ico) (cidet(j),j=0,nel)
      enddo

      if (allocated (iord)) then
        deallocate (iord,stat=ier)
        if (ier.ne.0) stop 'lmowfn.f: Cannot deallocate iord()'
      endif
      if (allocated (cj)) then
        deallocate (cj,stat=ier)
        if (ier.ne.0) stop 'lmowfn.f: Cannot deallocate cj()'
      endif
      if (allocated (cjcja)) then
        deallocate (cjcja,stat=ier)
        if (ier.ne.0) stop 'lmowfn.f: Cannot deallocate cjcja()'
      endif
      if (allocated (icj)) then
        deallocate (icj,stat=ier)
        if (ier.ne.0) stop 'lmowfn.f: Cannot deallocate icj()'
      endif

      write (stdout,113)
c
      call timer (4,ipid,'_lmowfn   ',-1)
 211  format (1x,'# Conf ',I5,2x,'Coef = ',F16.10,
     &        2x,'Localized MOs: ',100I4)
 112  format (1x,'#',/,1x,'# ',80('+'),/,
     &        1x,'# COMPUTING WFN IN TERMS OF LOCALIZED MOs',/,1x,'#')
 113  format (1x,'# DONE',/,1x,80('+'))
      return
      end
