c
c-----------------------------------------------------------------------
c
      subroutine twocendi (aom,epsbond,epsdet,probcut,ndets,ncent,nmo,
     &  nel,nval,ncore,nelact,moval,mav,mbv,ifilc,lw,lr,ngroup,ival,
     &  wfnf,largwr)
c
c.....Obtains the two-center delocalization indices between all pairs
c     of atoms of the current molecule described by a single-determi-
c     nant wavefunction.
c
c.......................................................................
c
      USE         space_for_cidet
      include    'implicit.inc'
      include    'param.inc'
      include    'constants.inc'
      include    'mline.inc'
c
      real(kind=8)  aom(ncent,nmo,nmo)
      real(kind=8), allocatable,dimension (:,:,:)  :: sg
      real(kind=8), allocatable,dimension (:,:)    :: pnew
      real(kind=8), allocatable,dimension (:)      :: p1a,p1b
      integer(kind=4), allocatable,dimension (:,:) :: indp2
      integer(kind=4), allocatable,dimension (:,:) :: resnca,resncb
      integer(kind=4), allocatable,dimension (:,:) :: conn
      real(kind=8), allocatable,dimension (:,:)    :: am,tpow,probs
      real(kind=8), allocatable,dimension (:,:)    :: oveaa
      real(kind=8), allocatable,dimension (:,:)    :: del,muts
      real(kind=8), allocatable,dimension (:,:,:)  :: ovea
      real(kind=8), allocatable,dimension (:)      :: w,ww
      integer(kind=4), allocatable,dimension (:)   :: ipvt,indx
      integer(kind=4), allocatable,dimension (:)   :: ordia,ordib
      integer(kind=4) atom1,atom2
c
      real(kind=8)     rcond,dumi,random,deter(2)
      integer     ival(nmo)
      integer     lw
      integer*4   idum
*     parameter  (mline=200)
      character*(mline) wfnf
      logical orderp,largwr
c
      if (icorr) then
        call ctwocen (aom,epsbond,epsdet,probcut,ndets,ncent,nmo,
     &    nel,nval,ncore,nelact,moval,mav,mbv,ifilc,lw,lr,ngroup,ival,
     &    wfnf,largwr)
          return
      endif
c
      call timer (2,itwodi,'_twocendi ',-1)

      if (ncent.le.1) return
      if (ncent.eq.2) ng=2
      if (ncent.gt.2) ng=3

      if (ng.gt.2) then
        nproba=(mav+1)*(mav+2)/2
        nprobb=(mbv+1)*(mbv+2)/2
      else
        nproba=(mav+1)
        nprobb=(mbv+1)
      endif

      if (.not.allocated(indp2)) then
        allocate (indp2(nproba*nprobb,2),stat=ier)
        if (ier.ne.0) stop 'twocendi.f: Cannot allocate indp2()'
      endif
      if (.not.allocated(pnew)) then
        allocate (pnew(0:nel,0:nel),stat=ier)
        if (ier.ne.0) stop 'twocendi.f: Cannot allocate pnew()'
      endif
      if (.not.allocated(p1a)) then
        allocate (p1a(0:nel),stat=ier)
        if (ier.ne.0) stop 'twocendi.f: Cannot allocate p1a()'
      endif
      if (.not.allocated(p1b)) then
        allocate (p1b(0:nel),stat=ier)
        if (ier.ne.0) stop 'twocendi.f: Cannot allocate p1b()'
      endif

      if (.not.allocated(ordia)) then
        allocate ( ordia(nmo), stat=ier ) 
        if (ier.ne.0) stop 'twocendi.f: Cannot allocate ordia()'
      endif
      if (.not.allocated(ordib)) then
        allocate ( ordib(nmo), stat=ier ) 
        if (ier.ne.0) stop 'twocendi.f: Cannot allocate ordib()'
      endif
      mela=0
      melb=0
      read (ifilc,rec=1) (cidet(m),m=0,nel)
      do m=1,nel
        mm=int(cidet(m))
        if (mm.gt.0) then
          mela=mela+1
          ordia(mela)=iabs(mm)
        else
          melb=melb+1
          ordib(melb)=iabs(mm)
        endif
      enddo
c
c-----Allocate resonance structures
c
      if (.not.allocated(resnca)) then
        allocate ( resnca(nproba,ng), stat=ier ) 
        if (ier.ne.0) stop 'twocendi.f: Cannot allocate resnca()'
      endif
      if (.not.allocated(resncb)) then
        allocate ( resncb(nprobb,ng), stat=ier ) 
        if (ier.ne.0) stop 'twocendi.f: Cannot allocate resncb()'
      endif
c
c-----Allocated connectivity matrix
c
      if (.not.allocated(conn)) then
        allocate ( conn(ncent,ncent), stat=ier ) 
        if (ier.ne.0) stop 'twocendi.f: Cannot allocate conn()'
      endif
      if (.not.allocated(del)) then
        allocate ( del(ncent,ncent), stat=ier ) 
        if (ier.ne.0) stop 'twocendi.f: Cannot allocate del()'
      endif
      if (.not.allocated(muts)) then
        allocate ( muts(ncent,ncent), stat=ier ) 
        if (ier.ne.0) stop 'twocendi.f: Cannot allocate muts()'
      endif
      conn(1:ncent,1:ncent)=0
      del(1:ncent,1:ncent)=0d0

      if (ng.gt.2) then
        ij=0
        do i=0,mav
          do j=0,mav-i
            ij=ij+1
            k=mav-i-j
            resnca(ij,1)=i
            resnca(ij,2)=j
            resnca(ij,3)=k
          enddo
        enddo      
        ij=0
        do i=0,mbv
          do j=0,mav-i
            ij=ij+1
            k=mbv-i-j
            resncb(ij,1)=i
            resncb(ij,2)=j
            resncb(ij,3)=k
          enddo
        enddo      
      else
        do i=0,mav
          resnca(i,1)=i
          resnca(i,2)=mav-i
        enddo      
        do i=0,mbv
          resncb(i,1)=i
          resncb(i,2)=mbv-i
        enddo      
      endif

      call semilla (idum)
      idum=-idum
      dumi=random(idum)
      npmax=max(nproba,nprobb)
      if (.not.allocated(sg)) then
        allocate ( sg(ng,nmo,nmo), stat=ier ) 
        if (ier.ne.0) stop 'twocendi.f: Cannot allocate sg()'
      endif
      if (.not.allocated(am)) then
        allocate ( am(npmax,npmax), stat=ier ) 
        if (ier.ne.0) stop 'twocendi.f: Cannot allocate am()'
      endif
      if (.not.allocated(tpow)) then
        allocate ( tpow(npmax,ng), stat=ier ) 
        if (ier.ne.0) stop 'twocendi.f: Cannot allocate tpow()'
      endif
      if (.not.allocated(ipvt)) then
        allocate ( ipvt(npmax), stat=ier ) 
        if (ier.ne.0) stop 'twocendi.f: Cannot allocate ipvt()'
      endif
      if (.not.allocated(w)) then
        allocate ( w(npmax), stat=ier ) 
        if (ier.ne.0) stop 'twocendi.f: Cannot allocate w()'
      endif
      if (.not.allocated(ovea)) then
        allocate ( ovea(ng,nmo,nmo), stat=ier ) 
        if (ier.ne.0) stop 'twocendi.f: Cannot allocate ovea()'
      endif
      if (.not.allocated(oveaa)) then
        allocate ( oveaa(nmo,nmo), stat=ier ) 
        if (ier.ne.0) stop 'twocendi.f: Cannot allocate oveaa()'
      endif
      if (.not.allocated(ww)) then
        allocate ( ww(nmo), stat=ier ) 
        if (ier.ne.0) stop 'twocendi.f: Cannot allocate ww()'
      endif
      if (.not.allocated(indx)) then
        allocate ( indx(nmo), stat=ier ) 
        if (ier.ne.0) stop 'twocendi.f: Cannot allocate indx()'
      endif
      if (.not.allocated(probs)) then
        allocate ( probs(npmax,2), stat=ier ) 
        if (ier.ne.0) stop 'twocendi.f: Cannot allocate probs()'
      endif

      write (lw,99) 
      do atom1=2,ncent
c
c-------Group overlap matrix for atom1
c
        do i=1,nmo
          do j=1,nmo
            sg1 = aom(atom1,i,j)
            sg(1,i,j)=sg1
          enddo
        enddo
        do atom2=1,atom1-1
          if (largwr) write (lw,881) atom1,atom2
c
c---------Group overlap matrix for atom2
c
          if (ng.eq.2) then
            do i=1,nmo
              do j=1,nmo
                deltaij=one
                if (i.ne.j) deltaij=zero
                sg(2,i,j) = deltaij - sg(1,i,j)
              enddo
            enddo
          else
            do i=1,nmo
              do j=1,nmo
                deltaij=one
                if (i.ne.j) deltaij=zero
                sg2 = aom(atom2,i,j)
                sg(2,i,j)=sg2
                sg(ng,i,j) = deltaij - sg(1,i,j) - sg2
              enddo
            enddo
          endif
          do ii=1,2
 2000       am=zero
            if (ii.eq.1) then
               nprobs=nproba
               malmbe=mav
            else
               nprobs=nprobb
               malmbe=mbv
            endif
            if (malmbe.gt.izero) then
              tpow=zero
              do i=1,nprobs
                do k=1,ng-1
                  tpow(i,k)=two*random(idum)-one
                enddo
              enddo
              do i=1,nprobs
                do j=1,nprobs
                  aco=one
                  do k=1,ng-1
                    if (ii.eq.1) nn=resnca(j,k)
                    if (ii.eq.2) nn=resncb(j,k)
                    aco=aco*tpow(i,k)**nn
                  enddo
                  am(i,j)=aco
                enddo
              enddo
              call dgeco (am,npmax,nprobs,ipvt,rcond,w)
              if (one + rcond .eq. one) goto 2000
              do m=1,malmbe
                if (ii.eq.1) then
                  ioma=ordia(m)
                else
                  ioma=ordib(m)
                endif
                do k=1,malmbe
                  if (ii.eq.1) then
                    ioka=ordia(k)
                  else
                    ioka=ordib(k)
                  endif
                  do igr=1,ng
                     ovea(igr,m,k)=sg(igr,ioma,ioka)
                  enddo
                enddo
              enddo
              do n=1,nprobs
                do m=1,malmbe
                  do k=1,malmbe
                    dumi=ovea(ng,m,k)
                    do igr=1,ng-1
                       dumi=dumi+tpow(n,igr)*ovea(igr,m,k)
                    enddo
                    oveaa(m,k)=dumi
                  enddo
                enddo
                call dgeco (oveaa,nmo,malmbe,indx,rcond,ww)
                job=10
                call dgedi (oveaa,nmo,malmbe,indx,deter,ww,job)
                probs(n,ii)=deter(1)*tenp**deter(2)
              enddo
              job=0
              call dgesl (am,npmax,nprobs,ipvt,probs(1,ii),job)
            else
              probs(1,ii)=one
            endif
          enddo
          p1   = zero
          p2   = zero
          p12  = zero
          pnew = zero
          p1a  = zero
          p1b  = zero
          ns   = 0
          do i=1,nproba
            do j=1,nprobb
              p=probs(i,1)*probs(j,2)
              if (p.gt.0d0) then
                n1is=resnca(i,1)+resncb(j,1)
                n2is=resnca(i,2)+resncb(j,2)
                n12=n1is+n2is
                p1=p1+n1is*p
                p2=p2+n2is*p
                p12=p12+n1is*n2is*p
                pnew(n1is,n2is)=pnew(n1is,n2is)+p
                p1a(n1is)=p1a(n1is)+p
                p1b(n2is)=p1b(n2is)+p

                if (ns.eq.0) then
                  ns=1
                  indp2(ns,1) = n1is
                  indp2(ns,2) = n2is
                else
                  do m=1,ns
                    if (n1is.eq.indp2(m,1).and.
     &                  n2is.eq.indp2(m,2)) goto 1
                  enddo
                  ns=ns+1
                  indp2(ns,1) = n1is
                  indp2(ns,2) = n2is
 1                continue
                endif
              endif
            enddo
          enddo

*         write (lw,*)'# Probabilities > 1D-4 for atom ',atom1
*         do i=0,nel
*           if (p1a(i).ge.1D-4) write (lw,683) i,p1a(i)
*         enddo
*         write (lw,*)'# Probabilities > 1D-4 for atom ',atom2
*         do i=0,nel
*           if (p1b(i).ge.1D-4) write (lw,683) i,p1b(i)
*         enddo
c
c-----------------------------------------------------------------------
c
c---------Compute and write entropy (S): 
c                                        |    p_2(na,nb)   |
c         S = SUM(na,nb) p_2(na,nb)  log | --------------- |
c                                        | p_1(na) p_1(nb) |
c-----------------------------------------------------------------------
c
          muts(atom1,atom2)=zero
          do m=1,ns
            j = indp2(m,1)
            k = indp2(m,2)
            p2g  = pnew(j,k)
            p1is = p1a(j)
            p2is = p1b(k)
            if (p2g.gt.0d0.and.p1is.gt.0d0.and.p2is.gt.0d0) then
              muts(atom1,atom2) = muts(atom1,atom2) + 
     &        p2g*log(p2g/p1is/p2is)
              if (p2g.gt.1D-4.and.largwr) then
                write (lw,681) j,k,p2g,p1is*p2is,p1is,p2is
              endif
            endif
          enddo
          muts(atom2,atom1)=muts(atom1,atom2)

          del(atom1,atom2)=-2d0*(p12-p1*p2)
          del(atom2,atom1)=delis
          delis=del(atom1,atom2)

          write (lw,96) atom1,atom2,delis,muts(atom1,atom2)
          if (del(atom1,atom2).lt.epsbond) then
            conn(atom1,atom2)=0
          else
            conn(atom1,atom2)=1
          endif
          conn(atom2,atom1)=conn(atom1,atom2)
        enddo
      enddo
c
c-----Write connectivity indices array
c
      if (largwr) then
        write (lw,*) '#'
        write (lw,*) '# SUMMARY'
        write (lw,98) 1,2,3,4
        do i=2,ncent
          mini=min(6,i-1)
          write (lw,'(" # ",I3,6F12.6)') I,(del(i,j),j=1,mini)
          if (mini+1.gt.6)
     &    write (lw,'(" #    ",6F12.6)') (del(i,j),j=mini+1,i-1)
        enddo
        write (lw,*)
        write (lw,980) 1,2,3,4
        do i=2,ncent
          mini=min(6,i-1)
          write (lw,'(" # ",I3,6F12.6)') I,(muts(i,j),j=1,mini)
          if (mini+1.gt.6)
     &    write (lw,'(" #    ",6F12.6)') (muts(i,j),j=mini+1,i-1)
        enddo
      endif
      write (lw,*)
      do ieps=1,10
        conn(1:ncent,1:ncent)=0D0
        do i=1,ncent
          do j=1,ncent
            if (del(i,j).lt.epsbond) then
              conn(i,j)=0
            else
              conn(i,j)=1
            endif
          enddo
        enddo
        write (lw,95) epsbond
        do i=1,ncent
          write (lw,'(" #   ",500(1x,I1))') (conn(i,j),j=1,ncent)
        enddo
        epsbond=epsbond*half
      enddo
c
c.....Deallocate arrays
c
      if (allocated(ordia)) then
        deallocate ( ordia, stat=ier ) 
        if (ier.ne.0) stop 'twocendi.f: Cannot deallocate ordia()'
      endif
      if (allocated(ordib)) then
        deallocate ( ordib, stat=ier ) 
        if (ier.ne.0) stop 'twocendi.f: Cannot deallocate ordib()'
      endif
      if (allocated(resnca)) then
        deallocate ( resnca, stat=ier ) 
        if (ier.ne.0) stop 'twocendi.f: Cannot deallocate resnca()'
      endif
      if (allocated(resncb)) then
        deallocate ( resncb, stat=ier ) 
        if (ier.ne.0) stop 'twocendi.f: Cannot deallocate resncb()'
      endif
      if (allocated(sg)) then
        deallocate ( sg, stat=ier ) 
        if (ier.ne.0) stop 'twocendi.f: Cannot deallocate sg()'
      endif
      if (allocated(am)) then
        deallocate ( am, stat=ier ) 
        if (ier.ne.0) stop 'twocendi.f: Cannot deallocate am()'
      endif
      if (allocated(tpow)) then
        deallocate ( tpow, stat=ier ) 
        if (ier.ne.0) stop 'twocendi.f: Cannot deallocate tpow()'
      endif
      if (allocated(ipvt)) then
        deallocate ( ipvt, stat=ier ) 
        if (ier.ne.0) stop 'twocendi.f: Cannot deallocate ipvt()'
      endif
      if (allocated(w)) then
        deallocate ( w, stat=ier ) 
        if (ier.ne.0) stop 'twocendi.f: Cannot deallocate w()'
      endif
      if (allocated(ovea)) then
        deallocate ( ovea, stat=ier ) 
        if (ier.ne.0) stop 'twocendi.f: Cannot deallocate ovea()'
      endif
      if (allocated(oveaa)) then
        deallocate ( oveaa, stat=ier ) 
        if (ier.ne.0) stop 'twocendi.f: Cannot deallocate oveaa()'
      endif
      if (allocated(ww)) then
        deallocate ( ww, stat=ier ) 
        if (ier.ne.0) stop 'twocendi.f: Cannot deallocate ww()'
      endif
      if (allocated(indx)) then
        deallocate ( indx, stat=ier ) 
        if (ier.ne.0) stop 'twocendi.f: Cannot deallocate indx()'
      endif
      if (allocated(probs)) then
        deallocate ( probs, stat=ier ) 
        if (ier.ne.0) stop 'twocendi.f: Cannot deallocate probs()'
      endif
      if (allocated(conn)) then
        deallocate ( conn, stat=ier ) 
        if (ier.ne.0) stop 'twocendi.f: Cannot deallocate conn()'
      endif
      if (allocated(del)) then
        deallocate ( del, stat=ier ) 
        if (ier.ne.0) stop 'twocendi.f: Cannot deallocate del()'
      endif
      if (allocated(indp2)) then
        deallocate (indp2,stat=ier)
        if (ier.ne.0) stop 'twocendi.f: Cannot deallocate indp2()'
      endif
      if (allocated(pnew)) then
        deallocate (pnew,stat=ier)
        if (ier.ne.0) stop 'twocendi.f: Cannot deallocate pnew()'
      endif
      if (allocated(p1a)) then
        deallocate (p1a,stat=ier)
        if (ier.ne.0) stop 'twocendi.f: Cannot deallocate p1a()'
      endif
      if (allocated(p1b)) then
        deallocate (p1b,stat=ier)
        if (ier.ne.0) stop 'twocendi.f: Cannot deallocate p1b()'
      endif

      call timer (4,itwodi,'_twocendi ',-1)
      return
 99   format (//' # Obtaining two-center DIs & Mutual entropies')
 95   format (' # Connectivity matrix with EPSBOND = ',F12.6)
 98   format (' # Delocalization indices',/,' #',4I12,' ...')
 980  format (' # Mutual Entropy indices',/,' #',4I12,' ...')
 96   format (' # DI & S for pair ( ',2I3,' ) = ',2(1x,E22.15))
 683  format (1x,'# p(',I3,') = ',F7.4)
 681  format (' # p(',2I3,' ) = ',F7.4,' ?= ',F7.4,' = ',F7.4,
     & ' * ', F7.4)
 682  format (' # Entropy (',2I3,' ) = ',E15.8)
 881  format (/' # PAIR ( ',2I3,' )')
c
c.....Formats.
c
      end
