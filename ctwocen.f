c
c-----------------------------------------------------------------------
c
      subroutine ctwocen (aom,epsbond,epsdet,probcut,ndets,ncent,nmo,
     &  nel,nval,ncore,nelact,moval,mav,mbv,ifilc,lw,lr,ngroup,ival,
     &  wfnf,largwr)
c
c.....Obtains the two-center delocalization indices between all pairs
c     of atoms of the current molecule described by a multi-determinant
c     wavefunction.
c
      USE         space_for_cidet
      include    'implicit.inc'
      include    'param.inc'
      include    'constants.inc'
      include    'mline.inc'
c
      real(kind=8)  aom(ncent,nmo,nmo)
      real(kind=8),    allocatable,dimension (:,:,:) :: sg
      real(kind=8),    allocatable,dimension (:)     :: pnew
      real(kind=8),    allocatable,dimension (:,:)   :: ptwo
      real(kind=8),    allocatable,dimension (:)     :: ponea,poneb
      integer(kind=4), allocatable,dimension (:,:)   :: indp2
      integer(kind=4), allocatable,dimension (:,:)   :: resnca,resncb
      integer(kind=4), allocatable,dimension (:,:)   :: resnc
      integer(kind=4), allocatable,dimension (:,:)   :: conn
      integer(kind=4), allocatable,dimension (:)     :: ind
      integer(kind=4), allocatable,dimension (:)     :: mogrp
      real(kind=8),    allocatable,dimension (:,:)   :: del,muts
      integer(kind=4)  atom1,atom2
      integer(kind=4)  ival(nmo),lw
*     parameter  (mline=200)
      character*(mline) wfnf
      logical orderp,inlist,largwr
c
      call timer (2,itwodi,'_ctwocen  ',-1)

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
c
c-----Allocate resonance structures
c
      if (.not.allocated(resnca)) then
        allocate ( resnca(nproba,ng), stat=ier ) 
        if (ier.ne.0) stop 'ctwocen.f: Cannot allocate resnca()'
      endif
      if (.not.allocated(resncb)) then
        allocate ( resncb(nprobb,ng), stat=ier ) 
        if (ier.ne.0) stop 'ctwocen.f: Cannot allocate resncb()'
      endif

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

      combi=one
      do i=0,nval-1
        combi=combi*dble(nval+ng-1-i)/dble(nval-i)
      enddo
      nprobt=int(combi+epsq10)
      if (.not.allocated(pnew)) then
        allocate (pnew(nprobt),stat=ier)
        if (ier.ne.0) stop 'ctwocen.f: Cannot allocate pnew()'
      endif
      if (.not.allocated(resnc)) then
        allocate (resnc(nprobt,ng),stat=ier)
        if (ier.ne.0) stop 'ctwocen.f: Cannot allocate resnc()'
      endif
      if (.not.allocated(ind)) then
        allocate (ind(ng),stat=ier)
        if (ier.ne.0) stop 'ctwocen.f: Cannot allocate ind()'
      endif

      n=0
      do ia=1,nproba
        do ib=1,nprobb
          n=n+1
          if (n.eq.1) then
            np=1
            do i=1,ng
              resnc(np,i)=resnca(ia,i)+resncb(ib,i)
            enddo
          else
            do i=1,ng
              ind(i)=resnca(ia,i)+resncb(ib,i)
            enddo
            do j=1,np
              inlist=.true.
              do k=1,ng
                inlist=inlist.and.(ind(k).eq.resnc(j,k))
              enddo
              if (inlist) goto 2
            enddo
            np=np+1
            do i=1,ng
              resnc(np,i)=resnca(ia,i)+resncb(ib,i)
            enddo
          endif
 2      enddo
      enddo
      if (.not.allocated(mogrp)) then
        allocate (mogrp(ng),stat=ier)
        if (ier.ne.0) stop 'ctwocen.f: Cannot allocate mogrp()'
      endif
      orderp=.false.
      mogrp(1:ng)=0
c
c-----Allocated connectivity matrix
c
      if (.not.allocated(conn)) then
        allocate ( conn(ncent,ncent), stat=ier ) 
        if (ier.ne.0) stop 'ctwocen.f: Cannot allocate conn()'
      endif
      if (.not.allocated(del)) then
        allocate ( del(ncent,ncent), stat=ier ) 
        if (ier.ne.0) stop 'ctwocen.f: Cannot allocate del()'
      endif
      if (.not.allocated(muts)) then
        allocate ( muts(ncent,ncent), stat=ier )
        if (ier.ne.0) stop 'ctwocen.f: Cannot allocate muts()'
      endif
      conn = 0
      del  = zero
      muts = zero

      if (.not.allocated(sg)) then
        allocate ( sg(ng,nmo,nmo), stat=ier ) 
        if (ier.ne.0) stop 'ctwocen.f: Cannot allocate sg()'
      endif
      if (.not.allocated(ptwo)) then
        allocate (ptwo(0:nel,0:nel),stat=ier)
        if (ier.ne.0) stop 'ctwocen.f: Cannot allocate ptwo()'
      endif
      if (.not.allocated(ponea)) then
        allocate (ponea(0:nel),stat=ier)
        if (ier.ne.0) stop 'ctwocen.f: Cannot allocate ponea()'
      endif
      if (.not.allocated(poneb)) then
        allocate (poneb(0:nel),stat=ier)
        if (ier.ne.0) stop 'ctwocen.f: Cannot allocate poneb()'
      endif
      if (.not.allocated(indp2)) then
        allocate (indp2(nprobt,2),stat=ier)
        if (ier.ne.0) stop 'ctwocen.f: Cannot allocate indp2()'
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
          call cindgr (epsdet,probcut,nproba,nprobb,ng,ndets,nmo,
     &      ncore,nelact,nel,moval,ival,mogrp,mav,mbv,ifilc,lw,wfnf,
     &      sg,resnca,resncb,orderp,.false.,pnew,nprobt)
          p1    = zero
          p2    = zero
          p12   = zero
          ptwo  = zero
          ponea = zero
          poneb = zero
          ns    = 0
          do i=1,nprobt
            p    = pnew(i)
            if (p.gt.0d0) then
              n1is = resnc(i,1)
              n2is = resnc(i,2)
              p1   = p1+n1is*p
              p2   = p2+n2is*p
              p12  = p12+n1is*n2is*p
              ptwo(n1is,n2is)=ptwo(n1is,n2is)+p
              ponea(n1is)=ponea(n1is)+p
              poneb(n2is)=poneb(n2is)+p

              if (ns.eq.0) then
                ns=1
                indp2(ns,1) = n1is
                indp2(ns,2) = n2is
              else
                do m=1,ns
                  if (n1is.eq.indp2(m,1).and.
     &                n2is.eq.indp2(m,2)) goto 1
                enddo
                ns=ns+1
                indp2(ns,1) = n1is
                indp2(ns,2) = n2is
 1              continue
              endif
            endif
          enddo
c
c-----------------------------------------------------------------------
c
c---------Compute and write entropy (S): 
c                                        |    p_2(na,nb)   |
c         S = SUM(na,nb) p_2(na,nb)  log | --------------- |
c                                        | p_1(na) p_1(nb) |
c-----------------------------------------------------------------------
c
          do m=1,ns
            j = indp2(m,1)
            k = indp2(m,2)
            p2g  = ptwo(j,k)
            p1is = ponea(j)
            p2is = poneb(k)
            if (p2g.gt.0d0.and.p1is.gt.0d0.and.p2is.gt.0d0) then
              muts(atom1,atom2) = muts(atom1,atom2) +
     &        p2g*log(p2g/p1is/p2is)
              if (p2g.gt.1D-4.and.largwr) then
                write (lw,681) j,k,p2g,p1is*p2is,p1is,p2is
              endif
            endif
          enddo
          muts(atom2,atom1)=muts(atom1,atom2)

          del(atom1,atom2)=-two*(p12-p1*p2)
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
c-----Write connectivity indices array
c
*     write (lw,98) 1,2,3,4
*     do i=2,ncent
*       mini=min(6,i-1)
*       write (lw,'(" # ",I3,6F12.6)') I,(del(i,j),j=1,mini)
*       if (mini+1.gt.6)
*    &  write (lw,'(" #    ",6F12.6)') (del(i,j),j=mini+1,i-1)
*     enddo
*     write (lw,*)

*     do ieps=1,10
*       conn(1:ncent,1:ncent)=0
*       do i=1,ncent
*         do j=1,ncent
*           if (del(i,j).lt.epsbond) then
*             conn(i,j)=0
*           else
*             conn(i,j)=1
*           endif
*         enddo
*       enddo
*       write (lw,95) epsbond
*       do i=1,ncent
*         write (lw,'(" #   ",500(1x,I1))') (conn(i,j),j=1,ncent)
*       enddo
*       epsbond=epsbond*half
*     enddo
c
c.....Deallocate arrays
c
      if (allocated(resnca)) then
        deallocate ( resnca, stat=ier ) 
        if (ier.ne.0) stop 'ctwocen.f: Cannot deallocate resnca()'
      endif
      if (allocated(resncb)) then
        deallocate ( resncb, stat=ier ) 
        if (ier.ne.0) stop 'ctwocen.f: Cannot deallocate resncb()'
      endif
      if (allocated(sg)) then
        deallocate ( sg, stat=ier ) 
        if (ier.ne.0) stop 'ctwocen.f: Cannot deallocate sg()'
      endif
      if (allocated(conn)) then
        deallocate ( conn, stat=ier ) 
        if (ier.ne.0) stop 'ctwocen.f: Cannot deallocate conn()'
      endif
      if (allocated(del)) then
        deallocate ( del, stat=ier ) 
        if (ier.ne.0) stop 'ctwocen.f: Cannot deallocate del()'
      endif
      if (allocated(pnew)) then
        deallocate (pnew,stat=ier)
        if (ier.ne.0) stop 'ctwocen.f: Cannot deallocate pnew()'
      endif
      if (allocated(resnc)) then
        deallocate (resnc,stat=ier)
        if (ier.ne.0) stop 'ctwocen.f: Cannot deallocate resnc()'
      endif
      if (allocated(ind)) then
        deallocate (ind,stat=ier)
        if (ier.ne.0) stop 'ctwocen.f: Cannot deallocate ind()'
      endif
      if (allocated(mogrp)) then
        deallocate (mogrp,stat=ier)
        if (ier.ne.0) stop 'ctwocen.f: Cannot deallocate mogrp()'
      endif
      if (allocated(ptwo)) then
        deallocate (ptwo,stat=ier)
        if (ier.ne.0) stop 'ctwocen.f: Cannot deallocate ptwo()'
      endif
      if (allocated(ponea)) then
        deallocate (ponea,stat=ier)
        if (ier.ne.0) stop 'ctwocen.f: Cannot deallocate ponea()'
      endif
      if (allocated(poneb)) then
        deallocate (poneb,stat=ier)
        if (ier.ne.0) stop 'ctwocen.f: Cannot deallocate poneb()'
      endif
      if (allocated(indp2)) then
        deallocate (indp2,stat=ier)
        if (ier.ne.0) stop 'ctwocen.f: Cannot deallocate indp2()'
      endif

      call timer (4,itwodi,'_ctwocen  ',-1)
      return
 99   format (//' # Obtaining two-center DIs & Mutual entropies')
 95   format (' # Connectivity matrix with EPSBOND = ',F12.6)
 98   format (' # Delocalization indices',/,' #',4I12,' ...')
 980  format (' # Mutual Entropy indices',/,' #',4I12,' ...')
 96   format (' # DI & S for pair ( ',2I3,' ) = ',2(1x,F22.15))
 683  format (1x,'# p(',I3,') = ',F7.4)
 681  format (' # p(',2I3,' ) = ',F7.4,' ?= ',F7.4,' = ',F7.4,
     & ' * ', F7.4)
 682  format (' # Entropy (',2I3,' ) = ',E15.8)
 881  format (/' # PAIR ( ',2I3,' )')
c
c.....Formats.
c
      end
