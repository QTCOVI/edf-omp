      subroutine sym (luw,largwr)
c
c.....computation of molecular symmetry from functions.
c
*     USE               mod_prec, only : dp, i4=>ip
      USE               space_for_wfnbasis
      USE               space_for_sym
      include           'implicit.inc'
      include           'param.inc'
      include           'integ.inc'
      include           'wfn.inc'
      include           'sym.inc'
      include           'datatm.inc'
c
c.....Parameters
c
      real(kind=8), parameter :: eps=1d-6
c
c.....Arguments
c
      logical, intent(in) :: largwr
c
c.....Local vars
c
      logical :: equiv
      character(len=6) :: symbol
      real(kind=8) :: xc(3) 
      real(kind=8), allocatable,dimension (:,:)  :: vector,vector2
      integer(kind=4) :: indx(2)
      real(kind=8) :: rplan(3,3), paxis(3), raux(3,3)
      real(kind=8) :: vec(3),vaux(3),aux(3)
      logical :: inversion, axis, ltemp, itis
      integer(kind=4) :: aiaxis
c
c.....Inline functions
c
      del(a,b)=abs(a-b)
      dis2(x,y,z,x1,y1,z1)=sqrt((x-x1)**2+(y-y1)**2+(z-z1)**2)
      vmod(x,y,z)=sqrt(x*x+y*y+z*z)
c
      if (.not.allocated(vector)) then
        allocate (vector(3,ncent),stat=ier)
        if (ier.ne.0) stop 'sym.f: Cannot allocate vector()'
      endif
      if (.not.allocated(vector2)) then
        allocate (vector2(3,ncent),stat=ier)
        if (ier.ne.0) stop 'sym.f: Cannot allocate vector2()'
      endif
c
c.....Write header 
c
      WRITE (LUW,112)
 112  FORMAT (' |',/,
     !' | ---------------------------------------------',/,' |',/,
     !' | B E G I N   S Y M M E T R Y   A N A L Y S I S',/,' |',/,
     !' | ---------------------------------------------',/,' |')
c
c.....Init data
c
      masksym(-1,1:maxis)=0
      masksym(+1,1:maxis)=0
      linear=.false.
      planar=.false.
c
c.....compute center of nuclear charge and displace.
c
      xc(:) = 0.0d0
      sumc  = 0.0d0
      do i=1,ncent
         ci=charge(i)
         sumc=sumc+ci
         xc(:)=xc(:)+xyz(i,:)*ci
      enddo
      xc(:)=xc(:)/sumc
      do i=1,3
         do j=1,ncent
            vector(i,j)=xyz(j,i)-xc(i)
         enddo
      enddo
c
c.....symmetry?
c
      if (nosym) then
         write (luw,*) '| Switching off symmetry module'
         nop=1
         call riden (oper(1,1,1),3)
         idict(1:ncent,1)=i
         idict(1:ncent,2)=i
         nmono=ncent
         goto 666
      endif
      if (lpointgroup) goto 667
c
c.....analyze dimensions of molecule.
c
      ndim=3
c
c     linear?
c
      s=0.0d0
      do i=1,ncent
        do j=i+1,ncent
          sa=vmod(vector(1,i),vector(2,i),vector(3,i))
          sa=sa*vmod(vector(1,j),vector(2,j),vector(3,j))
          if (sa.gt.eps) then
            sa=(vector(1,i)*vector(1,j)+vector(2,i)*vector(2,j)+
     &          vector(3,i)*vector(3,j))/sa
            s=s+abs(abs(sa)-1.0d0)
          endif
        enddo
      enddo
c
c.....linear case
c
      if (s.le.eps) then
         linear=.true.
         inversion=.true.
         ndim=1
         do 22 i=1,ncent
           ngo=0
           do j=1,ncent
             if (dis2(-vector(1,i),-vector(2,i),-vector(3,i),
     &                 vector(1,j), vector(2,j), vector(3,j)).le.eps) 
     &                 ngo=ngo+1
           enddo
           if (ngo.ne.1) then
              inversion=.false.
              goto 22
           endif    
 22      continue
         nop=1
         call riden (oper(1,1,1),3)
         if (inversion) then
           nop=2
           call riden (oper(1,1,2),3)
           oper(1,1,2)=-1.0d0
           oper(2,2,2)=-1.0d0
           oper(3,3,2)=-1.0d0
         endif
      endif
c
c     planar ?
c
      if (.not.linear) then
        s=0.0d0
        axis=.false.
        do i=1,ncent
           do j=i+1,ncent
              if (.not.axis) then
                call vecprod (vector(1,i),vector(1,j),paxis)
                xmod=vmod(paxis(1),paxis(2),paxis(3))
                if (xmod.gt.eps) then
                   axis=.true.
                   paxis(:)=paxis(:)/xmod
                endif
              endif
              do k=j+1,ncent
                 s=s+abs(r3volp(vector(1,i),vector(1,j),vector(1,k)))
              enddo
            enddo
        enddo
c
c.......planar case
c
        if (s.le.eps) then
           planar=.true.
           ndim=2
c
c          axis rotation
c
           call rotz (rplan,paxis)
c
c          rotate molecule
c
           do n=1,ncent
              do i=1,3
                 vaux(i)=0.0d0
                 do j=1,3
                    vaux(i)=vaux(i)+rplan(i,j)*vector(j,n)
                 enddo
              enddo
              vector2(1,n)=vaux(1)
              vector2(2,n)=vaux(2)
           enddo
        endif
      endif
c
c.....obtain the rotation matrices.
c
      if (ndim.gt.1) then
        indx(1)=0
        indx(2)=ncent
        if (ndim.eq.2) then
          call pntgrp (ndim,vector2,1,indx,tolsym,48,nop,oper2)
          nopold=nop
          nop=0
          do i=1,nopold
             nop=nop+1
             call riden (oper(1,1,nop),3)
             call riden (oper(1,1,nop+1),3)
             do j=1,2
                do k=1,2
                   oper(j,k,nop  )=oper2(j,k,i)
                   oper(j,k,nop+1)=oper2(j,k,i)
                enddo
             enddo
             oper(3,3,nop)  =+1.0d0
             oper(3,3,nop+1)=-1.0d0
             nop=nop+1
          enddo
          do i=1,nop
             call triprod (rplan,oper(1,1,i),raux,3)             
             call requal (oper(1,1,i),raux,3)
          enddo
        else
          call pntgrp (ndim,vector,1,indx,tolsym,48,nop,oper)
        endif
      endif
c
c.....Manually set point group, by now only most common groups
c     This help for example for benzene, as programs use d2h
c     for symmetry in correlated calcs as it is an abelian group
c     and no degenerancy is present
c
 667  if (lpointgroup) then
        oper(:,:,:)= 0.0d0
        if (cpointgroup(1:3).eq.'d2h') then
          ndim=3
          nop=8
          ! op1
          oper(1,1,1)= 1.0d0
          oper(2,2,1)= 1.0d0
          oper(3,3,1)= 1.0d0
          ! op2
          oper(1,1,2)= 1.0d0
          oper(2,2,2)= 1.0d0
          oper(3,3,2)=-1.0d0
          ! op3
          oper(1,1,3)=-1.0d0
          oper(2,2,3)=-1.0d0
          oper(3,3,3)= 1.0d0
          ! op4
          oper(1,1,4)=-1.0d0
          oper(2,2,4)=-1.0d0
          oper(3,3,4)=-1.0d0
          ! op5
          oper(1,1,5)=-1.0d0
          oper(2,2,5)= 1.0d0
          oper(3,3,5)= 1.0d0
          ! op6
          oper(1,1,6)=-1.0d0
          oper(2,2,6)= 1.0d0
          oper(3,3,6)=-1.0d0
          ! op7
          oper(1,1,7)= 1.0d0
          oper(2,2,7)=-1.0d0
          oper(3,3,7)= 1.0d0
          ! op8
          oper(1,1,8)= 1.0d0
          oper(2,2,8)=-1.0d0
          oper(3,3,8)=-1.0d0
        else if (cpointgroup(1:3).eq.'c2v') then
          ndim=3
          nop=4
          ! op1
          oper(1,1,1)= 1.0d0
          oper(2,2,1)= 1.0d0
          oper(3,3,1)= 1.0d0
          ! op2
          oper(1,1,2)=-1.0d0
          oper(2,2,2)= 1.0d0
          oper(3,3,2)= 1.0d0
          ! op3
          oper(1,1,3)= 1.0d0
          oper(2,2,3)=-1.0d0
          oper(3,3,3)= 1.0d0
          ! op4
          oper(1,1,4)=-1.0d0
          oper(2,2,4)=-1.0d0
          oper(3,3,4)= 1.0d0
        else if (cpointgroup(1:2).eq.'c1') then
          ndim=3
          nop=1
          ! op1
          oper(1,1,1)= 1.0d0
          oper(2,2,1)= 1.0d0
          oper(3,3,1)= 1.0d0
        else if (cpointgroup(1:2).eq.'ci') then
          ndim=3
          nop=2
          ! op1
          oper(1,1,1)= 1.0d0
          oper(2,2,1)= 1.0d0
          oper(3,3,1)= 1.0d0
          ! op2
          oper(1,1,1)=-1.0d0
          oper(2,2,1)=-1.0d0
          oper(3,3,1)=-1.0d0
        else if (cpointgroup(1:2).eq.'c2') then
          ndim=3
          nop=2
          ! op1
          oper(1,1,1)= 1.0d0
          oper(2,2,1)= 1.0d0
          oper(3,3,1)= 1.0d0
          ! op2
          oper(1,1,1)=-1.0d0
          oper(2,2,1)=-1.0d0
          oper(3,3,1)= 1.0d0
        else if (cpointgroup(1:2).eq.'d2') then
          ndim=3
          nop=4
          ! op1
          oper(1,1,1)= 1.0d0
          oper(2,2,1)= 1.0d0
          oper(3,3,1)= 1.0d0
          ! op2
          oper(1,1,1)= 1.0d0
          oper(2,2,1)=-1.0d0
          oper(3,3,1)=-1.0d0
          ! op3
          oper(1,1,1)=-1.0d0
          oper(2,2,1)= 1.0d0
          oper(3,3,1)=-1.0d0
          ! op4
          oper(1,1,1)=-1.0d0
          oper(2,2,1)=-1.0d0
          oper(3,3,1)= 1.0d0
        else if (cpointgroup(1:2).eq.'cs') then
          ndim=3
          nop=2
          ! op1
          oper(1,1,1)= 1.0d0
          oper(2,2,1)= 1.0d0
          oper(3,3,1)= 1.0d0
          ! op2
          oper(1,1,1)= 1.0d0
          oper(2,2,1)= 1.0d0
          oper(3,3,1)=-1.0d0
        else if (cpointgroup(1:3).eq.'c2h') then
          ndim=3
          nop=4
          ! op1
          oper(1,1,1)= 1.0d0
          oper(2,2,1)= 1.0d0
          oper(3,3,1)= 1.0d0
          ! op2
          oper(1,1,2)=-1.0d0
          oper(2,2,2)=-1.0d0
          oper(3,3,2)=-1.0d0
          ! op3
          oper(1,1,3)= 1.0d0
          oper(2,2,3)= 1.0d0
          oper(3,3,3)=-1.0d0
          ! op4
          oper(1,1,4)=-1.0d0
          oper(2,2,4)=-1.0d0
          oper(3,3,4)= 1.0d0
        end if
      end if
c
c.....clean operations. Write matrices. Obtain pointgroup
c
      if (largwr) write (luw,'(1x,a,i1)') "| Dimension : ", ndim
 666  do i=1,nop
         call typeop (oper(1,1,i),itipo,iaxis,iaxis1,irepet,
     &                eulang(1,i),vaxis(1,i))
c
         if (iaxis.gt.maxis) stop '# sym.f: Axis not supported'
         ioper(i,1)=itipo
         ioper(i,2)=iaxis        !rot-refl
         ioper(i,3)=iaxis1       !rot-inv
         aiaxis=abs(iaxis)
         masksym(itipo,aiaxis)=masksym(itipo,aiaxis)+1
         ! write
         if (largwr) then
           write (luw,15) i,itipo,iaxis,iaxis1,
     &           (vaxis(j,i),j=1,3),(eulang(j,i)*180.0d0/pi,j=1,3)
           do j=1,3
             write (luw,'(" | ",3f12.6)') (oper(j,k,i),k=1,3)
           enddo
         endif
      enddo
      call symgr (masksym,maxis,linear,symbol)
      write (luw,*)
      write (luw,'(1x,a,a)') '| Point Group is ', trim(symbol)
c
c.....obtain inverses
c
      do i=1,nop
        do j=1,nop
          ltemp=.true.
          do k=1,3
            do l=1,3
              ltemp=ltemp.and.abs(oper(k,l,i)-oper(l,k,j)).lt.1d-6
            enddo
          enddo
          if (ltemp) then
             invop(i)=j
             exit
          endif
        enddo
      enddo
c
c.....obtain nonequivalent atoms
c.....idx (i,2,2) is operation taking noneqiv to atom.
c
      mult(1:ncent)=0
      neq=0
      do i=1,ncent
         equiv=.false.
         zi=charge(i)
         do j=1,neq
            jj=ineq(j)
            zj=charge(jj)
            vec(:)=vector(:,jj)
            do n=2,nop
               do k=1,3
                  vaux(k)=0.0d0
                  do l=1,3
                    vaux(k)=vaux(k)+oper(k,l,n)*vec(l)
                  enddo
               enddo
               temp=dis(vector(1,i),vaux)
               if (temp.lt.eps .and. del(zi,zj).lt.eps) then
                  if (.not.equiv) then
                     equiv=.true.
                     mult(j)=mult(j)+1
                     idx(i,1)=j
                     idx(i,2)=n
                  endif
                  if (invop(n).eq.n) idx(i,2)=n ! to find two-orbits.
               endif
            enddo
         enddo
         if (.not.equiv) then
            neq=neq+1
            ineq(neq)=i
            mult(neq)=1
            idx(i,1)=neq
            idx(i,2)=1
         endif
      enddo
c
c.....Obtain integer multiplication table
c
      do i=1,ncent
         table(i,1)=i
         do n=2,nop
            do k=1,3
               aux(k)=0.0d0
               do l=1,3
                  aux(k)=aux(k)+oper(k,l,n)*vector(l,i)
               enddo
            enddo
            do j=1,ncent
               if (dis(vector(1,j),aux).lt.eps) table(i,n)=j
            enddo
         enddo
      enddo
c
c-----------------------------------------------------------------
c
      if (trim(symbol).ne.'C1') then
        write (luw,81) (invop(i),i=1,nop)
        write (luw,82) neq, (ineq(i),i=1,neq)
        write (luw,84) (mult(i),i=1,neq)
        write (luw,83) (idx(i,1),i=1,ncent)
        write (luw,85) (idx(i,2),i=1,ncent)
        write (luw,86) 
        do i=1,ncent
           write (luw,'(" | ",40i5)') (table(i,k),k=1,nop)
        enddo
      endif
c
c.....obtain nonequivalent pairs.
c
      indp(1:ncent,1:ncent)=0
      npair=0
      do i=1,neq
         ii=ineq(i)
         do j=1,ncent
            equiv=.false.
            if (j.ne.ii) then
               do k=1,npair
                  ip=ipair(k,1)
                  jp=ipair(k,2)
                  do l=1,nop
                     m=table(ip,l)
                     n=table(jp,l)
                     ltemp=(m.eq.ii.and.n.eq.j).or.(n.eq.ii.and.m.eq.j)
                     if (ltemp) then
                        equiv=.true.
                        indp(i,j)=k
                        if (n.eq.ii.and.m.eq.j) indp(i,j)=-k
                        goto 44
                     endif
                  enddo
               enddo
 44            if (.not.equiv) then
                  npair=npair+1
                  ipair(npair,1)=ii
                  ipair(npair,2)=j
                  indp(i,j)=npair
               endif
            endif
         enddo
      enddo
c
      write (luw,91) npair
      if (trim(symbol).ne.'C1') then
        write (luw,921) 
        write (luw,92) (ipair(i,1),ipair(i,2),i=1,npair)
        write (luw,93) 
        do i=1,neq
           write (luw,*) 
           write (luw,'(10i5)') (indp(i,j),j=1,ncent)
        enddo
      endif
cAQUI
c
c.....Transform orbitals ?
c
      if (.not. nosym) then 
         if (.not. icorr) then
            call transor (0,luw,largwr,symbol)       !natural orbitals 
         else
            call transor (1,luw,largwr,symbol)       !canonical orbitals
         endif
      else
         do i=1,nmo
           irrep(i,0)=1
           irrep(i,1)=i
           chartab(i,1,1)=1
         enddo
      endif
c
c.....List of atoms to compute nuclear attraction components.
c
      jpair(1:neq,0)=0
      do j=1,npair 
c
c........labels
c
         ip=ipair(j,1)
         jp=ipair(j,2)
         in=idx(ip,1)
         jn=idx(jp,1)
         iop=idx(jp,2)
c
c........direct 1,2
c
         jpair(in,0)=jpair(in,0)+1
         jpair(in,jpair(in,0))=jp
c
c........inverse 2,1: (i) Only for non two-orbits.
c
         if (invop(iop).ne.iop .or. in.ne.jn) then
            iop=invop(iop)
            icc=table(ip,iop)
            jpair(jn,0)=jpair(jn,0)+1
            jpair(jn,jpair(jn,0))=-icc
         endif
      enddo
c
c.....List of centers to be obtained/rotated
c
      do i=1,neq
         ic=ineq(i)
         idict(ic,1)=i
         idict(i,2)=ic
      enddo
      nmono=neq
      do i=1,neq
         do j=1,jpair(i,0)
            jc=abs(jpair(i,j))
            itis=.false.
            do k=1,nmono
               if (idict(k,2).eq.jc) then
                  itis=.true.
                  goto 444
               endif
            enddo
 444        if (.not. itis) then
               nmono=nmono+1
               idict(jc,1)=nmono
               idict(nmono,2)=jc
            endif
         enddo
      enddo
c
c----------------------------------------------------------------
c
c.....Maximum value of jpair(in,0), in=1,neq
c
      mnnuc = 0
      do i=1,neq
        mnnuc = max (jpair(i,0),mnnuc)
      enddo
c
c.....list of pairs of distances
c
      write (luw,110)
      do i=1,npair
         ic=ipair(i,1)
         jc=ipair(i,2)
         dist=dis(vector(1,ic),vector(1,jc))
         write(luw,1110) i,atnam(ic)(1:4),ic,atnam(jc)(1:4),jc,dist
      enddo
c
c-----Deallocate vector() and vector2() arrays.
c
      if (allocated(vector)) then
        deallocate (vector,stat=ier)
        if (ier.ne.0) stop 'sym.f: Cannot deallocate vector()'
      endif
      if (allocated(vector2)) then
        deallocate (vector2,stat=ier)
        if (ier.ne.0) stop 'sym.f: Cannot deallocate vector2()'
      endif

      WRITE (LUW,113)
 113  FORMAT (' |',/,
     !' | -----------------------------------------------',/,' |',/,
     !' | E N D   O F   S Y M M E T R Y   A N A L Y S I S',/,' |',/,
     !' | -----------------------------------------------',/,' |')
      return
c
c.....Formats
c
 15   format (1x,'| Symmetry Operation --> ',i3,/
     &        1x,'| Operation type (proper/improper):',
     &        ' order (rot-ref),(rot-inv)',  i3, ':', i3, ',', i3,/
     &        1x,'| axis vector:', 3f12.6,/,
     &        1x,'| Euler angles (rot-inv):  (alp,bet,gam):',
     &        3f12.6,/ 1x,'| Matrix:') 
 81   format (1x,'| Inverse Operations:',40i4)
 82   format (1x,'| Number of Inequivalent atoms: ',i4,/
     &        1x,'| Ineqs                    :', 1000(/,' |',40i4))
 84   format (1x,'| Multiplicity of each Ineq:',1000(/,' |',40i4))
 83   format (1x,
     & '| Atoms are equivalent to ineq (in order):',1000(/,' |',40i4))
 85   format(1x,
     & '| Operations transforming ineq to atom   :',1000(/,' |',40i4))
 86   format (1x,'| Multiplication table: m(i,j): op j takes atom i',
     &           ' into atom m(i,j)')
 91   format (1x,'| Number of inequivalent pairs:', i5)
 921  format (1x,'| Pairs:')
 92   format (' |',10(1x,'(',i3,',',i3,')'))
 93   format (' | Inequivalent pairs p(i,j): ineqp for pair (ineq i,j)')
 110  format (/1x,'| Distances between inequivalent atomic pairs:')   
 111  format (1x,2x,i5,': ',a2,x,i3,2x,a2,x,i3,x,f14.8)
 1110 format (1x,2x,i5,': ',a4,x,i3,2x,a4,x,i3,x,f14.8)
      end subroutine sym
c---------------------------------------------------------------
      function dis(a,b)
      include 'implicit.inc'
      real(kind=8) a(3),b(3)
      dis=0.0d0
      do i=1,3
         dis=dis+(a(i)-b(i))**2
      enddo
      dis=sqrt(dis)
      return
      end
