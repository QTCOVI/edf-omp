      subroutine transor (lcor,luw,largwr,symbol)
c
c.....computation of molecular symmetry from functions.
c
      USE                space_for_wfnbasis
      USE                space_for_sym
      include           'implicit.inc'
      include           'param.inc'
      include           'integ.inc'
      include           'wfn.inc'
      include           'sym.inc'
c
      parameter         (maxtrials=20)
c
      real(kind=8),allocatable,dimension (:,:) :: forbi, forbf, aux
      real(kind=8),allocatable,dimension (:,:) :: vec
      real(kind=8),allocatable,dimension (:) :: rho, rwork
      integer*4,allocatable,dimension (:) :: iwork
      real(kind=8),allocatable,dimension (:,:,:) :: orboper
      character (len=*) symbol
      dimension         xdet(2), vec1(3)
c
      integer           lcor
      logical           inversion, axis
      logical           ltemp,itis
      logical           largwr
      integer           aiaxis
c
      allocate(forbi(nmo,nmo),forbf(nmo,nmo),vec(3,nmo))
      allocate(rho(nmo),iwork(nmo),rwork(nmo),aux(nmo,nmo))
      allocate(orboper(nmo,nmo,48))
c
c.....construct idict  (list of centers that will be obtained/rotated)
c
      call rzero(aux,nmo)
c
c.....account for correlated case: lcor=1, ifirst=nmo
c     irrep and chartab contain two sectors. from 1 to nmo 
c     the symm. of natural orbitals. from nmo to nmo+nmo canonical orbs.
c     
*     if (lcor.eq.0) ifirst=0      !not necessary. only rotate canonical.
*     if (lcor.ne.0) ifirst=nmo
c
c.....Construction of a set of independent vectors by modular algebra.
c           We obtain nmo independent vectors
c
      ntrials=0
      r0=0.5d0
 123  n=4
      m=7
      r=r0+0.1234d0*ntrials
      do i=1,nmo
         th=mod(dble(n),pi)
         phi=mod(dble(m),2d0*pi)
         n=n+1
         m=m+1
         vec(1,i)=r*sin(th)*cos(phi)
         vec(2,i)=r*sin(th)*sin(phi)
         vec(3,i)=r*cos(th)
         r=r+0.1d0
      enddo
c
c.....Evaluate orbitals at the vec() coordinates.
c
      do i=1,nmo
        call orbfun(vec(1,i),rho,lcor)
        do j=1,nmo
           forbi(j,i)=rho(j)
        enddo
      enddo
c
c.....invert the matrix of basis vectors.
c
      call mydgeco(forbi,nmo,nmo,iwork,rcond,rwork)
      if (rcond.lt.epssym) then
         write (luw,*) '# problems in transor inversion !!'
         write (luw,*) '# epssym = ', epssym
         write (luw,*) '# condition number:', rcond
c
c........Guess a different set of basis vectors.
c
         ntrials=ntrials+1
         if (ntrials.lt.maxtrials) goto 123
      endif
      call mydgedi(forbi,nmo,nmo,iwork,xdet,rwork,1)
c
c.....Rotation of the vectors by the nop operations of the group
c
      do n=1,nop
        do i=1,nmo
          do j=1,3
            vec1(j)=0d0
            do k=1,3
              vec1(j)=vec1(j)+oper(j,k,n)*vec(k,i)
            enddo
          enddo
          call orbfun(vec1,rho,lcor)
          do j=1,nmo
            forbf(j,i)=rho(j)
          enddo
        enddo
c
c.......After inverting the matrix of basis vectors, multiply it 
c       by the rotated one
c
        call rmul (forbf,forbi,orboper(1,1,n),nmo)
c
c.......obtain the matrix of chi's squares. It determines the type
c       of irreducible representations.
c
        aaux=0d0
        do i=1,nmo
          do j=1,nmo
            aux(i,j)=aux(i,j)+(orboper(i,j,n))**2
            aaux=aaux+aux(i,j)
          enddo
        enddo
      enddo
c
c.....Analyze.
c
      if (abs(aaux-nmo*nop).gt.epsorb) then
         write (luw,*) aaux, nmo*nop
         write (luw,*) '# Orbitals break symmetry.'
         write (luw,*) '# Repeat with the NOSYMMETRY option '
         write (luw,*) '# OR increase the EPSORB parameter'
         stop
      endif
c
      ifirst=0
      do i=1,nmo
         ip=nint(nop/aux(i,i))
         irrep(i+ifirst,0)=ip
         l=0
         do j=1,nmo
            if (aux(i,j).gt.epsorb) then
               l=l+1
               irrep(i+ifirst,l)=j+ifirst
            endif
         enddo
         if (l.ne.ip) then
            write (luw,*) '# Problems with ',i,' orbital irrep.'
            write (luw,*) '# Repeat with the NOSYMMETRY option'
            write (luw,*) '# OR increase the EPSORB parameter'
            stop
         else if (ip.eq.1) then
            do n=1,nop
               chartab(i+ifirst,1,n)=sign(1d0,orboper(i,i,n))
            enddo
         else
            do n=1,nop
              do j=1,ip
                 jj=irrep(i+ifirst,j)-ifirst
                 chartab(i+ifirst,j,n)=orboper(i,jj,n)
              enddo
            enddo
         endif
      enddo
c
c.....write
c
      if (trim(symbol).ne.'C1') then
        if (lcor.eq.0) write (luw,11) 
        if (lcor.ne.0) write (luw,12) 
        do i=1,nmo
           write (luw,13) i,irrep(i+ifirst,0),
     &     (irrep(i+ifirst,j)-ifirst,j=1,irrep(i+ifirst,0))
        enddo
      endif
      if (largwr.and.trim(symbol).ne.'C1') then
        write (luw,14)
        do n=1,nop
           write (luw,15) n
           do j=1,nmo
              ip=irrep(j+ifirst,0)
              write (luw,16)j, (chartab(j+ifirst,k,n),k=1,ip)
           enddo
        enddo
      endif
c
      deallocate (forbi,forbf,vec,rho,iwork,rwork,aux)
      return
c
 11   format (/,1x,'| Analysis of the Symmetry of Standard  Orbitals:')
 12   format (/,1x,'| Analysis of the Symmetry of Canonical Orbitals:')
 13   format (1x,'|   Orbital ',i3,' belongs to irrep of order ',i3,
     &        ' and transforms with orbitals ',6i4)  
 14   format (/1x,'| Transformation matrices for orbitals:')
 15   format (/,1x,'|   Operation:',i3)
 16   format (' |',4x,' Orbital ',i3,':',6f12.6)
      end    
