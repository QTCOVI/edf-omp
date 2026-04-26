c
c-----------------------------------------------------------------------
c
      subroutine confrag (lw,covx,ncent,nfrag,atoms_in_frag,ifrag,maxc,
     &  largwr,nb,bonded,wh,madis)
c
c-----Given the cartesian coordinates of the n atoms of a molecule, this
c     routine determines the connectivity matrix, diagonalizes it, 
c     determines the characteristic polynomiun, the first n powers of the 
c     connectivity matrix, and the so-called distance matrix.
c
c     Evelio Francisco.
c     Universidad de Oviedo.
c     Nov, 2016.          
c
      USE  space_for_wfnbasis
      implicit real(kind=8) (a-h,o-z)
      parameter (zero = 0.0d+00)
      parameter (one  = 1.0d+00)
c
c.....error codes.
c
      parameter (infinity=-1)
      parameter (ams2bohr=1.889725989D0)
      integer(kind=4) atoms_in_frag(nfrag)
      integer(kind=4) ifrag(ncent,nfrag)
      integer(kind=4), allocatable,dimension (:)   :: catom,cdis
      integer(kind=4), allocatable,dimension (:)   :: iord
      integer(kind=4), allocatable,dimension (:)   :: iclus,iclaux
      real(kind=8),    allocatable,dimension (:)   :: d,c,ee
      real(kind=8),    allocatable,dimension (:,:) :: cnx,v
      real(kind=8),    allocatable,dimension (:)   :: work
      real(kind=8)     covx
      integer(kind=4)  coord(nfrag)
      integer(kind=4)  bonded(nfrag)
      integer(kind=4)  madis(nfrag,nfrag)

      character(len=1),allocatable,dimension (:,:) :: labdis
      integer(kind=4)  p(0:nfrag)
      integer(kind=4)  wh(nfrag,maxc)
      logical connected,largwr
      real(kind=8) covrad(85),xdis(3)
      character(len=1)  digs(-1:18)
      character(len=8)  this,oth(maxc)
      data digs / '-','-','1','2','3','4','5','6','7','8','9',
     &            'a','b','c','d','e','f','g','h','i'/

      data (covrad(i),i=1,85) 
     .   /0.32d0,0.93d0,1.23d0,0.90d0,0.82d0,0.77d0,0.75d0,0.73d0, 
     .    0.72d0,0.71d0,1.54d0,1.36d0,1.18d0,1.11d0,1.06d0,1.02d0, 
     .    0.99d0,0.98d0,2.03d0,1.74d0,1.44d0,1.32d0,1.22d0,1.18d0, 
     .    1.17d0,1.17d0,1.16d0,1.15d0,1.17d0,1.25d0,1.26d0,1.22d0, 
     .    1.20d0,1.16d0,1.14d0,1.12d0,2.16d0,1.91d0,1.62d0,1.45d0, 
     .    1.34d0,1.30d0,1.27d0,1.25d0,1.25d0,1.28d0,1.34d0,1.48d0, 
     .    1.44d0,1.41d0,1.40d0,1.36d0,1.33d0,1.31d0,2.35d0,1.98d0, 
     .    1.69d0,1.65d0,1.65d0,1.64d0,1.63d0,1.62d0,1.85d0,1.61d0, 
     .    1.59d0,1.59d0,1.58d0,1.57d0,1.56d0,1.56d0,1.56d0,1.44d0, 
     .    1.34d0,1.30d0,1.28d0,1.26d0,1.27d0,1.30d0,1.34d0,1.49d0, 
     .    1.48d0,1.47d0,1.46d0,1.46d0,1.45d0/
c
c-----Allocate some arrays
c
      n = nfrag
      allocate (catom(nfrag)   )
      allocate (cdis(nfrag)    )
      allocate (labdis(nfrag,nfrag))
      allocate (cnx(nfrag,nfrag)   )
      allocate (ee(nfrag)      )
      allocate (d(nfrag)       )
      allocate (v(nfrag,nfrag)     )
      allocate (c(0:nfrag)     )
      allocate (iord(nfrag)    )
      allocate (iclus(nfrag)   )
      allocate (iclaux(nfrag)  )
c
c
c-----Compute connectivity matrix. Two fragments are bonded if the 
c     distance between at least one of the pairs of atoms, one of them
c     in the first fragment and the other in the second fragment, is 
c     smaller the sum of their covalent radii multiplied by a factor COVX.
c
      cnx=zero 
      do k=1,nfrag-1
        do m=k+1,nfrag
          do l=1,atoms_in_frag(k)
            do n=1,atoms_in_frag(m)
              kis = ifrag(l,k)
              mis = ifrag(n,m)
              xdis(:) = xyz(k,:)-xyz(m,:)
              ichk = int(charge(kis))
              ichm = int(charge(mis))
              if (ichk.gt.85) then
                covk = 2.0d+00
              else
                covk = covrad(ichk)
              endif
              if (ichm.gt.85) then
                covm = 2.0d+00
              else
                covm = covrad(ichm)
              endif
              rbond = (covk+covm) * covx * ams2bohr
              if (dot_product(xdis,xdis).lt.rbond*rbond) then 
                cnx(k,m) = one
                cnx(m,k) = one
              endif
            enddo
          enddo
        enddo
      enddo
c
c-----Compute the coordinations of all the fragments.
c
      nb = 0
      do i=1,nfrag
        coord(i)=0
        do j=1,nfrag
          nb = nb+nint(cnx(i,j))
          coord(i) = coord(i)+nint(cnx(i,j))
        enddo
      enddo
      nb = nb/2
c
c-----Write coordination indices and connectivity matrix.
c
      write (lw,132) ' # Bonding Criterion = ',covx
      write (lw,1) '# Coordination indices and Connectivity Matrix'
      write (lw,1) '# --------------------------------------------'
      do k=1,nfrag
        nwh=0
        do m=1,nfrag
          if (nint(cnx(k,m)).eq.1) then
            nwh=nwh+1
            if (nwh.gt.maxc) then
              stop 'confrag.f: Increase the value of MAXCOORD'
            endif
            wh(k,nwh)=m
          endif
        enddo
       
*       this(1:2)=atnam(k)(3:4)
*       forall(m=1:nwh) oth(m)(1:2)=atnam(wh(k,m))(3:4)
        this(1:8)='FRAGMENT'
        forall(m=1:nwh) oth(m)(1:8)='FRAGMENT'
        if (nfrag.lt.10) then
          write (lw,2201) this,k,coord(k),(wh(k,m),m=1,nwh)
        elseif (nfrag.lt.100) then
          write (lw,2202) this,k,coord(k),(wh(k,m),m=1,nwh)
        elseif (nfrag.lt.1000) then
          write (lw,2202) this,k,coord(k),(wh(k,m),m=1,nwh)
        else
          write (lw,2204) this,k,coord(k),(wh(k,m),m=1,nwh)
        endif
      enddo
      bonded = coord
      write (lw,10) nb
c
c-----Order the coordinations of all the fragments.
c
      forall (i=1:nfrag) iord(i)=i
      call iqcksort (coord,iord,nfrag,1,nfrag)
      forall (i=1:nfrag) catom(i)=coord(iord(i))
c
c-----Diagonalize the connectivity matrix.
c
      v = cnx

      n3 = 3*nfrag-1
      allocate (work(n3))
      call dsyev ('V','L',nfrag,v,nfrag,d,work,n3,info)
      deallocate (work)
c
c-----Detemine the characteristic polynomium.
c
      call polich (d,c,nfrag)
      forall (i=0:nfrag) p(i) = nint(c(i))
c
c-----Compute the distance matrix. Algorithm: See Chemical Graph Theory.
c     Chapter 2 by O. E. Polansky, Section 2.9, Eq. (51).
c
      connected = .true.
      do i=1,nfrag-1
        madis(i,i) = 0
        cyclej: do j=i+1,nfrag
          do nu=1,nfrag
            ankl=0.0d+00
            do k=1,nfrag
              ankl = ankl+v(i,k)*v(j,k)*d(k)**nu
            enddo
            if (nint(ankl).gt.0) then
              madis(i,j) = nu
              madis(j,i) = nu
              cycle cyclej
            endif
          enddo
c
c---------This cluster is a disconnected graph.
c
          connected = .false.
          madis(i,j) = infinity
          madis(j,i) = infinity
 4        continue
        enddo cyclej
      enddo
      madis(nfrag,nfrag) = 0
c
c-----Compute the sum of distance indices for all the fragments and order them.
c
      forall (i=1:nfrag) coord(i)=sum(madis(i,:))
c
c-----Write the sum of indices of distances and the distance matrix.
c
      if (     connected) write (lw,1)     '# Connected graph'
      if (.not.connected) write (lw,1) '# Non-Connected graph'
      write (lw,1) "# Distance Matrix: '-' ==> non connected fragments"
      write (lw,1) "# Distance Matrix: Only values < 10 are printed"
      if (nfrag.le.100) then
        labdis(1:nfrag,1:nfrag) = digs(0)
        do i=1,nfrag
          do j=1,nfrag
            labdis(i,j) = digs(madis(i,j))
          enddo
          write (lw,'(1000(1x,"# ",100a1))') (labdis(i,j),j=1,nfrag)
        enddo
      endif
c
c-----Largest number of bonds between two 'connected' fragments
c
      longest_chain = maxval(madis)
      write (lw,21) longest_chain
c
c-----Order the sums of indices of distances.
c
      forall (i=1:nfrag) iord(i)=i
      call iqcksort (coord,iord,nfrag,1,nfrag)
      forall (i=1:nfrag) cdis(i)=coord(iord(i))
c
c-----Determine and identify the different non-connected clusters.
c
      forall (i=1:nfrag) iclus(i)=0
      nclus = 0
      do i=1,nfrag
        if (iclus(i).eq.0) then
          nclus = nclus+1
          iclus(i) = nclus
          do j=i+1,nfrag
            if (madis(i,j).ne.infinity) iclus(j)=nclus
          enddo
        endif
      enddo
      if (nclus.gt.1) then
        write (lw,7) nclus
        do i=1,nclus
          nclaux = 0
          do j=1,nfrag
            if (iclus(j).eq.i) then 
              nclaux = nclaux+1
              iclaux(nclaux) = j
            endif
          enddo
          write (lw,8) i, nclaux
          write (lw,9) (iclaux(j),j=1,nclaux)
        enddo
      endif
c
c-----Deallocate some arrays
c
      deallocate (catom )
      deallocate (cdis  )
      deallocate (labdis)
      deallocate (cnx   )
      deallocate (ee    )
      deallocate (d     )
      deallocate (v     )
      deallocate (c     )
      deallocate (iord  )
      deallocate (iclus )
      deallocate (iclaux)
      return
c
 1    format (1x,a)
 7    format (1x,'Molecule made of ',I2,' non-connected fragmets')
 8    format (1x,'Fragment ',I4,' contains ', I4,' fragments :')
 9    format (20(1x,I4))
 10   format (I4,' bonds')
 21   format (' #',/,' # Largest distance index is ',I2,/,' #')
 132  format (a,F12.5,' * Sum of covalent radii')
 2201 format (' # (',a8,1x,I3,'): Coord = ',I3,' -->',30(1x,'(',I3,')'))
 2202 format (' # (',a8,1x,I3,'): Coord  =',I3,' -->',30(1x,'(',I3,')'))
 2204 format (' # (',a8,1x,I3,'): Coord  =',I3,' -->',30(1x,'(',I3,')'))
      end
