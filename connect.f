c
c-----------------------------------------------------------------------
c
      subroutine connect (lw,covx,ncent,maxc,largwr,nb,bonded,wh,madis)
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
      parameter (zero = 0d0)
      parameter (one  = 1d0)
c
c.....error codes.
c
      parameter (infinity=-1)
      parameter (a2b=1.889725989D0)
      integer(kind=4), allocatable,dimension (:)   :: catom,cdis
      integer(kind=4), allocatable,dimension (:)   :: iord
      integer(kind=4), allocatable,dimension (:)   :: iclus,iclaux
      real(kind=8),    allocatable,dimension (:)   :: d,c,ee
      real(kind=8),    allocatable,dimension (:,:) :: cnx,v
      real(kind=8),    allocatable,dimension (:)   :: work
      real(kind=8)     covx
      integer(kind=4)  coord(ncent)
      integer(kind=4)  bonded(ncent)
      integer(kind=4)  madis(ncent,ncent)

      character(len=1),allocatable,dimension (:,:) :: labdis
      integer(kind=4)  p(0:ncent)
      integer(kind=4)  wh(ncent,maxc)
      logical connected,largwr
      real(kind=8) covrad(85),xdis(3)
      character(len=1)  digs(-1:18)
      character(len=2)  this,oth(maxc)
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
      n = ncent
      allocate (catom(n)   )
      allocate (cdis(n)    )
      allocate (labdis(n,n))
      allocate (cnx(n,n)   )
      allocate (ee(n)      )
      allocate (d(n)       )
      allocate (v(n,n)     )
      allocate (c(0:n)     )
      allocate (iord(n)    )
      allocate (iclus(n)   )
      allocate (iclaux(n)  )
c
c
c-----Compute connectivity matrix. Two atoms are bonded if the dis-
c     tance between them is smaller the sum of their covalent radius 
c     multiplied by a factor 1.2.
c
      cnx=zero 
      do k=1,n-1
        do m=k+1,n
          xdis(:)=xyz(k,:)-xyz(m,:)
          dis2=xdis(1)**2+xdis(2)**2+xdis(3)**2
          ichk=int(charge(k))
          ichm=int(charge(m))
          if (ichk.gt.85) then
            covk=2D0
          else
            covk=covrad(ichk)
          endif
          if (ichm.gt.85) then
            covm=2D0
          else
            covm=covrad(ichm)
          endif
          rbond=(covk+covm) * COVX * a2b
          if (dis2.lt.rbond*rbond) then 
            cnx(k,m)=one
            cnx(m,k)=one
          endif
        enddo
      enddo
c
c-----Compute the coordinations of all the atoms.
c
      nb=0
      do i=1,n
        coord(i)=0
        do j=1,n
          nb=nb+nint(cnx(i,j))
          coord(i)=coord(i)+nint(cnx(i,j))
        enddo
      enddo
      nb=nb/2
c
c-----Write coordination indices and connectivity matrix.
c
      write (lw,132) covx
      write (lw,1) 'Coordination indices and Connectivity Matrix'
      write (lw,1) '--------------------------------------------'
      do k=1,n
        nwh=0
        do m=1,n
          if (nint(cnx(k,m)).eq.1) then
            nwh=nwh+1
            if (nwh.gt.maxc) then
              stop 'connect.f: Increase the value of MAXCOORD'
            endif
            wh(k,nwh)=m
          endif
        enddo
        this(1:2)=atnam(k)(3:4)
        forall(m=1:nwh) oth(m)(1:2)=atnam(wh(k,m))(3:4)
        if (n.lt.10) then
          write (lw,2201) this,k,coord(k),(oth(m)(1:2),wh(k,m),m=1,nwh)
        elseif (n.lt.100) then
          write (lw,2202) this,k,coord(k),(oth(m)(1:2),wh(k,m),m=1,nwh)
        elseif (n.lt.1000) then
          write (lw,2202) this,k,coord(k),(oth(m)(1:2),wh(k,m),m=1,nwh)
        else
          write (lw,2204) this,k,coord(k),(oth(m)(1:2),wh(k,m),m=1,nwh)
        endif
      enddo
      bonded = coord
      write (lw,10) nb
 1    format (1x,a)
 10   format (I4,' bonds')
 2201 format (' (',a2,1x,I1,') Coor = ',I1,' -->',
     &  30(1x,'(',a2,1x,I1,')'))
 2202 format (' (',a2,1x,I2,') Coor = ',I2,' -->',
     &  30(1x,'(',a2,1x,I2,')'))
 2204 format (' (',a2,1x,I4,') Coor = ',I4,' -->',
     &  30(1x,'(',a2,1x,I4,')'))
 132  format (/
     & ' # Bonding Criterion = ',F12.5,' * Sum of covalent radii')
c
c-----Order the coordinations of all the atoms.
c
      forall (i=1:n) iord(i)=i
      call iqcksort (coord,iord,n,1,n)
      forall (i=1:n) catom(i)=coord(iord(i))
c
c-----Diagonalize the connectivity matrix.
c
      v=cnx

      n3=3*n-1
      allocate (work(n3))
      call dsyev ('V','L',n,v,n,d,work,n3,info)
      deallocate (work)
c
c-----The eigenvalues d(i) are returned by dsyev in ascending order
c
*     call tred2  (v,n,n,d,ee,.true.)
*     call tqli   (d,ee,n,n,v,.true.)
***   call eigsrt (cnx,d,v,n,n)
c
c-----Detemine the characteristic polynomium.
c
      call polich (d,c,n)
      forall (i=0:n) p(i)=nint(c(i))
c
c-----Compute the distance matrix. Algorithm: See Chemical Graph Theory.
c     Chapter 2 by O. E. Polansky, Section 2.9, Eq. (51).
c
      connected = .true.
      do i=1,n-1
        madis(i,i)=0
        do j=i+1,n
          do nu=1,n
            ankl=0d0
            do k=1,n
              ankl=ankl+v(i,k)*v(j,k)*d(k)**nu
            enddo
            if (nint(ankl).gt.0) then
              madis(i,j)=nu
              madis(j,i)=nu
              go to 4
            endif
          enddo
c
c---------This cluster is a disconnected graph.
c
          connected = .false.
          madis(i,j)=infinity
          madis(j,i)=infinity
 4        continue
        enddo
      enddo
      madis(n,n)=0
c
c-----Compute the sum of distance indices for all the atoms and order them.
c
      forall (i=1:n) coord(i)=sum(madis(i,:))
c
c-----Write the sum of indices of distances and the distance matrix.
c
      if (connected) then
        write (lw,1) 'Connected graph'
      else
        write (lw,1) 'Non-connected graph'
      endif
      write (lw,1) "Distance Matrix: '-' means non connected atoms"
      write (lw,1) "Distance Matrix: Only values < 10 are printed"
      if (n.le.100) then
        labdis(1:n,1:n)=digs(0)
        do i=1,n
          do j=1,n
            labdis(i,j)=digs(madis(i,j))
          enddo
          write (lw,'(1000(1x,100a1))') (labdis(i,j),j=1,n)
        enddo
      endif
c
c-----Largest number of bonds between two 'connected' atoms
c
      longest_chain=maxval(madis)
      write (lw,21) longest_chain
 21   format (' #',/,' # Largest distance index is ',I2,/,' #')
c
c-----Order the sums of indices of distances.
c
      forall (i=1:n) iord(i)=i
      call iqcksort (coord,iord,n,1,n)
      forall (i=1:n) cdis(i)=coord(iord(i))
c
c-----Determine and identify the different non-connected clusters.
c
      forall (i=1:n) iclus(i)=0
      nclus=0
      do i=1,n
        if (iclus(i).eq.0) then
          nclus=nclus+1
          iclus(i)=nclus
          do j=i+1,n
            if (madis(i,j).ne.infinity) iclus(j)=nclus
          enddo
        endif
      enddo
      if (nclus.gt.1) then
*       write (lw,11)
        write (lw,7) nclus
        do i=1,nclus
          nclaux=0
          do j=1,n
            if (iclus(j).eq.i) then 
              nclaux=nclaux+1
              iclaux(nclaux)=j
            endif
          enddo
          write (lw,8) i, nclaux
          write (lw,9) (iclaux(j),j=1,nclaux)
        enddo
      endif
c
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
 7    format (1x,'Molecule made of ',I2,' non-connected fragmets')
 8    format (1x,'Fragment ',I4,' contains ', I4,' atoms :')
 9    format (20(1x,I4))
c
      end
c
c-----------------------------------------------------------------------
c
      subroutine polich (x,c,n)
c
c-----Let x(i), i=1,n, the n eigenvalues of a real symmetric matrix A
c     of logical size n x n; the n solutions of the secular equation of A
c           
c     det | x I - A | = SUM (k=0,N) c(k) x^k = 0
c             -   -
c     are x(1), x(2),..., x(n). c(k) are the cofficients of the charac-
c     teristic polynomial and are the objective of this routine.
c
      implicit real(kind=8) (a-h,o-z)
      parameter (zero = 0d0)
      parameter (one  = 1d0)
      dimension x(n), c(0:n)
c
      if (n.gt.0) then
        c(0) = -x(1)
        c(1) = one
        do i=2,n
          c(i)=c(i-1)
          do k=i-1,1,-1
            c(k)=c(k-1)-x(i)*c(k)
          enddo
          c(0) = -x(i)*c(0)
        enddo
      else
         stop ' # _Polich_: Improper dimension'
      endif
      return
      end
