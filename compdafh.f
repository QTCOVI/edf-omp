c
c------read the rho^xc matrix as computed by dorho and obtains
c      the DAFH in each atom of the molecule.
c
      subroutine compdafh (aom,dafh,nmo,ncent,lerr,wfnfile)
c
      implicit none
      include   'mline.inc'
      real(kind=8)  aom(ncent,nmo,nmo)
      real(kind=8) dafh(ncent,nmo,nmo)
      real(kind=8) val1,val2,val3,val4,factor
      integer(kind=4) nmo,ncent,lu2,ios,ier,nmopair,m,n,i,j,ij,k,l,kl
      integer(kind=4) lerr
*     parameter (mline   = 200)  
      character*(mline) word
      character(len=*) wfnfile
c
      real(kind=8), allocatable,dimension (:,:)  :: c2exc
c
      lu2=55
      word=trim(wfnfile)//".2rdm"
      open (lu2,file=trim(word),form='unformatted',
     &  iostat=ios,status='old')
      if (ios.ne.0) write (lerr,100) trim(word)
      if (ios.ne.0) stop
      read (lu2) nmopair
      if (nmopair.ne.nmo*(nmo+1)/2) then
        write (lerr,1123) trim(word),nmopair,nmo*(nmo+1)/2
        stop
      endif
      if (.not.allocated(c2exc)) then
        allocate (c2exc(nmopair,nmopair),stat=ier)
        if (ier.ne.0) stop '# compdafh.f: Cannot allocate c2exc()'
      endif
      c2exc = 0d0
      do
        read (lu2,iostat=ios) m,n,val1,val2,val3,val4
        if (ios.ne.0) exit
CCCCC   c2e  (m,n) = val1
CCCCC   c2e  (n,m) = val1
        c2exc(m,n) = val2
        c2exc(n,m) = val2
CCCCC   c2ex (m,n) = val3
CCCCC   c2ex (n,m) = val3
      enddo
      dafh=0d0
      ij=0
      do i=1,nmo
        do j=1,i
          ij=ij+1
          kl=0
          do k=1,nmo
            do l=1,k
              kl=kl+1
              dafh(:,i,j)=dafh(:,i,j)+c2exc(ij,kl)*aom(:,k,l)
            enddo
          enddo
          if (i.eq.j) factor = 1.0d0
          if (i.ne.j) factor = 0.5d0
          dafh(:,i,j) = - dafh(:,i,j) * 0.5d0 * factor
          dafh(:,j,i) =   dafh(:,i,j)
        enddo
      enddo
      return
 100  format (' # compdafh.f: File ',a,' NOT FOUND')
 1123 format (' # NMOPAIR in file ',a,'(',I4,') # ',I4,'=Current value')
      end
