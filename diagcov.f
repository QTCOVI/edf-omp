      subroutine diagcov (ngroup,nprob,lw,navg,xli,d1,pnew,resnc)
c
c.......................................................................
c
      implicit real(kind=8) (a-h,o-z)
      include 'param.inc'

      real   (kind=8) :: covmat(ngroup,ngroup)
      real   (kind=8) :: covvec(ngroup,ngroup)
      real   (kind=8) :: coveig(ngroup)
      real   (kind=8) :: wcov(ngroup*3-1)
      real   (kind=8) :: navg(ngroup)
      real   (kind=8) :: xli(ngroup)
      real   (kind=8) :: d1(ngroup*(ngroup-1)/2)
      real   (kind=8) :: pnew(nprob)
      integer(kind=4) :: resnc(nprob,ngroup)

      forall(i=1:ngroup) covmat(i,i)=navg(i)-xli(i)
      ipair=0
      do i=2,ngroup
        do j=1,i-1
          ipair=ipair+1
          covmat(i,j) = -0.5d0*d1(ipair)
          covmat(j,i) = covmat(i,j)
        enddo
      enddo
      write (lw,23)
      do i=1,ngroup
        write (lw,'(8(1x,F15.8))')(covmat(i,j),j=1,ngroup)
      enddo
!     call jacobi (covmat,ngroup,ngroup,coveig,covvec,nrot)
!     write (lw,24) (coveig(i),i=1,ngroup)
      n3=ngroup+ngroup+ngroup-1
      call dsyev ('V','U',ngroup,covmat,ngroup,coveig,wcov,n3,info)
      call errdsyev (info,'binedf.f')
      write (lw,24) (coveig(i),i=1,ngroup)
      write (lw,25)
      do i=1,ngroup
        write (lw,'(8(1x,F15.8))') (covmat(i,j),j=1,ngroup)
      enddo
      return
 23   format (' # Covariance Matrix')
 24   format (' # Covariance Eigenvalues',/,100(8(1x,F15.8),/))
 25   format (' # Covariance Eigenvectors')
 26   format (8(1x,F15.8))
      end
