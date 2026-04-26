      subroutine variance (ngroup,nprob,lw,navg,pnew,resnc,core)
c
c.......................................................................
c
      implicit real(kind=8) (a-h,o-z)
      include 'param.inc'
      real   (kind=8) :: delta0(ngroup)
      real   (kind=8) :: delta (ngroup)
      real   (kind=8) :: navg  (ngroup)
      real   (kind=8) :: ddisquare(nprob)
      real   (kind=8) :: pnew(nprob)
      integer(kind=4) :: ioprob(nprob)
      integer(kind=4) :: resnc(nprob,ngroup)
      integer(kind=4) :: core(ngroup)
c
c-----Variance of the groups
c
      write (lw,230)
      totvar=0d0
      do i=1,ngroup
        delta0(i) = 0d0
        do k=1,nprob
          delta0(i) = delta0(i) + pnew(k) * (navg(i)-resnc(k,i))**2
        enddo
        totvar=totvar+delta0(i)
        write (lw,23) i,delta0(i)
      enddo
      write (lw,24) totvar
      forall (k=1:nprob) ioprob(k)=k
      do k=1,nprob
        ddi = 0d0
        ddisquare(k) = totvar
        do i=1,ngroup
          ddi = navg(i) - resnc(k,i)
          ddisquare(k) = ddisquare(k) + ddi * ddi
        enddo
      enddo
      call qqsort (ddisquare,ioprob,1,nprob,nstack)
      do k=1,nprob
        io=ioprob(k)
        this = ddisquare(io)
        if (pnew(io).gt.1D-4) then
          write (lw,25,advance='no')
     &       this,(resnc(io,i)+2*core(i),i=1,ngroup)
          write (lw,26) pnew(io)
        endif
      enddo

 230  format (' # Variance of Group i: SUM_{S} p{S} [n_i-<n_i>]^2')
 23   format (' # Variance of Group ',I4,' = ',F16.9)
 24   format (' # Total variance         = ',F16.9)
 25   format (' # Variance = ',F16.9,' for S = (',30I4,')')
 26   format (4x,'with PROB = ',F16.9)
 
      return
      end
