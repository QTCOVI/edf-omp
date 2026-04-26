      subroutine aomnorma (line,aom,tol,lw,lr,lu18,aomnorm,wfnf)
      USE space_for_wfnbasis
      USE space_for_aomspin
      include 'implicit.inc'
      include 'wfn.inc'
      include 'param.inc'
      include 'corr.inc'
      real    (kind=8)   aom(ncent,nmo,nmo)
      real    (kind=8)   tol,sumaom,valueaom,diffaom
      integer (kind=4)   lw,lr,i,m,j,iren,lp
      character*(*)      line,wfnf
      character (len=80) word
      logical            setint,ok,aomnorm
c
c-----In UHF WFNs the AOM between alpha and beta MOs is separately
c     renormalized.
c
      call aomspin (ncent,nmo,aom,nalpha,nbeta)

      lp=1
      ok=setint(iren,line,lp)
      if (ok) then
        if (iren.eq.0) iren=1
        iren=min(iabs(iren),ncent)
      else
        iren=ncent
      endif
      if (uhf) then
        call uhfnorm (iren,ncent,nmo,nalpha,tol,lw,lr,ialpha,aom)
        call uhfnorm (iren,ncent,nmo,nbeta ,tol,lw,lr,ibeta ,aom)
        do i=1,ncent
          aomalpha(i,:,:) = aom(i,ialpha(1:nalpha),ialpha(1:nalpha))
          aombeta (i,:,:) = aom(i,ibeta (1:nbeta ),ibeta (1:nbeta ))
        enddo
      else
c
c.......Test consistency of AOM elements
c
        do m=1,nmo
          do j=1,m
            sumaom=0d0
            do ii=1,ncent
              sumaom=sumaom+aom(ii,j,m)
            enddo
            valueaom=1d0
            if (j.ne.m) valueaom=0d0
            diffaom=abs(sumaom-valueaom)
            if (diffaom.gt.tol) write (lr,544) m,j,sumaom
          enddo
        enddo
c     
c.......Modify the AOM elements of center IREN such that 
c       SUM (i=1,ncent) aom(i,m,j) will be exactly 0 or 1.
c
        write (lw,335) iren
        do j=1,nmo
          do m=1,nmo
            aom2=aom(iren,m,j)
            aom1=0d0
            do i=1,ncent
              if (i.ne.iren) aom1=aom1+aom(i,m,j)
            enddo
            aomtot=aom1+aom2
            if (m.eq.j) then
              if (aomtot.lt.0d0) then
                aom(iren,m,j)=aom2-(1d0-abs(aomtot))
              endif
              if (aomtot.ge.0d0) then
                aom(iren,m,j)=aom2+(1d0-abs(aomtot))
              endif
            else
              aom(iren,m,j)=-aom1
            endif
          enddo
        enddo
      endif
      aomnorm   = .true.
      return

 544  format (' # !! AOM is not Ok !! : SUM_n AOM(n,',2I3,') = ',E15.8)
 335  format (1x,'# AOM IS RENORMALIZED, THE REFERENCE ATOM IS ',I5)
 80   format (6(1x,e16.10))
      end

