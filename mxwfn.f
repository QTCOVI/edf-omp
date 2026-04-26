c
c-----------------------------------------------------------------------
c
*     subroutine mxwfn (epsdet,nproba,nprobb,ngroup,ndets,nmo,
*    &  ncore,nelact,nel,moval,ival,naval,nbval,ifilc,lw,wfnfile,sg)
      subroutine mxwfn (nvar,xvar,psi,gpsi) 
c
      USE        space_for_cidet
      USE        space_for_conf

      include   'implicit.inc'
      include   'param.inc'
      include   'constants.inc'
      include   'opt.inc'
      include   'corr.inc'
      include   'wfn.inc'
c
      real(kind=8)   xvar(nvar)
      real(kind=8)   gpsi(nvar)
      integer, allocatable,dimension (:)     :: iorda,iordb
      integer, allocatable,dimension (:)     :: ordia,ordib,ordja,ordjb
      integer, allocatable,dimension (:)     :: indxa,indxb
      integer, allocatable,dimension (:)     :: nordia,nordib
      real(kind=8),  allocatable,dimension (:,:)   :: dia,dib,diad,dibd
      real(kind=8),  allocatable,dimension (:)     :: detdia,detdib
      real(kind=8),  allocatable,dimension (:,:)   :: detdiad,detdibd
      real(kind=8),  allocatable,dimension (:)     :: wa,wb
      real(kind=8)   xmo(nmo,nel),xmo1(nmo,3,nel),p(3)
      real(kind=8)   deter(2)
      logical  canmos
c
c-----Update current electron coordinates.
c
      k=0
      do i=1,nel
        do j=1,3
          if (iopt(i,j).eq.1) then
            k=k+1
            xyzel(i,j)=xvar(k)
          endif
        enddo
      enddo
c
c-----Obtain MOs and their gradients at each xyzel() point
c
      canmos=.true.
      do i=1,nel
        p(1:3)=xyzel(i,1:3)
        call molorb1 (p,xmo(:,i),xmo1(:,:,i),canmos)
      enddo 
c
c-----run over all pairs of different alpha configs
c
      allocate (ordia(nmo),ordja(nmo))
      allocate (ordib(nmo),ordjb(nmo))
      allocate (nordia(nmo),nordib(nmo))
      mal=ncore+nelal
      mbe=ncore+nelbe
      do m=1,ncore
        ordia(m)=m
        ordja(m)=m
        ordib(m)=m
        ordjb(m)=m
      enddo
c
c-----Obtain all the D_i\alpha determinants and their gradients
c
      allocate (indxa(mal))
      allocate (wa(mal))
      allocate (detdia(npa))
      allocate (detdiad(npa,3*mal))
      allocate (dia (mal,mal))
      allocate (diad(mal,mal))
      do k=1,npa
        ordia(ncore+1:mal)=nconfa(k,1:nelal)+ncore
        do j=1,mal
          jc=ordia(j)
          do i=1,mal
            dia(i,j)=xmo(jc,i)
          enddo
        enddo
c
c-------Gradients
c
        kis=0
        do i=1,mal
          do j=1,3
            if (iopt(i,j).eq.1) then
              kis=kis+1
              diad=dia
              do kk=1,mal
                jc=ordia(kk)
                diad(i,kk)=xmo1(jc,i,j)
              enddo
              call dgeco (diad,mal,mal,indxa,rcond,wa)
              if (rcond.lt.1d-20) then
                stop 'mxwfn.f: !! The DIAD() matrix is singular !!'
              endif
              job=10
              deter=zero
              call dgedi (diad,mal,mal,indxa,deter,wa,job)
              detdiad(k,kis)=deter(1)*tenp**deter(2) 
            endif
          enddo
        enddo
        call dgeco (dia,mal,mal,indxa,rcond,wa)
        if (rcond.lt.1d-20) then
          stop 'mxwfn.f: !! The DIA() matrix is singular !!'
        endif
        job=10
        deter=zero
        call dgedi (dia,mal,mal,indxa,deter,wa,job)
        detdia(k)=deter(1)*tenp**deter(2)  ! D_i\alpha
      enddo
c
c-----Obtain all the D_i\beta determinants and their gradients
c
      allocate (indxb(mbe))
      allocate (wb(mbe))
      allocate (detdib(npb))
      allocate (dib (mbe,mbe))
      allocate (dibd(mbe,mbe))
      do k=1,npb
        ordib(ncore+1:mbe)=nconfb(k,1:nelbe)+ncore
        do j=1,mbe
          jc=ordib(j)
          do i=1,mbe
            dib(i,j)=xmo(jc,i+mal)
          enddo
        enddo
c
c-------Gradients
c
        kis=0
        do i=1,mbe
          do j=1,3
            if (iopt(i+mal,j).eq.1) then
              kis=kis+1
              dibd=dib
              do kk=1,mbe
                jc=ordib(kk)
                dibd(i,kk)=xmo1(jc,i+mal,j)
              enddo
              call dgeco (dibd,mbe,mbe,indxb,rcond,wb)
              if (rcond.lt.1d-20) then
                stop 'mxwfn.f: !! The DIBD() matrix is singular !!'
              endif
              job=10
              deter=zero
              call dgedi (dibd,mbe,mbe,indxb,deter,wb,job)
              detdibd(k,kis)=deter(1)*tenp**deter(2) 
            endif
          enddo
        enddo
        call dgeco (dib,mbe,mbe,indxb,rcond,wb)
        if (rcond.lt.1d-20) then
          stop 'mxwfn.f: !! The DIB() matrix is singular !!'
        endif
        job=10
        deter=zero
        call dgedi (dib,mbe,mbe,indxb,deter,wb,job)
        detdib(k)=deter(1)*tenp**deter(2)
      enddo

      psi=0d0
      do i=1,ndets
         cd=cdet(i)*nsiga(i)*nsigb(i)
         if (abs(cd).gt.abs(epsdet)) then
           psi=psi+cd*detdia(kwa(i))*detdib(kwb(i))
         endif
      enddo
      psi=-psi*psi/xnorm
      return
      end
