c
c-----------------------------------------------------------------------
c
      subroutine locrohf (aoma,aomb,coef,udat,ialpha,ibeta,
     &    ovc,ncent,nmo,nprims,na,nb,lw,aomfil,wfnfil,okw,larg)
c
c.......................................................................
c
c-----This routine localizes separately the ALPHA and BETA MOs of a ROHF
c     wavefunction. It takes the AOM between the ALPHA MOs, calls the 
c     localization routine (ruedmis), and finds a transformed set of lo-
c     calized ALPHA MOs (LAMOs). The same task is performed with the 
c     BETA MOs, obtaining a transformed set of localized BETA MOs (LBMOs). 
c     After these two jobs are done, LAMOs and LBMOs localized more than 
c     an 'ovc' percentaje in any of the atoms of the molecule is discarded. 
c     The remaining MOs that are not so-localized are written in the file 
c     wfnf//'loc'. The AOM integrals between these not so-localized MOs 
c     are also written in another file called wfnf//'loc.aom'
c
c-----------------------------------------------------------------------
c
c
      implicit none
      integer(kind=4) :: udat,na,nb,nprims,ncent,lw,nmo
      real   (kind=8) :: aoma(ncent,na,na)
      real   (kind=8) :: aomb(ncent,nb,nb)
      real   (kind=8) :: coef(nmo+nmo,nprims)
      integer(kind=4) :: ialpha(nmo)
      integer(kind=4) :: ibeta (nmo)
      real   (kind=8) :: ovc
      logical         :: okw,larg
      character*(*)   :: aomfil,wfnfil
c
c-----Local variables
c
      real   (kind=8),allocatable,dimension (:,:,:) :: aomxa,aomxb
      real   (kind=8),allocatable,dimension (:,:,:) :: aomt
      real   (kind=8),allocatable,dimension (:,:)   :: vrot
      real   (kind=8),allocatable,dimension (:,:)   :: coefa,coefb,cc
      real   (kind=8),allocatable,dimension (:)     :: oc
      real   (kind=8),allocatable,dimension (:)     :: ocup
      integer(kind=4),allocatable,dimension (:,:)   :: ncore
      integer(kind=4),allocatable,dimension (:)     :: ivala
      integer(kind=4),allocatable,dimension (:)     :: ivalb
      integer(kind=4)  :: i,j,k,ktarget,nt,ios,nva,nvb,inpr,luw,ier
      real   (kind=8)  :: temp
      logical          :: exfil,iscore
      character*(1000) :: wfnloc,fileloc

      allocate (aomxa(ncent,na,na))
      allocate (aomxb(ncent,nb,nb))
      allocate (ncore(max(na,nb),2))
      allocate (ivala(na))
      allocate (ivalb(nb))
      allocate (vrot(na,na))
      allocate (oc(na))      
      forall (i=1:na) oc(i)=1d0

      write (6,82) 'ALPHA',na,'BETA ',nb
      vrot=0d0
      aomxa=aoma
      write (lw,83) 'ALPHA',ovc
      call ruedmis (aomxa,vrot,oc,na,ncent,lw,larg,okw)
c
c.....Expand the LAMOs in terms of primitive Gaussians.
c        
      allocate (coefa(na,nprims))
      do i=1,na
        do j=1,nprims
          temp = 0d0
          do k=1,na
            temp = temp + vrot(i,k) * coef(ialpha(k),j)
          enddo
          coefa(i,j) = temp
        enddo
      enddo
      write (lw,65) 'ALPHA'
      nva = 0
      ivala = 0
      ncore = 0
      do i=1,na
        iscore=.false.
        do k=1,ncent
          iscore=iscore.or.(abs(aomxa(k,i,i)).gt.ovc)
          if (iscore) then
            ktarget=k
            write (lw,*) '# ALPHA MSO',i,' is CORE of atom',k
            exit
          endif
        enddo
        if (iscore) then
          ncore(ktarget,1)=ncore(ktarget,1)+1
        else
          nva=nva+1
          ivala(nva)=i
        endif
      enddo 
      deallocate (vrot,oc)      
      do k=1,ncent
        write (lw,81) k,ncore(k,1),'ALPHA'
      enddo
      write (lw,11) nva,(ivala(k),k=1,nva)

      allocate (vrot(nb,nb))
      allocate (oc(nb))      
      forall (i=1:nb) oc(i)=1d0
      vrot=0d0
      write (lw,83) 'BETA ',ovc
      aomxb=aomb
      call ruedmis (aomxb,vrot,oc,nb,ncent,lw,larg,okw)
c
c.....Expand the LBMOs in terms of primitive Gaussians.
c        
      allocate (coefb(nb,nprims))
      do i=1,nb
        do j=1,nprims
          temp = 0d0
          do k=1,nb
            temp = temp + vrot(i,k) * coef(ibeta(k),j)
          enddo
          coefb(i,j) = temp
        enddo
      enddo
      write (lw,65) 'BETA '
      nvb = 0
      ivalb = 0
      do i=1,nb
        iscore=.false.
        do k=1,ncent
          iscore=iscore.or.(aomxb(k,i,i).gt.ovc)
          if (iscore) then
            write (lw,*) '# BETA  MSO',i,' is CORE of atom',k
            ktarget=k
            exit
          endif
        enddo
        if (iscore) then
          ncore(ktarget,2)=ncore(ktarget,2)+1
        else
          nvb=nvb+1
          ivalb(nvb)=i
        endif
      enddo 
      deallocate (vrot,oc)      
      do k=1,ncent
        write (lw,81) k,ncore(k,2),'BETA '
      enddo
      write (lw,12) nvb,(ivalb(k),k=1,nvb)

      write (6,84) 'ALPHA',nva,'BETA ',nvb
c
c.....Write a WFN file with the LAMOs and LBMOs
c
      nt = nva+nvb
      allocate (ocup(nt),stat=ier)
      if (ier.ne.0) stop '# locrohf.f: Cannot allocate ocup()'
      forall (i=1:nva) ocup(i)     = +1d0
      forall (i=1:nvb) ocup(i+nva) = -1d0
      allocate (cc(nt,nprims))
      forall (i=1:nva) cc(i,    :)=coefa(ivala(i),:)
      forall (i=1:nvb) cc(i+nva,:)=coefb(ivalb(i),:)

      inpr = nprims
      wfnloc = trim(wfnfil)//"loc"
      luw  =  19
      exfil = .true.
      do while (exfil)
        inquire (unit=luw,opened=exfil)
        if (.not.exfil) then
          open (unit=luw,file=trim(wfnloc),iostat=ios)
          if (ios.ne.0) then
            write (0,*) ' # locrohf.f: Error openning '//trim(wfnloc)
            stop
          endif
        else
          luw=luw+1
        endif
      enddo
      write (lw,351) trim(wfnloc)
      call cdafhmos (udat,luw,0,ocup,cc,nt,inpr)
      deallocate (ocup,stat=ier)
      if (ier.ne.0) stop '# locrohf.f: Cannot deallocate ocup()'
      deallocate (cc,stat=ier)
      close (luw)

      allocate (aomt(ncent,nt,nt))
      aomt=0d0
      aomt(:,1:nva,1:nva)=aomxa(:,ivala(1:nva),ivala(1:nva))
      aomt(:,nva+1:nt,nva+1:nt)=aomxb(:,ivalb(1:nvb),ivalb(1:nvb))
c
c-----Fill in AOM file with localized MOs
c
      luw  =  20
      exfil = .true.
      fileloc = trim(wfnfil)//'loc.aom'
      do while (exfil)
         inquire (unit=luw,opened=exfil)
         if (.not.exfil) then
           open (unit=luw,file=trim(fileloc),iostat=ios)
           if (ios.ne.0) then
             write (0,*) ' # locrohf.f: Error openning '//trim(fileloc)
             stop
           endif
         else
           luw=luw+1
         endif
      enddo
      write (lw,3510) trim(fileloc)
      do i=1,ncent
        write (luw,'(I4,a)') i,' <--- AOMLOC within this center'
        write (luw,80) ((aomt(i,k,j),k=1,j),j=1,nt)
      enddo
      close (luw)

 80   format (6(1x,e16.10))
 3510 format (" #",/,1x,"# Writing file ","'",a,"' with AOMs over LMOs")
 65   format (//,1x,'# Automatic analysis of ',a,' VALENCE MOs')
 351  format (1x,"#",/,1x,"# Writing file ","'",a,"' with LMOs")
 81   format (1x,'# Atom ',I4,' has ',I4,1x,a,1x,'CORE spin-orbitals')
 84   format (1x,'# Number of valence ',a,' spin-orbitals is ',I4,/,
     &        1x,'# Number of valence ',a,' spin-orbitals is ',I4)
 83   format (//1x,'# Calling the localizaiton routine with ',
     & a,' SMOs',/,1x,'# The value of critical overlap is ',F15.8)
 82   format (1x,'# Number of ',a,' spin-orbitals is ',I4,/,
     &        1x,'# Number of ',a,' spin-orbitals is ',I4)
 11   format (1x,'# ',I4,' ALPHA MSOs',/,1000(1x,'#',10I4,/))
 12   format (1x,'# ',I4,' BETA  MSOs',/,1000(1x,'#',10I4,/))

      return
      end
