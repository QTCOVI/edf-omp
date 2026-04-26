c
c-----------------------------------------------------------------------
c
      subroutine suppress (aom,iwfn,filedat,wfnfile,line,lw)
c 
c.......................................................................
c
      USE      space_for_wfnbasis
      USE      space_for_wfncoef
      include  'implicit.inc'
      include  'param.inc'
      include  'constants.inc'
      include  'wfn.inc'
      include  'corr.inc'
      include  'stderr.inc'
      include  'mline.inc'
      integer(kind=4), allocatable,dimension (:)      :: isup,isupx,iord
      character*(*) line,filedat,wfnfile
      character*(mline) filenew,word
      real(kind=8)       aom(ncent,nmo,nmo)
      logical       setint,ok,todelete,notinlist,setword,exfil
      integer       lp,nsup,nosup
      character*(4) fourchar,mode
      character*80  wfnttl,label
c
      lp = 1
      ok = setword(word,line,lp)
      if (.not.ok) return
      if (trim(word).eq.'AOM'.or.trim(word).eq.'WFN') then
c
c.......Determine the MOs to be elliminated.
c
        allocate (isup(nmo))
        allocate (isupx(nmo))
        line=line(lp:)
        lp=1
        nsup=0
        do while (setint(n,line,lp))
          if (n.le.nmo.and.n.gt.0) then
            if (nsup.eq.0) then
              nsup = 1
              isup(1) = n
            else
              i=1
              notinlist=.true.
              do while (i.le.nsup.and.notinlist) 
                if (n.eq.isup(i)) notinlist=.false.
                i=i+1
              enddo
              if (notinlist) then
                nsup=nsup+1
                isup(nsup)=n
              endif
            endif
          endif
        enddo
        if (nsup.eq.0) return ! None MO is elliminated
c
c-------Ordering of the MOs to be elliminated
c
        allocate (iord(nsup))
        forall (i=1:nsup) iord(i)=i
        call iqcksort (isup,iord,nsup,1,nsup)
c
c.......Determine the MOs that ARE NOT elliminated.
c
        nosup=0
        do i=1,nmo
          todelete=.false.
          do k=1,nsup
            if (i.eq.isup(iord(k))) todelete=.true.
          enddo
          if (.not.todelete) then
            nosup=nosup+1
            isupx(nosup)=i
          endif
        enddo

        lu19     =  19
        exfil    = .true.
        if (trim(word).eq.'AOM') then
          filenew = trim(filedat)//'-supp'
        else
          filenew = trim(wfnfile)//'-supp'
        endif
        do while (exfil)
          inquire (unit=lu19,opened=exfil)
          if (.not.exfil) then
            open (unit=lu19,file=trim(filenew),iostat=ios)
            if (ios.ne.0) then
              write (0,*) '# suppress.f: Cannot open '//trim(filenew)
              stop
            endif
          else
            lu19=lu19+1
          endif
        enddo
        if (trim(word).eq.'AOM') then
          write (lw,3511) 'AOM',trim(filenew),trim(filedat),
     &                    (isup(iord(m)),m=1,nsup)
          write (lw,3512) 
          do i=1,ncent
            write (lu19,'(I4,a)') i,' <--- AOM in this center'
            write (lu19,80) ((aom(i,isupx(m),isupx(j)),m=1,j),j=1,nosup)
          enddo
        else
          write (lw,3511) 'WFN',trim(filenew),trim(wfnfile),
     &                    (isup(iord(m)),m=1,nsup)
          write (lw,3512) 
          iwfnew=lu19
          read (iwfn,101) wfnttl
          write (iwfnew,'(a)') trim(wfnttl)//' - suppressing MOs'
          read (iwfn,102) mode,ndummy1,ndummy2,ndummy3
          write (iwfnew,1022) nmo-nsup,nprims,ncent
          do i=1,ncent
            read  (iwfn,101)  wfnttl
            write (iwfnew,'(a)') trim(wfnttl)
          enddo
          read (iwfn,104) (ndummy,i=1,nprims)
          nrec = nprims/20
          nres = mod (nprims,20)
          inic = 1
          ifin = 20
          do i=1,nrec
            write (iwfnew,1041) (icen(j),j=inic,ifin)
            inic=inic+20
            ifin=ifin+20
          enddo
          if (nres.gt.0) then
            write (iwfnew,1041) (icen(i),i=20*nrec+1,nprims)
          endif
          nlines=nrec
          if (nres.gt.0) nlines=nlines+1
          read (iwfn,'(a)') word
          inic = 1
          ifin = 20
          do i=1,nrec
            write (iwfnew,1042) (ityp(j),j=inic,ifin)
            inic=inic+20
            ifin=ifin+20
          enddo
          if (nres.gt.0) then
            write (iwfnew,1042) (ityp(i),i=20*nrec+1,nprims)
          endif
          read (iwfn,105) (dummy,i=1,nprimsold)
          nrec = nprims/5
          nres = mod (nprims,5)
          inic = 1
          ifin = 5
          do i=1,nrec
            write (iwfnew,1051) (oexp(j),j=inic,ifin)
            inic=inic+5
            ifin=ifin+5
          enddo
          if (nres.gt.0) then
            write (iwfnew,1051) (oexp(i),i=5*nrec+1,nprims)
          endif
          do i = 1,nosup
            write (iwfnew,1061) i,occ(isupx(i)),eorb(isupx(i))
            write (iwfnew,1071) (coef(isupx(i),j),j=1,nprims)
          enddo
          write (iwfnew,108) 'END DATA'
          write (iwfnew,109) 0d0,0d0
        endif
      else
        return
      endif
      deallocate (isup,isupx,iord)
      return
 3511 format (1x,"#",/,' # ',90('='),/," # Writing ",a," file ",
     & "'",a,"'",/," # obtained from file ","'",a,
     & "' by deleting the MOs",/,1000(' #',30I3,/))
 3512 format (' # ',90('='))
 80   format (6(1x,e16.10))
 1022 format ('GAUSSIAN',10X,I5,15X,I5,15X,I5,' NUCLEI')
 109  format (' TOTAL ENERGY =',F22.12,' THE VIRIAL(-V/T)=',F13.8)
 101  FORMAT (A80)
 102  FORMAT (4X,A4,10X,3(I5,15X))
 103  FORMAT (A8,11X,I3,2X,3F12.8,10X,F5.1)
 104  FORMAT (20X,20I3)
 1041 FORMAT ('CENTRE ASSIGNMENTS  ',20I3)
 1042 FORMAT ('TYPE ASSIGNMENTS    ',20I3)
 105  FORMAT (10X,5E14.7)
 1051 FORMAT ('EXPONENTS ',5(1PE14.7))
 106  FORMAT (35X,F12.8,15X,F12.8)
 1061 FORMAT ('MO',I5,5X,'MO 0.0',8x,'OCC NO = ',F12.7,
     &  '  ORB. ENERGY =',F12.6)
 107  FORMAT (5E16.8)
 1071 FORMAT (5(1PE16.8))
 108  FORMAT (A8)
      end
