c
c-----------------------------------------------------------------------
c
      SUBROUTINE cdafhmos (iwfn,iwfnew,stderr,ocup,cc,nloc,inpr)
c
c.....Write WFN file with the transformed orbitals.
c
      USE            space_for_wfnbasis
      USE            space_for_cidet
      include       'implicit.inc'
      include       'param.inc'
      include       'wfn.inc'
      include       'corr.inc'
      include       'mline.inc'
      real(kind=8)   ocup(nloc)
      real(kind=8)   cc(nloc,inpr)
      integer        leng,stderr
      character*(mline)  cicoef
      character*80   line,rot0,rot1,rot2
c
c.....local variables
c
      character*80 wfnttl,jobttl
      character*4 mode
      character*8 check
      character*80 label
      data zero/0d0/
c
      rot0 = ' - Localized MO basis'
      rot1 = 'CANONICAL ORBITALS'
      rot2 = 'CANMO'
c
      if (nprims.ne.inpr) then
        stop '# cdafhmos.f: Incompatible nprims and inpr values'
      endif
      rewind (iwfn)
      read (iwfn,101) wfnttl
      write (iwfnew,'(a)') trim(wfnttl)//trim(rot0)
      read (iwfn,102) mode,nmoold,nprimsold,ncentold
 1022 format ('GAUSSIAN',10x,i5,15x,i5,15x,i5,' NUCLEI')
      write (iwfnew,1022) nloc,inpr,ncent
C
      do i=1,ncent
        read (iwfn,101)   wfnttl
        write (iwfnew,'(a)') trim(wfnttl)
      enddo
      read (iwfn,104) (ndummy,i=1,nprimsold)
      nrec = inpr/20
      nres = mod (inpr,20)
      inic = 1
      ifin = 20
      do i=1,nrec
        write (iwfnew,1041) (icen(j),j=inic,ifin)
        inic=inic+20
        ifin=ifin+20
      enddo
      if (nres.gt.0) then
        write (iwfnew,1041) (icen(i),i=20*nrec+1,inpr)
      endif
c
      read (iwfn,104) (ndummy,i=1,nprimsold)
      inic = 1
      ifin = 20
      do i=1,nrec
        write (iwfnew,1042) (ityp(j),j=inic,ifin)
        inic=inic+20
        ifin=ifin+20
      enddo
      if (nres.gt.0) write (iwfnew,1042) (ityp(i),i=20*nrec+1,inpr)
c
      read (iwfn,105) (dummy,i=1,nprimsold)
      nrec = inpr/5
      nres = mod (inpr,5)
      inic = 1
      ifin = 5
      do i=1,nrec
        write (iwfnew,1051) (oexp(j),j=inic,ifin)
        inic=inic+5
        ifin=ifin+5
      enddo
c
c-----Here, we read the Electron Density Function (EDF) associated to the
c     core electron with the same format used in the AIMALL code
c
      call readedf (iwfn,thereisEDF)
c
      if (nres.gt.0) write (iwfnew,1051) (oexp(i),i=5*nrec+1,inpr)
c
      do i=1,nmoold
        read (iwfn,106) occdummy,eorbdummy 
        read (iwfn,107) (coefdummy,j=1,nprimsold)
      enddo
      do i=1,nloc
        write (iwfnew,1061) i,ocup(i),zero
        write (iwfnew,1071) (cc(i,j),j=1,inpr)
      enddo

      read (iwfn,108) check
      write (iwfnew,108) 'END DATA'
      read (iwfn,'(a80)') line
      label(1:5)=line(1:5)
      write (iwfnew,'(a)') 'RHF      ENERGY =        0.0000000000'
      close (iwfnew)
      rewind (iwfn)
      return
c
 101  FORMAT (A80)
 102  FORMAT (4X,A4,10X,3(I5,15X))
 103  FORMAT (A8,11X,I3,2X,3F12.8,10X,F5.1)
 104  FORMAT (20X,20I3)
 1041 FORMAT ('CENTRE ASSIGNMENTS  ',20I3)
 1042 FORMAT ('TYPE ASSIGNMENTS    ',20I3)
 105  FORMAT (10X,5E14.7)
 1051 FORMAT ('EXPONENTS ',5(1PE14.7))
 106  FORMAT (35X,F12.8,15X,F12.8)
 1061 FORMAT ('Localized MO',I3,11X,'OCC NO = ',F12.8,
     &        '  ORB. ENERGY =',F12.8)
 107  FORMAT (5E16.8)
 1071 FORMAT (5(1PE16.8))
 108  FORMAT (A8)
 109  FORMAT (a17,F20.12,18X,F13.8)
 1091 FORMAT (a17,F20.12,18X,F13.8)
 110  format (E22.16,1x,100I4)
 2000 format (' rdwfn: !!! WFN file is not of RHF, MCSCF, ALDET 
     &  or GENCI type type   !!!')
      END
