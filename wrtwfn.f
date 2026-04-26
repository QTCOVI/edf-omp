c
c-----------------------------------------------------------------------
c
      SUBROUTINE wrtwfn (iwfn,iwfnew,stderr,ifiln,cc,inmo,inpr)
c
c.....Write WFN file with the transformed orbitals.
c
      USE            space_for_wfnbasis
      USE            space_for_cidet
      USE            space_for_wfxedf
      include       'implicit.inc'
      include       'param.inc'
      include       'wfn.inc'
      include       'corr.inc'
      include       'mline.inc'
      real(kind=8)   cc(inmo,inpr)
      integer        leng,stderr
*     parameter      (mline=200)
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
      rewind (iwfn)
      read (iwfn,101) wfnttl
      write (iwfnew,'(a)') trim(wfnttl)//trim(rot0)
      read (iwfn,102) mode,nmoold,nprimsold,ncentold
      write (iwfnew,1022) inmo,nprims,ncent
C
      do 100 i = 1,ncent
        read  (iwfn,101)   wfnttl
        write (iwfnew,'(a)') trim(wfnttl)
100   continue
      read (iwfn,104) (ndummy,i=1,nprimsold)
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
c
      read (iwfn,104) (ndummy,i=1,nprimsold)
      inic = 1
      ifin = 20
      do i=1,nrec
        write (iwfnew,1042) (ityp(j),j=inic,ifin)
        inic=inic+20
        ifin=ifin+20
      enddo
      if (nres.gt.0) write (iwfnew,1042) (ityp(i),i=20*nrec+1,nprims)
c
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
      if (nres.gt.0) write (iwfnew,1051) (oexp(i),i=5*nrec+1,nprims)
c
c-----Here, we read the Electron Density Function (EDF) associated to the
c     core electron with the same format used in the AIMALL code
c
      if (thereisEDF) then
        do i=1,nrrec
          read  (iwfn,'(a)') line
        enddo
        write (iwfnew,131)
        write (iwfnew,'(I3)') nedf
        write (iwfnew,132)
        write (iwfnew,'(5I20)')(icenedf(i),i=1,nedf)
        write (iwfnew,133)
        write (iwfnew,'(5I20)')(itypedf(i),i=1,nedf)
        write (iwfnew,134)
        write (iwfnew,'(5E20.12)')(expedf(i),i=1,nedf)
        write (iwfnew,135)
        write (iwfnew,'(5E20.12)')(cedf(i),i=1,nedf)
        write (iwfnew,136)
      endif
c
      do i = 1,nmo
        read (iwfn,106) occdummy,eorbdummy 
        read (iwfn,107) (coefdummy,j=1,nprimsold)
      enddo
      do i = 1,inmo
        if (rhf) then
          write (iwfnew,1061) i,occ(i),zero
        else
          write (iwfnew,1061) i,zero,zero
        endif
        write (iwfnew,1071) (cc(i,j),j=1,nprims)
      enddo
      read (iwfn,108) check
      write (iwfnew,108) 'END DATA'
      if (.not.icorr.and.ndets.eq.1) then
        read (iwfn,'(A80)') line
        write (iwfnew,'(a)') trim(line)
        close (iwfnew)
        rewind (iwfn)
        return
      endif
c
      read (iwfn,'(A80)') line
      label(1:5)=line(1:5)
      write (iwfnew,'(a)') 'ALDET'//line(6:leng(line))
c
      if (label(1:5).eq.'MCSCF' .or.
     &    label(1:5).eq.'ALDET' .or.
     &    label(1:5).eq.'GENCI') then
c
c.......correlated case. Equal number of NMO
c
        read (iwfn,'(a)',end=9999) label
        write (iwfnew,'(a)') rot1(1:leng(rot1))
        do i=nmo+1,nmo+nmo
           read (iwfn,*) label
           write (iwfnew,'(a,I3)') rot2(1:leng(rot2)),i-nmo
           read (iwfn,107) (coefdummy,j=1,nprimsold)
           write (iwfnew,1071) (cc(i,j),j=1,nprims)
        enddo
        read (iwfn,108) check
        write (iwfnew,108) 'END DATA'
        if (check .ne. 'END DATA') stop ' RDWFN : END CARD NOT FOUND '
        read (iwfn,*) label
        write (iwfnew,'(a)') 'NELACTIVE, NDETS, NORBCORE, NORBACTIVE'
c
c.......NELACTIVE, NDETS, NORBCORE, NORBACTIVE
c
        read (iwfn,*) nelact, ndets, ncore, nact
        read (iwfn,*) label
        write (iwfnew,*) nelact, ndets, ncore, nact
c
        write (iwfnew,'(a)')'COEFFICIENT/ OCCUPIED ACTIVE SPIN ORBITALS'
        do i=1,ndets
           read (ifiln,rec=i) (cidet(m),m=0,nelact)
           write (iwfnew,110) cidet(0),(int(cidet(j)),j=1,nelact)
        enddo
      else
        write (iwfnew,'(a)') rot1(1:leng(rot1))
        do i=nmo+1,nmo+nmo
           write (iwfnew,'(a,I3)') rot2(1:leng(rot2)),i-nmo
           do j=1,nprims
             if (abs(cc(i,j)).lt.1.0D-15) cc(i,j)=0D0
           enddo
           write (iwfnew,1071) (cc(i,j),j=1,nprims)
        enddo
        write (iwfnew,108) 'END DATA'
        write (iwfnew,'(a)') 'NELACTIVE, NDETS, NORBCORE, NORBACTIVE'
        write (iwfnew,*) nelact, ndets, ncore, nact
        write (iwfnew,'(a)')'COEFFICIENT/ OCCUPIED ACTIVE SPIN ORBITALS'
        do i=1,ndets
           read (ifiln,rec=i) (cidet(m),m=0,nelact)
           write (iwfnew,110) cidet(0),(int(cidet(j)),j=1,nelact)
        enddo
      endif
 9999 continue
      close (iwfnew)
      rewind (iwfn)
      return
c
 1022 FORMAT ('GAUSSIAN',10X,I5,15X,I5,15X,I5,' NUCLEI')
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
 110  format (E22.16,1x,1000I4)
 2000 format (' rdwfn: !!! WFN file is not of RHF, MCSCF, ALDET 
     &  or GENCI type type   !!!')
 131  format ('<Additional Electron Density Function (EDF)>',/,
     &        '<Number of EDF Primitives>')
 132  format ('</Number of EDF Primitives>',/,
     &        '<EDF Primitive Centers>')
 133  format ('</EDF Primitive Centers>',/,'<EDF Primitive Types>')
 134  format ('</EDF Primitive Types>',/,'<EDF Primitive Exponents>')
 135  format ('</EDF Primitive Exponents>',/,
     &        '<EDF Primitive Coefficients>')
 136  format ('</EDF Primitive Coefficients>',/,
     &        '</Additional Electron Density Function (EDF)>')
      END
