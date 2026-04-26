c 
c.......................................................................
c
      subroutine coremos (sg,ng,nmo,ifilc,ndets,ncore,nelcorea,
     & nelcoreb,nact,
     & nelact,nel,nalpha,nbeta,ovc,mocore,core,corea,coreb,icore,
     & icorea,icoreb,moval,ival,lw,ok,largwr,icorr,line,wfnfile)
      USE        space_for_cidet
      include     'implicit.inc'
      include     'mline.inc'
      include     'lengrec.inc'
      real(kind=8) sg(ng,nmo,nmo)
      integer(kind=4) icore(nmo),ival(nmo)
      integer(kind=4) icorea(nmo)
      integer(kind=4) icoreb(nmo)
      integer(kind=4) core(ng),corea(ng),coreb(ng)
      integer(kind=4) icogrp(ng,nmo)
      integer(kind=4) n2(nmo)
      logical iscore,ok,largwr,icorr
      parameter (lline  = 2000)  
      integer(kind=4) totcore
      character(len=4) fourchar
      character*(lline) line
      character*(mline) cicoef
      character*(*) wfnfile
c
c-----------------------------------------------------------------------
c
c
c-----In multideterminant wave funtions, all MOs very localized in a
c     given fragment must be present in all the determinants and must
c     be doubly occupied. If this is not satisfied, the MO IS NOT taken 
c     as a core orbital. The ith MO will be present in all the dete-
c     rminant with both spins if, after the following loop over over k
c     n2(i) = 2 * ndets
c
      if (ndets.gt.1) then
        n2(1:nmo) = 0
        cicoef = trim(wfnfile)//".CIcoef"
        open (ifilc,file=cicoef,access='direct',
     &      recl=RecLength*(nelact+1),form='unformatted')
        n2(1:ncore) = ndets * 2
        do k=1,ndets
          read (ifilc,rec=k) (cidet(j),j=0,nelact)
          do m=1,nelact
            imois = ncore+abs(nint(cidet(m)))
            n2(imois) = n2(imois) + 1
          enddo
        enddo
      endif

      write (lw,65) 
      mocore = 0
      core   = 0
      corea  = 0
      coreb  = 0
      nelcorea = 0
      nelcoreb = 0
      do ig=1,ng
        do i=1,nmo
          if (largwr) write (lw,66) ig,i,sg(ig,i,i)
          if (sg(ig,i,i).gt.ovc) then
            if (ndets.eq.1) then
              mocore=mocore+1
              core(ig)=core(ig)+1
              icore(mocore)=i
              icogrp(ig,core(ig))=i
              if (i.le.ncore) then
                nelcorea = nelcorea + 1
                nelcoreb = nelcoreb + 1
                corea(ig) = corea(ig) + 1
                coreb(ig) = coreb(ig) + 1
                icorea(nelcorea) = i   
                icoreb(nelcoreb) = i   
              else
c
c---------------The following lines in between 'do k=1,nelact' and enddo
c               replace the commented lines below. The latter are WRONG
c
                do k=1,nelact
                  if (abs(int(cidet(k))).eq.i) then
                    if (int(cidet(k)).gt.0) then
                      nelcorea = nelcorea + 1
                      corea(ig) = corea(ig) + 1
                      icorea(nelcorea) = i   
                    else
                      nelcoreb = nelcoreb + 1
                      coreb(ig) = coreb(ig) + 1
                      icoreb(nelcoreb) = i   
                    endif
                  endif
                enddo 
CCC           if (int(cidet(i)).gt.0) then
CCC             nelcorea = nelcorea + 1
CCC             corea(ig) = corea(ig) + 1
CCC             icorea(nelcorea) = i   
CCC           else
CCC             nelcoreb = nelcoreb + 1
CCC             coreb(ig) = coreb(ig) + 1
CCC             icoreb(nelcoreb) = i   
CCC           endif
              endif
            else
c
c-------------Check that the MO number i is in all determinants
c
              if (n2(i).eq.2*ndets) then
                mocore = mocore+1
                core(ig) = core(ig)+1
                icore(mocore) = i
                icogrp(ig,core(ig)) = i
                nelcorea = nelcorea + 1
                nelcoreb = nelcoreb + 1
                corea(ig) = corea(ig) + 1
                coreb(ig) = coreb(ig) + 1
                icorea(nelcorea) = i   
                icoreb(nelcoreb) = i   
              endif
            endif
          endif
        enddo
      enddo
c
c.....Re-compute the valence orbitals
c
      moaux = nmo-mocore
      moval = 0
      do i=1,nmo
        j=1
        iscore=.false.
        do while (j.le.mocore.and.(.not.iscore))
          if (i.eq.icore(j)) iscore=.true.
          j=j+1
        enddo
        if (.not.iscore) then
          moval=moval+1
          ival(moval)=i
        endif
      enddo
      if (moval.ne.moaux) then
        write (0,*)
        write (0,*) '# coremos.f: Wrong COREMO order !!!'
        write (0,*)
        stop
      endif
c
      write (lw,223) ovc
      totcore=0
      do ig=1,ng
        write (lw,224) ig,core(ig),(icogrp(ig,k),k=1,core(ig))
        totcore=totcore+core(ig)
        do k=1,core(ig)
           ick=icogrp(ig,k)
           write (lw,226) ick,ick,ig,sg(ig,ick,ick)
        enddo
      enddo
      ok=.false.
      if (totcore.le.lline/4-1) then
        ok=.true.
        line(1:4) = '    '
        iw=5
        do ig=1,ng
          if (core(ig).gt.0) then
            do k=1,core(ig)
             line(iw:iw+4) = fourchar(icogrp(ig,k))
             do j=iw,iw+4
               if (line(j:j).eq.'0') then
                 line(j:j)=' '
               else
                 exit
               endif
             enddo
             iw=iw+4
            enddo
          endif
        enddo
      endif
c
c.....Write list of CORE and VALENCE orbitals.
c
      if (mocore.gt.0) then
        write (lw,29) mocore,(icore(i),i=1,mocore)
      else
        write (lw,29) mocore
      endif
      if (moval.gt.0) then
         write (lw,31) moval,(ival(i),i=1,moval)
      else
         write (lw,31) moval
      endif
      return
      stop
 223  format (1x,'#',/,
     & 1x,'#',/,1x,'# Automatic analysis of CORE MOs',/,
     & 1x,'# Overlap criterium of Localization = ',F16.10)
 224  format (1x,'# FRAGMENT ',I3,' HAS ',I3,
     & ' MOs ALMOST FULLY LOCALIZED ON IT: ',/,100(1x,'# ---> ',20I3,/))
 226  format (1x,'# <MO_',I3,'|MO_',I3,'>_',I2,' = ',F16.10)
 29   format (//,' # THERE ARE ',I3,' CORE    ORBITALS :',/,
     &   1000(' # ',20(1x,I3),/))
 31   format (//,' # THERE ARE ',I3,' VALENCE ORBITALS :',/,
     &   1000(' # ',20(1x,I3),/))
 66   format (' # Group ',I2,'  MO ',I4,'  Self-Overlap = ',F15.8)
 65   format (//,'# Automatic analysis of CORE MOs')
      end
